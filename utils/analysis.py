### Perform data analysis ###

## Using objects created previously: metadata, df_protAbundance
import os
import subprocess
import json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy.spatial.distance import pdist, squareform
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from gprofiler import GProfiler
from utils.check_env import get_repo_root
import itertools
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from glob import glob
from PIL import Image
from pathlib import Path


##### there is a warning to suppress when fitting models
warnings.filterwarnings("ignore", message="Negative binomial dispersion parameter alpha not set. Using default value alpha=1.0.")
# Suppress only ConvergenceWarning
warnings.simplefilter("ignore", ConvergenceWarning)

def generate_pca(df: pd.DataFrame,
                 output_dir: str,
                 metadata: pd.DataFrame = None,
                 components: int = 2) -> pd.DataFrame:
    """
    Perform PCA on the given DataFrame and save the results as a plot.

    Args:
        df (pd.DataFrame): protein abundance data.
        output_dir (str): Directory to save the PCA plot and data.
        metadata (pd.DF): sample metadata
        components (int): Number of principal components to compute.

    Returns:
        pd.DataFrame: DataFrame with principal components and metadata.
    """
    # Perform PCA
    pca = PCA(n_components=components)
    principal_components = pca.fit_transform(df)
    n_prot = df.shape[1]
    plot_title = 'PCA Plot (n proteins = ' + str(n_prot) + ')'
    # Create DataFrame with PCA results
    df_PCA = pd.DataFrame(
        data=principal_components,
        columns=[f'principal component {i + 1}' for i in range(components)],
        index = pd.DataFrame(df).index
    )
    # pull over treatment from metadata
    # df_PCA . look up by index . and map to
    # metadata . temporarily set index to sample_id [ and find treatment ]
    df_PCA['treatment'] = df_PCA.index.map(metadata.set_index('sample_rep')['treatment'])
    # Calculate variance explained by each component
    explained_variance = pca.explained_variance_ratio_
    pc_labels = [
        f'PC{i + 1} ({round(var * 100, 2)}%)'
        for i, var in enumerate(explained_variance[:components])
    ]
    # Plot PCA
    plt.figure(figsize=(4, 3))
    plot_pca = sns.scatterplot(
        data=df_PCA,
        x='principal component 1',
        y='principal component 2',
        hue='treatment'
    )
    # Annotate points with sample IDs
    for i in range(df_PCA.shape[0]):
        plt.text(
            x=df_PCA['principal component 1'][i] + 0.1, 
            y=df_PCA['principal component 2'][i] + 0.1,
            s=df_PCA.index[i],  # Sample ID
            fontsize=9,
            ha='center',
            va='bottom'
        )
    # Set labels and title
    plot_pca.set(xlabel=pc_labels[0], ylabel=pc_labels[1])
    plt.title(plot_title)
    plt.tight_layout()
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, 'plots', 'pca_plot.png')
    data_path = os.path.join(output_dir, 'data', 'pca_data.csv')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    df_PCA.to_csv(data_path, index=False)
    print(f"PCA plot saved to {plot_path}")
    print(f"PCA data saved to {data_path}")
    return df_PCA

def generate_MDS(df: pd.DataFrame,
                 output_dir: str,
                 metadata: pd.DataFrame = None) -> pd.DataFrame:
    """
    Perform MDS on the given DataFrame and save the results as a plot.

    Args:
        df (pd.DataFrame): protein abundance data.
        output_dir (str): Directory to save the plot and data.
        metadata (pd.DF): sample metadata

    Returns:
        pd.DataFrame: DataFrame with MDS and metadata.
    """
    # Perform MDS
    # pdist(x) computes the Euclidean distances between each pair of points in an array
    dissimilarity_array = pdist(df, metric='euclidean')
    n_prot = df.shape[1]
    plot_title = 'MDS Dissimilarity based on Protein Abundance (n proteins = ' + str(n_prot) + ')'
    # suqareform() returns the matrix
    dissimilarity_matrix = squareform(dissimilarity_array)
    # Perform NMDS
    mds = MDS(n_components=2, # n_components is number of dimensions
            dissimilarity="precomputed", # we calculate above
            random_state=42,
            metric = True ) 
    # Fit the data  and return the embedded coordinates.
    mds_coords = mds.fit_transform(dissimilarity_matrix) # If dissimilarity=='precomputed', the input should be the dissimilarity matrix.
    # Convert the NMDS coordinates (nmds_coords) to a pandas DataFrame
    mds_coords_df = pd.DataFrame(mds_coords,
                                 columns=['MDS1', 'MDS2'],
                                 index = pd.DataFrame(df).index)
    # pull over treatment from metadata
    # df_PCA . look up by index . and map to
    # metadata . temporarily set index to sample_id [ and find treatment ]
    mds_coords_df['treatment'] = mds_coords_df.index.map(metadata.set_index('sample_rep')['treatment'])
    #### Generate NMDS plot
    plot_nmds = sns.scatterplot(data = mds_coords_df, x = 'MDS1', y = 'MDS2',
                    hue = 'treatment')
    ## seaborn returns an axis-object rather than a figure, so you can still alter features of it. E.g. axes names:
    plot_nmds.set(xlabel = 'MDS 1', ylabel = 'MDS 2')
    # Add sample IDs from the 'target' column with an offset
    for i in range(mds_coords_df.shape[0]):
        plt.text(
            x=mds_coords_df['MDS1'][i] + 1,  # Add a small x-offset
            y=mds_coords_df['MDS2'][i] + 1,  # Add a small y-offset
            s=mds_coords_df.index[i],
            fontsize=9,
            ha='center',  # Horizontal alignment
            va='bottom'   # Vertical alignment
        )
    # Set labels and title
    plt.title(plot_title)
    plt.tight_layout()
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, 'plots', 'mds_plot.png')
    data_path = os.path.join(output_dir, 'data', 'mds_data.csv')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    mds_coords_df.to_csv(data_path, index=False)
    print(f"MDS plot saved to {plot_path}")
    print(f"MDS data saved to {data_path}")
    return mds_coords_df


def generate_heatmap(df: pd.DataFrame,
                 output_dir: str) -> pd.DataFrame:
    """
    Generate a heatmap using the abundance data

    Args:
        df (pd.DataFrame): protein abundance data.
        output_dir (str): Directory to save the PCA plot and data.

    Returns:
        
    """
    #### need the other orientation for heatmap
    # transpose and add sample ids as colnames
    df_heat = pd.DataFrame(df.T)
    df_heat.columns = df.index
    n_prot = df_heat.shape[0]
    plot_title = 'heatmap of euclidean distance (n protein = ' + str(n_prot) + ')'
    # Create clustermap
    heatmap = sns.clustermap(df_heat, 
                    #annot=True, 
                    fmt=".2f",
                    cmap="viridis")
    # Move x-axis labels (column names) to the top
    heatmap.ax_heatmap.xaxis.set_ticks_position("top")  # Move ticks to the top
    heatmap.ax_heatmap.xaxis.set_label_position("top")  # Move axis label to the top
    # Remove dendrogram tick labels
    heatmap.ax_heatmap.set_xticks(np.arange(len(df_heat.columns)))  # Set tick positions for columns
    heatmap.ax_heatmap.set_xticklabels(df_heat.columns, rotation=45, ha='left')  # Use only column names
    # Explicitly disable extra tick labels from the dendrogram
    heatmap.ax_heatmap.tick_params(axis='x', which='both', bottom=False, top=True)
    # Set labels and title
    plt.title(plot_title)
    # plt.tight_layout()
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, 'plots', 'heatmap_plot.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"heatmap saved to {plot_path}")
    return df_heat

def make_volcano(df_pair: pd.DataFrame,
                 output_dir: str,
                 pair_name: str,
                 config: dict,
                 metadata_pair: pd.DataFrame = None
                 ) -> pd.DataFrame:
    """
    fit model for each protein. Calls limma script

    Args:
        df (pd.DataFrame): protein abundance data.
        output_dir (str): Directory to save the output
        pair_name (str): name of treatments being compared (pairwise)
        metadata (pd.DF): sample metadata
        pair_name (str): treatment groups to be compared, separated by an underscore

    Returns:
        df_model_out (pd.DataFrame): gene name, F statistic, p value, FDR corrected p value.
    """
    #### calculate DE using limma (R package) ####
    pair_data_path = os.path.join(output_dir, 'data', pair_name, 'prots.csv')
    df_pair.to_csv(pair_data_path, index = True)
    # Save sample metadata
    pair_metadata_path = os.path.join(output_dir, 'data', pair_name, 'metadata.csv')
    metadata_pair.to_csv(pair_metadata_path, index=False)
    # Define where to save the limma results
    pair_result_path = os.path.join(output_dir, 'data', pair_name, 'limmaOut.csv')
    # run R script - note: r-limma-env conda env required
    subprocess.run([
    "conda", "run", "-n", "r-limma-env",
    "Rscript", "utils/DE-limma.R",
    pair_data_path.replace("\\", "/"),
    pair_metadata_path.replace("\\", "/"),
    pair_result_path.replace("\\", "/")
    ], check=True)
    # read results back in
    diffExpr_df = pd.read_csv(pair_result_path, index_col=0)
    n_prot = diffExpr_df.shape[0]
    diffExpr_df['Log10_FDR_P_Value'] = -np.log10(diffExpr_df['adj.P.Val'])
    diffExpr_df['Log10_unadjusted_p_Value'] = -np.log10(diffExpr_df['P.Value'])
    ### whether to plot the -log10(p_value) i.e. unadjusted or -log10(FDR_p_value) is specified in json field "LFC_plot_p_or_FDRp"
    Volcano_y_axis = config.get("LFC_plot_p_or_FDRp")
    Volcano_y_data = diffExpr_df[Volcano_y_axis]
    # Add the Colour column based on LOG2FC and p_values_FDR
    diffExpr_df['Colour'] = diffExpr_df.apply(
        lambda row: 'blue' if (abs(row['logFC']) > 2 and row['adj.P.Val'] < 0.05) else 'gray', axis=1
    )
    diffExpr_path = os.path.join(output_dir, 'data', pair_name, 'limma_output.csv')
    diffExpr_df.to_csv(diffExpr_path, index=True)
    ### Create volcano plot
    plot_title = 'Protein Abundance Log Fold Change for treatments \n' + pair_name + '(n = ' + str(n_prot) + ')'    
    sns.scatterplot(
        data=diffExpr_df, 
        x='logFC',  # Log fold change
        y=Volcano_y_data,  # -log10(FDR-corrected p-value)
        hue='Colour',  # Color based on significance
        palette={'blue': 'blue', 'gray': 'gray'},
        legend=False,  # No legend for this plot
        alpha=0.7  # Transparency for points
    )
    # Customize the plot
    plt.axvline(x=2, color='red', linestyle='--', linewidth=1)
    plt.axvline(x=-2, color='red', linestyle='--', linewidth=1) 
    plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=1) 
    plt.title(plot_title, fontsize=16)
    plt.xlabel('Log2 Fold Change (LFC)', fontsize=12)
    plt.ylabel(Volcano_y_axis, fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    # Set symmetrical x-axis
    # Save plot and PCA data
    #plt.xlim(-4, 4)
    plot_path = os.path.join(output_dir, 'plots', pair_name, 'volcano_plot.png')
    if os.name == 'nt':
        plot_path = '\\\\?\\' + os.path.abspath(plot_path)
    plt.savefig(plot_path, dpi=300)
    plt.close()
    # save the top 20 rows to csv for display in final report
    # Sort by a relevant column (modify column name as needed)
    sorted_df = diffExpr_df.sort_values(by="logFC", ascending=False, key = abs)
    # Select the top 20 rows
    top_20_df = sorted_df.head(20).round(decimals = 2)
    top_20_df = top_20_df.round({"P.Value":4, "adj.P.Value":4})
    # Save to CSV
    top_prot_path = os.path.join(output_dir, 'data', pair_name, 'top_20_by_LFC.csv')
    if os.name == 'nt':
        top_prot_path = '\\\\?\\' + os.path.abspath(top_prot_path)
    top_20_df.to_csv(top_prot_path, index=True)
    return diffExpr_df, top_20_df

def plot_abundance(df: pd.DataFrame,
                   top_20_df: pd.DataFrame,
                   metadata: pd.DataFrame,
                 output_dir: str) :
    """
    generate box plots of abundances of most differentially abundant proteins

    Args:
        df (pd.DataFrame): df with raw abundance
        top_20_df (pd.DataFrame):  abundance data for top 20 proteins by LFC
        metadata (pd.DataFrame): project metadata
        output_dir (str): Directory to save the output.

    Returns:
    """
    top_10_df = top_20_df.head(10)
    top_10_df = top_10_df.merge(
        df,
        left_on='Gene',
        right_index=True,
        how = 'left'
    ).set_index( 'Gene' )
    top_10_df = top_10_df.loc[:, top_10_df.columns.isin(metadata['sample_rep'])].T
    top_10_df['treatment'] = top_10_df.index.map(metadata.set_index('sample_rep')['treatment'])
    # rename index
    top_10_df.index.name = 'sample_rep'
    if top_10_df.isna().any().any():
        raise ValueError(f"Error: Missing data in protein abundances. Expected only complete data")
    # Create the panel of 10 plots with gene names added
    fig, axes = plt.subplots(2, 5, figsize=(20, 10), sharey=False)  # Set sharey=False for varying scales
    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    treatment_levels = top_10_df['treatment'].unique()

    ##### Generate plots for each protein ######
    for i, protein in enumerate(top_10_df.columns.drop('treatment')):
        # Create the data for this protein
        plot_data = pd.DataFrame({
            'Abundance': top_10_df[protein],
            'Treatment': top_10_df['treatment']
        })
        # Create the boxplot with points overlaid
        sns.boxplot(
            data=plot_data,
            x='Treatment',
            y='Abundance',
            ax=axes[i],
            hue='Treatment',
            legend=False,
            palette={treatment_levels[0]: 'blue', treatment_levels[1]: 'orange'},
            showmeans=True,
            meanline=True
        )
        sns.stripplot(
            data=plot_data,
            x='Treatment',
            y='Abundance',
            color='black',
            alpha=0.7,
            jitter=True,
            ax=axes[i]
        )
        # Customize each subplot
        axes[i].set_title(f'{protein} Abundance', fontsize=12)
        axes[i].set_xlabel('Treatment', fontsize=10)
        axes[i].set_ylabel('Abundance', fontsize=10)
    # Adjust layout to avoid overlap and show the plot
    plt.tight_layout()
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, 'plots', 'abundance_top10_plot.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()

def plot_venn(anova_lm_df: pd.DataFrame,
                        output_dir: str) :
    """
    generate box plots of abundances of most differentially abundant proteins

    Args:
        df (pd.DataFrame): df with raw abundance
        top_20_df (pd.DataFrame):  abundance data for top 20 proteins by LFC
        metadata (pd.DataFrame): project metadata
        output_dir (str): Directory to save the output.

    Returns:
    """

def enrichment_analysis(anova_lm_df: pd.DataFrame,
                        pair_name: str,
                        output_dir: str) :
    """
    run enrichment analysis to identify overrepresented genes, functions and pathways

    Args:
        anova_lm_df (pd.DataFrame): df with anova outputs including LFC, p value and FDR adjusted p value
        pair_name (str): name of treatments being compared (pairwise)
        output_dir (str): Directory to save the output.

    Returns:
        enrichment data (pd.DataFrame)
    """    
    ##### Calculate enrichment #####
    gp = GProfiler(return_dataframe=True)
    ##### G Profiler options #####
    # for ORA, just need a list of genes
    pathway_query_genes = anova_lm_df.loc[
    (anova_lm_df['adj.P.Val'] < 0.05)
    ].index
    # in the case of phosphoproteomic data, gene names have a double __ with phosphorylation state added,
    # for now, we remove the phospho data from this set of genes
    # may want to look at separately later
    if any(isinstance(gene, str) and '__' in gene for gene in pathway_query_genes):
        pathway_query_genes = [str(gene).split("__")[0] for gene in pathway_query_genes]
    # pathway database can be REAC, GO or KEGG. Also less common but available: CORUM, HPA, TF and MIRNA
    # defaults to REAC
    source=['GO']
    plot_title = 'Pathway Enrichment (' + source[0] + ')' 
    # p value threshold defaults to 0.05
    p_threshold=0.05
    # all results returns all results, not just those below p threshold
    all_results = False
    # a background set can be specified e.g. background=["BRCA1", "TP53", "AKT1", "MTOR", "EGFR", "MYC"]
    # multiple testing correction can be g_SCS (default, Set Counts and Sizes), bonferroni, or fdr
    # from quick look, g_SCS seems to be similar to bonferroni, which are both less strict than fdr.
    significance_threshold_method='g_SCS'
    ##### Running pathway enrichment #####
    if len(pathway_query_genes) > 0:
        pathway_result = gp.profile(organism='hsapiens',
            query=list(pathway_query_genes),
            sources=source,
            user_threshold = p_threshold,
            significance_threshold_method = significance_threshold_method,
            all_results = all_results
        )  # REAC for Reactome
        ### save results to file
        enrichment_path = os.path.join(output_dir, 'data', pair_name, 'pathway_enrichment.csv')
        ## windows sometimes rejects long paths. Workaround:
        if os.name == 'nt':
            enrichment_path = '\\\\?\\' + os.path.abspath(enrichment_path)
        pathway_result.round({'precision': 2, 'recall':2}).to_csv(enrichment_path, index=False)
        ##### Plot enrichment #####
        pathway_plot_df = pathway_result.sort_values('p_value', ascending = True).head(20)
        pathway_plot_df['-log10(p_value)'] = -np.log10(pathway_plot_df['p_value'])
        # Setting up the plot using Seaborn and Matplotlib adjustments
        # plt.figure(figsize=(15, 20))
        #ax = plt.gca()
        sns.relplot(
            data = pathway_plot_df,
            x = "-log10(p_value)",
            y = "name",
            size = "recall", # the proportion of query genes associated with the term
            color = "green",
            height=15,       # Adjust figure height for better fit
            aspect=0.5       # Maintain a suitable aspect ratio)
        )
        # Adjustments for axes padding and limits
        # Adjustments for labels and titles
        plt.xlabel('-log10(p-value)', fontsize=14)
        plt.ylabel('Pathway Names', fontsize=14)
        plt.title(plot_title, fontsize=16)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)  # Ensure y-axis labels are readable
        # Adjust layout to ensure labels are fully visible
        plt.tight_layout()
        # Save plot data
        plot_path = os.path.join(output_dir, 'plots', pair_name, 'pathway_enrichment_plot.png')
        ## windows sometimes rejects long paths. Workaround:
        if os.name == 'nt':
            plot_path = '\\\\?\\' + os.path.abspath(plot_path)
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")  # bbox_inches ensures labels aren't cut off
        plt.close()
        return pathway_result
    
#### combine plots made for different treatment groups
def combine_plots(search_term,
                  search_path,
                  output_dir, 
                  output_filename = None,
                  img_size=(800, 600),
                  max_cols=3):
    """
    Finds all images matching the search_term in subdirectories of search_path,
    arranges them in a grid with up to max_cols columns, and saves a single PNG.

    Parameters:
    - search_term (str): The filename pattern to search for (e.g., "volcano_plot.png").
    - search_path (str): The base directory where images are stored.
    - output_filename (str): The filename for the combined image (auto-generated if None).
    - img_size (tuple): (width, height) to resize images.
    - max_cols (int): Maximum number of columns in the grid.

    Returns:
    - str: Path to the saved combined image, or None if no images found.
    """
    if os.name == 'nt':
        search_path = '\\\\?\\' + os.path.abspath(search_path)
    # Find all matching images
    image_paths = []
    for root, dirs, files in os.walk(search_path):
        for file in files:
            if search_term in file:
                image_paths.append(os.path.join(root, file))
    image_paths = [img for img in image_paths if 'combined_volcano_plot.png' not in img]
    # image_paths = sorted(glob(os.path.join(search_path, "**", search_term), recursive=True))
    if not image_paths:
        print(f"No plots found for '{search_term}'.")
        return None
    # Load and resize images
    images = [Image.open(img).resize(img_size, Image.LANCZOS) for img in image_paths]
    # Determine grid layout
    cols = min(max_cols, len(images))
    rows = len(images) // cols
    rows = (len(images) + cols - 1) // cols  # Round up to fit all images by adding cols - 1
    # Create a blank canvas
    combined_width = cols * img_size[0]
    combined_height = rows * img_size[1]
    combined_image = Image.new("RGB", (combined_width, combined_height), (255, 255, 255)) # 255,255,255 specifies background colour = white
    # Paste images into grid
    for idx, img in enumerate(images):
        x_offset = (idx % cols) * img_size[0]
        y_offset = (idx // cols) * img_size[1]
        combined_image.paste(img, (x_offset, y_offset))
    # Generate output filename if not provided
    if output_filename is None:
        output_filename = os.path.join(output_dir, f"plots/combined_{search_term.replace('.png', '')}.png")
        ## windows sometimes rejects long paths. Workaround:
        if os.name == 'nt':
            output_filename = '\\\\?\\' + os.path.abspath(output_filename)
    # Save the final image
    combined_image.save(output_filename)
    print(f"Combined plot saved to: {output_filename}")

### for combining data from different treatments for display in the report ###
def combine_csv_files(filename,
                      output_dir,
                      output_filename=None,
                      top_n=10,
                      new_column="treatment_pair"):
    """
    General function to combine CSV files from subdirectories into a single file.

    Parameters:
    - filename (str): The name of the CSV file to search for (e.g., "top_20_by_LFC.csv").
    - output_dir (str): The root directory where data folders are stored.
    - output_filename (str or None): The output filename for the combined CSV.
                                     If None, it's auto-generated based on `filename`.
    - top_n (int): Number of rows to take from each CSV file.
    - new_column (str): Column name to store the extracted folder name (e.g., "treatment_pair").

    Returns:
    - pd.DataFrame: The combined DataFrame.
    - str: The path where the final CSV is saved.
    """
    # Find all matching images
    csv_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file == filename:
                csv_files.append(os.path.join(root, file))
    # # Search for matching CSV files in subdirectories
    # search_pattern = os.path.join(output_dir, "data", "**", filename)
    # csv_files = glob(os.path.join(output_dir, "data", "*", filename)) + glob(os.path.join(output_dir, "data", "*", "*", filename))
    # csv_files = sorted(glob(search_pattern, recursive=True))
    # check if csv files exist
    if not csv_files:
        print(f"No files found matching '{filename}'.")
        return None, None
    combined_data = []
    # Process each found CSV file
    for file in csv_files:
        if os.name == 'nt':
            file = '\\\\?\\' + os.path.abspath(file)
        # Extract the folder name (used as the category column)
        folder_name = os.path.basename(os.path.dirname(file))  
        # Read CSV and select top `n` rows
        df = pd.read_csv(file).head(top_n)
        # Add the extracted folder name as a new column
        df[new_column] = folder_name  
        # Append to the list
        combined_data.append(df)
    # Merge all data
    combined_df = pd.concat(combined_data, ignore_index=True)
    # Generate output filename if not provided
    if output_filename is None:
        output_filename = os.path.join(output_dir, f"data/combined_{filename}")
        ## windows sometimes rejects long paths. Workaround:
        if os.name == 'nt':
            output_filename = '\\\\?\\' + os.path.abspath(output_filename)
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    # Save combined CSV
    combined_df.to_csv(output_filename, index=False)
    print(f"Combined file saved at: {output_filename}")


##########################################################################
#### This is the main function for analysing data, combining the above ###
##########################################################################

def run_analysis(df: pd.DataFrame,
                 metadata: pd.DataFrame,
                 output_dir: str,
                 config: dict, 
                 json_out: str) -> dict:
    """
    Full analysis pipeline: performs PCA,

    Parameters:
        df (pd.DataFrame): Raw protein abundance data.
        metadata (pd.DataFrame): Metadata containing sample information.
        output_dir (str): Directory to save analysis outputs.
        json_out (str): File for saving information to go into the final report

    Returns:
        dict: Dictionary containing results from all analyses.
    """

    # Initialize results dictionary
    results = {}

    ##### Analyses for all treatment groups #####

    # Perform PCA and save results
    print("Performing PCA...")
    pca_results = generate_pca(df.T, output_dir, metadata=metadata)
    results['pca'] = pca_results

    # Perform MDS and save results
    print("Performing MDS...")
    mds_coords_df = generate_MDS(df.T, output_dir, metadata=metadata)
    results['mds'] = mds_coords_df

    # Generate and save heatmap
    print("Generating heatmap...")
    df_heatmap = generate_heatmap(df.T, output_dir)
    results['heatmap'] = df_heatmap

    # Look at overlaps of differnential abundance

    #### if 2 or 3 groups, Venn

    #### if 4 or more, upset

    ###### Pairwise Analyses #####
    # if there are > 2 treatment groups, pairwise analyses will have to be run separately for each pair of treatments
    treatment_pairs = list(itertools.combinations(metadata['treatment'].unique(), 2))

    for pair in treatment_pairs:
        print("starting analysis for pair ", pair)
        metadata_pair = metadata[metadata['treatment'].isin(pair)]
        df_pair = df[metadata_pair['sample_rep'].tolist()]
        pair_name = "_".join(map(str, pair) )
        if not os.path.exists(os.path.join(output_dir, 'plots', pair_name)):
            os.mkdir(os.path.join(output_dir, 'plots', pair_name))
        if not os.path.exists(os.path.join(output_dir, 'data', pair_name)):
            os.mkdir(os.path.join(output_dir, 'data', pair_name))
        # Generate and save volcano plot
        print("Generating volcano plot for pair ", pair_name, "...")
        anova_lm_df, top_20_df = make_volcano(df_pair,
                                              output_dir,
                                              pair_name = pair_name,
                                              config = config,
                                              metadata_pair = metadata_pair)
        results_name = 'df_lm_' + pair_name
        results['results_name'] = anova_lm_df
        # Find overrepresented pathways and save output
        print("Running enrichment analysis for pair ", pair_name, "...")
        enrichment_analysis(anova_lm_df,
                            pair_name,
                            output_dir
                            )

    # combine plots from different pairs
    combine_plots(search_path = output_dir,
                  search_term = "volcano_plot.png",
                  output_dir=output_dir) 
    
    combine_plots(search_path = output_dir,
                  search_term = "pathway_enrichment_plot.png",
                  output_dir=output_dir) 
    
    # Combine the top 10 most differentially abundant proteins
    combine_csv_files(filename="top_20_by_LFC.csv",
                      output_dir=output_dir)
    # Combine pathway enrichment data
    combine_csv_files(filename="pathway_enrichment.csv",
                      output_dir=output_dir)

    ### write to file the version of this script
    REPO_ROOT = get_repo_root()
    analysis_version =  subprocess.check_output(["git", "log", "-n", "1", "--format=%H", "--", 
                                                 os.path.join(REPO_ROOT, 'utils/analysis.py'),
                                                 ]).strip().decode("utf-8")

    # read data from json file
    with open(json_out, "r") as f:
        existing_data = json.load(f)

    analysis_meta = {"ANALYSIS_VERSION": analysis_version}

    # Append new data
    existing_data.update(analysis_meta)

    # Write back to JSON file
    with open(json_out, "w") as f:
        json.dump(existing_data, f, indent=4)
    
    print("Analysis pipeline completed.")
    return results
