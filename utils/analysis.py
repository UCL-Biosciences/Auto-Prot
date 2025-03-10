### Perform data analysis ###

## Using objects created previously: metadata, df_protAbundance, df_protAbundance_standardised

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
from matplotlib import patches



##### there is a warning to suppress when fitting models
warnings.filterwarnings("ignore", message="Negative binomial dispersion parameter alpha not set. Using default value alpha=1.0.")
# Suppress only ConvergenceWarning
warnings.simplefilter("ignore", ConvergenceWarning)

def generate_histogram(df: pd.DataFrame,
                        output_dir: str,
                        metadata: pd.DataFrame,
                        json_out: dict): 
    """
    Generate histograms showing distributions for the different treatments

    Args:
        df (pd.DataFrame): protein abundance data.
        output_dir (str): Directory to save the output
        metadata (pd.DF): sample metadata
        config (dict): configuration info for the pipeline. we add mean and sd per treatment, which can go into report
        json_out (str): File for saving information to go into the final repo

    Returns:
        
    """
    ### to create histograms, need to re format
    # need cols: sample replicate, treatment, abundance
    df_long = df.melt(var_name='Sample', value_name='Abundance', ignore_index=False)
    df_long = df_long.reset_index().rename(columns={'index': 'Protein'})  # Ensure 'Protein' column exists
    df_long = df_long.merge(metadata[['sample_rep', 'treatment']], left_on='Sample', right_on='sample_rep')
    # Log-transform the abundance values to make it easier to see
    df_long['Log_Abundance'] = np.log1p(df_long['Abundance'])
    # Create the figure
    plt.figure(figsize=(10, 6))
    ax = sns.kdeplot(data=df_long, x='Log_Abundance', hue='treatment', common_norm=False, linewidth=2)
    # Force legend display
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles=handles, labels=labels, title="Treatment")
    # Set labels and title
    plt.xlabel("Log(Protein Abundance + 1)")
    plt.ylabel("Density")
    plt.title("Density of Log-Transformed Protein Abundance by Treatment")
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, 'plots', 'histogram_all_treatments_plot.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    ### for the full dataset, we want to record the mean and SD for abundance for each treatment ###
    if "full_dataset" in output_dir:
        # Compute mean Abundance per treatment
        mean_abundance = df_long.groupby('treatment')['Abundance'].mean().round().astype('int')
        mean_for_json = ", ".join([f"{mean}: {int(n):,}" for mean, n in mean_abundance.items()])
        sd_abundance = df_long.groupby('treatment')['Abundance'].std().round().astype('int')
        sd_for_json = ", ".join([f"{sd}: {int(n):,}" for sd, n in sd_abundance.items()])
        # read data from json file
        with open(json_out, "r") as f:
            existing_data = json.load(f)
        abundance_stats = {
                "MEAN_ABUNDANCE": mean_for_json,
                "SD_ABUNDANCE": sd_for_json
        }
        # Append new data
        existing_data.update(abundance_stats)
        # Write back to JSON file
        with open(json_out, "w") as f:
            json.dump(existing_data, f, indent=4)


def generate_pca(df_standardised: pd.DataFrame,
                 output_dir: str,
                 metadata: pd.DataFrame = None,
                 components: int = 2) -> pd.DataFrame:
    """
    Perform PCA on the given DataFrame and save the results as a plot.

    Args:
        df_standardised (pd.DataFrame): Standardised protein abundance data.
        output_dir (str): Directory to save the PCA plot and data.
        metadata (pd.DF): sample metadata
        components (int): Number of principal components to compute.

    Returns:
        pd.DataFrame: DataFrame with principal components and metadata.
    """
    # Perform PCA
    pca = PCA(n_components=components)
    principal_components = pca.fit_transform(df_standardised)
    n_prot = df_standardised.shape[1]
    plot_title = 'PCA Plot (n proteins = ' + str(n_prot) + ')'
    # Create DataFrame with PCA results
    df_PCA = pd.DataFrame(
        data=principal_components,
        columns=[f'principal component {i + 1}' for i in range(components)],
        index = pd.DataFrame(df_standardised).index
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

def generate_MDS(df_standardised: pd.DataFrame,
                 output_dir: str,
                 metadata: pd.DataFrame = None) -> pd.DataFrame:
    """
    Perform MDS on the given DataFrame and save the results as a plot.

    Args:
        df_standardised (pd.DataFrame): Standardised protein abundance data.
        output_dir (str): Directory to save the plot and data.
        metadata (pd.DF): sample metadata

    Returns:
        pd.DataFrame: DataFrame with MDS and metadata.
    """
    # Perform MDS
    # pdist(x) computes the Euclidean distances between each pair of points in an array
    dissimilarity_array = pdist(df_standardised, metric='euclidean')
    n_prot = df_standardised.shape[1]
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
                                 index = pd.DataFrame(df_standardised).index)
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


def generate_heatmap(df_standardised: pd.DataFrame,
                 output_dir: str) -> pd.DataFrame:
    """
    Generate a heatmap using the standardised data

    Args:
        df_standardised (pd.DataFrame): Standardised protein abundance data.
        output_dir (str): Directory to save the PCA plot and data.

    Returns:
        
    """
    #### need the other orientation for heatmap
    # transpose and add sample ids as colnames
    df_heat = pd.DataFrame(df_standardised.T)
    df_heat.columns = df_standardised.index
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

# Function to run ANOVA for each protein
def run_anova(row, metadata):
    """
    fit linear models and ANOVA models to assess differences between treatment groups. Called by make_volcano

    Args:
        row : iterating through rows with df.apply()
        metadata (pd.DF): sample metadata
        
    Returns:
        
    """
    gene_name = row.name  # Use row index as gene name
    data = pd.DataFrame(row)  # Convert row to DataFrame
    data.columns=['Abundance'] # rename column so can be used in model
    # Map sample_id to treatment
    data["treatment"] = data.index.map(metadata.set_index("sample_rep")["treatment"])
    ### count data are usually modelled with a generalised linear model with a negative binomial error distribution
    ### However, abundance data often violate the assumption of poisson GLMs that mean = variance
    ### To deal with this "overdispersion", we use a generalised linear model with a negative binomial error distribution 
    ### this won't work for cases where all values are 0 as we can't estimate deviance (or deviance returns NaN)
    ### exclude these for now.
    if (data["Abundance"] == 0).all():
        print("All values are zero. NB GLM cannot be fitted. For gene ", gene_name)
    if (data["Abundance"] != 0).any():
        model = smf.negativebinomial('Abundance ~ C(treatment)', data=data).fit(disp=0)
        ## some times these models fail to estimate parameters. Only continue when models have been fitted correctly
        if not np.any(np.isnan(model.bse)):
            # For generalISED linear models, we can't run anova tests to look at how much variance is explained by the different parameters.
            # instead, we use the Wald statistic. Wald stat is the coefficient divided by the standard error.
            # This is good because the wald statistic will be higher - indicating stronger effect - when there is less variation (less noise) and/or a larger sample size - both contribute to stronger evidence
            # larger wald stat = stronger evidence. Wald test stat is the chi2 value
            wald_test = model.wald_test_terms(scalar = True)
            wald_stat, p_value = wald_test.table.loc["C(treatment)", "statistic"], wald_test.table.loc["C(treatment)", "pvalue"]
            #### Extract group means and LFC
            # Compute group means (sorted alphabetically)
            group_means = data.groupby("treatment")["Abundance"].mean().sort_index()
            grp1_name = group_means.index.tolist()[0]
            grp1_mean = group_means.iloc[0] + 1 ## add small constant for LFC
            grp2_name = group_means.index.tolist()[1]
            grp2_mean = group_means.iloc[1] + 1 ## add small constant for LFC
            fc = ( grp1_mean ) / ( grp2_mean )
            log2fc = np.log2( fc )
            #return outputs
            results_anova = pd.Series({"Gene": gene_name,
                            "Group_1": grp1_name,
                            "Group_1_mean": grp1_mean,
                            "Group_2": grp2_name,
                            "Group_2_mean": grp2_mean,
                            "Raw_Fold_Change": fc,
                            "Log2_Fold_Change": log2fc,
                            "Wald_stat": wald_stat,
                            "p_value": p_value})
            return results_anova


def volcano_plot(anova_lm_df,
                  config,
                  plot_title
                  ):
    """
    Adds an inset zoomed-in plot to a volcano plot if necessary.
    """
    # Define y-axis cutoff threshold
    fdr_cut_off = config.get("volcano_y_cutoff") # note this is in terms of FDR for readability, converted to -log10() here:
    y_cutoff = -np.log10(fdr_cut_off)
    # Determine if an inset is needed (i.e., if any points exceed the y-axis limit)
    max_y_value = np.max(anova_lm_df['Log10_FDR_P_Value'])
    truncate = max_y_value > y_cutoff
    ### whether to plot the -log10(p_value) i.e. unadjusted or -log10(FDR_p_value) is specified in json field "LFC_plot_p_or_FDRp"
    Volcano_y_axis = config.get("LFC_plot_p_or_FDRp")
    Volcano_y_data = anova_lm_df[Volcano_y_axis]
    LFC_threshold = config.get("LFC_threshold")
    # if max y NOT above the cutoff, do normal plot
    if not truncate:
        # Create the figure and axis for plotting
        fig, ax = plt.subplots(figsize=(8, 6))
        ### Create volcano plot
        sns.scatterplot(
            data=anova_lm_df, 
            x='Log2_Fold_Change',  # Log fold change
            y=Volcano_y_data,  # -log10(FDR-corrected p-value)
            hue='Colour',  # Color based on significance
            palette={'blue': 'blue', 'gray': 'gray'},
            legend=False,  # No legend for this plot
            alpha=0.7  # Transparency for points
        )
        # Customize the plot
        plt.axvline(x=LFC_threshold, color='red', linestyle='--', linewidth=1)
        plt.axvline(x=-LFC_threshold, color='red', linestyle='--', linewidth=1)
        ax.axhline(y=threshold_value, color='red', linestyle='--', linewidth=1)
        plt.title(plot_title, fontsize=16)
        plt.xlabel('Log2 Fold Change (LFC)', fontsize=12)
        plt.ylabel(Volcano_y_axis, fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
    # if max y above cutoff, make inset.
    if truncate:
        threshold_value = config.get("FDR_threshold")
        # Identify points exceeding the cutoff
        high_values = anova_lm_df['Log10_FDR_P_Value'] > y_cutoff
        # Create the figure and axis for plotting
        fig, ax = plt.subplots(figsize=(8, 6))
        # plot outliers, with prominent black rings to distinguish from normal points
        sns.scatterplot(
            data=anova_lm_df[high_values], 
            x='Log2_Fold_Change',  # Log fold change
            y=np.full(sum(high_values), y_cutoff),  # -log10(FDR-corrected p-value)
            hue='Colour',  # Color based on significance
            palette={'blue': 'blue', 'gray': 'gray'},
            legend=False,  # No legend for this plot
            alpha=1, # Transparency for points
            edgecolor='black',
            linewidth=1.5
        )
        # plot non-outliers
        sns.scatterplot(
            data=anova_lm_df[~high_values], 
            x='Log2_Fold_Change',  # Log fold change
            y=Volcano_y_data[~high_values],  # -log10(FDR-corrected p-value)
            hue='Colour',  # Color based on significance
            palette={'blue': 'blue', 'gray': 'gray'},
            legend=False,  # No legend for this plot
            alpha=0.7 # Transparency for points
        )
        # Customize the plot
        plt.axvline(x=LFC_threshold, color='red', linestyle='--', linewidth=1)
        plt.axvline(x=-LFC_threshold, color='red', linestyle='--', linewidth=1)
        plt.axhline(y=threshold_value, color='red', linestyle='--', linewidth=1)
        plt.title(plot_title, fontsize=16)
        plt.xlabel('Log2 Fold Change (LFC)', fontsize=12)
        plt.ylabel(Volcano_y_axis, fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        ### return the plot
        return fig, ax

def make_volcano(df_pair: pd.DataFrame,
                 output_dir: str,
                 pair_name: str,
                 config: dict,
                 metadata: pd.DataFrame = None
                 ) -> pd.DataFrame:
    """
    fit linear model for each protein. Calls run_anova, defined above

    Args:
        df (pd.DataFrame): protein abundance data.
        output_dir (str): Directory to save the output
        pair_name (str): name of treatments being compared (pairwise)
        metadata (pd.DF): sample metadata
        pair_name (str): treatment groups to be compared, separated by an underscore

    Returns:
        df_model_out (pd.DataFrame): gene name, F statistic, p value, FDR corrected p value.
    """
    anova_lm_df = df_pair.apply(run_anova, axis=1, metadata=metadata).dropna()
    n_prot = anova_lm_df.shape[0]
    plot_title = 'Protein Abundance Log Fold Change for treatments \n' + pair_name + '(n = ' + str(n_prot) + ')'
    # Apply FDR correction (Benjamini-Hochberg)
    _, fdr_corrected_pvals, _, _ = multipletests(anova_lm_df["p_value"].values , method="fdr_bh")
    # Add FDR-adjusted p-values to DataFrame
    anova_lm_df["FDR_p_value"] = fdr_corrected_pvals
    # and -log10(FDR) for plot
    anova_lm_df['Log10_FDR_P_Value'] = -np.log10(anova_lm_df['FDR_p_value'])
    anova_lm_df['Log10_unadjusted_p_Value'] = -np.log10(anova_lm_df['p_value'])
    # find thresholds
    LFC_threshold = config.get("LFC_threshold")
    FDR_threshold = config.get("FDR_threshold")
    # Add the Colour column based on LOG2FC and p_values_FDR
    anova_lm_df['Colour'] = anova_lm_df.apply(
        lambda row: 'blue' if (abs(row['Log2_Fold_Change']) >= LFC_threshold and row['FDR_p_value'] <= FDR_threshold) else 'gray', axis=1
    )
    lm_path = os.path.join(output_dir, 'data', pair_name, 'lm_output.csv')
    anova_lm_df.to_csv(lm_path, index=False)
    # make plot
    fig, ax = volcano_plot( anova_lm_df, config, plot_title )
    plot_path = os.path.join(output_dir, 'plots', pair_name, 'volcano_plot.png')
    fig.savefig(plot_path, dpi=300)
    plt.close()
    return anova_lm_df

def enrichment_analysis(anova_lm_df: pd.DataFrame,
                        pair_name: str,
                        config: dict,
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
    # threshold to define genes of interest
    LFC_threshold = config.get("LFC_threshold")
    ##### Calculate enrichment #####
    gp = GProfiler(return_dataframe=True)
    ##### G Profiler options #####
    # for ORA, just need a list of genes
    pathway_query_genes = anova_lm_df.loc[
    ( (anova_lm_df['FDR_p_value'] < 0.05) & (abs(anova_lm_df['Log2_Fold_Change']) >= LFC_threshold) ) ##################################### change back to FDR_p_value
    ]['Gene']
    # in the case of phosphoproteomic data, gene names have a double __ with phosphorylation state added,
    # for now, we remove the phospho data from this set of genes
    # may want to look at separately later
    if any(isinstance(gene, str) and '__' in gene for gene in pathway_query_genes):
        pathway_query_genes = [str(gene).split("__")[0] for gene in pathway_query_genes]
    # for GSEA, query genes can be weighted e.g. genes_weighted = {'BRCA1': 2.3, 'TP53': 1.8, 'AKT1': -1.2, 'MTOR': -2.1}
    genes_weighted_dict = anova_lm_df[ (anova_lm_df[[ 'FDR_p_value' ]]<0.05).all(axis=1) ] \
    .set_index('Gene')['Log2_Fold_Change'].to_dict()
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
    # Find all matching images
    image_paths = sorted(glob(os.path.join(search_path, "**", search_term), recursive=True))
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
    # Save the final image
    combined_image.save(output_filename)
    print(f"Combined plot saved to: {output_filename}")
    return output_filename  # Return the path for reference

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
    # Search for matching CSV files in subdirectories
    search_pattern = os.path.join(output_dir, "data", "**", filename)
    csv_files = sorted(glob(search_pattern, recursive=True))
    # check if csv files exist
    if not csv_files:
        print(f"No files found matching '{filename}'.")
        return None, None
    combined_data = []
    # Process each found CSV file
    for file in csv_files:
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
    # Auto-generate output filename if not provided
    if output_filename is None:
        output_filename = os.path.join(output_dir, f"data/combined_{filename}")
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    # Save combined CSV
    combined_df.to_csv(output_filename, index=False)
    print(f"Combined file saved at: {output_filename}")
    return combined_df, output_filename


##########################################################################
#### This is the main function for analysing data, combining the above ###
##########################################################################

def run_analysis(df: pd.DataFrame,
                 df_standardised: pd.DataFrame,
                 metadata: pd.DataFrame,
                 output_dir: str,
                 config: dict, 
                 json_out: str) -> dict:
    """
    Full analysis pipeline: performs PCA,

    Parameters:
        df (pd.DataFrame): Raw protein abundance data.
        df_standardised (pd.DataFrame): Standardised protein abundance data.
        metadata (pd.DataFrame): Metadata containing sample information.
        output_dir (str): Directory to save analysis outputs.
        json_out (str): File for saving information to go into the final report

    Returns:
        dict: Dictionary containing results from all analyses.
    """

    # Initialize results dictionary
    results = {}

    ##### Analyses for all treatment groups #####
    # Generate histogram
    print("Generating histogram")
    generate_histogram(df,
                        output_dir,
                        metadata,
                        json_out) 

    # Perform PCA and save results
    print("Performing PCA...")
    pca_results = generate_pca(df_standardised, output_dir, metadata=metadata)
    results['pca'] = pca_results

    # Perform MDS and save results
    print("Performing MDS...")
    mds_coords_df = generate_MDS(df_standardised, output_dir, metadata=metadata)
    results['mds'] = mds_coords_df

    # Generate and save heatmap
    print("Generating heatmap...")
    df_heatmap = generate_heatmap(df_standardised, output_dir)
    results['heatmap'] = df_heatmap

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
        anova_lm_df = make_volcano(df_pair,
                                   output_dir,
                                   metadata=metadata_pair,
                                   pair_name = pair_name,
                                   config = config)
        results_name = 'df_lm_' + pair_name
        results[results_name] = anova_lm_df
        # Find overrepresented pathways and save output
        print("Running enrichment analysis for pair ", pair_name, "...")
        enrichment_analysis(anova_lm_df,
                            pair_name,
                            config,
                            output_dir
                            )

    # combine plots from different pairs
    combine_plots(search_path = output_dir,
                  search_term = "volcano_plot.png",
                  output_dir=output_dir) 
    
    combine_plots(search_path = output_dir,
                  search_term = "pathway_enrichment_plot.png",
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
