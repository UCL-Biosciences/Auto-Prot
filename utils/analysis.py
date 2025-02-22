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


def generate_pca(df_standardised: pd.DataFrame,
                 output_dir: str,
                 metadata: pd.DataFrame = None,
                 components: int = 2) -> pd.DataFrame:
    """
    Perform PCA on the given DataFrame and save the results as a plot.

    Args:
        df_standardised (pd.DataFrame): Standardised protein abundance data.
        output_dir (str): Directory to save the PCA plot and data.
        components (int): Number of principal components to compute.

    Returns:
        pd.DataFrame: DataFrame with principal components and metadata.
    """

    # Perform PCA
    pca = PCA(n_components=components)
    principal_components = pca.fit_transform(df_standardised)

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
    plt.title('PCA Plot')
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
                 metadata: pd.DataFrame = None,
                 components: int = 2) -> pd.DataFrame:
    """
    Perform MDS on the given DataFrame and save the results as a plot.

    Args:
        df_standardised (pd.DataFrame): Standardised protein abundance data.
        output_dir (str): Directory to save the PCA plot and data.

    Returns:
        pd.DataFrame: DataFrame with MDS and metadata.
    """

    # Perform MDS
    # pdist(x) computes the Euclidean distances between each pair of points in an array
    dissimilarity_array = pdist(df_standardised, metric='euclidean')
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
    plt.title('MDS Dissimilarity based on Protein Abundance')
    plt.tight_layout()

    # Save plot and PCA data
    plot_path = os.path.join(output_dir, 'plots', 'mds_plot.png')
    data_path = os.path.join(output_dir, 'data', 'mds_data.csv')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    mds_coords_df.to_csv(data_path, index=False)

    print(f"PCA plot saved to {plot_path}")
    print(f"PCA data saved to {data_path}")

    return mds_coords_df


def generate_heatmap(df_standardised: pd.DataFrame,
                 output_dir: str,
                 metadata: pd.DataFrame = None) -> pd.DataFrame:
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
    plt.title('Heatmap of Sample Dissimilarity based on Protein Abundance')
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
        
    Returns:
        
    """
    gene_name = row.name  # Use row index as gene name
    data = pd.DataFrame(row)  # Convert row to DataFrame
    data.columns=['Abundance'] # rename column so can be used in model
    # Map sample_id to treatment
    data["treatment"] = data.index.map(metadata.set_index("sample_rep")["treatment"])
    # Fit OLS model
    model = smf.ols("Abundance ~ C(treatment)", data=data).fit()
    # extract F stat and p value for the model
    # anova for looks at within and between group variance
    # between group variance is based on differences of group means from overall means
    # within-group variance is residuals of observations from group mean
    # Mean Squares stats are SS / df for treatment and residuals
    # F stat is MS(treatment) / MS(residual)
    # a large ratio implies more difference between groups
    # a large ratio implies more difference within groups
    anova_table = sm.stats.anova_lm(model, typ=2)  # type 2 output doesn't depend on order of predictors
    f_stat, p_value = anova_table.iloc[0]["F"], anova_table.iloc[0]["PR(>F)"]
    #### Extract group means and LFC
    # Compute group means (sorted alphabetically)
    group_means = data.groupby("treatment")["Abundance"].mean().sort_index()
    grp1_name = group_means.index.tolist()[0]
    grp1_mean = group_means.iloc[0]
    grp2_name = group_means.index.tolist()[1]
    grp2_mean = group_means.iloc[1]
    fc = ( grp2_mean + 1) / ( grp1_mean + 1 )
    log2fc = np.log2( fc )
    #return outputs
    results_anova = pd.Series({"Gene": gene_name,
                      "Group_1": grp1_name,
                      "Group_1_mean": grp1_mean,
                      "Group_2": grp2_name,
                      "Group_2_mean": grp2_mean,
                      "Raw_Fold_Change": fc,
                      "Log2_Fold_Change": log2fc,
                      "F_Stat": f_stat,
                      "p_value": p_value})
    return results_anova


def make_volcano(df_pair: pd.DataFrame,
                 output_dir: str,
                 metadata: pd.DataFrame = None,
                 pair_name: str) -> pd.DataFrame:
    """
    fit linear model for each protein. Calls run_anova, defined above

    Args:
        df (pd.DataFrame): protein abundance data.
        output_dir (str): Directory to save the output
        metadata (pd.DataFrame): metadata
        pair_name (str): treatment groups to be compared, separated by an underscore

    Returns:
        df_model_out (pd.DataFrame): gene name, F statistic, p value, FDR corrected p value.
    """
    anova_lm_df = df_pair.apply(run_anova, axis=1, metadata=metadata)
    # Apply FDR correction (Benjamini-Hochberg)
    _, fdr_corrected_pvals, _, _ = multipletests(anova_lm_df["p_value"].values , method="fdr_bh")
    # Add FDR-adjusted p-values to DataFrame
    anova_lm_df["FDR_p_value"] = fdr_corrected_pvals
    # and -log10(FDR) for plot
    anova_lm_df['Log10_FDR_P_Value'] = -np.log10(anova_lm_df['FDR_p_value'])
    # Add the Colour column based on LOG2FC and p_values_FDR
    anova_lm_df['Colour'] = anova_lm_df.apply(
        lambda row: 'blue' if (abs(row['Log2_Fold_Change']) > 2 and row['FDR_p_value'] < 0.05) else 'gray', axis=1
    )
    ### Create volcano plot
    sns.scatterplot(
        data=anova_lm_df, 
        x='Log2_Fold_Change',  # Log fold change
        y='Log10_FDR_P_Value',  # -log10(FDR-corrected p-value)
        hue='Colour',  # Color based on significance
        palette={'blue': 'blue', 'gray': 'gray'},
        legend=False,  # No legend for this plot
        alpha=0.7  # Transparency for points
    )
    # Customize the plot
    plt.axvline(x=2, color='red', linestyle='--', linewidth=1)  # Threshold for LFC > 1
    plt.axvline(x=-2, color='red', linestyle='--', linewidth=1)  # Threshold for LFC < -1
    plt.title('Protein Abundance Log Fold Change', fontsize=16)
    plt.xlabel('Log2 Fold Change (LFC)', fontsize=12)
    plt.ylabel('-log10(FDR-corrected P-value)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    # Set symmetrical x-axis
    # Save plot and PCA data
    #plt.xlim(-4, 4)
    plot_path = os.path.join(output_dir, 'plots', pair_name, 'volcano_plot.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    # save the top 20 rows to csv for display in final report
    # Sort by a relevant column (modify column name as needed)
    sorted_df = anova_lm_df.sort_values(by="Log2_Fold_Change", ascending=False, key = abs)
    # Select the top 20 rows
    top_20_df = sorted_df.head(20).round(decimals = 2)
    top_20_df = top_20_df.round({"p_value":4, "FDR_p_value":4})
    # Save to CSV
    top_prot_path = os.path.join(output_dir, 'data', pair_name, 'top_20_by_LFC.csv')
    top_20_df.to_csv(top_prot_path, index=False)
    return anova_lm_df, top_20_df

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
                        output_dir: str) :
    """
    run enrichment analysis to identify overrepresented genes, functions and pathways

    Args:
        anova_lm_df (pd.DataFrame): df with anova outputs including LFC, p value and FDR adjusted p value
        output_dir (str): Directory to save the output.

    Returns:
        enrichment data (pd.DataFrame)
    """    

    ##### Calculate enrichment #####
    gp = GProfiler(return_dataframe=True)

    ##### G Profiler options #####
    # for ORA, just need a list of genes
    pathway_query_genes = anova_lm_df.loc[
    (anova_lm_df['FDR_p_value'] < 0.05)
    ]['Gene']

    # for GSEA, query genes can be weighted e.g. genes_weighted = {'BRCA1': 2.3, 'TP53': 1.8, 'AKT1': -1.2, 'MTOR': -2.1}
    genes_weighted_dict = anova_lm_df[ (anova_lm_df[[ 'FDR_p_value' ]]<0.05).all(axis=1) ] \
    .set_index('Gene')['Log2_Fold_Change'].to_dict()

    # pathway database can be REAC, GO or KEGG. Also less common but available: CORUM, HPA, TF and MIRNA
    # defaults to REAC
    source=['REAC']

    # p value threshold defaults to 0.05
    p_threshold=0.05

    # all results returns all results, not just those below p threshold
    all_results = False

    # a background set can pe specified e.g. background=["BRCA1", "TP53", "AKT1", "MTOR", "EGFR", "MYC"]

    # multiple testing correction can be g_SCS (default, Set Counts and Sizes), bonferroni, or fdr
    # from quick look, g_SCS seems to be similar to bonferroni, which are both less strict than fdr.
    significance_threshold_method='g_SCS'

    ##### Running pathway enrichment #####
    pathway_result = gp.profile(organism='hsapiens',
        query=pathway_query_genes.to_list(),
        sources=source,
        user_threshold = p_threshold,
        significance_threshold_method = significance_threshold_method,
        all_results = all_results
    )  # REAC for Reactome

    ### save results to file
    enrichment_path = os.path.join(output_dir, 'data', 'pathway_enrichment.csv')
    pathway_result.round(decimals = 2).to_csv(enrichment_path, index=False)

    ##### Plot enrichment #####
    pathway_plot_df = pathway_result.sort_values('p_value', ascending = True).head(20)

    pathway_plot_df['-log10(p_value)'] = -np.log10(pathway_plot_df['p_value'])

    # Setting up the plot using Seaborn and Matplotlib adjustments
    plt.figure(figsize=(15, 20))
    #ax = plt.gca()

    sns.relplot(
        data = pathway_plot_df,
        x = "-log10(p_value)",
        y = "name",
        size = "recall", # the proportion of query genes associated with the term
        color = "green",
        #height=20,       # Adjust figure height for better fit
        aspect=0.6       # Maintain a suitable aspect ratio)
    )
    # Adjustments for axes padding and limits
    # Adjustments for labels and titles
    plt.xlabel('-log10(p-value)', fontsize=14)
    plt.ylabel('Pathway Names', fontsize=14)
    plt.title('Pathway Enrichment', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)  # Ensure y-axis labels are readable

    # Adjust layout to ensure labels are fully visible
    plt.tight_layout()

    # Save plot data
    plot_path = os.path.join(output_dir, 'plots', 'pathway_enrichment_plot.png')
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")  # bbox_inches ensures labels aren't cut off
    plt.close()
    return pathway_result
    

##########################################################################
#### This is the main function for analysing data, combining the above ###
##########################################################################

def run_analysis(df: pd.DataFrame,
                 df_standardised: pd.DataFrame,
                 metadata: pd.DataFrame,
                 output_dir: str,
                 json_out: str) -> dict:
    """
    Full analysis pipeline: performs PCA,

    Parameters:
        df (pd.DataFrame): Raw protein abundance data.
        df_standardised (pd.DataFrame): Standardised protein abundance data.
        metadata (pd.DataFrame): Metadata containing sample information.
        output_dir (str): Directory to save analysis outputs.
        json_out (str): File for saving information to go into the final repo

    Returns:
        dict: Dictionary containing results from all analyses.
    """

    # Initialize results dictionary
    results = {}

    ##### Analyses for all treatment groups #####

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
    df_heatmap = generate_heatmap(df_standardised, output_dir, metadata=metadata)
    results['heatmap'] = df_heatmap

    # Look at overlaps of differnential abundance

    #### if 2 or 3 groups, Venn

    #### if 4 or more, upset

    ###### Pairwise Analyses #####
    # if there are > 2 treatment groups, pairwise analyses will have to be run separately for each pair of treatments
    treatment_pairs = list(itertools.combinations(metadata['treatment'].unique(), 2))

    for pair in treatment_pairs:
        metadata_pair = metadata[metadata['treatment'].isin(pair)]
        df_pair = df[metadata_pair['sample_rep'].tolist()]
        pair_name = "_".join(map(str, pair) )
        if not os.path.exists(os.path.join(output_dir, 'plots', pair_name)):
            os.mkdir(os.path.join(output_dir, 'plots', pair_name))
        if not os.path.exists(os.path.join(output_dir, 'data', pair_name)):
            os.mkdir(os.path.join(output_dir, 'data', pair_name))

        # Generate and save volcano plot
        print("Generating volcano plot...")
        anova_lm_df, top_20_df = make_volcano(df_pair,
                                              output_dir,
                                              metadata=metadata_pair,
                                              pair_name = pair_name)
        results_name = 'volcano_' + pair_name
        results['results_name'] = anova_lm_df




    # Generate and save volcano plot
    print("Generating volcano plot...")
    anova_lm_df, top_20_df = make_volcano(df, output_dir, metadata=metadata)
    results['volcano'] = anova_lm_df

    # Generate and abundance of top proteins by LFC
    print("Plotting abundance...")
    plot_abundance(df=df,
                   top_20_df=top_20_df,
                   output_dir=output_dir,
                   metadata=metadata)
    
    # Find overrepresented pathways and save output
    print("Running enrichment analysis...")
    enrichment_analysis(anova_lm_df,
                   output_dir
    )
    
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



