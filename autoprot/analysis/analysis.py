### Perform data analysis ###

## Using objects created previously: metadata, df_protAbundance
import itertools
import json
import os
import subprocess
import warnings

import pandas as pd
import yaml

from sklearn.preprocessing import StandardScaler

from autoprot.analysis.all_samples import generate_heatmap, generate_MDS, generate_pca, run_clustering_analysis
from autoprot.analysis.pairwise import enrichment_analysis, make_volcano
from autoprot.utils.check_env import get_repo_root
from autoprot.utils.data_utils import combine_csv_files, combine_plots

##### there is a warning to suppress when fitting models
warnings.filterwarnings(
    "ignore",
    message="Negative binomial dispersion parameter alpha not set. Using default value alpha=1.0.",
)


##########################################################################
#### This is the main function for analysing data, combining the above ###
##########################################################################
def run_analysis(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: str,
    config: dict,
    formula: str,
) -> dict:
    """
    Runs the full analysis pipeline for a protein abundance dataset, including clustering (PCA, MDS, heatmap)
    sample-level visualisations and pairwise comparisons (differential abundance, volcano plot, pathway enrichment).

    The pipeline includes:
        - Optionally z-score transforming the data for clustering
        - PCA, MDS, and heatmap generation for all samples
        - Pairwise differential abundance analysis using limma (via an R script)
        - Optional functional enrichment analysis using g:Profiler. Default is GO
        - Aggregation of plots and results into summary files

    Args:
        df (pd.DataFrame): Protein abundance data (samples as rows, proteins as columns).
        metadata (pd.DataFrame): Sample metadata, including 'sample_rep' and 'treatment' columns.
        output_dir (str): Path to directory where outputs (plots, data files) will be written.
        config (dict): Configuration parameters for differential analysis and plotting, including:
            - "LFC_threshold" (float): Suggested log fold change threshold for classification and plotting.
            - "FDR_threshold" (float): Suggested p-value threshold for volcano plot annotation.
            - "LFC_plot_p_or_FDRp" (str): Column to use for y-axis in volcano plot.
        formula (str): the formula passted to the DE calculation. May need to be different for full dataset and subsets.

    Returns:
        dict: Dictionary containing result DataFrames from PCA, MDS, heatmap, and each pairwise analysis.
              Keys include "pca", "mds", "heatmap", and "df_lm_<pair_name>" for each treatment comparison.
    """
    # Initialize results dictionary
    results = {}

    ##### Analyses for all treatment groups #####
    print(df.head())

    # Perform clustering (PCA, MDS, heatmap) and save results
    if config.get("z_score_for_clustering"):
        scaler = StandardScaler().set_output(transform="pandas")
        df_scaled = scaler.fit_transform(df)
        results = run_clustering_analysis(df = df_scaled, metadata = metadata, output_dir = output_dir)
    else:
        results = run_clustering_analysis(df = df, metadata = metadata, output_dir = output_dir)

    print(df.head())
    ###### Pairwise Analyses #####
    # if there are > 2 treatment groups, pairwise analyses will have to be run separately for each pair of treatments
    treatment_pairs = list(itertools.combinations(metadata["treatment"].unique(), 2))

    for pair in treatment_pairs:
        print("starting analysis for pair ", pair)

        metadata_pair = metadata[metadata["treatment"].isin(pair)]
        df_pair = df[metadata_pair["sample_rep"].tolist()]
        pair_name = "_".join(map(str, pair))

        os.makedirs(os.path.join(output_dir, "plots", pair_name), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "data", pair_name), exist_ok=True)

        # Generate and save volcano plot
        print("Generating volcano plot for pair ", pair_name, "...")
        lm_results_df = make_volcano(
            df_pair,
            output_dir,
            metadata_pair=metadata_pair,
            pair_name=pair_name,
            config=config,
            formula=formula,
        )
        results_name = "df_lm_" + pair_name
        results[results_name] = lm_results_df

        # Find overrepresented pathways and save output
        print("Running enrichment analysis for pair ", pair_name, "...")
        enrichment_analysis(lm_results_df, pair_name, config, output_dir)

    # combine plots from different pairs
    combine_plots(
        search_path=output_dir, search_term="volcano", output_dir=output_dir
    )

    combine_plots(
        search_path=output_dir,
        search_term="pathway_enrichment",
        output_dir=output_dir,
    )

    # Combine most DE proteins for each treatment pair
    combine_csv_files(
        filename="limma_output.csv",
        output_dir=output_dir,
        output_filename=os.path.join(output_dir, "data/combined_topLFC.csv"),
        sort_by_logfc=True,
    )

    # Combine pathway enrichment data
    combine_csv_files(
        filename="*pathway_enrichment.csv",
        output_dir=output_dir,
        output_filename=os.path.join(
            output_dir, "data/combined_top_pathway_enrichment.csv"
        ),
    )

    print("Analysis pipeline completed.")
    return results
