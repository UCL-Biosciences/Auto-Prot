### Analysis for pairwise datasets
## for when comparing two groups or treatments
## e.g. calculate diff abundance, make volcanoes, run pathway enrichment analysis
import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from gprofiler import GProfiler


def make_volcano(
    df_pair: pd.DataFrame,
    output_dir: str,
    pair_name: str,
    config: dict,
    formula: str,
    metadata_pair: pd.DataFrame = None,
) -> pd.DataFrame:
    """
    Performs differential abundance analysis for a pairwise comparison using an external R script (limma),
    generates a volcano plot, and saves the top proteins to CSV.

    This function exports the relevant data, runs an R script to perform the differential expression analysis
    (via the limma package), reads the results back in, and makes a volcano plot based on user-specified thresholds.

    Args:
        df_pair (pd.DataFrame): Protein abundance data for two groups (samples as rows, proteins as columns).
        output_dir (str): Root directory where results (plots and data) will be saved.
        pair_name (str): Label for this treatment comparison (used in filenames and output folders).
        config (dict): Configuration containing:
            - "LFC_threshold" (float): Minimum absolute log2 fold change to consider a protein differentially expressed.
            - "FDR_threshold" (float): Maximum adjusted p-value (FDR) to be differentially expressed.
            - "LFC_plot_p_or_FDRp" (str): Column to use for y-axis in volcano plot ("Log10_FDR_P_Value" or "Log10_unadjusted_p_Value").
        metadata_pair (pd.DataFrame, optional): Metadata for the subset of samples in this comparison.
        formula (str): the formula passted to the DE calculation. May need to be different for full dataset and subsets.

    Returns:
        pd.DataFrame: DataFrame with differential expression results, including logFC, p-values, adjusted p-values,
                      and plot colour classification.
    """
    #### calculate DE using limma (R package) ####
    pair_data_path = os.path.join(output_dir, "data", pair_name, "prots.csv")
    df_pair.to_csv(pair_data_path, index=True)
    # Save sample metadata
    pair_metadata_path = os.path.join(output_dir, "data", pair_name, "metadata.csv")
    metadata_pair.to_csv(pair_metadata_path, index=False)
    # Define where to save the limma results
    pair_result_path = os.path.join(output_dir, "data", pair_name, "limma_output.csv")
    # run R script - note: r-limma-env conda env required
    subprocess.run(
        [
            "conda",
            "run",
            "-n",
            "r-limma-env",
            "Rscript",
            "autoprot/r_scripts/DE-limma.R",
            pair_data_path.replace("\\", "/"),
            pair_metadata_path.replace("\\", "/"),
            pair_result_path.replace("\\", "/"),
            formula,  ## formula used in DE analysis
        ],
        check=True,
    )
    # read results back in
    diffExpr_df = pd.read_csv(pair_result_path, index_col=0)
    n_prot = diffExpr_df.shape[0]
    diffExpr_df["Log10_FDR_P_Value"] = -np.log10(diffExpr_df["adj.P.Val"])
    diffExpr_df["Log10_unadjusted_p_Value"] = -np.log10(diffExpr_df["P.Value"])
    ### whether to plot the -log10(p_value) i.e. unadjusted or -log10(FDR_p_value) is specified in json field "LFC_plot_p_or_FDRp" ("Log10_FDR_P_Value" or "Log10_unadjusted_p_Value")
    LFC_threshold = config["LFC_threshold"]
    FDR_threshold = config["FDR_threshold"]
    Volcano_y_axis = config["LFC_plot_p_or_FDRp"]
    Volcano_y_data = diffExpr_df[Volcano_y_axis]
    # Add the Colour column based on LOG2FC and p_values_FDR
    ### colour blue when p value < threshold defined in config (LFC_threshold)
    ### by default, use FDR adjusted p value
    ### but might want to use unadjusted p value to compare with others
    ### defined by LFC_plot_p_or_FDRp in config:
    if Volcano_y_axis == "Log10_FDR_P_Value":
        p_cutoff_column = "adj.P.Val"
    elif Volcano_y_axis == "Log10_unadjusted_p_Value":
        p_cutoff_column = "P.Value"
    diffExpr_df["Colour"] = diffExpr_df.apply(
        lambda row: (
            "blue"
            if (
                abs(row["logFC"]) > LFC_threshold
                and row[p_cutoff_column] < FDR_threshold
            )
            else "gray"
        ),
        axis=1,
    )
    diffExpr_path = os.path.join(output_dir, "data", pair_name, f"{pair_name}_limma_output.csv")
    diffExpr_df.to_csv(diffExpr_path, index=True)
    ### Create volcano plot
    plot_title = (
        "Protein Abundance Log Fold Change for treatments \n"
        + pair_name
        + "(n = "
        + str(n_prot)
        + ")"
    )
    sns.scatterplot(
        data=diffExpr_df,
        x="logFC",  # Log fold change
        y=Volcano_y_data,  # -log10(FDR-corrected p-value)
        hue="Colour",  # Color based on significance
        palette={"blue": "blue", "gray": "gray"},
        legend=False,  # No legend for this plot
        alpha=0.7,  # Transparency for points
    )
    ### label top 10 blue proteins
    # Sort by adjusted p-value and take top 10 (or all if < 10)
    top_proteins = diffExpr_df[diffExpr_df["Colour"] == "blue"].nsmallest(
        10, "adj.P.Val"
    )
    # Prepare text objects
    texts = []
    for _, row in top_proteins.iterrows():
        texts.append(
            plt.text(
                row["logFC"],
                row[Volcano_y_axis],
                row.name,  # Assumes protein names are in the index
                fontsize=12,
            )
        )

    # Automatically adjust text positions to minimise overlap
    adjust_text(
        texts,
        arrowprops=dict(arrowstyle="->", color="gray", lw=0.5),
        force_text=0.5,
        force_points=0.3,
        expand_points=(1.2, 1.4),
        expand_text=(1.2, 1.4),
        only_move={"points": "y", "text": "xy"},
    )
    # Customize the plot
    plt.axvline(x=LFC_threshold, color="red", linestyle="--", linewidth=1)
    plt.axvline(x=-LFC_threshold, color="red", linestyle="--", linewidth=1)
    plt.axhline(y=-np.log10(FDR_threshold), color="red", linestyle="--", linewidth=1)
    plt.title(plot_title, fontsize=16)
    plt.xlabel("Log2 Fold Change (LFC)", fontsize=12)
    plt.ylabel(Volcano_y_axis, fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.6)
    # Set symmetrical x-axis
    # Save plot and PCA data
    # plt.xlim(-4, 4)
    plot_path = os.path.join(output_dir, "plots", pair_name, "volcano_plot.png")
    if os.name == "nt":
        plot_path = "\\\\?\\" + os.path.abspath(plot_path)
    plt.savefig(plot_path, dpi=300)
    plt.close()
    # save the top 20 rows to csv for display in final report
    # Sort by a relevant column (modify column name as needed)
    sorted_df = diffExpr_df.sort_values(by="logFC", ascending=False, key=abs)
    # Select the top 20 rows
    top_20_df = sorted_df.head(20).round(decimals=2)
    top_20_df = top_20_df.round({"P.Value": 4, "adj.P.Value": 4})
    # Save to CSV
    top_prot_path = os.path.join(output_dir, "data", pair_name, "top_20_by_LFC.csv")
    if os.name == "nt":
        top_prot_path = "\\\\?\\" + os.path.abspath(top_prot_path)
    top_20_df.to_csv(top_prot_path, index=True)
    return diffExpr_df


def enrichment_analysis(
    lm_results_df: pd.DataFrame, pair_name: str, config: dict, output_dir: str
):
    """
    Performs pathway enrichment analysis using g:Profiler for differentially abundant proteins.

    This function applies over-representation analysis (ORA) via g:Profiler to identify enriched biological pathways
    (e.g., KEGG) among significantly changing proteins. It selects genes based on log fold change and adjusted p-value
    thresholds specified in the `config`, doesn't currently remove post-translational modification labels (e.g., "__phospho") if present,
    and performs enrichment analysis for the configured species.

    The results are saved as a CSV file, and a bubble plot of the top 20 pathways (by p-value) is saved as a PNG image.

    Args:
        lm_results_df (pd.DataFrame): Differential expression results, with 'logFC', 'adj.P.Val', and 'P.Value'.
        pair_name (str): Label for the treatment comparison, used in output filenames and titles.
        config (dict): Configuration dictionary containing:
            - "LFC_threshold" (float): Minimum |log2 fold change| to define differentially expressed genes.
            - "FDR_threshold" (float): Maximum adjusted p-value to define significance.
            - "species" (str): g:Profiler organism code (e.g., 'hsapiens').
        output_dir (str): Root output directory where results and plots will be saved.

    Returns:
        pd.DataFrame or None: g:Profiler enrichment results (may be empty or None if no significant genes were found).
    """

    # threshold to define genes of interest
    LFC_threshold = config["LFC_threshold"]
    FDR_threshold = config["FDR_threshold"]
    ##### Calculate enrichment #####
    gp = GProfiler(return_dataframe=True)
    ##### G Profiler options #####
    # for ORA, just need a list of genes
    pathway_query_genes = lm_results_df.loc[
        (
            (lm_results_df["adj.P.Val"] < FDR_threshold)
            & (abs(lm_results_df["logFC"]) >= LFC_threshold)
        )
    ].index
    # in the case of phosphoproteomic data, gene names have a double __ with phosphorylation state added,
    # for now, we remove the phospho data from this set of genes
    # may want to look at separately later
    # if any(isinstance(gene, str) and "__" in gene for gene in pathway_query_genes):
    #     pathway_query_genes = [str(gene).split("__")[0] for gene in pathway_query_genes]
    # pathway database can be REAC, GO or KEGG. Also less common but available: CORUM, HPA, TF and MIRNA
    # defaults to REAC
    source = ["KEGG"]
    plot_title = "Pathway Enrichment (" + source[0] + ") for treatments \n " + pair_name
    # p value threshold defaults to 0.05
    p_threshold = 0.05
    # all results returns all results, not just those below p threshold
    all_results = False
    # a background set can be specified e.g. background=["BRCA1", "TP53", "AKT1", "MTOR", "EGFR", "MYC"]
    # multiple testing correction can be g_SCS (default, Set Counts and Sizes), bonferroni, or fdr
    # from quick look, g_SCS seems to be similar to bonferroni
    significance_threshold_method = "g_SCS"
    ##### Running pathway enrichment #####
    if len(pathway_query_genes) > 0:
        pathway_result = gp.profile(
            organism=config["species"],
            query=list(pathway_query_genes),
            sources=source,
            user_threshold=p_threshold,
            significance_threshold_method=significance_threshold_method,
            all_results=all_results,
        )  # REAC for Reactome
        ### gprofiler can return empty results that write to file. In that case, return nothing
        if pathway_result.empty:
            print(
                f"⚠️ No enriched pathways found for {pair_name}. Skipping save and plot."
            )
            return None
        ### save results to file
        enrichment_path = os.path.join(
            output_dir, "data", pair_name, (pair_name + "_pathway_enrichment.csv")
        )
        ## windows sometimes rejects long paths. Workaround:
        if os.name == "nt":
            enrichment_path = "\\\\?\\" + os.path.abspath(enrichment_path)
        pathway_result.round({"precision": 2, "recall": 2}).to_csv(
            enrichment_path, index=False
        )
        ##### Plot enrichment #####
        pathway_plot_df = pathway_result.sort_values("p_value", ascending=True).head(20)
        pathway_plot_df["-log10(p_value)"] = -np.log10(pathway_plot_df["p_value"])
        # Setting up the plot using Seaborn and Matplotlib adjustments
        # plt.figure(figsize=(15, 20))
        # ax = plt.gca()
        sns.relplot(
            data=pathway_plot_df,
            x="-log10(p_value)",
            y="name",
            size="recall",  # the proportion of query genes associated with the term
            color="green",
            height=15,  # Adjust figure height for better fit
            aspect=0.5,  # Maintain a suitable aspect ratio)
        )
        # Adjustments for axes padding and limits
        # Adjustments for labels and titles
        plt.xlabel("-log10(p-value)", fontsize=14)
        plt.ylabel("Pathway Names", fontsize=14)
        plt.title(plot_title, fontsize=16)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)  # Ensure y-axis labels are readable
        # Adjust layout to ensure labels are fully visible
        plt.tight_layout()
        # Save plot data
        plot_path = os.path.join(
            output_dir, "plots", pair_name, (pair_name + "_pathway_enrichment_plot.png")
        )
        ## windows sometimes rejects long paths. Workaround:
        if os.name == "nt":
            plot_path = "\\\\?\\" + os.path.abspath(plot_path)
        plt.savefig(
            plot_path, dpi=300, bbox_inches="tight"
        )  # bbox_inches ensures labels aren't cut off
        plt.close()
        return pathway_result
