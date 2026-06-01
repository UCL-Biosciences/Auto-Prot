import sys
from pathlib import Path

# Add parent of `autoprot/` to sys.path to import code.utils.check_env
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import json
import os
import subprocess

import markdown2  # conda env info in configs/auto-prot-env-markdown-macOS.yml
import pandas as pd
import yaml

from autoprot.reporting.image_conversion import inline_base64_images
from autoprot.utils.check_env import get_repo_root


def generate_report_html(config = None):
    """
    Generate an HTML report by populating a markdown template with outputs from the analysis.

    This function fills placeholders in a markdown report template using data and metadata from differential
    expression and pathway enrichment analysis. It then converts the markdown to HTML, embeds base64-encoded
    images (if present), and writes a self-contained HTML report to file.

    Git commit hashes of the template and report-generation script are also recorded in a JSON metadata file.

    Args:
        config (dict): Configuration dictionary with keys including:
            - "outPath" (str): Path to output directory used for resolving image paths.

    Returns:
        None. Writes HTML report to file and updates metadata JSON.
    """
    ### find repo root
    REPO_ROOT = get_repo_root()

    #### set location of report template and where html will be stored
    report_MD = os.path.join(REPO_ROOT, "./report/report-template.md")
    report_html = os.path.join(REPO_ROOT, config["outPath"] + "/report-out.html")
    top_LFC_prots_path = os.path.join(
        REPO_ROOT, config["outPath"] + "/full_dataset/data/combined_topLFC.csv"
    )
    enrichment_path = os.path.join(
        REPO_ROOT,
        config["outPath"] + "/full_dataset/data/combined_top_pathway_enrichment.csv",
    )
    enrichment_plot_path = os.path.join(
        REPO_ROOT,
        config["outPath"] + "/full_dataset/plots/combined_pathway_enrichment_plot.png",
    )
    json_out = os.path.join(REPO_ROOT, config["json_outPath"])

    # read data from json file
    with open(json_out) as f:
        values = json.load(f)

    script_meta = {
        "PERCENT_MISSING" : round( 100 * config.get("missing_threshold") ),
        "DF_USED" : config.get("df_to_use"),
        "NORM_METHOD": config.get("normalise_method"),
        "IMP_METHOD" : config.get("imputation_method"),
        "FDR_THRESHOLD" : config.get("FDR_threshold"),
        "IQR_THRESHOLD" : config.get("IQR_threshold") *100,
    }

    # Append new data
    values.update(script_meta)

    # Write back to JSON file
    with open(json_out, "w") as f:
        json.dump(values, f, indent=4)

    # Load data from csv
    # read straight to html so it can slot into the template
    # Default message if no enrichment results are found
    no_enrichment_html = (
        "<p><em>No pathway enrichment was found for the full dataset. "
        "This may be due to a lack of differentially expressed proteins.</em></p>"
    )

    if os.path.exists(enrichment_path):
        enrichment_df = pd.read_csv(enrichment_path)[
            [
                "source",
                "native",
                "name",
                "p_value",
                "term_size",
                "query_size",
                "intersection_size",
                "precision",
                "recall",
                "treatment_pair",
            ]
        ].to_html(index=False, border=1)
    else:
        enrichment_df = no_enrichment_html

    # similarly for the enrichment plot
    if os.path.exists(enrichment_plot_path):
        enrichment_plot_md = f'<img src="{enrichment_plot_path}">'
    else:
        enrichment_plot_md = (
            "<p><em>No pathway enrichment plot available. "
            "This may be due to lack of enriched terms.</em></p>"
        )
    top_LFC_df = pd.read_csv(top_LFC_prots_path)

    # P.Value and adj.P.Val round to 4dp
    top_LFC_df["P.Value"] = pd.to_numeric(top_LFC_df["P.Value"], errors="coerce").round(4)
    top_LFC_df["adj.P.Val"] = pd.to_numeric(top_LFC_df["adj.P.Val"], errors="coerce").round(4)

    # logFC, AveExpr, t, B all to 2dp
    for col in [ "logFC", "AveExpr", "t", "B" ]:
        top_LFC_df[ col ] = top_LFC_df[ col ].round(2)

    # convert to html for template
    top_LFC_df = top_LFC_df.to_html(index=False, border=1)

    # Open the file for reading and read the input to a temp variable
    with open(report_MD) as f:
        tempMd = f.read()

    # Replace placeholders
    for key, val in values.items():
        # first escape underscores so they are not interpreted as italics by markdown
        val = str(val).replace("_", r"\_")
        tempMd = tempMd.replace(f"${key}", str(val))

    # Replace placeholders with tables
    table_replacements = {
        "{{ENRICHMENT_DF}}": enrichment_df,  ### pathway/functional enrichment
        "{{ENRICHMENT_PLOT}}": enrichment_plot_md,
        "{{TOP_LFC_PROTS}}": top_LFC_df,  ### top 20 protein abundance data
    }

    for key, value in table_replacements.items():
        tempMd = tempMd.replace(key, str(value))

    # also replacing the {outPath} placeholder with real out dir
    tempMd = tempMd.replace("{outPath}", config.get("outPath"))

    # Convert the input to HTML
    tempHtml = markdown2.markdown(tempMd)

    # Inline images as base64 so HTML is portable
    tempHtml = inline_base64_images(
        html=tempHtml )

    # If necessary, could print or edit the results at this point.
    # Open the HTML file and write the output.
    with open(report_html, "w") as f:
        f.write(tempHtml)
    print(f"Report successfully generated: {report_html}")


if __name__ == "__main__":
    REPO_ROOT = get_repo_root()
    config_path = os.path.join(REPO_ROOT, "configs/auto-prot-config.yaml")
    ### Read in configuration data, stored in a json
    with open(config_path) as f:
        config = yaml.safe_load(f)
    generate_report_html(config)
