import sys
from pathlib import Path

# Add parent of `src/` to sys.path to import code.utils.check_env
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import json
import yaml
import os
import subprocess

import markdown2  # conda env info in configs/auto-prot-env-markdown-macOS.yml
import pandas as pd

from src.utils.check_env import get_repo_root
from src.reporting.image_conversion import inline_base64_images


def generate_report_html(
    report_MD: str,
    report_html: str,
    top_LFC_prots_path: str,
    enrichment_path: str,
    enrichment_plot_path: str,
    json_out: str,
    config: dict
):
    """
    Generate an HTML report by populating a markdown template with outputs from the analysis.

    This function fills placeholders in a markdown report template using data and metadata from differential
    expression and pathway enrichment analysis. It then converts the markdown to HTML, embeds base64-encoded
    images (if present), and writes a self-contained HTML report to file.

    Git commit hashes of the template and report-generation script are also recorded in a JSON metadata file.

    Args:
        report_MD (str): Path to the markdown template file.
        report_html (str): Path where the final HTML report will be written.
        top_LFC_prots_path (str): Path to CSV containing top proteins by log fold change.
        enrichment_path (str): Path to CSV containing pathway enrichment results.
        enrichment_plot_path (str): Path to image file of the enrichment plot (e.g., PNG).
        json_out (str): Path to a JSON file where report metadata and script/template versions are saved.
        config (dict): Configuration dictionary with keys including:
            - "outPath" (str): Path to output directory used for resolving image paths.

    Returns:
        None. Writes HTML report to file and updates metadata JSON.
    """
    ### write to file the version of this script
    REPO_ROOT = get_repo_root()
    generate_version = (
        subprocess.check_output(
            [
                "git",
                "log",
                "-n",
                "1",
                "--format=%H",
                "--",
                os.path.join(REPO_ROOT, "utils/generate_report.py"),
            ]
        )
        .strip()
        .decode("utf-8")
    )
    template_version = (
        subprocess.check_output(
            [
                "git",
                "log",
                "-n",
                "1",
                "--format=%H",
                "--",
                os.path.join(REPO_ROOT, "report/report-template.md"),
            ]
        )
        .strip()
        .decode("utf-8")
    )

    # read data from json file
    with open(json_out) as f:
        existing_data = yaml.safe_load(f)

    script_meta = {
        "GENERATE_REPORT_VERSION": generate_version,
        "TEMPLATE_VERSION": template_version,
    }

    # Append new data
    existing_data.update(script_meta)

    # Write back to JSON file
    with open(json_out, "w") as f:
        json.dump(existing_data, f, indent=4)

    # Load metadata values
    with open(json_out) as f:
        values = yaml.safe_load(f)

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
        enrichment_plot_md = (
            f'<img src="{enrichment_plot_path}">'
        )
    else:
        enrichment_plot_md = (
            "<p><em>No pathway enrichment plot available. "
            "This may be due to lack of enriched terms.</em></p>"
        )
    top_LFC_df = pd.read_csv(top_LFC_prots_path).to_html(index=False, border=1)

    # Open the file for reading and read the input to a temp variable
    with open(report_MD) as f:
        tempMd = f.read()

    # Replace placeholders
    for key, val in values.items():
        tempMd = tempMd.replace(f"${key}", str(val))

    # Replace placeholders with tables
    table_replacements = {
        "{{ENRICHMENT_DF}}": enrichment_df,  ### pathway/functional enrichment
        "{{ENRICHMENT_PLOT}}": enrichment_plot_md,
        "{{TOP_LFC_PROTS}}": top_LFC_df  ### top 20 protein abundance data
    }

    for key, value in table_replacements.items():
        tempMd = tempMd.replace(key, str(value))

    # also replacing the {outPath} placeholder with real out dir
    tempMd = tempMd.replace("{outPath}", config.get("outPath"))

    # Convert the input to HTML
    tempHtml = markdown2.markdown(tempMd)

    # Inline images as base64 so HTML is portable
    tempHtml = inline_base64_images(html = tempHtml, base_dir=os.path.dirname(config["outPath"]))

    # If necessary, could print or edit the results at this point.
    # Open the HTML file and write the output.
    with open(report_html, "w") as f:
        f.write(tempHtml)
    print(f"Report successfully generated: {report_html}")


if __name__ == "__main__":
    ### find repo root
    REPO_ROOT = get_repo_root()
    ### path to config file containing key info - must be present!
    config_path = os.path.join(REPO_ROOT, "configs/auto-prot-config.yaml")
    ### Read in configuration data, stored in a json
    with open(config_path) as f:
        config = yaml.safe_load(f)

    #### set location of report template and where html will be stored
    report_MD = os.path.join(REPO_ROOT, "./report/report-template.md")
    report_html = os.path.join(REPO_ROOT, config["outPath"] + "/report-out.html")
    top_LFC_prots_path = os.path.join(
        REPO_ROOT, config["outPath"] + "/full_dataset/data/combined_topLFC.csv"
    )
    enrichment_path = os.path.join(
        REPO_ROOT, config["outPath"] + "/full_dataset/data/combined_top_pathway_enrichment.csv"
    )
    enrichment_plot_path = os.path.join(
        REPO_ROOT, config["outPath"] + "/full_dataset/plots/combined_pathway_enrichment_plot.png"
    )
    json_out = os.path.join(REPO_ROOT, config["json_outPath"])

    generate_report_html(
        report_MD,
        report_html,
        top_LFC_prots_path,
        enrichment_path,
        enrichment_plot_path,
        json_out,
        config
    )
