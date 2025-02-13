import markdown2 # conda env info in configs/auto-prot-env-markdown-macOS.yml
import os
import subprocess
import json
import pandas as pd
from check_env import get_repo_root

def generate_report_html(report_MD: str,
                         report_html: str,
                         top_20_prots_path: str,
                         enrichment_path: str,
                         json_out: str):
    """
    Generate html report from outputs and markdown template

    Parameters:
    - report_MD (str): Path to the markdown template file.
    - report_html (str): Path to the html outputfile.
    - top_20_prots_path (str): path to csv containing data on top 20 proteins by LFC
    - enrichment_path (str): path to csv containing functional/pathway enrichment
    - json_out (str): File for saving information to go into the final repo


    Returns:
    
    """
    ### write to file the version of this script
    REPO_ROOT = get_repo_root()
    generate_version =  subprocess.check_output(["git", "log", "-n", "1", "--format=%H", "--", 
                                                 os.path.join(REPO_ROOT, 'utils/generate_report.py'),
                                                 ]).strip().decode("utf-8")
    template_version =  subprocess.check_output(["git", "log", "-n", "1", "--format=%H", "--", 
                                                 os.path.join(REPO_ROOT, 'report/report-template.md'),
                                                 ]).strip().decode("utf-8")

    # read data from json file
    with open(json_out, "r") as f:
        existing_data = json.load(f)

    script_meta = {"GENERATE_REPORT_VERSION": generate_version,
                   "TEMPLATE_VERSION": template_version}

    # Append new data
    existing_data.update(script_meta)

    # Write back to JSON file
    with open(json_out, "w") as f:
        json.dump(existing_data, f, indent=4)


    # Load metadata values
    with open(json_out, "r") as f:
        values = json.load(f)

    # Load data from csv
    # read straight to html so it can slot into the template
    top_20_df = pd.read_csv(top_20_prots_path).to_html(index=False, border =1 )

    # read straight to html so it can slot into the template
    enrichment_df = pd.read_csv(enrichment_path)[['source', 'native', 'name', 'p_value',
                                                  'term_size', 'query_size', 'intersection_size',
                                                  'precision', 'recall']] \
        .head(20).to_html(index=False, border =1 )
    
    # Open the file for reading and read the input to a temp variable
    with open(report_MD, 'r') as f:
        tempMd= f.read()
    
    # Replace placeholders
    for key, val in values.items():
        tempMd = tempMd.replace(f"${key}", str(val))

    # Replace placeholders with tables
    table_replacements = {
    "{{TOP_20_PROTS}}": top_20_df, ### top 20 protein abundance data
    "{{ENRICHMENT_DF}}": enrichment_df ### pathway/functional enrichment
    }

    for key, value in table_replacements.items():
        tempMd = tempMd.replace(key, str(value))

    # Convert the input to HTML
    tempHtml = markdown2.markdown(tempMd)
    # If necessary, could print or edit the results at this point.
    # Open the HTML file and write the output.
    with open(report_html, 'w') as f:
        f.write(tempHtml)
    print(f"Report successfully generated: {report_html}")

if __name__ == "__main__":
    ### find repo root
    REPO_ROOT = get_repo_root()

    #### set location of report template and where html will be stored
    report_MD=os.path.join(REPO_ROOT, "./report/report-template.md")
    report_html =os.path.join(REPO_ROOT, "./output/report-out.html")
    top_20_prots_path =os.path.join(REPO_ROOT, "./output/data/top_20_by_LFC.csv")
    enrichment_path =os.path.join(REPO_ROOT, "./output/data/pathway_enrichment.csv")
    json_out = os.path.join(REPO_ROOT, 'output/data/data_for_report.json')

    generate_report_html(report_MD,
                         report_html,
                         top_20_prots_path,
                         enrichment_path,
                         json_out)
