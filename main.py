# Generate Proteomics Report

## Generate common proteomics outputs, save to file and combine in a report.
## Inputs must be protein abundance file and metadata, as specified below and in the github README.

## Import libraries
import os
import subprocess

import yaml

import autoprot.analysis.analysis as an
import autoprot.processing.data_processing as dp

## Import functions
# Functions are saved in separate files and imported here
# Separated by module
# More detail in relevant files and on github.
import autoprot.utils.check_env as env
from autoprot.utils.data_io import make_outdir
from autoprot.utils.data_utils import get_subset, tidy_up_files
from autoprot.reporting.generate_report import generate_report_html


##### Define main function for creating outputs
def main():
    """
    Main script to run the complete data analysis pipeline.
    """
    # Check environment
    # Check if running in VSCode, if not, compare environments
    if "VSCODE_PID" not in os.environ:
        env.compare_envs()
    # File paths
    # getting repo dir automatically is useful as should mean we don't need to specify
    # machine-specific paths and code should run on different users' machines
    REPO_ROOT = env.get_repo_root()
    ### path to config file containing key info - must be present!
    config_path = os.path.join(REPO_ROOT, "configs/auto-prot-config.yaml")
    ### Read in configuration data, stored in a json
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # make sure imputation method is not pimms if running on a mac - see docs
    env.check_pimms_support(config)

    # File paths for input and output, defined in config file
    proteinDataPath = os.path.join(REPO_ROOT, config["protPath"])
    metadataPath = os.path.join(REPO_ROOT, config["metaPath"])
    outPath = os.path.join(REPO_ROOT, config["outPath"])
    json_out = os.path.join(REPO_ROOT, config["json_outPath"])
    full_formula = config["DE_full_formula"]
    subset_formula = config["DE_subset_formula"]
    # Create the output directory
    make_outdir(outPath, make_subdirs=True)
    # Data processing
    print("Loading and processing data...")
    # metadata
    metadata = dp.process_data(
        file_path=metadataPath, json_out=json_out, outPath=outPath, config=config
    )
    # protein abundance data
    df_protAbundance = dp.process_data(
        file_path=proteinDataPath,
        metadata=metadata,
        json_out=json_out,
        outPath=outPath,
        config=config,
    )
    # if you have already processed the data and want to read in a previously saved version
    # import pandas as pd; df_protAbundance = pd.read_csv(os.path.join(outPath, "data/proteinAbundance.csv"))

    print("Data loaded and processed...")
    # Analysis
    print("Running analysis...")

    # if subsetting not required, go through with full datasets
    if config["analyse_full_dataset"] is True:
        full_outPath = os.path.join(outPath, "full_dataset")
        make_outdir(full_outPath)
        an.run_analysis(
            df=df_protAbundance,
            metadata=metadata,
            json_out=json_out,
            output_dir=full_outPath,
            config=config,
            formula=full_formula,
        )
        print("Full analysis complete.")
    if config["analyse_subsets"] is True:
        # Read subsets from config
        if not config["subsets"]:
            subset_terms = metadata[config["subset_variable"]].unique()
        if config["subsets"]:
            subset_terms = list(config["subsets"])
        # Loop through subsets
        for subset in subset_terms:
            print(f"Processing subset: {subset}")
            subset_variable = config["subset_variable"]
            ### find the rows in metadata that match the subset term
            subset_metadata = metadata.loc[
                metadata[subset_variable].astype(str) == str(subset),
            ]
            # Subset data based on index search term
            subset_df = get_subset(
                df=df_protAbundance,
                subset_term=subset,
                metadata=subset_metadata,
                subset_variable=subset_variable,
            )
            # Create a new output directory for the subset
            # if subset has a space in it, replace with _
            if " " in str(subset):
                subset = str(subset).replace(" ", "_")
            subset_outPath = os.path.join(
                outPath, "subsets", subset_variable + "_" + str(subset)
            )
            make_outdir(subset_outPath, make_subdirs=True)
            # Run analysis for the subset
            print(f"Running analysis for {subset}...")
            an.run_analysis(
                df=subset_df,
                metadata=subset_metadata,
                json_out=json_out,
                output_dir=subset_outPath,
                config=config,
                formula=subset_formula,
            )
        print("All subsets processed successfully.")

    print("generating html report...")
    generate_report_html()

    print("Analysis complete. Outputs saved to output directory. Tidying up intermediate files")
    
    # remove intermediate files
    tidy_up_files(outPath)


    
#### If executed in main script, run the function to produce the output
if __name__ == "__main__":
    main()
