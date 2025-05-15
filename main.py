# Generate Proteomics Report

## Generate common proteomics outputs, save to file and combine in a report.
## Inputs must be protein abundance file and metadata, as specified below and in the github README.

## Import libraries
import json
import os

import src.analysis.analysis as an
import src.processing.data_processing as dp

## Import functions
# Functions are saved in separate files and imported here
# Separated by module
# More detail in relevant files and on github.
import src.utils.check_env as env
from src.utils.data_io import make_outdir
from src.utils.data_utils import get_subset


##### Define main function for creating outputs
def main():
    """
    Main script to run the complete data analysis pipeline.
    """
    # Check environment - needs fixing
    if "VSCODE_PID" not in os.environ:
        env.compare_envs()
    # File paths
    # getting repo dir automatically is useful as should mean we don't need to specify
    # machine-specific paths and code should run on different users' machines
    REPO_ROOT = env.get_repo_root()
    
    ### path to config file containing key info - must be present!
    config_path = os.path.join(REPO_ROOT, "configs/auto-prot-config.json")
    ### Read in configuration data, stored in a json
    with open(config_path) as f:
        config = json.load(f)

    # File paths for input and output, defined in config file
    proteinDataPath = os.path.join(REPO_ROOT, config["protPath"])
    metadataPath = os.path.join(REPO_ROOT, config["metaPath"])
    outPath = os.path.join(REPO_ROOT, config["outPath"])
    json_out = os.path.join(REPO_ROOT, config["json_outPath"] )
    
    # Create the output directory
    make_outdir(outPath, make_subdirs=True)
    # Data processing
    print("Loading and processing data...")

    # metadata
    metadata = dp.process_data(
        file_path=metadataPath, json_out=json_out, outPath=outPath, config = config
    )
    # protein abundance data
    df_protAbundance = dp.process_data(
        file_path=proteinDataPath,
        metadata=metadata,
        json_out=json_out,
        outPath=outPath,
        config=config,
    )

    print("Data loaded and processed...")
    # Analysis
    print("Running analysis...")

    # if subsetting not required, go through with full datasets
    if config.get("analyse_full_dataset") is True:
        full_outPath = os.path.join(outPath, "full_dataset")
        make_outdir(full_outPath)
        an.run_analysis(
            df=df_protAbundance,
            metadata=metadata,
            json_out=json_out,
            output_dir=full_outPath,
            config=config,
        )
        print("Analysis complete.")
    if config.get("analyse_subsets") is True:
        # Read subsets from config
        subset_terms = config.get("subsets", [])
        # Loop through subsets
        for subset in subset_terms:
            print(f"Processing subset: {subset}")
            # Subset data based on index search term
            subset_df = get_subset(df_protAbundance, subset)
            # Create a new output directory for the subset
            subset_outPath = os.path.join(outPath, "subsets", subset.replace(" ", "_"))
            make_outdir(subset_outPath, make_subdirs=True)
            # Run analysis for the subset
            print(f"Running analysis for {subset}...")
            an.run_analysis(
                df=subset_df,
                metadata=metadata,
                json_out=json_out,
                output_dir=subset_outPath,
                config=config,
            )
        print("All subsets processed successfully.")


#### If executed in main script, run the function to produce the output
if __name__ == "__main__":
    main()
