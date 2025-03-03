# Generate Proteomics Report

## Generate common proteomics outputs, save to file and combine in a report.
## Inputs must be protein abundance file and metadata, as specified below and in the github README.

## Import libraries
import os
import json

## Import functions
# Functions are saved in separate files and imported here
# Separated by module
# More detail in relevant files and on github.
import utils.check_env as env
import utils.data_processing as dp
import utils.analysis as an

##### Define main function for creating outputs
def main():
    """
    Main script to run the complete data analysis pipeline.
    """
    # Check environment
    if not "VSCODE_PID" in os.environ:
        env.compare_envs()
    # File paths
    # getting repo dir automatically is useful as should mean we don't need to specify
    # machine-specific paths and code should run on different users' machines
    REPO_ROOT = env.get_repo_root()
    # File path to the dataset
    proteinDataPath = os.path.join(REPO_ROOT, 'input/data/proteindata.csv')
    metadataPath = os.path.join(REPO_ROOT, 'input/data/metadata.csv')
    outPath = os.path.join(REPO_ROOT, 'output')
    json_out = os.path.join(REPO_ROOT, 'output/data/data_for_report.json')
    config_path = os.path.join(REPO_ROOT, 'configs/auto-prot.json')
    # Create the output directory
    dp.make_outdir(outPath, make_subdirs = True)
    # Data processing
    print("Loading and processing data...")

    # metadata
    metadata = dp.preprocess_data(file_path=metadataPath,
                                  json_out=json_out,
                                  outPath = outPath)
    # protein abundance data
    df_protAbundance, df_protAbundance_standardised = dp.preprocess_data(file_path=proteinDataPath,
                                                                         metadata=metadata,
                                                                         json_out=json_out,
                                                                         outPath = outPath)    
    ### Read in configuration data, stored in a json
    with open(config_path, "r") as f:
        config = json.load(f)
    print("Data loaded and processed...")
    # Analysis
    print("Running analysis...")

    # if subsetting not required, go through with full datasets
    if config.get("analyse_full_dataset") is True:
        full_outPath=os.path.join(outPath, "full_dataset")
        dp.make_outdir(full_outPath)
        analysis_results = an.run_analysis(df = df_protAbundance,
                                            df_standardised = df_protAbundance_standardised,
                                            metadata = metadata,
                                            json_out=json_out,
                                            output_dir = full_outPath,
                                            config = config)
        print("Analysis complete.")
    if config.get("analyse_subsets") is True:
        # Read subsets from config
        subset_terms = config.get("subsets", [])
    # Loop through subsets
        for subset in subset_terms:
            print(f"Processing subset: {subset}")
            # Subset data based on index search term
            subset_df = df_protAbundance[df_protAbundance.index.str.contains(subset, regex=False)]
            subset_df_standardised = df_protAbundance_standardised.loc[:, df_protAbundance_standardised.columns.str.contains(subset, regex=False)]
            # Raise an error if no matching rows are found
            if subset_df.empty:
                raise ValueError(f"No matches found for subset: {subset}")
            # Create a new output directory for the subset
            subset_outPath = os.path.join(outPath, 'subsets', subset.replace(" ", "_"))
            dp.make_outdir(subset_outPath)
            # Run analysis for the subset
            print(f"Running analysis for {subset}...")
            analysis_results = an.run_analysis(
                df=subset_df,
                df_standardised=subset_df_standardised,
                metadata=metadata,
                json_out=json_out,
                output_dir=subset_outPath,
                config = config
            )
        print("All subsets processed successfully.")

#### If executed in main script, run the function to produce the output
if __name__ == "__main__":
    main()