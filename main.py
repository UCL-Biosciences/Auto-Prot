# Generate Proteomics Report

## Generate common proteomics outputs, save to file and combine in a report.
## Inputs must be protein abundance file and metadata, as specified below and in the github README.

## Import libraries
import os

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
    
    # Create the output directory
    dp.make_outdir(outPath)

    # Data processing
    print("Loading and processing data...")
    # metadata
    metadata = dp.preprocess_data(file_path=metadataPath,
                                  json_out=json_out,
                                  outPath = outPath)
    #dp.validate_metadata(metadata)

    # protein abundance data
    df_protAbundance, df_protAbundance_standardised = dp.preprocess_data(file_path=proteinDataPath,
                                                                         metadata=metadata,
                                                                         json_out=json_out,
                                                                         outPath = outPath)    
    #dp.validate_proteindata(df_protAbundance, metadata)

    print("Data loaded and processed...")

    # Analysis
    print("Running analysis...")
    analysis_results = an.run_analysis(df = df_protAbundance,
                                        df_standardised = df_protAbundance_standardised,
                                        metadata = metadata,
                                        json_out=json_out,
                                        output_dir = outPath)
    print("Analysis complete.")


#### If executed in main script, run the function to produce the output
if __name__ == "__main__":
    main()