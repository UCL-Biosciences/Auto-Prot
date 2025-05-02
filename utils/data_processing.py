### Read in and pre-process data ###

import json
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn
import sklearn.preprocessing


from utils.data_io import load_data
from utils.data_utils import normalise_column_names, validate_metadata, validate_proteindata
import utils.data_preprocess as dpp


# specify location of errors to standard output
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s - %(levelname)s - %(message)s"
)

#### a few steps for cleaning metadata
## rename protein abundance vars, make a unique sample_replicate id, create a colour for each treatment group
def clean_meta(df, json_out):
    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    df["treatment"] = df["treatment"].astype(str).str.strip()
    df["protein_abundance_name"] = (
            df["protein_abundance_name"].str.lower().str.replace(" ", "_")
        )
    df["sample_rep"] = (
        df["sample_id"] + "_" + df["replicate"].astype(str)
    )  # for a unique ID for each sample_replicate
    ### for plotting, assign a colour to each unique treatment group
    colours = sns.color_palette(
        "colorblind", len(df["treatment"].unique())
    ).as_hex()  # list of colour-blind-friendly colours as long as number of unique treatment levels
    # colours = px.colors.qualitative.Alphabet[:len(df['treatment'].unique())]
    col_map = dict(
        zip(df["treatment"].unique(), colours)
    )  # zip (combine) the treatments to the colours, store as a dictionary
    df["colours"] = df["treatment"].map(
        col_map
    )  # create a new var, mapping treatment in metadata to the colours in the dictionary
    # generate some info to go into html report
    # Compute treatment-wise sample counts
    treatment_counts = df["treatment"].value_counts().to_dict()
    # Format as "group 1: n1, group 2: n2"
    treatment_summary = ", ".join(
        [f"{treatment}: {n}" for treatment, n in treatment_counts.items()]
    )
    metadata_values = {
        "NUM_SAMPLES": df["sample_rep"].nunique(),
        "NUM_TREATMENTS": df["treatment"].nunique(),
        "TREATMENTS": treatment_summary,
    }
    # Save to a JSON file
    with open(json_out, "w") as f:
        json.dump(metadata_values, f)
    return df

### clean prot data
def clean_prot(df, metadata):
    # Filter columns based on metadata['protein_abundance_name']
    # extract column names to keep
    valid_columns = metadata["protein_abundance_name"].tolist()
    # filter df to keep prot abundance columns
    df = df.loc[:, df.columns.isin(valid_columns)]
    # protein columns should have only numeric data
    # convert non-numeric values to NaN and print warning message
    if not df.equals(df.select_dtypes(include=[np.number])):
        print(
            "Warning: DataFrame contains non-numeric values! Converting to NaN: these proteins will be imputed!"
        )
        df = df.apply(pd.to_numeric, errors="coerce")
    ### protein columns can have long names - better to have just sample name
    # Create a mapping of old column names to new column names
    rename_mapping = dict(
        zip(metadata["protein_abundance_name"], metadata["sample_rep"])
    )
    # Rename columns in data_in based on the mapping
    df = df.rename(columns=rename_mapping)
    nrow_original = len(df.index)
    return df, nrow_original

def prot_summary(df, nrow_original, json_out):
    #### protein summary
    NUM_PROTS = len(df.index)
    NUM_PROTS_REMOVED = nrow_original - NUM_PROTS

    # Compute mean value per row (per protein across all samples)
    mean_abundance = df.mean(axis=1)

    # Compute min, max, and mean of these mean abundances
    abundance_stats = {
        "NUM_PROTS_OG": nrow_original,
        "NUM_PROTS_REMOVED": NUM_PROTS_REMOVED,
        "NUM_PROTS": f"{NUM_PROTS:,.0f}",
        "MIN_AVERAGE_ABUNDANCE": f"{mean_abundance.min():,.0f}",
        "MAX_AVERAGE_ABUNDANCE": f"{mean_abundance.max():,.0f}",
        "MEDIAN_AVERAGE_ABUNDANCE": f"{mean_abundance.median():,.0f}",
    }

    # read data from json file
    with open(json_out) as f:
        existing_data = json.load(f)

    # Append new data
    existing_data.update(abundance_stats)

    # Write back to JSON file
    with open(json_out, "w") as f:
        json.dump(existing_data, f, indent=4)     


def clean_data(df, file_path=None, metadata=None, outPath=None, config=None, json_out=None):
    """
    Perform basic cleaning on the data: first filter protein data to include only protein abundance cols
    and make sure they are all present.
    Non-numeric protein abundances are converted to NaNs, then duplicate rows or rows with NAs are removed.

    Parameters:
    - df (pd.DataFrame): The dataframe to clean.
    - metadata (pd.DataFrame): Optional metadata dataframe containing 'protein_abundance_name'.

    Returns:
    - pd.DataFrame: Cleaned dataframe.

    Raises:
    - ValueError: If metadata is not provided or does not contain 'protein_abundance_name'.

    """
    if "metadata" in file_path:
        df = clean_meta(df = df, json_out = json_out)
    if "proteindata" in file_path:
        if metadata is None:
            raise ValueError("Error: Metadata is required but not provided.")
        # Check if the required column is in metadata
        elif "protein_abundance_name" not in metadata.columns:
            raise ValueError(
                "Error: 'protein_abundance_name' column is missing in the metadata."
            )
        else:
            df, nrow_original = clean_prot(df, metadata)
            ### pre process protein abundance data
            ## replace 0 with NA, remove prots with lots of missing data
            ## log2 transform, normalise using each sample's median, and impute using random forest
            dfs = dpp.process_prot_data(df, metadata, config)
            ## to be shown when plotting distributions
            plot_titles = [
                "Raw Intensities",
                "Log₂ Transformed",
                "Sample-Median Normalised",
                "Normalised then Imputed",
            ]
            ## note dfs is name of df and data itself. .values() access the data
            dpp.view_prot_distributions(dfs.values(), plot_titles, metadata, outPath)
            ## which df to use?
            df_to_use = config.get("df_to_use")
            df = dfs[df_to_use]
            # some proteins do not produce any associated genes. these values are left blank in the index
            # we replace the NaNs with Unknown-Gene-X, where X is a unique number for each unknown gene.\
            # Convert index to a Series to manipulate NaNs
            index_series = df.index.to_series()
            # Find NaN values in index
            nan_mask = index_series.isna()
            # Replace NaNs with "Unknown-gene-N"
            index_series[nan_mask] = [
                f"Unknown-gene-{i+1}" for i in range(nan_mask.sum())
            ]
            # Set updated index
            df.index = index_series
            # drop rows with duplicated index values <<<< This is a weird quirk of the phosphoproteomics data. Need to find a way of uniquely identifying rows. Added to github issue
            df = df[~df.index.duplicated(keep="first")]
            ### save some summary info to file for report
            prot_summary(df, nrow_original, json_out)

    df = df.drop_duplicates()
    if "proteindata" in file_path:
        return df, nrow_original
    else:
        return df


###########################################################################
#### This is the main function for processing data, combining the above ###
###########################################################################

def process_data(file_path, metadata=None, json_out=None, outPath=None, config=None):
    """
    Full preprocessing pipeline: load, clean, and normalize column names.

    Parameters:
    - file_path (str): Path to the data file.

    Returns:
    - pd.DataFrame: Preprocessed dataframe.
    """
    df = load_data(file_path)
    if df is not None:
        df = normalise_column_names(df, file_path=file_path)
    if "metadata" in file_path:
        ### clean metadata
        df = clean_data(df, file_path=file_path, config=config, json_out=json_out)
        ### run function to validate metadata
        validate_metadata(df)

    if "proteindata" in file_path:
        df, nrow_original = clean_data(
            df, file_path=file_path, metadata=metadata, outPath=outPath, config=config, json_out=json_out
        )
        ### run function to validate protein abundance data
        validate_proteindata(data=df, metadata=metadata)
    return df


