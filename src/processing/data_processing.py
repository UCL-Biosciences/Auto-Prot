### Read in and pre-process data ###
import json
import logging
import os

import numpy as np
import pandas as pd
import seaborn as sns
import yaml

import src.processing.data_preprocess as dpp
from src.utils.data_io import load_data
from src.utils.data_utils import (
    normalise_column_names,
    validate_metadata,
    validate_proteindata,
)

# specify location of errors to standard output
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s - %(levelname)s - %(message)s"
)


#### a few steps for cleaning metadata
## rename protein abundance vars, make a unique sample_replicate id, create a colour for each treatment group
def clean_meta(df, json_out):
    """
    Clean metadata for downstream use.

    Standardises string fields, generates a unique sample + replicate ID (`sample_rep`), assigns
    colour labels to treatment groups, and writes summary info to a JSON file.

    Args:
        df (pd.DataFrame): Metadata containing sample and treatment info.
        json_out (str): Path to JSON file where summary info will be saved.

    Returns:
        pd.DataFrame: Cleaned and augmented metadata.
    """
    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    df["treatment"] = df["treatment"].astype(str).str.strip()
    df["protein_abundance_name"] = (
        df["protein_abundance_name"].str.lower().str.replace(" ", "_")
    )
    df["sample_rep"] = (
        df["sample_id"] + "_r" + df["replicate"].astype(str)
    )  # for a unique ID for each sample_replicate
    # if time point present - include that to ensure sample_rep stays unique
    if "timepoint" in df.columns:
        df["sample_rep"] += "_t" + df["timepoint"].astype(str)
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
    """
    Filter and rename protein abundance columns using metadata.

    Ensures only valid sample columns are kept, converts non-numeric entries to NaN,
    and simplifies long sample names using metadata.

    Args:
        df (pd.DataFrame): Raw protein abundance DataFrame.
        metadata (pd.DataFrame): Metadata containing valid sample names and mappings.

    Returns:
        Tuple[pd.DataFrame, int]: Cleaned DataFrame and original number of proteins (rows).
    """
    # Filter columns based on metadata['protein_abundance_name']
    # extract column names to keep
    valid_columns = metadata["protein_abundance_name"].tolist()
    # filter df to keep prot abundance columns
    df = df.loc[:, df.columns.isin(valid_columns)]
    # At this point, remove rows fully duplicated including the index
    df = df[~df.reset_index().duplicated(keep="first").values]
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
    """
    Append summary statistics about the protein data to a JSON file.

    Calculates number of proteins retained/removed and min/max/median mean abundance per protein.

    Args:
        df (pd.DataFrame): Cleaned protein abundance DataFrame.
        nrow_original (int): Number of proteins before filtering.
        json_out (str): Path to JSON file to update.
    """
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
        existing_data = yaml.safe_load(f)

    # Append new data
    existing_data.update(abundance_stats)

    # Write back to JSON file
    with open(json_out, "w") as f:
        json.dump(existing_data, f, indent=4)


def clean_data(
    df, file_path=None, metadata=None, outPath=None, config=None, json_out=None
):
    """
    Main cleaning function for either metadata or protein abundance data.

    Cleans metadata (if file_path includes 'metadata') or processes protein data
    with transformation, normalisation, imputation, and summary export.
    Renames unknown genes from NA to "Unknown-gene-1" and removes duplicate rows.
    Produces plots of distributions after different processing steps.

    Args:
        df (pd.DataFrame): Raw data to clean.
        file_path (str): Path or identifier indicating file type.
        metadata (pd.DataFrame, optional): Metadata required for protein data cleaning.
        outPath (str, optional): Directory where plots will be saved.
        config (dict, optional): Configuration dict, must include 'df_to_use'.
        json_out (str, optional): Path to JSON file to write summary info.

    Returns:
        pd.DataFrame | Tuple[pd.DataFrame, int]: Cleaned data, plus original row count if protein data.

    Raises:
        ValueError: If metadata is missing or incomplete when required.
    """
    if "metadata" in file_path:
        df = clean_meta(df=df, json_out=json_out)
    if "protein" in file_path:
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
            dfs = dpp.process_prot_data(df, config, metadata=metadata, outPath=outPath)
            ## to be shown when plotting distributions
            plot_titles = [
                "Raw Intensities",
                "Log₂ Transformed",
                "Normalised",
                "Normalised then Imputed",
            ]
            ## note dfs is name of df and data itself. .values() access the data
            dpp.view_prot_distributions(dfs.values(), plot_titles, metadata, outPath)
            ## which df to use?
            df_to_use = config["df_to_use"]
            df = dfs[df_to_use].dropna()
            ### save some summary info to file for report
            prot_summary(df, nrow_original, json_out)
    df = df.drop_duplicates()
    return df


###########################################################################
#### This is the main function for processing data, combining the above ###
###########################################################################


def process_data(file_path, metadata=None, json_out=None, outPath=None, config=None):
    """
    End-to-end data preprocessing pipeline.

    Loads and standardises data, cleans based on file type,
    and validates metadata or protein abundance input.

    Args:
        file_path (str): Path to CSV or TSV file to load.
        metadata (pd.DataFrame, optional): Metadata to use for proteindata cleaning.
        json_out (str, optional): Path to write summary info.
        outPath (str, optional): Directory to save plots and output.
        config (dict, optional): Configuration dict for processing logic.

    Returns:
        pd.DataFrame: Cleaned and validated data.
    """
    df_in = load_data(file_path)
    if df_in is not None:
        df_renamed = normalise_column_names(
            df=df_in, file_path=file_path, config=config
        )
    if "metadata" in file_path:
        ### clean metadata
        df = clean_data(
            df=df_renamed, file_path=file_path, config=config, json_out=json_out
        )
        ### run function to validate metadata
        validate_metadata(df)

    if "protein" in file_path:
        df = clean_data(
            df=df_renamed,
            file_path=file_path,
            metadata=metadata,
            outPath=outPath,
            config=config,
            json_out=json_out,
        )
        ### run function to validate protein abundance data
        validate_proteindata(data=df, metadata=metadata)

        # protein data takes a while to produce and you might want to read or look without processing every time
        # write to file for easy access
        os.makedirs(os.path.join(outPath, "data"), exist_ok=True)
        df.to_csv(os.path.join(outPath, "data/proteinAbundance.csv"), index=True)

    return df
