
### Data utils
## general functions for data handling

import pandas as pd

def normalise_column_names(df, file_path=None, metadata=None):
    """
    Standardise column names by making them lowercase and replacing spaces with underscores.
    Note, need to do this in metadata cols too - see preprocess_data

    Parameters:
    - df (pd.DataFrame): The dataframe to modify.

    Returns:
    - pd.DataFrame: Dataframe with normalized column names.
    """
    df.columns = [col.lower().replace(" ", "_") for col in df.columns]
    df.columns = df.columns.astype(str)
    ### For phosphoproteomic data, there are abundances for phosphorylated proteins
    ### Each protein can be present multiple times - once per phosphorylation state
    ### For the analysis to proceed, we need a unique ID for the protein-phosphorylation state combination
    ## we append the phosporylation state (in column PTM.ModificationTitle and PTM.SiteAA) to the pg.genes column
    if (
        ("ptm.modificationtitle" in df)
        and ("ptm.siteaa" in df)
        and ("ptm.sitelocation" in df)
        and ("pg.genes" in df)
    ):
        print(
            "Gene names and phosphorylation state present. Combining to make unique gene names"
        )
        df["pg.genes"] = (
            df["pg.genes"]
            + "__"
            + df["ptm.modificationtitle"]
            + "_"
            + df["ptm.siteaa"]
            + "_"
            + df["ptm.sitelocation"].astype(str)
        )
    # If 'proteindata' is in the file path, set a column containing 'genes' as the index
    if "proteindata" in file_path:
        genes_columns = [col for col in df.columns if "genes" in col.lower()]
        if genes_columns:  # If any column contains 'genes'
            df = df.set_index(genes_columns[0])
    return df

## helper function for subsetting based on phosphoproteomics. data to subset for must be indicated in protein abundance index, and must match subset term in config.
def get_subset(df, subset_term):
    """
    Subsets a dataframe based on a search term in the index.
    Raises ValueError if no matches are found.
    """
    subset_df = df[df.index.str.contains(subset_term, regex=False)]
    if subset_df.empty:
        raise ValueError(f"No matches found for subset: {subset_term}")
    return subset_df



### Validation metadata function ###
# checks column names, sample id duplicates, replicate is numeric, sample data columns are unique
def validate_metadata(metadata):
    required_columns = ["sample_id", "treatment", "replicate", "protein_abundance_name"]
    # Check if all required columns are present
    missing_columns = [col for col in required_columns if col not in metadata.columns]
    if missing_columns:
        raise ValueError(f"Error: Missing required columns: {missing_columns}")
    # Check for missing values (NaNs) in the entire DataFrame
    if metadata.isna().any().any():
        raise ValueError("Error: The metadata contains missing (NaN) values!")
    # Check if the combination of 'sample_id' and 'replicate' is unique
    if metadata[["sample_id", "replicate"]].duplicated().any():
        raise ValueError("Error: Each (sample_id, replicate) pair must be unique.")
    # Check if 'replicate' is numeric
    if not pd.api.types.is_numeric_dtype(metadata["replicate"]):
        raise ValueError("Error: 'replicate' must be numeric.")
    # Check for duplicate protein_abundance_name values
    if metadata["protein_abundance_name"].duplicated().any():
        raise ValueError(
            "Error: 'protein_abundance_name' must not contain duplicate values."
        )


### Validation protein abundance data ###
def validate_proteindata(data, metadata):
    ##### checks for raw protein abundance data #####
    # Check if data is empty
    if data.empty:
        raise ValueError("Error: The protein abundance data is empty.")
    # Check if row indices (proteins) are unique
    if not data.index.str.lower().is_unique:
        raise ValueError("Error: Protein identifiers (row indices) must be unique.")
    # Check if column names (samples) are unique
    if not data.columns.str.lower().is_unique:
        raise ValueError("Error: Sample identifiers (columns) must be unique.")
    # Validate protein abundance columns from metadata
    # names of abundance columns:
    abundance_columns = metadata[
        "sample_rep"
    ]  # Columns to check. Note this depends on whether cols have been renamed
    # Ensure one-to-one mapping: all sample_id values in columns, and all columns are sample_ids
    if set(abundance_columns) != set(data.columns):
        raise ValueError(
            "Error: Mismatch between metadata samples and protein abundance columns. Ensure a strict one-to-one mapping."
        )
    # If all protein_abundance_name columns are present,
    # Check numeric values only in specified abundance columns
    non_numeric_columns = [
        col for col in abundance_columns if not pd.api.types.is_numeric_dtype(data[col])
    ]
    if non_numeric_columns:
        raise ValueError(
            f"Error: The following columns contain non-numeric values: {non_numeric_columns}"
        )
    # Check for missing values (NaNs) in the entire DataFrame
    if data.isna().any().any():
        raise ValueError(
            "Error: The protein abundance df contains missing (NaN) values!"
        )

