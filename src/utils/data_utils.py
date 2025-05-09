
### Data utils
## general functions for data handling

import pandas as pd

def normalise_column_names(df, file_path=None):
    """
    Standardise column names and optionally combine phosphorylation information.

    Converts all column names to lowercase with underscores, and if phosphoproteomic
    columns are present, combines them to create unique protein identifiers. Also sets
    index to the 'pg.genes' column for proteindata.

    Args:
        df (pd.DataFrame): Input dataframe to clean.
        file_path (str, optional): Path used to detect 'proteindata' context.

    Returns:
        pd.DataFrame: DataFrame with cleaned column names and possibly updated index.
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
            + "__" # having double __ after genes allows easier splitting later on
            + df["ptm.modificationtitle"]
            + "_"
            + df["ptm.siteaa"]
            + "_"
            + df["ptm.sitelocation"].astype(str)
        )
    # If 'proteindata' is in the file path, set a column containing 'genes' as the index
    if "proteindata" in file_path:
        genes_columns = [col for col in df.columns if "pg.genes" in col.lower()]
        if genes_columns:  # If any column contains 'pg.genes'
            df = df.set_index(genes_columns[0])
    return df

## helper function for subsetting based on phosphoproteomics. data to subset for must be indicated in protein abundance index, and must match subset term in config.
def get_subset(df, subset_term):
    """
    Subset a DataFrame by searching for a term in the index.

    Args:
        df (pd.DataFrame): Protein abundance dataframe.
        subset_term (str): Substring to match in the index.

    Returns:
        pd.DataFrame: Subset of the original DataFrame.

    Raises:
        ValueError: If no index values contain the subset term.
    """

    subset_df = df[df.index.str.contains(subset_term, regex=False)]
    if subset_df.empty:
        raise ValueError(f"No matches found for subset: {subset_term}")
    return subset_df



### Validation metadata function ###
# checks column names, sample id duplicates, replicate is numeric, sample data columns are unique
def validate_metadata(metadata):
    """
    Validate that metadata meets all required structure and constraints.

    Checks for:
    - Required columns present
    - No missing values
    - Unique sample_id + replicate combinations
    - Numeric replicate values
    - Unique protein_abundance_name entries

    Args:
        metadata (pd.DataFrame): The metadata to validate.

    Raises:
        ValueError: If any validation checks fail.
    """
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
    """
    Validate structure and contents of a protein abundance dataframe.

    Checks include:
    - Data is not empty
    - Unique, case-insensitive row and column identifiers
    - Columns match sample_rep from metadata
    - Abundance values are numeric
    - No missing values

    Args:
        data (pd.DataFrame): Cleaned protein abundance data.
        metadata (pd.DataFrame): Corresponding metadata containing sample_rep.

    Raises:
        ValueError: If any check fails.
    """
    ##### checks for raw protein abundance data #####
    # Check if data is empty
    if data.empty:
        raise ValueError("Error: The protein abundance data is empty.")
    # Check if row indices (proteins) are unique
    index_str = pd.Index(data.index.map(str))
    if not index_str.str.lower().is_unique:
        raise ValueError("Error: Protein identifiers (row indices) must be unique.")
    # Check if column names (samples) are unique
    col_str = pd.Index(data.columns.map(str))
    if not col_str.str.lower().is_unique:
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

