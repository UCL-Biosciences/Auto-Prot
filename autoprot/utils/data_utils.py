### Data utils
## general functions for data handling

import fnmatch
import glob
import os
import warnings

import numpy as np
import pandas as pd
from pdf2image import convert_from_path
from PIL import Image
import math


def apply_row_id_config(df, config):
    if config["data_type"] != "phospho":
        return df  # Skip modification

    row_id_cfg = config.get("phospho_row_id", {})
    fields = [f.lower().replace(" ", ".") for f in row_id_cfg.get("fields", [])]
    missing = row_id_cfg.get("missing_value", "NA")

    if not fields:
        print("No fields specified for row ID construction; skipping.")
        return df

    print("Combining index with fields to make unique row IDs:", fields)

    # Safely get values from each field, fill missing if needed
    unique_phos_name = []
    for field in fields:
        print(field)
        if field in df.columns:
            part = df[field].fillna(missing).astype(str)
            unique_phos_name.append(part)
        else:
            print(
                f"Warning: field '{field}' not in DataFrame; not adding to unique row IDs"
            )  # if the field is missing

    # Combine extra parts with underscores
    extra_str = unique_phos_name[0]
    for part in unique_phos_name[1:]:
        extra_str += "_" + part

    # Combine with existing index
    df.index = df.index.astype(str) + "__" + extra_str

    return df


def normalise_column_names(df, file_path=None, outPath = None, config=None):
    """
    Standardise column names and optionally combine phosphorylation information.

    Converts all column names to lowercase with underscores, and if phosphoproteomic
    columns are present, combines them to create unique protein identifiers. Also sets
    index to the 'pg.genes' column for proteindata.

    Args:
        df (pd.DataFrame): Input dataframe to clean.
        file_path (str): Path used to detect 'proteindata' context.
        outPath (str): Output path for saving intermediate files.
        config (dict): Configuration dictionary with data type and row ID settings.

    Returns:
        pd.DataFrame: DataFrame with cleaned column names and possibly updated index.
    """
    ### standardise column names
    df.columns = [col.lower().replace(" ", "_") for col in df.columns] ## convert to lowercase and replace spaces with underscores
    df.columns = df.columns.astype(str) # ensure all columns are strings

    # If 'protein' is in the file path, set a column containing 'gene' as the index
    if "protein" in os.path.basename(file_path):

        genes_columns = [col for col in df.columns if "gene" in col.lower()]
        
        if genes_columns:  # If any column contains 'genes'
               
            # Set initial index
            ## genes can have multiple names separated by ; - we take the first non-empty one (sometimes the first is empty)
            # Keep NaN as NaN, only split for valid strings
            base_index = (
                df[genes_columns[0]]
                .apply(lambda x: next((v.strip() for v in str(x).split(";") if v.strip()), None)
                                if pd.notna(x) else np.nan)
            )

            # Count how many times each base value appears
            counts = base_index.value_counts()
            
            # Identify which ones are duplicated
            duplicated = base_index.isin(counts[counts > 1].index)
            
            # Create a suffix only for duplicates
            suffixes = (
                base_index[duplicated].groupby(base_index[duplicated]).cumcount() + 1
            ).astype(str)
            
            # Create the final index
            new_index = base_index.copy() # start with base index
            new_index[duplicated] = new_index[duplicated] + "-" + suffixes # then add suffixes to duplicates

            # Set the new index
            df.index = new_index
           
        ### Genes with no name in original data ###
        ## if there is no gene name, we need to replace the NaN with Unknown-Gene-X, where X is a unique number for each unknown gene.

        # Convert index to a Series to manipulate NaNs
        index_series_original = df.index.to_series().astype("object")  ## 'object' allows strings and NAs

        # Strip whitespace and handle empty strings
        index_series = index_series_original.astype(str).str.strip()

        # Find NaN values in index
        nan_mask = index_series.isna() | (index_series == "") | (index_series == "nan")

        # Replace NaNs with "Unknown-gene-N"
        index_series[nan_mask] = [f"Unknown-gene-{i+1}" for i in range(nan_mask.sum())]

        # Set updated index
        df.index = index_series

        ### mapping original index to new index ###
        # we want to keep track of which original index maps to which new index
        # this is useful for debugging and for tracking proteins through the analysis
        # at this point in the pipeline, the old names are in the column with genes in the name (this will soon be removed)
        # and the new names are the index. so it is a good time to save
        (df.loc[ df.index != df[genes_columns[0]], ]
         .rename(columns={'pg.genes': 'pg.genes_original'})
         .to_csv(os.path.join(outPath, "data/prots_name_mapping.csv"),
                                                           index = True)
        )

        ### For phosphoproteomic data, there are abundances for phosphorylated proteins
        ### Each protein can be present multiple times - once per phosphorylation state
        ### For the analysis to proceed, we need a unique ID for the protein-phosphorylation state combination
        ## we append the phosporylation state (in column PTM.ModificationTitle and PTM.SiteAA) to the pg.genes column
        if config["data_type"] == "phospho":
            print(
                "Gene names and phosphorylation state present. Combining to make unique gene names"
            )
            df = apply_row_id_config(df, config)

    return df


## helper function for subsetting based on phosphoproteomics. data to subset for must be indicated in protein abundance index, and must match subset term in config.
def get_subset(df, subset_term, metadata, subset_variable):
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
    ### find the rows that match the subset term, then take the sample_rep column (protein intensities are named after sample_rep)
    matching_sample_ids = metadata.loc[
        metadata[subset_variable].astype(str) == str(subset_term), "sample_rep"
    ].tolist()

    subset_df = df[matching_sample_ids]
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
    if metadata[["sample_rep"]].duplicated().any():
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
        warnings.warn(
            "Warning: The protein abundance df contains missing (NaN) values!"
        )


def combine_plots(search_path, search_term, output_dir):
    """
    Combines plots matching a search term into a single image.

    Args:
    - search_path (str): Path to search for images.
    - search_term (str): Term to search for in image filenames.
    - max_cols (int): Maximum number of columns in the grid.
    - max_combined_width (int): Maximum width of the combined image.
    - max_combined_height (int): Maximum height of the combined image.

    Returns:
    - str: Path to the saved combined image, or None if no images found.
    """

    if os.name == "nt":
        # Windows long path workaround
        search_path = "\\\\?\\" + os.path.abspath(search_path)

    # Find matching image paths
    image_paths = []
    # Use os.walk to find all files matching the search term
    # extract the root directory and file names
    for root, _, files in os.walk(search_path):
        for file in files:
            # Check if the file matches the search term and is not a combined plot (don't want combined volcano plot inception)
            if search_term in file and not fnmatch.fnmatch(file, "combined_*_plot.png") and not fnmatch.fnmatch(file, "*.csv"):
                # add the path to the image_paths list
                image_paths.append(os.path.join(root, file))

    if not image_paths:
        print(f"No plots found for '{search_term}'.")
        return None

    # Load images at original size
    images = []
    for path in image_paths:
        if path.endswith(".pdf"):
            img = convert_from_path(path, dpi=300)[0]  # Use first page
        else:
            img = Image.open(path)
        images.append(img)

    # Width is the widest image 
    total_width = max(img.width for img in images)
    # Height is the sum of heights
    total_height = sum(img.height for img in images)

    # Determine grid layout
    combined_image = Image.new("RGB", (total_width, total_height), (255, 255, 255))

    y_offset = 0
    for img in images:
        combined_image.paste(img, (0, y_offset))
        y_offset += img.height
    
    # reduce colour depth to reduce file size
    combined_image = combined_image.convert("P", palette=Image.ADAPTIVE)

    # Save the combined image
    combined_image_path = os.path.join(output_dir, f"plots/combined_{search_term}_plot.png")
    combined_image.save(combined_image_path, optimize = True)
    print(f"Combined image saved to {combined_image_path}")
    return combined_image_path



### for combining data from different treatments for display in the report ###
def combine_csv_files(
    filename,
    output_dir,
    output_filename=None,
    top_n=10,
    new_column="treatment_pair",
    sort_by_logfc=False,
):
    """
    General function to combine CSV files from subdirectories into a single file.

    Parameters:
    - filename (str): The name of the CSV file to search for (e.g., "top_20_by_LFC.csv").
    - output_dir (str): The root directory where data folders are stored.
    - output_filename (str or None): The output filename for the combined CSV.
                                     If None, it's auto-generated based on `filename`.
    - top_n (int): Number of rows to take from each CSV file.
    - new_column (str): Column name to store the extracted folder name (e.g., "treatment_pair").
    - sort_by_logfc (bool): Whether to sort each file by absolute logFC (default: False).


    Returns:
    - pd.DataFrame: The combined DataFrame.
    - str: The path where the final CSV is saved.
    """
    # Find all matching images
    csv_files = []
    for root, _, files in os.walk(output_dir):
        for file in files:
            if fnmatch.fnmatch(file, filename):  # supports wildcards like * and ?
                csv_files.append(os.path.join(root, file))
    csv_files = [
        csv
        for csv in csv_files
        if not fnmatch.fnmatch(os.path.basename(csv), "combined*.csv")
    ]
    # # Search for matching CSV files in subdirectories
    # search_pattern = os.path.join(output_dir, "data", "**", filename)
    # csv_files = glob(os.path.join(output_dir, "data", "*", filename)) + glob(os.path.join(output_dir, "data", "*", "*", filename))
    # csv_files = sorted(glob(search_pattern, recursive=True))
    # check if csv files exist
    if not csv_files:
        print(f"No files found matching '{filename}'.")
        return None, None
    combined_data = []
    # Process each found CSV file
    for file in csv_files:
        if os.name == "nt":
            file = "\\\\?\\" + os.path.abspath(file)
        # Extract the folder name (used as the category column)
        folder_name = os.path.basename(os.path.dirname(file))
        # Read CSV and select top `n` rows
        df = pd.read_csv(file)
        if sort_by_logfc and "logFC" in df.columns:
            df = df.sort_values(by="logFC", ascending=False, key=abs)
        df = df.head(top_n)
        # Add the extracted folder name as a new column
        df[new_column] = folder_name
        # Append to the list
        combined_data.append(df)
    # Merge all data
    combined_df = pd.concat(combined_data, ignore_index=True)
    # Generate output filename if not provided
    if output_filename is None:
        output_filename = os.path.join(output_dir, f"data/combined_{filename}")
        ## windows sometimes rejects long paths. Workaround:
        if os.name == "nt":
            output_filename = "\\\\?\\" + os.path.abspath(output_filename)
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    # Save combined CSV

    combined_df.to_csv(output_filename, index=False)
    print(f"Combined file saved at: {output_filename}")


def tidy_up_files(outPath) -> list:
    """
    Remove temporary files created during processing.

    Specifically removes:
    - 'prots_name_mapping.csv' from the 'data' directory.

    Returns:
        List of relative paths (relative to outPath) of files that were deleted.
    """
    files_to_remove = [
        glob.glob(os.path.join(outPath, "full_dataset/data/*/limma_output.csv")),
        glob.glob(os.path.join(outPath, "full_dataset/data/*/prots.csv")),
        glob.glob(os.path.join(outPath, "full_dataset/data/*/top_20_by_LFC.csv")),
    ]

    # glob returns a list of files so files_to_remove is currently a list of lists
    # flatten the list of lists
    flattened_files_to_remove = []
    for sublist in files_to_remove: # loop through list of lists
        for file in sublist: # take each file
            flattened_files_to_remove.append(file)

    # remove each file if it exists and track what was deleted
    deleted = []
    for file in flattened_files_to_remove:
        if os.path.exists(file):
            os.remove(file)
            deleted.append(os.path.relpath(file, outPath))

    print("Temporary files removed.")
    return deleted