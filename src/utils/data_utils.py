### Data utils
## general functions for data handling

import fnmatch
import os

import pandas as pd

from PIL import Image

def apply_row_id_config(df, config):
    if config.get("data_type") != "phospho":
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
        if field in df.columns:
            part = df[field].fillna(missing).astype(str)
        else:
            print(f"Warning: field '{field}' not in DataFrame; not adding to unique row IDs") # if the field is missing
        unique_phos_name.append(part)

    # Combine extra parts with underscores
    extra_str = unique_phos_name[0]
    for part in unique_phos_name[1:]:
        extra_str += "_" + part

    # Combine with existing index
    df.index = df.index.astype(str) + "__" + extra_str

    return df


def normalise_column_names(df, file_path=None, config = None):
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
    # If 'protein' is in the file path, set a column containing 'gene' as the index
    if "protein" in file_path:
        genes_columns = [col for col in df.columns if "gene" in col.lower()]
        if genes_columns:  # If any column contains 'genes'
            df = df.set_index(df[genes_columns[0]].str.split(";").str[0])
        # some proteins do not produce any associated genes. these values are left blank in the index
        # we replace the NaNs with Unknown-Gene-X, where X is a unique number for each unknown gene.
        # Convert index to a Series to manipulate NaNs
        index_series = df.index.to_series().astype(
            "object"
        )  ## 'object' allows strings and NAs
        # Find NaN values in index
        nan_mask = index_series.isna()
        # Replace NaNs with "Unknown-gene-N"
        index_series[nan_mask] = [
            f"Unknown-gene-{i+1}" for i in range(nan_mask.sum())
        ]
        # Set updated index
        df.index = index_series
        ### For phosphoproteomic data, there are abundances for phosphorylated proteins
        ### Each protein can be present multiple times - once per phosphorylation state
        ### For the analysis to proceed, we need a unique ID for the protein-phosphorylation state combination
        ## we append the phosporylation state (in column PTM.ModificationTitle and PTM.SiteAA) to the pg.genes column
        if config.get("data_type") == "phospho":
            print("Gene names and phosphorylation state present. Combining to make unique gene names" )
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
    matching_sample_ids = metadata.loc[metadata[subset_variable].astype(str) == str(subset_term), "sample_rep"].tolist()
    ## subset target columns
    print(f"subset_variable: {subset_variable}")
    print(f"subset_term: {subset_term}")
    print("metadata[subset_variable].unique():", metadata[subset_variable].unique())
    print("metadata.shape:", metadata.shape)
    print("metadata.head():\n", metadata.head())
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
        raise ValueError(
            "Error: The protein abundance df contains missing (NaN) values!"
        )

#### combine plots made for different treatment groups
def combine_plots(
    search_term,
    search_path,
    output_dir,
    output_filename=None,
    img_size=(800, 600),
    max_cols=3,
):
    """
    Finds all images matching the search_term in subdirectories of search_path,
    arranges them in a grid with up to max_cols columns, and saves a single PNG.

    Parameters:
    - search_term (str): The filename pattern to search for (e.g., "volcano_plot.png").
    - search_path (str): The base directory where images are stored.
    - output_filename (str): The filename for the combined image (auto-generated if None).
    - img_size (tuple): (width, height) to resize images.
    - max_cols (int): Maximum number of columns in the grid.

    Returns:
    - str: Path to the saved combined image, or None if no images found.
    """
    if os.name == "nt":
        search_path = "\\\\?\\" + os.path.abspath(search_path)
    # Find all matching images
    image_paths = []
    for root, _, files in os.walk(search_path):
        for file in files:
            if search_term in file:
                image_paths.append(os.path.join(root, file))
    image_paths = [img for img in image_paths if not fnmatch.fnmatch(os.path.basename(img), "combined_*_plot.png")]
    # Determine resize size
    if len(image_paths) == 1:
        resized_size = (int(img_size[0] * 0.6), int(img_size[1] * 0.6))  # shrink solo image
    else:
        resized_size = img_size
    if not image_paths:
        print(f"No plots found for '{search_term}'.")
        return None
    # Load and resize images
    images = [Image.open(img).resize(resized_size, Image.LANCZOS) for img in image_paths]
    
    # images = [Image.open(img).resize(img_size, Image.LANCZOS) for img in image_paths]
    # Determine grid layout
    cols = min(max_cols, len(images))
    rows = len(images) // cols
    rows = (
        len(images) + cols - 1
    ) // cols  # Round up to fit all images by adding cols - 1
    # Create a blank canvas
    combined_width = cols * resized_size[0]
    combined_height = rows * resized_size[1]
    combined_image = Image.new(
        "RGB", (combined_width, combined_height), (255, 255, 255)
    )  # 255,255,255 specifies background colour = white
    # Paste images into grid
    for idx, img in enumerate(images):
        x_offset = (idx % cols) * img_size[0]
        y_offset = (idx // cols) * img_size[1]
        combined_image.paste(img, (x_offset, y_offset))
    # Generate output filename if not provided
    if output_filename is None:
        output_filename = os.path.join(
            output_dir, f"plots/combined_{search_term.replace('.png', '')}.png"
        )
        ## windows sometimes rejects long paths. Workaround:
        if os.name == "nt":
            output_filename = "\\\\?\\" + os.path.abspath(output_filename)
    # Save the final image
    combined_image.save(output_filename)
    print(f"Combined plot saved to: {output_filename}")


### for combining data from different treatments for display in the report ###
def combine_csv_files(
    filename, output_dir, output_filename=None, top_n=10, new_column="treatment_pair",
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