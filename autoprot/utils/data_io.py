### Data io
## functions for reading / writing data

import os

import pandas as pd


def load_data(file_path):
    """
    Load a CSV or TSV file into a pandas DataFrame and strip whitespace from strings.

    This function:
    - Supports `.csv` and `.tsv` file extensions.
    - Strips leading/trailing whitespace from column names.
    - Strips whitespace from string values in the DataFrame.

    Args:
        file_path (str): Path to the input file. Must end with `.csv` or `.tsv`.

    Returns:
        pd.DataFrame: Loaded data with whitespace stripped from column names and cells.

    Raises:
        ValueError: If the file extension is not `.csv` or `.tsv`.
    """
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path)
    elif file_path.endswith(".tsv"):
        df = pd.read_csv(file_path, sep="\t")
    else:
        raise ValueError("Unsupported file format. Please use CSV or TSV.")

    # Drop rows where all values are NA or empty string
    # We replace empty strings with NA first, then drop rows that are all NA
    df = df.replace(r'^\s*$', pd.NA, regex=True)  # treat empty strings as NA
    df = df.dropna(how="all")

    # Strip whitespace from column names
    df.columns = df.columns.str.strip()

    # Strip whitespace from string entries (cell values)
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)

    if "proteindata" in file_path:
        # Set the first column as the index
        df.set_index(df.columns[0], inplace=True)

    return df


### make outdir
def make_outdir(out_path, make_subdirs=True):
    """
    Create an output directory if it doesn't exist, with optional subdirectories.

    Args:
        out_path (str): Path to the main output directory.
        make_subdirs (bool): If True, also creates 'data' and 'plots' subdirectories.

    Side effects:
        Creates directories on the file system. Prints status messages.
    """
    ### some previous outputs are brought into the report unintentionally
    ### remove output dir if it exists
    if os.path.exists(out_path):
        
        import shutil
        shutil.rmtree(out_path) # remove directory and all its contents

    if not os.path.exists(out_path):

        try:
            os.makedirs(out_path)
            print(f"Directory '{out_path}' created successfully.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{out_path}'.")
        except Exception as e:
            print(f"An error occurred: {e}")

    if make_subdirs:
        out_subDirs = ["data", "plots"]
        for subDir in out_subDirs:
            path = os.path.join(out_path, subDir)
            if not os.path.exists(path):
                os.mkdir(path)
