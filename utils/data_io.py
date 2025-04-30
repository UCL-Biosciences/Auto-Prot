### Data io
## functions for reading / writing data

import os
import pandas as pd

def load_data(file_path):
    """
    Load a CSV or TSV into a DataFrame.
    Raises ValueError if format unsupported.
    """
    if file_path.endswith(".csv"):
        return pd.read_csv(file_path)
    elif file_path.endswith(".tsv"):
        return pd.read_csv(file_path, sep="\t")
    else:
        raise ValueError("Unsupported file format. Please use CSV or TSV.")

### make outdir
def make_outdir(out_path, make_subdirs=True):
    """
    Make output dir, including checks

    Parameters:
    - out_path (str): Path to where project output should be stored.

    Returns:

    """
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
