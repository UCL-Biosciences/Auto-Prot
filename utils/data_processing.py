### Read in and pre-process data ###

import json
import logging
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn
import sklearn.ensemble
import sklearn.preprocessing
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer  # noqa: F401

# specify location of errors to standard output
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s - %(levelname)s - %(message)s"
)

# to write to log file, add e.g. filename='data_processing.log'. filemode = 'a' would append instead of overwriting


def load_data(file_path):
    """
    Load data from a CSV or TSV file into a pandas DataFrame.

    Parameters:
    - file_path (str): Path to the file.

    Returns:
    - pd.DataFrame: Loaded data.
    """
    try:
        if file_path.endswith(".csv"):
            df = pd.read_csv(file_path)
        elif file_path.endswith(".tsv"):
            df = pd.read_csv(file_path, sep="\t")
        else:
            raise ValueError("Unsupported file format. Please use CSV or TSV.")
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        return None


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


def clean_data(df, file_path=None, metadata=None, outPath=None, config=None):
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
    if "proteindata" in file_path:
        if metadata is None:
            raise ValueError("Error: Metadata is required but not provided.")
        # Check if the required column is in metadata
        elif "protein_abundance_name" not in metadata.columns:
            raise ValueError(
                "Error: 'protein_abundance_name' column is missing in the metadata."
            )
        else:
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
            df = df.replace(0, np.nan)
            df = df[df.isnull().mean(axis=1) <= 0.2]
            #### impute and normalise from AlphaPepStats functions #####
            ### log2 all vals
            df_log2 = np.log2(df + 1)
            # Note It has been shown that normalizing the data first and then imputing the data performs better, than the other way around (Karpievitch et al. 2012).
            # This preprocessing order is also acquired in AlphaPepStats (unless preprocessing is done in several steps).
            #### normalise ####
            # Subtract the median of each sample (column) from each value
            df_median_norm = df_log2.sub(df_log2.median(axis=0), axis=1)
            # Transpose: sklearn expects features in columns, samples in rows
            df_median_t = df_median_norm.T  # shape: samples × proteins
            df_median_t.columns = df_median_t.columns.astype(str)
            #### impute ####
            # Setup imputer using random forest approach (see docs for refs)
            # df_test = df_median_t.iloc[:, 0:2000]
            estimator = sklearn.ensemble.HistGradientBoostingRegressor(
                max_iter=30, max_depth=4, min_samples_leaf=5, random_state=0
            )
            imputer = sklearn.impute.IterativeImputer(
                max_iter=3,
                tol=0.01,
                random_state=0,
                estimator=estimator,
                verbose=2,
                n_nearest_features=30,
            )
            # Time the imputation
            start = time.time()
            imputation_array = imputer.fit_transform(df_median_t)
            end = time.time()
            print(f"Imputation took {end - start:.2f} seconds")
            # imputation_array = imp.fit_transform(df_log2_T)
            df_imp = pd.DataFrame(
                imputation_array.transpose(), index=df.index, columns=df.columns
            )
            ### view distributions after different steps #####
            #### By Sample ######
            # List of dataframes and titles for the subplots
            dfs = [df, df_log2, df_median_norm, df_imp]
            titles = [
                "Raw Intensities",
                "Log₂ Transformed",
                "Sample-Median Normalised",
                "Normalised then Imputed",
            ]
            # Set up the figure with a 2x2 grid
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            axes = axes.flatten()
            for ax, data, title in zip(axes, dfs, titles):
                # Convert the DataFrame to long format: one row per measurement with sample id and intensity.
                long_df = data.melt(var_name="sample_rep", value_name="intensity")
                # Merge with metadata to obtain treatment information for each sample.
                # Make sure metadata has 'sample_id' and 'treatment' columns.
                long_df = long_df.merge(
                    metadata[["sample_rep", "treatment"]], on="sample_rep", how="left"
                )
                # Plot boxplots: each sample's distribution is shown on the x-axis.
                # The boxes are colored by treatment.
                sns.boxplot(
                    x="sample_rep", y="intensity", data=long_df, hue="treatment", ax=ax
                )
                ax.set_title(title)
                ax.set_xlabel("Sample ID")
                ax.set_ylabel("Intensity")
                # Rotate x tick labels for better readability
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
                # Remove redundant legends in subplots (optional)
                if ax != axes[0]:
                    ax.get_legend().remove()
            # Add a single legend for the entire figure
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, title="Treatment", loc="upper right")
            plt.tight_layout()
            plot_path = os.path.join(
                outPath, "plots", "boxplots_preProcessing_all_samples_plot.png"
            )
            plt.savefig(plot_path, dpi=300)
            plt.close()
            #### for treatments ####
            # Set up the figure with a 2x2 grid for KDE plots
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            axes = axes.flatten()
            for ax, data, title in zip(axes, dfs, titles):
                # Convert the DataFrame to long format: one row per measurement with sample id and intensity.
                long_df = data.melt(var_name="sample_rep", value_name="intensity")
                # Merge with metadata to obtain treatment information for each sample.
                # Make sure metadata has 'sample_rep' and 'treatment' columns.
                long_df = long_df.merge(
                    metadata[["sample_rep", "treatment"]], on="sample_rep", how="left"
                )
                # Plot KDE for each treatment on the same subplot
                for treatment, group in long_df.groupby("treatment"):
                    sns.kdeplot(data=group, x="intensity", ax=ax, label=treatment)
                ax.set_title(title)
                ax.set_xlabel("Intensity")
                ax.set_ylabel("Density")
                ax.legend(title="Treatment")
            plt.tight_layout()
            plot_path = os.path.join(
                outPath, "plots", "KDE_preProcessing_all_treatments_plot.png"
            )
            plt.savefig(plot_path, dpi=300)
            plt.close()
            ## which df to use?
            df_to_use = config.get("df_to_use")
            df = locals()[df_to_use]
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
            #### We also want to remove low abundance proteins. Explained in the docs
    df = df.drop_duplicates()
    if "proteindata" in file_path:
        return df, nrow_original
    else:
        return df


###########################################################################
#### This is the main function for processing data, combining the above ###
###########################################################################


def preprocess_data(file_path, metadata=None, json_out=None, outPath=None, config=None):
    """
    Full preprocessing pipeline: load, clean, and normalize column names.

    Parameters:
    - file_path (str): Path to the data file.

    Returns:
    - pd.DataFrame: Preprocessed dataframe.
    """
    df = load_data(file_path)
    if df is not None:
        df = normalise_column_names(df, file_path=file_path, metadata=metadata)
    #     # Clean data with or without metadata
    #     if metadata is None:
    #         df = clean_data(df, file_path = file_path)
    #     else:
    #         df, row_nan_count = clean_data(df, file_path = file_path, metadata=metadata)

    if "metadata" in file_path:
        df = clean_data(df, file_path=file_path, config=config)
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

        ### run function to validate metadata
        validate_metadata(df)

        return df  # if just metadata

    if "proteindata" in file_path:
        df, nrow_original = clean_data(
            df, file_path=file_path, metadata=metadata, outPath=outPath, config=config
        )

        #### protein summary
        NUM_PROTS = len(df.index)
        NUM_PROTS_REMOVED = nrow_original - NUM_PROTS

        # Compute mean value per row (per protein across all samples)
        mean_abundance = df.iloc[:, 1:].mean(axis=1)  # Exclude protein column

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

        ### run function to validate protein abundance data
        validate_proteindata(data=df, metadata=metadata)

        return df


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
    if not data.index.is_unique:
        raise ValueError("Error: Protein identifiers (row indices) must be unique.")
    # Check if column names (samples) are unique
    if not data.columns.is_unique:
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
