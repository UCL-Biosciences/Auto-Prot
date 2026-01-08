### Data pre-processing
## log, normalise and impute protein abundance data

import glob
import os
import subprocess
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.ensemble
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer  # noqa: F401


## impute function
def impute_prot_data_histgradboost(df_filtered, df_norm_t):
    """
    Impute missing values in a protein abundance matrix using a tree-based model.

    Uses a gradient-boosted regressor with sklearn's IterativeImputer on a transposed
    median-normalised DataFrame. Warns if imputation shifts sample means by >0.2 SD.

    Args:
        df_filtered (pd.DataFrame): Filtered (pre-imputation) DataFrame, used for index/column names.
        df_norm_t (pd.DataFrame): Transposed, normalised DataFrame for imputation.

    Returns:
        pd.DataFrame: Imputed protein abundance DataFrame with original shape and labels.
    """
    # Setup imputer using random forest approach (see docs for refs)
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
    imputation_array = imputer.fit_transform(df_norm_t)
    end = time.time()
    print(f"Imputation took {end - start:.2f} seconds")
    # convert back to df
    df_imp = pd.DataFrame(
        imputation_array.transpose(),
        index=df_filtered.index,
        columns=df_filtered.columns,
    )
    # After creating df_imp, check whether any sample has undergone large change in mean value
    # Means per sample BEFORE (in transposed matrix → axis=1 = sample axis)
    mean_before = df_norm_t.mean(axis=1)  # rows = samples
    # Means per sample AFTER (transpose df_imp to match)
    mean_after = df_imp.T.mean(
        axis=1
    )  # rows = samples. note df_imp needs transposing to match df_norm_t
    # SDs per sample BEFORE
    sd_before = df_norm_t.std(axis=1)
    # Compute per-sample shift in SDs
    shift_in_sds = (mean_before - mean_after).abs() / sd_before
    # Warn if any column changed too much
    for col, shift in shift_in_sds.items():
        if shift > 0.2:
            print(
                f"⚠️ Warning: Imputation shifted sample '{col}' mean by {shift:.2f} SDs"
            )
    return df_imp


def normalise_vsn(file_path_in, file_path_normalised_out, meanSdPlot_path):
    """
    Normalises proteomics data using Variance Stabilising Normalisation (VSN)
    via an R script, and generates a mean–SD diagnostic plot.

    Parameters
    ----------
    file_path_in : str
        Path to the input CSV/TSV file containing raw protein abundance data.

    file_path_normalised_out : str
        Path where the normalised output file should be saved.

    meanSdPlot_path : str
        Path where the mean–SD diagnostic plot (e.g. PNG or PDF) will be saved.

    Raises
    ------
    subprocess.CalledProcessError
        If the R script fails or returns a non-zero exit code.

    Notes
    -----
    - This function runs an R script (`normalise-vsn.R`) that must accept
      three arguments (input file, output file, and plot path).
    - The script is executed within a Conda environment called `r-limma-env`.
    - `check=True` ensures Python raises an error if the R script fails.
    """
    subprocess.run(
        [
            "conda",
            "run",
            "-n",
            "r-limma-env",
            "Rscript",
            "autoprot/r_scripts/normalise-vsn.R",
            file_path_in,
            file_path_normalised_out,
            meanSdPlot_path,
        ],
        check=True,
    )


def filter_proteins_by_group_missingness(
    df, metadata, sample_col="sample_rep", group_col="treatment", threshold=None, config=None
):
    """
    Filter proteins based on missingness within each group.

    Parameters:
        df (pd.DataFrame): Protein abundance DataFrame (rows = proteins, columns = samples).
        metadata (pd.DataFrame): Metadata with sample and group columns.
        sample_col (str): Column name in metadata with sample names matching df columns.
        group_col (str): Column name in metadata to group by (e.g., treatment).
        config (dict): Configuration dict, must include 'missing_threshold'.


    Returns:
        pd.DataFrame: Filtered protein DataFrame.
    """
    threshold = config.get("missing_threshold")
    # create empty list to store valid protein sets for each group
    valid_sets = []
    # .groupby groups the metadata by the specified group column
    # outputs group name and DataFrame for each group
    for group, group_df in metadata.groupby(group_col):
        # get the sample names for this group
        samples = group_df[sample_col]
        # filter the original DataFrame to only include these samples
        sub_df = df[samples]
        # calculate the proportion of non-missing values for each protein
        # and keep those that are present in at least `threshold` proportion of samples
        valid = sub_df.notna().mean(axis=1) >= threshold
        # add the indices of valid proteins to the list
        # this will be a set of protein indices that are valid for this group
        valid_sets.append(set(df.index[valid]))
    # find intersection of all valid sets to get proteins present > threshold in all groups
    keep_proteins = set.intersection(*valid_sets)
    print(
        "found ",
        len(keep_proteins),
        " proteins in ",
        (100 * threshold),
        "% of each treatment group",
    )
    return df.loc[list(keep_proteins)]


def impute_pimms_cf(
    df: pd.DataFrame,
    n_factors: int = 30,
    batch_size: int = 4096,
    epochs_max: int = 20,
    cuda: bool = False,
    target_column: str = "intensity",
) -> pd.DataFrame:
    """
    Impute missing values in a proteomics intensity matrix using PIMMS's collaborative-filtering.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-form DataFrame of shape (n_samples, n_proteins), where rows are sample IDs,
        columns are protein IDs, and values are log2-transformed intensities (floats). Missing values should
        be represented as NaN.
    n_factors : int, optional
        Number of latent factors to learn for both samples and proteins (default=30).
    batch_size : int, optional
        Batch size for training the collaborative-filtering model (default=4096).
    epochs_max : int, optional
        Maximum number of training epochs (default=20).
    cuda : bool, optional
        Whether to use GPU acceleration (if available). Defaults to False.
    target_column : str, optional
        Name to assign to the intensity column when stacking into a long-form Series.
        Defaults to "intensity".

    Returns
    -------
    pd.DataFrame
        DataFrame of the same shape and indexing as `df`, with all NaNs imputed.
        Imputed values are the model’s predictions, combined with the original observed values.

    Example
    -------
    >>> df_imputed = impute_pimms_cf(raw_df, n_factors=50, epochs_max=30, cuda=True)
    """
    # have to import here because pimms not available on mac OS
    from pimmslearn.sklearn.cf_transformer import CollaborativeFilteringTransformer
    # 1. Stack to long form
    df.index.name = "sample_id"
    df.columns.name = "protein_id"
    series = df.stack()
    series.name = target_column
    # series.index.names = ["sample_id", "protein_id"]
    # 2. Initialize transformer
    # controlling batch size
    cf = CollaborativeFilteringTransformer(
        target_column=target_column,
        sample_column="sample_id",
        item_column="protein_id",
        n_factors=n_factors,
        batch_size=(int(len(series) / 10)),
    )
    # 3. Fit and transform
    start_time = time.time()
    cf.fit(series, cuda=cuda, epochs_max=epochs_max)
    imputed_long = cf.transform(series)
    elapsed = time.time() - start_time
    print(f"[PIMMS CF] Imputation completed in {elapsed:.1f} seconds.")
    plt.close("all")
    ## clean up files produced by pimms cf
    for file in glob.glob("collab_training*"):
        os.remove(file)
    for file in glob.glob("model_params*"):
        os.remove(file)
    # 4. Unstack back to wide form
    df_imputed = pd.DataFrame(imputed_long.unstack(level="protein_id").transpose())
    return df_imputed


def process_prot_data(df, config, outPath, metadata):
    """
    Preprocess protein abundance data by filtering, transforming, normalising, and imputing.

    Applies:
    - Filtering of high-missingness proteins
    - Log2 transform
    - Sample-wise median normalisation
    - Imputation using a tree-based model

    Args:
        df (pd.DataFrame): Raw protein abundance data (proteins in rows, samples in columns).

    Returns:
        dict: Dictionary of intermediate DataFrames:
            - 'df': Filtered
            - 'df_log2': Log2-transformed
            - 'df_norm_norm': Normalised
            - 'df_imp': Final imputed data
    """
    df = df.replace(0, np.nan)
    #### Duplicates ####
    # sometimes there are identical values for multiple rows (proteins or PTMs)
    # We can't say for all cases what the cause is
    # most likely the same peptide or PTM in a different context i.e. peptide has been cleaved at different sites.
    # which means the same protein or PTMP is represented by multiple peptides
    # generally we think better to keep only one of them
    # count rows before dropping
    n_before = len(df)
    #
    # drop duplicates
    df = df.drop_duplicates()
    #
    # count rows after
    n_after = len(df)
    #
    print("Rows with duplicated values removed:", n_before - n_after)
    # filter for proteins found in XX% per treatment group
    threshold = config["missing_threshold"]
    df_filtered = filter_proteins_by_group_missingness(
        df, metadata, threshold=threshold, config = config
    )
    # df = df[df.isnull().mean(axis=1) <= 0.2]
    ### log2 all vals
    df_log2 = np.log2(df_filtered)
    #### normalise and impute #####
    # Note It has been shown that normalizing the data first and then imputing the data performs better than the other way around (Karpievitch et al. 2012).
    # see e.g. AlphaPepStats
    #### normalise ####
    ## currently two options supported: vsn (recommended) and sample-median
    ## vsn requires raw positive intensities, sample median works on log2-transformed data
    normalise_method = config["normalise_method"]
    if normalise_method == "vsn":
        prot_path = os.path.join(outPath, "data/prots_no_zero_values.csv").replace(
            "\\", "/"
        )
        normalised_path = os.path.join(
            outPath, "data/prots_vsn_normalised.csv"
        ).replace("\\", "/")
        meanSdPlot_path = os.path.join(outPath, "plots/vsn_meanSDplot.png").replace(
            "\\", "/"
        )
        df_filtered.to_csv(prot_path, index=True)
        normalise_vsn(
            file_path_in=prot_path,
            file_path_normalised_out=normalised_path,
            meanSdPlot_path=meanSdPlot_path,
        )
        df_norm = pd.read_csv(normalised_path, index_col=0)
    ### for normalise by sample median, normalise log2 transformed data
    if normalise_method == "sample-median":
        # Subtract the median of each sample (column) from each value
        df_norm = df_log2.sub(df_log2.median(axis=0), axis=1)
    # Transpose: sklearn expects features in columns, samples in rows
    df_norm_t = df_norm.T  # shape: samples × proteins
    df_norm_t.columns = df_norm_t.columns.astype(str)
    #### impute ####
    if config["imputation_method"] == "hist_grad_boost":
        
        print("imputing with histogram gradient booster")
        df_imp = impute_prot_data_histgradboost(df_filtered, df_norm_t)

    elif config["imputation_method"] == "pimms_collabfilter":
        
        print("imputing with pimms: collaborative filtering")
        df_imp = impute_pimms_cf(df=df_log2.T)
        
    return {
        "df": df,
        "df_log2": df_log2,
        "df_norm": df_norm,
        "df_imp": df_imp,
    }


def view_prot_distributions(dfs_values, plot_titles, metadata, outPath):
    """
    Visualise protein abundance distributions before and after processing.

    Creates a subplot for each df (clean raw, log2 transofrmed, sample-median normalisted, imputed):
    - Per-sample boxplots (grouped by treatment)
    - Kernel density estimates (KDEs) grouped by treatment

    Args:
        dfs_values (Iterable[pd.DataFrame]): A sequence of processed DataFrames (e.g. raw, log2, normalised, imputed).
        plot_titles (list[str]): Titles for each corresponding DataFrame.
        metadata (pd.DataFrame): Metadata containing 'sample_rep' and 'treatment'.
        outPath (str): Output directory for saving plots.

    Side effects:
        Saves boxplot and KDE plots as PNG files in the output directory.
    """
    ### view distributions after different steps #####
    ### for each sample as boxplots
    # Set up the figure with a 2x2 grid
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()
    for ax, data, title in zip(axes, dfs_values, plot_titles):
        # Convert the DataFrame to long format: one row per measurement with sample id and intensity.
        long_df = data.melt(var_name="sample_rep", value_name="intensity")
        # Merge with metadata to obtain treatment information for each sample.
        # Make sure metadata has 'sample_id' and 'treatment' columns.
        long_df = long_df.merge(
            metadata[["sample_rep", "treatment"]], on="sample_rep", how="left"
        )
        # Sort samples by treatment, then sample_id, then sample_rep
        sample_order = (
            long_df.sort_values(["treatment", "sample_rep"])["sample_rep"]
            .drop_duplicates()
            .tolist()
        )
        # Set plot_label as a categorical with the desired order
        long_df["sample_rep"] = pd.Categorical(
            long_df["sample_rep"], categories=sample_order, ordered=True
        )
        # Plot boxplots: each sample's distribution is shown on the x-axis.
        # The boxes are colored by treatment.
        sns.boxplot(x="sample_rep", y="intensity", data=long_df, hue="treatment", ax=ax)
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
    if os.name == "nt":
        plot_path = "\\\\?\\" + os.path.abspath(plot_path)
    plt.savefig(plot_path, dpi=300)
    plt.close()
    #### for treatments as kernel densities ####
    # Set up the figure with a 2x2 grid for KDE plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()
    for ax, data, title in zip(axes, dfs_values, plot_titles):
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
    if os.name == "nt":
        plot_path = "\\\\?\\" + os.path.abspath(plot_path)
    plt.savefig(plot_path, dpi=300)
    plt.close()
