### Data pre-processing
## log, normalise and impute protein abundance data

import os
import time
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.ensemble
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer  # noqa: F401



## impute function
def impute_prot_data(df, df_median_t):
    """
    Impute missing values in a protein abundance matrix using a tree-based model.

    Uses a gradient-boosted regressor with sklearn's IterativeImputer on a transposed
    median-normalised DataFrame. Warns if imputation shifts sample means by >0.2 SD.

    Args:
        df (pd.DataFrame): Original (pre-imputation) DataFrame, used for index/column names.
        df_median_t (pd.DataFrame): Transposed, normalised DataFrame for imputation.

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
    imputation_array = imputer.fit_transform(df_median_t)
    end = time.time()
    print(f"Imputation took {end - start:.2f} seconds")
    # convert back to df
    df_imp = pd.DataFrame(
        imputation_array.transpose(), index=df.index, columns=df.columns
    )
    # After creating df_imp, check whether any sample has undergone large change in mean value
    # Means per sample BEFORE (in transposed matrix → axis=1 = sample axis)
    mean_before = df_median_t.mean(axis=1)  # rows = samples
    # Means per sample AFTER (transpose df_imp to match)
    mean_after = df_imp.T.mean(
        axis=1
    )  # rows = samples. note df_imp needs transposing to match df_median_t
    # SDs per sample BEFORE
    sd_before = df_median_t.std(axis=1)
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
    subprocess.run(
    [
        "conda",
        "run",
        "-n",
        "r-limma-env",
        "Rscript",
        "src/r_scripts/normalise-vsn.R",
        file_path_in,
        file_path_normalised_out,
        meanSdPlot_path,
    ],
    check=True,
)
    
def process_prot_data(df, config, outPath):
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
            - 'df_median_norm': Normalised
            - 'df_imp': Final imputed data
    """
    df = df.replace(0, np.nan)
    df = df[df.isnull().mean(axis=1) <= 0.2]
    ### log2 all vals
    df_log2 = np.log2(df)
    #### normalise and impute #####
    # Note It has been shown that normalizing the data first and then imputing the data performs better than the other way around (Karpievitch et al. 2012).
    # see e.g. AlphaPepStats
    #### normalise ####
    ## currently two options supported: vsn (recommended) and sample-median
    ## vsn requires raw positive intensities, sample median works on log2-transformed data
    normalise_method = config.get("normalise_method")
    if normalise_method == "vsn":
        prot_path = os.path.join(outPath, "data/prots_no_zero_values.csv").replace("\\", "/")
        normalised_path = os.path.join(outPath, "data/prots_vsn_normalised.csv").replace("\\", "/")
        meanSdPlot_path = os.path.join(outPath, "plots/vsn_meanSDplot.png").replace("\\", "/")
        df.to_csv(prot_path, index=True)
        normalise_vsn(file_path_in = prot_path,
                      file_path_normalised_out = normalised_path,
                      meanSdPlot_path = meanSdPlot_path)
        df_median_norm = pd.read_csv(normalised_path, index_col=0)
    ### for normalise by sample median, normalise log2 transformed data
    if normalise_method == "sample-median":
        # Subtract the median of each sample (column) from each value
        df_median_norm = df_log2.sub(df_log2.median(axis=0), axis=1)
    # Transpose: sklearn expects features in columns, samples in rows
    df_median_t = df_median_norm.T  # shape: samples × proteins
    df_median_t.columns = df_median_t.columns.astype(str)
    #### impute ####
    df_imp = impute_prot_data(df, df_median_t)
    return {
        "df": df,
        "df_log2": df_log2,
        "df_median_norm": df_median_norm,
        "df_imp": df_imp,
    }


def view_prot_distributions(dfs, plot_titles, metadata, outPath):
    """
    Visualise protein abundance distributions before and after processing.

    Creates a subplot for each df (clean raw, log2 transofrmed, sample-median normalisted, imputed):
    - Per-sample boxplots (grouped by treatment)
    - Kernel density estimates (KDEs) grouped by treatment

    Args:
        dfs (Iterable[pd.DataFrame]): A sequence of processed DataFrames (e.g. raw, log2, normalised, imputed).
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
    for ax, data, title in zip(axes, dfs, plot_titles):
        # Convert the DataFrame to long format: one row per measurement with sample id and intensity.
        long_df = data.melt(var_name="sample_rep", value_name="intensity")
        # Merge with metadata to obtain treatment information for each sample.
        # Make sure metadata has 'sample_id' and 'treatment' columns.
        long_df = long_df.merge(
            metadata[["sample_rep", "treatment"]], on="sample_rep", how="left"
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
    for ax, data, title in zip(axes, dfs, plot_titles):
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
