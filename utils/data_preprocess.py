### Data pre-processing
## log, normalise and impute protein abundance data

import os 
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer  # noqa: F401
import sklearn.ensemble

## impute function
def impute_prot_data(df, df_median_t):
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
    # imputation_array = imp.fit_transform(df_log2_T)
    df_imp = pd.DataFrame(
        imputation_array.transpose(), index=df.index, columns=df.columns
    )
    return df_imp

def process_prot_data(df, metadata, config):
    df = df.replace(0, np.nan)
    df = df[df.isnull().mean(axis=1) <= 0.2]
    #### impute and normalise #####
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
    df_imp = impute_prot_data(df, df_median_t)
    return {
        "df": df,
        "df_log2": df_log2,
        "df_median_norm": df_median_norm,
        "df_imp": df_imp,
    }

def view_prot_distributions(dfs, plot_titles, metadata, outPath):
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
    plt.savefig(plot_path, dpi=300)
    plt.close()