#### analysis of all_samples
## e.g. pca, mds, heatmap
## as opposed to pairwise treatment analysis
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.patches import Ellipse
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.manifold import MDS


### function for drawing centroids ellipses
def add_group_ellipses_and_centroids(
    ax, df, x_col, y_col, group_col, palette=None, scale=4
):
    """
    Add ellipses and centroids to a scatterplot grouped by a categorical column.

    Args:
        ax (matplotlib.axes.Axes): The axes to draw on.
        df (pd.DataFrame): DataFrame containing data to group and plot.
        x_col (str): Column name for X-axis values.
        y_col (str): Column name for Y-axis values.
        group_col (str): Column name used for grouping (e.g. treatment).
        palette (dict, optional): Dict mapping group labels to colours.
        scale (float): Multiplier on std dev for ellipse size (default 4 = ~95% CI).
    """
    for name, group in df.groupby(group_col):
        ellipse_color = palette[name] if palette and name in palette else None
        # Centroid
        centroid_x = group[x_col].mean()
        centroid_y = group[y_col].mean()
        ax.plot(
            centroid_x,
            centroid_y,
            marker="o",
            markersize=8,
            markerfacecolor="white",
            markeredgewidth=1.5,
            markeredgecolor=ellipse_color,
            zorder=5,
        )
        # Ellipse
        cov = np.cov(group[[x_col, y_col]].values.T)
        lambda_, v = np.linalg.eig(cov)
        lambda_ = np.sqrt(lambda_)
        angle = np.rad2deg(np.arccos(v[0, 0]))
        ell = Ellipse(
            xy=(centroid_x, centroid_y),
            width=lambda_[0] * scale,
            height=lambda_[1] * scale,
            angle=angle,
            edgecolor=ellipse_color,
            facecolor="none",
            linewidth=1,
        )
        ell.set_zorder(3)  # Ensure it's drawn above grid
        ax.add_patch(ell)


def generate_pca(
    df: pd.DataFrame,
    output_dir: str,
    metadata: pd.DataFrame = None,
    variance_threshold: float = 0.95,
    plot_components: tuple[int, int] = (1, 2),
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Performs PCA on protein abundance data, plots selected components, and saves results.

    The number of principal components is automatically chosen to capture at least
    `variance_threshold` proportion of total variance. A 2D plot of the specified components
    (default PC1 vs PC2) is saved.

    Args:
        df (pd.DataFrame): Protein abundance data (samples as rows, proteins as columns).
        output_dir (str): Directory to save PCA outputs.
        metadata (pd.DataFrame, optional): Metadata containing 'sample_rep' and 'treatment' columns.
        variance_threshold (float): Minimum cumulative variance to retain in PCA (default is 0.95).
        plot_components (Tuple[int, int]): Principal components to plot (1-based index, e.g., (1, 2)).

    Returns:
        Tuple[pd.DataFrame, np.ndarray]:
            - DataFrame with PCA coordinates and treatment.
            - 1D array of explained variance ratios.
    """
    # Drop duplicates to get one colour per treatment
    treatment_palette = (
        metadata[["treatment", "colours"]]
        .drop_duplicates("treatment")
        .set_index("treatment")["colours"]
        .to_dict()
    )
    df = df.dropna(axis="columns")  # note PCA can't handle missing
    n_prot = df.shape[1]

    ### PC calculated twice - first to identify optimal number of components
    # i.e. the miminum number of components required to explain a given proportion of variance
    # First PCA to get full spectrum of variance
    full_pca = PCA()
    full_pca.fit(df)
    explained_variance = full_pca.explained_variance_ratio_

    # Determine how many components to keep
    cumulative = np.cumsum(explained_variance)
    n_components = np.searchsorted(cumulative, variance_threshold) + 1

    # Perform PCA with optimal number of components
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(df)

    plot_title = "PCA Plot (n proteins = " + str(n_prot) + ")"
    # Create DataFrame with PCA results
    df_PCA = pd.DataFrame(
        data=principal_components,
        columns=[f"PC{i+1}" for i in range(principal_components.shape[1])],
        index=df.index,
    )

    # pull over treatment from metadata
    # df_PCA . look up by index . and map to
    # metadata . temporarily set index to sample_id [ and find treatment ]
    df_PCA["treatment"] = df_PCA.index.map(
        metadata.set_index("sample_rep")["treatment"]
    )

    ## check mapping has worked
    if df_PCA["treatment"].isnull().any():
        missing = df_PCA[df_PCA["treatment"].isnull()].index.tolist()
        raise ValueError(f"PCA df is missing treatment labels for samples: {missing}")

    # Calculate variance explained by each component
    explained_variance = pca.explained_variance_ratio_

    # Plot selected components
    pc_x, pc_y = plot_components
    pc_labels = [
        f"PC{pc_x} ({explained_variance[pc_x-1]*100:.2f}%)",
        f"PC{pc_y} ({explained_variance[pc_y-1]*100:.2f}%)",
    ]

    # Plot PCA
    plt.figure(figsize=(6, 4.5))
    plot_pca = sns.scatterplot(
        data=df_PCA,
        x=f"PC{pc_x}",
        y=f"PC{pc_y}",
        hue="treatment",
        palette=treatment_palette,
    )
    ax = plt.gca()
    add_group_ellipses_and_centroids(
        ax,
        df_PCA,  # or mds_coords_df
        x_col=f"PC{pc_x}",
        y_col=f"PC{pc_y}",
        group_col="treatment",
        palette=treatment_palette,
    )
    # Annotate points
    # Add non-overlapping text labels
    texts = []
    for i in range(df_PCA.shape[0]):
        texts.append(
            plt.text(
                x=df_PCA[f"PC{pc_x}"].iloc[i],
                y=df_PCA[f"PC{pc_y}"].iloc[i],
                s=df_PCA.index[i],
                fontsize=8,
                ha="center",
                va="bottom",
            )
        )

    # Smart adjustment to avoid overlaps
    adjust_text(
        texts,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
        expand_points=(1.1, 1.1),
        expand_text=(1.2, 1.2),
        force_text=(0.5, 0.5),
    )

    plt.xlabel(pc_labels[0])
    plt.ylabel(pc_labels[1])

    # Set labels and title
    plot_pca.set(xlabel=pc_labels[0], ylabel=pc_labels[1])
    plt.title(plot_title)
    plt.tight_layout()
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, "plots", "pca_plot.png")
    data_path = os.path.join(output_dir, "data", "pca_data.csv")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    df_PCA.to_csv(data_path, index=False)
    print(f"PCA plot saved to {plot_path}")
    print(f"PCA data saved to {data_path}")
    return df_PCA


def generate_MDS(
    df: pd.DataFrame, output_dir: str, metadata: pd.DataFrame = None
) -> pd.DataFrame:
    """
    Performs metric multidimensional scaling (MDS) using Euclidean distance on protein abundance data,
    and saves a scatter plot of the resulting coordinates.

    This method visualises pairwise sample dissimilarities in a reduced 2D space. Sample points are
    optionally coloured by treatment group using metadata.

    Args:
        df (pd.DataFrame): Protein abundance data (samples as rows, proteins as columns).
        output_dir (str): Directory to save the MDS plot and coordinate data.
        metadata (pd.DataFrame, optional): Sample metadata containing 'sample_rep' and 'treatment' columns.

    Returns:
        pd.DataFrame: DataFrame containing MDS coordinates and associated metadata.
    """
    # Drop duplicates to get one colour per treatment
    treatment_palette = (
        metadata[["treatment", "colours"]]
        .drop_duplicates("treatment")
        .set_index("treatment")["colours"]
        .to_dict()
    )
    # Perform MDS
    df = df.dropna(axis="columns")  # note MDS can't handle missing
    # pdist(x) computes the Euclidean distances between each pair of points in an array
    dissimilarity_array = pdist(df, metric="euclidean")
    n_prot = df.shape[1]
    plot_title = (
        "MDS Dissimilarity based on Protein Abundance (n proteins = "
        + str(n_prot)
        + ")"
    )
    # suqareform() returns the matrix
    dissimilarity_matrix = squareform(dissimilarity_array)
    # Perform NMDS
    mds = MDS(
        n_components=2,  # n_components is number of dimensions
        dissimilarity="precomputed",  # we calculate above
        random_state=42,
        metric=True,
    )
    # Fit the data  and return the embedded coordinates.
    mds_coords = mds.fit_transform(
        dissimilarity_matrix
    )  # If dissimilarity=='precomputed', the input should be the dissimilarity matrix.
    # Convert the NMDS coordinates (nmds_coords) to a pandas DataFrame
    mds_coords_df = pd.DataFrame(mds_coords, columns=["MDS1", "MDS2"], index=df.index)
    # pull over treatment from metadata
    # df_PCA . look up by index . and map to
    # metadata . temporarily set index to sample_id [ and find treatment ]
    mds_coords_df["treatment"] = mds_coords_df.index.map(
        metadata.set_index("sample_rep")["treatment"]
    )
    #### Generate NMDS plot
    plot_nmds = sns.scatterplot(
        data=mds_coords_df,
        x="MDS1",
        y="MDS2",
        hue="treatment",
        palette=treatment_palette,
    )
    add_group_ellipses_and_centroids(
        plt.gca(),
        mds_coords_df,
        x_col="MDS1",
        y_col="MDS2",
        group_col="treatment",
        palette=treatment_palette,
    )
    ## seaborn returns an axis-object rather than a figure, so you can still alter features of it. E.g. axes names:
    plot_nmds.set(xlabel="MDS 1", ylabel="MDS 2")
    # Add sample IDs from the 'target' column with an offset
    # Create list of text objects to adjust
    texts = []
    for i in range(mds_coords_df.shape[0]):
        texts.append(
            plt.text(
                x=mds_coords_df["MDS1"].iloc[i],
                y=mds_coords_df["MDS2"].iloc[i],
                s=mds_coords_df.index[i],
                fontsize=9,
                ha="center",
                va="bottom",
            )
        )

    # Adjust the labels to avoid overlap, with optional arrows
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
    # Set labels and title
    plt.title(plot_title)
    plt.tight_layout()
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, "plots", "mds_plot.png")
    data_path = os.path.join(output_dir, "data", "mds_data.csv")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    mds_coords_df.to_csv(data_path, index=False)
    print(f"MDS plot saved to {plot_path}")
    print(f"MDS data saved to {data_path}")
    return mds_coords_df


def generate_heatmap(
    df: pd.DataFrame, output_dir: str, metadata: pd.DataFrame = None
) -> pd.DataFrame:
    """
    Creates a clustered heatmap of protein abundance data using Euclidean distance.

    Proteins are clustered along rows and samples along columns. Dendrograms and clustering
    are based on the transposed abundance matrix.

    Args:
        df (pd.DataFrame): Protein abundance data (samples as rows, proteins as columns).
        output_dir (str): Directory to save the heatmap plot.
        metadata (pd.DataFrame): Metadata with sample_rep, treatment, and colours columns.

    Returns:
        pd.DataFrame: The (unchanged) abundance DataFrame, returned for chaining if needed.
    """
    # 1) Build the colour strip from metadata in the SAME order as df.columns
    # Map each sample_rep to its colour
    colour_map = metadata.set_index("sample_rep")["colours"]
    # Create a list/Series of colours aligned to df.columns
    col_colors = df.columns.map(colour_map)

    # 2) Generate the clustermap
    #    - df: your data matrix
    #    - cmap: the heat-map palette
    #    - col_colors: tells Seaborn to draw a full-width strip under the dendrogram
    #    - cbar_pos=None: suppress the tiny default legend in the corner
    
    df = df.dropna(axis="rows")  # note clustermap can't handle missing
    cg = sns.clustermap(
        df,
        fmt=".2f",
        cmap="viridis",
        col_colors=col_colors,
        xticklabels=True,
        yticklabels=False,
        # cbar_pos=None
    )

    # 3) Tidy up the x-axis labels to appear on top, matching the clustered order
    #    - reordered_ind gives you the new column ordering
    ordered_cols = df.columns[cg.dendrogram_col.reordered_ind]
    cg.ax_heatmap.xaxis.set_ticks_position("top")
    cg.ax_heatmap.xaxis.set_label_position("top")
    cg.ax_heatmap.set_xticks(np.arange(len(ordered_cols)))
    cg.ax_heatmap.set_xticklabels(ordered_cols, rotation=45, ha="left")

    # 4) Save the figure at high resolution and close it to free memory
    plot_path = os.path.join(output_dir, "plots", "heatmap_plot.png")
    cg.savefig(plot_path, dpi=300)
    plt.close(cg.fig)
    print(f"Heatmap saved to {plot_path}")

    # Return the original DataFrame (in case the caller wants to chain further operations)
    return df

def run_clustering_analysis( df: pd.DataFrame, metadata: pd.DataFrame, output_dir: str) -> dict:
    """
    Runs clustering analyses (PCA, MDS, heatmap) on the full dataset and saves results.
    Args:
        df (pd.DataFrame): Protein abundance data (columns = samples, rows = proteins or features).
        metadata (pd.DataFrame): Sample metadata, including 'sample_rep' and 'treatment' columns.
        output_dir (str): Path to directory where outputs (plots, data files) will be written.
    Returns:
        dict: Dictionary containing result DataFrames from PCA, MDS, and heatmap.
              Keys include "pca", "mds", and "heatmap".
    """
    results = {}

    # Perform PCA and save results
    print("Performing PCA...")
    pca_results = generate_pca(df.T, output_dir, metadata=metadata)
    results["pca"] = pca_results

    # Perform MDS and save results
    print("Performing MDS...")
    mds_coords_df = generate_MDS(df.T, output_dir, metadata=metadata)
    results["mds"] = mds_coords_df

    # Generate and save heatmap
    print("Generating heatmap...")
    df_heatmap = generate_heatmap(df, output_dir, metadata=metadata)
    results["heatmap"] = df_heatmap

    return results