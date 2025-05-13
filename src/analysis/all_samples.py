#### analysis of all_samples 
## e.g. pca, mds, heatmap
## as opposed to pairwise treatment analysis
import os
import pandas as pd
import numpy as np

from sklearn.decomposition import PCA
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import seaborn as sns

def generate_pca(
    df: pd.DataFrame,
    output_dir: str,
    metadata: pd.DataFrame = None,
    variance_threshold: float = 0.95,
    plot_components: tuple[int, int] = (1, 2)
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
    n_prot = df.shape[1]

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
    plt.figure(figsize=(4, 3))
    plot_pca = sns.scatterplot(
        data=df_PCA,
        x=f"PC{pc_x}",
        y=f"PC{pc_y}",
        hue="treatment",
    )
    # Annotate points
    for i in range(df_PCA.shape[0]):
        plt.text(
            x=df_PCA[f"PC{pc_x}"].iloc[i] + 0.1,
            y=df_PCA[f"PC{pc_y}"].iloc[i] + 0.1,
            s=df_PCA.index[i],
            fontsize=9,
            ha="center",
            va="bottom",
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
    # Perform MDS
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
    mds_coords_df = pd.DataFrame(
        mds_coords, columns=["MDS1", "MDS2"], index=df.index
    )
    # pull over treatment from metadata
    # df_PCA . look up by index . and map to
    # metadata . temporarily set index to sample_id [ and find treatment ]
    mds_coords_df["treatment"] = mds_coords_df.index.map(
        metadata.set_index("sample_rep")["treatment"]
    )
    #### Generate NMDS plot
    plot_nmds = sns.scatterplot(data=mds_coords_df, x="MDS1", y="MDS2", hue="treatment")
    ## seaborn returns an axis-object rather than a figure, so you can still alter features of it. E.g. axes names:
    plot_nmds.set(xlabel="MDS 1", ylabel="MDS 2")
    # Add sample IDs from the 'target' column with an offset
    for i in range(mds_coords_df.shape[0]):
        plt.text(
            x=mds_coords_df["MDS1"].iloc[i] + 1,  # Add a small x-offset
            y=mds_coords_df["MDS2"].iloc[i] + 1,  # Add a small y-offset
            s=mds_coords_df.index[i],
            fontsize=9,
            ha="center",  # Horizontal alignment
            va="bottom",  # Vertical alignment
        )
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


def generate_heatmap(df: pd.DataFrame, output_dir: str,
                     metadata: pd.DataFrame = None) -> pd.DataFrame:
    """
    Creates a clustered heatmap of protein abundance data using Euclidean distance.

    Proteins are clustered along rows and samples along columns. Dendrograms and clustering
    are based on the transposed abundance matrix.

    Args:
        df (pd.DataFrame): Protein abundance data (samples as rows, proteins as columns).
        output_dir (str): Directory to save the heatmap plot.

    Returns:
        pd.DataFrame: Transposed abundance data used to generate the heatmap.
    """
    #### need the other orientation for heatmap
    # transpose and add sample ids as colnames
    df_heat = df.T
    df_heat.columns = df.index
    n_prot = df_heat.shape[0]
    plot_title = "heatmap of euclidean distance (n protein = " + str(n_prot) + ")"
    
    # Create col_colors using metadata['Colour'] if provided
    col_colors = None
    if metadata is not None:
        required_cols = {"sample_rep", "treatment", "colours"}
        if not required_cols.issubset(metadata.columns):
            raise ValueError(f"Metadata must contain columns: {required_cols}")

        colour_map = metadata.set_index("sample_rep")["colours"]

        missing = df_heat.columns.difference(colour_map.index)
        if not missing.empty:
            raise ValueError(f"Missing colour info for samples: {list(missing)}")

        col_colors = df_heat.columns.map(colour_map)

    # Create clustermap
    heatmap = sns.clustermap(
        df_heat,
        fmt=".2f",
        cmap="viridis",
        col_colors=col_colors,

    )

    # Move x-axis labels (column names) to the top
    heatmap.ax_heatmap.xaxis.set_ticks_position("top")  # Move ticks to the top
    heatmap.ax_heatmap.xaxis.set_label_position("top")  # Move axis label to the top
    # Remove dendrogram tick labels
    heatmap.ax_heatmap.set_xticks(
        np.arange(len(df_heat.columns))
    )  # Set tick positions for columns
    heatmap.ax_heatmap.set_xticklabels(
        df_heat.columns, rotation=45, ha="left"
    )  # Use only column names
    # Explicitly disable extra tick labels from the dendrogram
    heatmap.ax_heatmap.tick_params(axis="x", which="both", bottom=False, top=True)
    
    # Set labels and title
    plt.title(plot_title)
    # plt.tight_layout()
    # Save plot and PCA data
    plot_path = os.path.join(output_dir, "plots", "heatmap_plot.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"heatmap saved to {plot_path}")
    return df_heat