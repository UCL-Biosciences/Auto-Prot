import pandas as pd
import numpy as np
import tempfile
import os
from pathlib import Path
from src.analysis.all_samples import generate_pca, generate_MDS, generate_heatmap

def gen_spiked_data(effect_size=2.0, n_prots=100, n_samples=10):
    """
    Generate synthetic protein abundance data with controlled signal spiking.

    Simulates a dataset where a subset of proteins (25%) have increased abundance
    in the "treated" group. Returns the abundance DataFrame, metadata DataFrame,
    and a list of spiked protein names.

    Args:
        effect_size (float): The fold-change to add to spiked proteins in treated samples.
        n_prots (int): Total number of proteins (columns).
        n_samples (int): Total number of samples (rows).

    Returns:
        df (pd.DataFrame): Protein abundance values (samples × proteins).
        metadata (pd.DataFrame): Sample metadata including treatment and colour.
        spiked_proteins (list): Names of proteins with increased abundance.
    """
    np.random.seed(0)
    df = pd.DataFrame(
        np.random.rand(n_samples, n_prots),
        index=[f"S{i}" for i in range(n_samples)],
        columns=[f"Prot{i}" for i in range(n_prots)]
    )
    metadata = pd.DataFrame({
        "sample_rep": [f"S{i}" for i in range(n_samples)],
        "treatment": ["ctrl"] * int(n_samples / 2) + ["treated"] * int(n_samples / 2),
        "colours": ["#1f77b4"] * 5 + ["#ff7f0e"] * 5
    })

    # Spike signal in the first 25 proteins for treated samples
    treated_samples = metadata[metadata["treatment"] == "treated"]["sample_rep"]
    spiked_proteins = [f"Prot{i}" for i in range(25)]
    df.loc[treated_samples, spiked_proteins] += effect_size

    return df, metadata, spiked_proteins


# Shared test fixture data
df, metadata, spiked_proteins = gen_spiked_data()


def test_generate_pca_basic(df=df, metadata=metadata):
    """
    Test that PCA runs and produces expected outputs including:
    - Correct number of samples
    - Presence of PC1 and PC2
    - Group separation based on treatment
    - Saved plot and data
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "plots"))
        os.makedirs(os.path.join(tmpdir, "data"))

        df_pca = generate_pca(df, tmpdir, metadata)

        # Basic output checks
        assert "PC1" in df_pca.columns
        assert "PC2" in df_pca.columns
        assert "treatment" in df_pca.columns
        assert df_pca.shape[0] == df.shape[0]

        # Files written
        assert Path(tmpdir, "plots", "pca_plot.png").exists()
        assert Path(tmpdir, "data", "pca_data.csv").exists()

        # Dimensionality check
        assert df_pca.filter(like="PC").shape[1] >= 2

        # Ensure separation of treated vs control along PC1
        treated = metadata[metadata["treatment"] == "treated"]["sample_rep"]
        control = metadata[metadata["treatment"] == "ctrl"]["sample_rep"]
        pc1 = df_pca["PC1"]
        if (pc1.loc[treated] > 0).all() and (pc1.loc[control] < 0).all():
            pass
        elif (pc1.loc[treated] < 0).all() and (pc1.loc[control] > 0).all():
            pass
        else:
            raise AssertionError("function did not cleanly separate groups on axis 1")


def test_generate_mds_basic(df=df, metadata=metadata):
    """
    Test that MDS runs and outputs expected structure:
    - Includes MDS1, MDS2 and treatment
    - Correct number of samples
    - Saved plot and data
    - Treated and control means differ
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "plots"))
        os.makedirs(os.path.join(tmpdir, "data"))

        df_mds = generate_MDS(df, tmpdir, metadata)

        # Output structure
        assert "MDS1" in df_mds.columns
        assert "MDS2" in df_mds.columns
        assert "treatment" in df_mds.columns
        assert df_mds.shape[0] == df.shape[0]

        # Output files
        assert Path(tmpdir, "plots", "mds_plot.png").exists()
        assert Path(tmpdir, "data", "mds_data.csv").exists()

        # No missing treatment values
        assert df_mds["treatment"].isnull().sum() == 0

        # Group separation check
        treated = df_mds[df_mds["treatment"] == "treated"]
        control = df_mds[df_mds["treatment"] == "ctrl"]
        mean_treated = treated.iloc[:, :2].mean().mean()
        mean_control = control.iloc[:, :2].mean().mean()
        assert abs(mean_treated - mean_control) > 1.0, "MDS failed to separate spiked groups"


def test_generate_heatmap_basic(df=df.T, metadata=metadata):
    """
    Test that the heatmap function runs:
    - Output shape matches input
    - Sample columns are preserved
    - Heatmap image is saved
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "plots"))

        df_out = generate_heatmap(df, tmpdir, metadata)

        assert df_out.shape == df.shape
        assert all(df_out.columns == [f"S{i}" for i in range(10)])
        assert Path(tmpdir, "plots", "heatmap_plot.png").exists()
