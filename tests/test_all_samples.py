import pandas as pd
import numpy as np
import tempfile
import os
from pathlib import Path
from src.analysis.all_samples import generate_pca, generate_MDS, generate_heatmap

#########
# test pca
#########


def test_generate_pca_basic():
    # Simulate input data
    np.random.seed(0)
    df = pd.DataFrame(
        np.random.rand(5, 10),  # 5 samples, 10 proteins
        index=[f"S{i}" for i in range(5)],
        columns=[f"Prot{i}" for i in range(10)]
    )
    metadata = pd.DataFrame({
        "sample_rep": [f"S{i}" for i in range(5)],
        "treatment": ["ctrl", "ctrl", "treated", "treated", "ctrl"]
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create expected subdirs
        os.makedirs(os.path.join(tmpdir, "plots"))
        os.makedirs(os.path.join(tmpdir, "data"))

        df_pca = generate_pca(df, tmpdir, metadata)

        # Check output shape
        assert "PC1" in df_pca.columns
        assert "PC2" in df_pca.columns
        assert "treatment" in df_pca.columns
        assert df_pca.shape[0] == df.shape[0]

        # Check that plot and data files were created
        plot_path = Path(tmpdir) / "plots" / "pca_plot.png"
        data_path = Path(tmpdir) / "data" / "pca_data.csv"
        assert plot_path.exists()
        assert data_path.exists()

        # Check explained variance is reasonable
        assert df_pca.filter(like="PC").shape[1] >= 2

        print("Test passed: PCA function works with minimal data.")

#########
# test mds
#########

def test_generate_mds_basic():
    # Create dummy data
    np.random.seed(1)
    df = pd.DataFrame(
        np.random.rand(5, 8),  # 5 samples, 8 proteins
        index=[f"S{i}" for i in range(5)],
        columns=[f"Prot{i}" for i in range(8)]
    )
    metadata = pd.DataFrame({
        "sample_rep": [f"S{i}" for i in range(5)],
        "treatment": ["A", "A", "B", "B", "A"]
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        # Ensure plots and data folders exist
        os.makedirs(os.path.join(tmpdir, "plots"))
        os.makedirs(os.path.join(tmpdir, "data"))

        df_mds = generate_MDS(df, tmpdir, metadata)

        # Check structure of result
        assert "MDS1" in df_mds.columns
        assert "MDS2" in df_mds.columns
        assert "treatment" in df_mds.columns
        assert df_mds.shape[0] == df.shape[0]

        # Check plot and data file existence
        plot_path = Path(tmpdir) / "plots" / "mds_plot.png"
        data_path = Path(tmpdir) / "data" / "mds_data.csv"
        assert plot_path.exists()
        assert data_path.exists()

        # Check for missing treatment
        assert df_mds["treatment"].isnull().sum() == 0

        print("Test passed: MDS function runs and outputs valid files.")
    

def test_generate_heatmap_basic():
    # Dummy protein abundance data
    np.random.seed(42)
    df = pd.DataFrame(
        np.random.rand(4, 6),  # 4 samples, 6 proteins
        index=["S1", "S2", "S3", "S4"],
        columns=[f"P{i}" for i in range(6)]
    )

    # Dummy metadata with treatment and colour
    metadata = pd.DataFrame({
        "sample_rep": ["S1", "S2", "S3", "S4"],
        "treatment": ["control", "control", "treated", "treated"],
        "colours": ["#1f77b4", "#1f77b4", "#ff7f0e", "#ff7f0e"]
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "plots"))

        df_out = generate_heatmap(df, tmpdir, metadata)

        # Check output DataFrame structure
        assert df_out.shape == (6, 4)  # Transposed: 6 proteins, 4 samples
        assert all(df_out.columns == ["S1", "S2", "S3", "S4"])

        # Check that heatmap image is saved
        plot_path = Path(tmpdir) / "plots" / "heatmap_plot.png"
        assert plot_path.exists(), "Heatmap plot was not saved."

        print("Test passed: heatmap with metadata colour bar generated successfully.")