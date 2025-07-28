import pandas as pd
import numpy as np
import tempfile
import os
from pathlib import Path
from src.analysis.all_samples import generate_pca, generate_MDS, generate_heatmap

## for generating mini datasets
def gen_spiked_data(effect_size=2.0, n_prots = 100, n_samples = 10):
    np.random.seed(0)
    df = pd.DataFrame(
        np.random.rand(n_samples, n_prots),
        index=[f"S{i}" for i in range(n_samples)],
        columns=[f"Prot{i}" for i in range(n_prots)]
    )
    # Add metadata
    metadata = pd.DataFrame({
        "sample_rep": [f"S{i}" for i in range(n_samples)],
        "treatment": ["ctrl"] * int((n_samples / 2)) + ["treated"] * int((n_samples/2)),
        "colours": ["#1f77b4"] * 5 + ["#ff7f0e"] * 5
    })
    #
    # Spike in signal: increase abundance for first 5 proteins in treated samples
    treated_samples = metadata[metadata["treatment"] == "treated"]["sample_rep"]
    spiked_proteins = [f"Prot{i}" for i in range(25)]
    df.loc[treated_samples, spiked_proteins] += effect_size
    #
    return df, metadata, spiked_proteins

# Generate a small dataset for testing
df, metadata, spiked_proteins = gen_spiked_data()

#########
# test pca
#########
def test_generate_pca_basic(df = df, metadata = metadata):

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

        ### the above makes a big difference in 25% of proteins.
        ### It adds the same amount so the variance should be easy to capture in one PC
        ### therefore, we should see big difference in control and treated
        # Check that samples separate by treatment
        treated = metadata[metadata["treatment"] == "treated"]["sample_rep"]
        control = metadata[metadata["treatment"] == "ctrl"]["sample_rep"]
        #
        # all samples within trt should be either positive or negative on PC1
        pc1 = df_pca.iloc[:,0]
        if (pc1.loc[treated] > 0).all() and (pc1.loc[control] < 0).all():
            pass  # OK
        elif (pc1.loc[treated] < 0).all() and (pc1.loc[control] > 0).all():
            pass  # Also OK
        else:
            raise AssertionError("function did not cleanly separate groups on axis 1")

#########
# test mds
#########

def test_generate_mds_basic(df = df, metadata = metadata):

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

        # Check that MDS separates treated and control groups
        treated = df_mds[df_mds["treatment"] == "treated"]
        control = df_mds[df_mds["treatment"] == "ctrl"]

        # Basic separation check on MDS1
        mean_treated = treated.iloc[:, :2].mean().mean()
        mean_control = control.iloc[:, :2].mean().mean()

        assert abs(mean_treated - mean_control) > 1.0, "MDS failed to separate spiked groups"
    

def test_generate_heatmap_basic(df = df.T, # note heatmap expects samples in cols
                                metadata = metadata):

    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "plots"))

        df_out = generate_heatmap(df, tmpdir, metadata)

        # Check output DataFrame structure
        assert df_out.shape == df.shape
        assert all(df_out.columns == [f"S{i}" for i in range(10)])

        # Check that heatmap image is saved
        plot_path = Path(tmpdir) / "plots" / "heatmap_plot.png"
        assert plot_path.exists(), "Heatmap plot was not saved."
