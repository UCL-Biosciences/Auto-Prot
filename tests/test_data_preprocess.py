import numpy as np
import pandas as pd
import tempfile
import os
import subprocess

from src.processing.data_preprocess import process_prot_data, view_prot_distributions, filter_proteins_by_group_missingness, impute_pimms_cf, impute_prot_data_histgradboost, normalise_vsn
from tests.test_all_samples import gen_spiked_data

import pytest

## generate a small dataset for testing
df, metadata, spiked_proteins = gen_spiked_data(n_prots=1000)

df = df.T # Transpose to have proteins as rows, samples as columns

# ######### test for filtering based on missingness
def test_filter_proteins_by_group_missingness_with_spike(df = df,
                                                         metadata = metadata, 
                                                         spiked_proteins = spiked_proteins):
    
    # Simulate missingness: remove values randomly in some proteins
    # Prot90–99 will have high missingness in 'ctrl' group
    for prot in [f"Prot{i}" for i in range(90, 100)]:
        df.loc[ prot, ["S0", "S1", "S2"] ] = np.nan

    # Run filter at 0.6 threshold (must be present in ≥60% of samples per group)
    filtered_df = filter_proteins_by_group_missingness(df, metadata, threshold=0.6)

    # Make sure spiked proteins are kept (they should be complete)
    assert all(p in filtered_df.index for p in spiked_proteins[:5])

    # Make sure high-missingness proteins are removed
    for prot in [f"Prot{i}" for i in range(90, 100)]:
        assert prot not in filtered_df.index


def test_normalise_vsn_runs():
    df, _, _ = gen_spiked_data()
    df = df.T
    df = df.replace(0, np.nan)  # VSN doesn't like 0s

    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = os.path.join(tmpdir, "input.csv")
        out_path = os.path.join(tmpdir, "normalised.csv")
        plot_path = os.path.join(tmpdir, "meanSdPlot.png")

        df.to_csv(in_path)

        try:
            normalise_vsn(in_path, out_path, plot_path)
        except FileNotFoundError:
            pytest.skip("Rscript or conda environment not found — skipping VSN test")
        except subprocess.CalledProcessError as e:
            pytest.fail(f"VSN R script failed: {e}")

        # Check outputs
        assert os.path.exists(out_path), "Expected normalised output file not found"
        assert os.path.exists(plot_path), "Expected plot file not found"

        # Optionally check structure of output
        df_norm = pd.read_csv(out_path, index_col=0)
        assert df_norm.shape == df.shape



@pytest.mark.parametrize("method", ["hist_grad_boost", "pimms_collabfilter"])
def test_imputation_methods_run(method):
    df, metadata, _ = gen_spiked_data()

    # Log2 transform and randomly mask 10% of the values
    df_log2 = np.log2(df)
    df_masked = df_log2.mask(np.random.rand(*df_log2.shape) < 0.1)

    # Filter out sparse proteins (rows)
    df_filtered = df_masked.loc[df_masked.notna().mean(axis=1) > 0.5]

    # Transpose for imputation (samples × proteins)
    df_filtered_T = df_filtered.T

    # Run selected method
    if method == "hist_grad_boost":
        # Normalise by median before using sklearn-based imputer
        df_norm = df_filtered.sub(df_filtered.median(axis=0), axis=1)
        df_norm_t = df_norm.T
        df_imp = impute_prot_data_histgradboost(df_filtered, df_norm_t)
    elif method == "pimms_collabfilter":
        df_imp = impute_pimms_cf(df_filtered_T)
    else:
        raise ValueError(f"Unknown imputation method: {method}")

    # Structural checks
    # df_imp comes back with same shape as df_filtered
    # ML imputers require features (prots) in cols, but the function returns them in original shape
    assert df_imp.shape == df_filtered.shape
    assert df_imp.isna().sum().sum() == 0
    assert all(df_imp.columns == df_filtered.columns)
    assert all(df_imp.index == df_filtered.index)


def test_process_prot_data_runs():
    df, metadata, _ = gen_spiked_data()
    df = df.T  # Transpose to have proteins as rows, samples as columns
    df = df.replace(0, np.nan)  # Ensure no zeros before log2

    config = {
        "missing_threshold": 0.6,
        "normalise_method": "sample-median",          # avoid VSN dependency
        "imputation_method": "hist_grad_boost",       # faster and local
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "plots"), exist_ok=True)
        os.makedirs(os.path.join(tmpdir, "data"), exist_ok=True)

        outputs = process_prot_data(df, config, tmpdir, metadata)

        # Check expected keys
        expected_keys = {"df", "df_log2", "df_norm", "df_imp"}
        assert expected_keys.issubset(outputs.keys())

        # Check shapes and no NaNs in imputed data
        df_imp = outputs["df_imp"]
        assert isinstance(df_imp, pd.DataFrame)
        assert df_imp.isna().sum().sum() == 0
        assert df_imp.shape[0] > 0 and df_imp.shape[1] > 0

        # Check that log2 transformation was applied
        assert (np.log2(df).sort_index() == outputs["df_log2"].sort_index()).all().all()

        # Optional: check that filtering/log2/normalisation reduced values
        assert (outputs["df_log2"] <= np.log2(df.max().max())).all().all()


