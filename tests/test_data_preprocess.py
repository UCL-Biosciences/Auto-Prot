import os
import subprocess
import tempfile

import numpy as np
import pandas as pd
import pytest

from autoprot.processing.data_preprocess import (
    filter_proteins_by_group_missingness,
    impute_pimms_cf,
    impute_prot_data_histgradboost,
    normalise_vsn,
    process_prot_data,
    view_prot_distributions,
)
from tests.test_all_samples import gen_spiked_data

## generate a small dataset for testing
df, metadata, spiked_proteins = gen_spiked_data(n_prots=1000)

df = df.T  # Transpose to have proteins as rows, samples as columns


# ######### test for filtering based on missingness
def test_filter_proteins_by_group_missingness_with_spike(
    df=df, metadata=metadata, spiked_proteins=spiked_proteins
):
    """Test filtering proteins based on missingness across groups.
    This checks that proteins with high missingness in a group are removed,
    while ensuring that spiked proteins (which should be complete) are retained.
    """

    # Simulate missingness: remove values randomly in some proteins
    # Prot90–99 will have high missingness in 'ctrl' group
    for prot in [f"Prot{i}" for i in range(90, 100)]:
        df.loc[prot, ["S0", "S1", "S2"]] = np.nan

    # filter_proteins_by_group_missingness expects config file with field 'missing_threshold'
    config = {"missing_threshold": 0.6}

    # Run filter at 0.6 threshold (must be present in ≥60% of samples per group)
    filtered_df = filter_proteins_by_group_missingness(df, metadata, config = config)

    # Make sure spiked proteins are kept (they should be complete)
    assert all(p in filtered_df.index for p in spiked_proteins[:5])

    # Make sure high-missingness proteins are removed
    for prot in [f"Prot{i}" for i in range(90, 100)]:
        assert prot not in filtered_df.index


def test_normalise_vsn_runs():
    """Test that the VSN normalisation runs without errors and produces expected outputs.
    This requires an R environment with the VSN package installed.
    """
    # Generate a small dataset with spikes
    df, _, _ = gen_spiked_data()
    df = df.T
    df = df.replace(0, np.nan)  # VSN doesn't like 0s

    # Create temporary files for input and output
    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = os.path.join(tmpdir, "input.csv")
        out_path = os.path.join(tmpdir, "normalised.csv")
        plot_path = os.path.join(tmpdir, "meanSdPlot.png")

        # Save the DataFrame to CSV
        df.to_csv(in_path)

        # Run the VSN normalisation script
        # note the normalise_vsn function calls a separate R script and requires r-limma-env environment
        # full info is given in the docs
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
        # Check that the normalised DataFrame has the same shape as the input
        assert df_norm.shape == df.shape


@pytest.mark.parametrize("method", ["hist_grad_boost", "pimms_collabfilter"])
def test_imputation_methods_run(method):
    """Test that the specified imputation methods run without errors and produce expected outputs.
    This tests both histogram gradient boosting and collaborative filtering methods.
    """
    # Generate a small dataset with spikes
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
        # note no normalisation is needed for PIMMS
        df_imp = impute_pimms_cf(df_filtered_T)
    else:
        raise ValueError(f"Unknown imputation method: {method}")

    # Structural checks
    # df_imp comes back with same shape as df_filtered
    # ML imputers require features (prots) in cols, but the function returns them in original shape
    assert df_imp.shape == df_filtered.shape
    # Check that no NaNs remain after imputation
    assert df_imp.isna().sum().sum() == 0
    # Check that the columns and index match the original filtered DataFrame
    assert all(df_imp.columns == df_filtered.columns)
    assert all(df_imp.index == df_filtered.index)


def test_process_prot_data_runs():
    """Test the full processing pipeline for proteomics data.
    This includes filtering, log2 transformation, normalisation, and imputation.
    It checks that the outputs are as expected and that no NaNs remain after imputation.
    """
    df, metadata, _ = gen_spiked_data()
    df = df.T  # Transpose to have proteins as rows, samples as columns
    df = df.replace(0, np.nan)  # Ensure no zeros before log2

    # Define configuration for processing
    # Note: this requires the VSN R package to be installed in the r-limma-env conda environment
    config = {
        "missing_threshold": 0.6,
        "normalise_method": "vsn",  # test VSN dependency
        "imputation_method": "hist_grad_boost",  # faster and local
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "plots"), exist_ok=True)
        os.makedirs(os.path.join(tmpdir, "data"), exist_ok=True)

        # Run the processing function
        # This function will filter, log2 transform, normalise, and impute the data
        # `process_prot_data` returns a dictionary with processed DataFrames
        outputs = process_prot_data(df, config, tmpdir, metadata)

        # Check expected keys
        expected_keys = {"df", "df_log2", "df_norm", "df_imp"}
        # Ensure all expected keys are present in the outputs
        assert expected_keys.issubset(outputs.keys())

        # Check shapes and no NaNs in imputed data
        df_imp = outputs["df_imp"]
        assert isinstance(df_imp, pd.DataFrame)
        assert df_imp.isna().sum().sum() == 0
        assert df_imp.shape[0] > 0 and df_imp.shape[1] > 0

        # Check that log2 transformation was applied
        # We expect df_log2 to be the log2 of the original df
        assert (np.log2(df).sort_index() == outputs["df_log2"].sort_index()).all().all()


def test_view_prot_distributions_creates_plots():
    """Test that the distribution viewing function creates the expected plots.
    This checks that boxplots and KDE plots are generated for processed data stages.
    """
    # Generate a small dataset with 4 sample replicates
    sample_reps = ["S1_1", "S1_2", "S2_1", "S2_2"]
    dfs = []
    for _ in range(4):
        df = pd.DataFrame(
            np.random.rand(10, len(sample_reps)) * 100,
            columns=sample_reps,
            index=[f"Gene{i}" for i in range(10)],
        )
        dfs.append(df)

    # Create titles for each DataFrame
    # Assuming these are the stages of processing
    titles = ["Raw", "Log2", "Norm", "Imputed"]

    #  Create metadata for the samples
    metadata = pd.DataFrame(
        {"sample_rep": sample_reps, "treatment": ["A", "A", "B", "B"]}
    )

    # tempfile is used to create a temporary directory for plots
    with tempfile.TemporaryDirectory() as tmpdir:
        plots_dir = os.path.join(tmpdir, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        # Call the function to view distributions
        # This will generate boxplots and KDE plots for the provided DataFrames
        # the output of the function is not captured, but it generates plots in the specified directory
        view_prot_distributions(dfs, titles, metadata, tmpdir)

        # Check that the plots were created
        boxplot_path = os.path.join(
            plots_dir, "boxplots_preProcessing_all_samples_plot.png"
        )
        kdeplot_path = os.path.join(
            plots_dir, "KDE_preProcessing_all_treatments_plot.png"
        )

        assert os.path.exists(boxplot_path), "Boxplot not generated"
        assert os.path.exists(kdeplot_path), "KDE plot not generated"
