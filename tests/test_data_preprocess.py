import numpy as np
import pandas as pd

from src.processing.data_preprocess import process_prot_data, view_prot_distributions

#########
# test process_prot_data (log, median normalise, impute)
#########


def test_process_prot_data_full_pipeline():
    df = pd.DataFrame(
        {
            "s1": [0, 1, 2, 3],
            "s2": [4, 5, 6, 0],
            "s3": [0, 8, 0, 10],
            "s4": [0, 11, 12, 13],
            "s5": [0, 14, 15, 16],
        },
        index=["drop_me", "p1", "p2", "p3"],
    )  # drop_me has too many zeros/NaNs

    result = process_prot_data(df.copy())

    # Expect all keys returned
    assert set(result.keys()) == {"df", "df_log2", "df_median_norm", "df_imp"}

    # 1. Zeros replaced with NaN, drop_me removed
    df_clean = result["df"]
    assert "drop_me" not in df_clean.index
    assert (df_clean == 0).sum().sum() == 0

    # 2. Log2 transformation is applied correctly
    expected_log2 = np.log2(df_clean)
    pd.testing.assert_frame_equal(result["df_log2"], expected_log2)

    # 3. Median-normalisation is applied
    expected_norm = expected_log2.sub(expected_log2.median(axis=0), axis=1)
    pd.testing.assert_frame_equal(result["df_median_norm"], expected_norm)

    # 4. Imputation output has same shape and no NaNs
    df_imp = result["df_imp"]
    assert df_imp.shape == df_clean.shape
    assert df_imp.isna().sum().sum() == 0
    assert df_imp.index.equals(df_clean.index)
    assert df_imp.columns.equals(df_clean.columns)


#########
# test plot creating prot distributions
#########


def test_view_prot_distributions_creates_plots(tmp_path):

    # 1. Create dummy DataFrames
    dfs = []
    for _i in range(4):  # one for each subplot
        dfs.append(pd.DataFrame({"s1": [1, 2, 3], "s2": [4, 5, 6]}))
    plot_titles = [f"Step {i}" for i in range(4)]

    # 2. Create metadata
    metadata = pd.DataFrame({"sample_rep": ["s1", "s2"], "treatment": ["Ctl", "Drug"]})

    # 3. Create output structure
    out_dir = tmp_path / "output"
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True)

    # 4. Run the plotting function
    view_prot_distributions(dfs, plot_titles, metadata, str(out_dir))

    # 5. Assert that output files exist
    boxplot_path = plot_dir / "boxplots_preProcessing_all_samples_plot.png"
    kdeplot_path = plot_dir / "KDE_preProcessing_all_treatments_plot.png"

    assert boxplot_path.exists()
    assert kdeplot_path.exists()
