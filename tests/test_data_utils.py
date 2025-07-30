import sys
from pathlib import Path

import fnmatch
import pandas as pd
import pytest
import os
import shutil
import tempfile
import pytest

from PIL import Image


from src.utils.data_utils import (apply_row_id_config, get_subset,
                                    validate_metadata, validate_proteindata,
                                    combine_plots, combine_csv_files)

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))


# ————————————————————————————————————————————
# test phospho row ID function
# ————————————————————————————————————————————

def test_apply_row_id_config_phospho_fields_present():
    df = pd.DataFrame({
        "site": ["S10", "T20", None],
        "modification": ["phospho", "phospho", "acetyl"]
    }, index=["P1", "P2", "P3"])

    config = {
        "data_type": "phospho",
        "phospho_row_id": {
            "fields": ["site", "modification"],
            "missing_value": "NA"
        }
    }

    result = apply_row_id_config(df.copy(), config)
    expected_index = [
        "P1__S10_phospho",
        "P2__T20_phospho",
        "P3__NA_acetyl"
    ]
    assert list(result.index) == expected_index

def test_apply_row_id_config_non_phospho():
    df = pd.DataFrame({"A": [1, 2]}, index=["X", "Y"])
    config = {"data_type": "other"}
    result = apply_row_id_config(df.copy(), config)
    assert result.equals(df)

def test_apply_row_id_config_missing_field():
    df = pd.DataFrame({
        "site": ["S10", "T20"]
    }, index=["P1", "P2"])

    config = {
            "data_type": "phospho",
            "phospho_row_id": {
                "fields": ["site", "modification"],  # list of strings
                "missing_value": "NA"
            }
        }

    result = apply_row_id_config(df.copy(), config)
    # Only 'site' is available
    expected_index = [
        "P1__S10",
        "P2__T20"
    ]
    assert list(result.index) == expected_index

def test_apply_row_id_config_no_fields():
    df = pd.DataFrame({"A": [1, 2]}, index=["X", "Y"])
    config = {
        "data_type": "phospho",
        "phospho_row_id": {
            "fields": []
        }
    }
    result = apply_row_id_config(df.copy(), config)
    assert result.equals(df)


# ————————————————————————————————————————————
# test get_subset()
# ————————————————————————————————————————————
def test_get_subset_success():
    # Mock protein abundance dataframe
    df = pd.DataFrame({
        "A_rep1": [1.0, 2.0],
        "A_rep2": [1.5, 2.5],
        "B_rep1": [3.0, 4.0],
    }, index=["Prot1", "Prot2"])

    # Metadata with subset variable
    metadata = pd.DataFrame({
        "sample_rep": ["A_rep1", "A_rep2", "B_rep1"],
        "treatment": ["A", "A", "B"],
    })

    # Expected subset term and variable
    subset_term = "A"
    subset_variable = "treatment"

    result = get_subset(df, subset_term, metadata, subset_variable)

    assert list(result.columns) == ["A_rep1", "A_rep2"]
    assert result.shape == (2, 2)

def test_get_subset_no_match():
    df = pd.DataFrame({
        "A_rep1": [1.0, 2.0],
        "B_rep1": [3.0, 4.0],
    }, index=["Prot1", "Prot2"])

    metadata = pd.DataFrame({
        "sample_rep": ["A_rep1", "B_rep1"],
        "treatment": ["A", "B"],
    })

    subset_term = "C"  # Doesn't exist in metadata
    subset_variable = "treatment"

    with pytest.raises(ValueError, match="No matches found for subset: C"):
        get_subset(df, subset_term, metadata, subset_variable)


# ————————————————————————————————————————————
#   test validate_metadata()
# ————————————————————————————————————————————
def good_metadata():
    return pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "treatment": ["Ctl", "Drug"],
            "replicate": [1, 2],
            "protein_abundance_name": ["col1", "col2"],
            "sample_rep": ["S1_1", "S2_1"]
        }
    )


def test_validate_metadata_passes():
    df = good_metadata()
    validate_metadata(df)  # should not raise


def test_missing_column_raises():
    df = good_metadata().drop(columns="treatment")
    with pytest.raises(ValueError, match="Missing required columns"):
        validate_metadata(df)


def test_nan_raises():
    df = good_metadata()
    df.loc[1, "sample_id"] = None
    with pytest.raises(ValueError, match="contains missing"):
        validate_metadata(df)

def test_duplicate_sample_replicate_pair_raises():
    df = good_metadata()
    df.loc[1, "sample_id"] = "S1"
    df.loc[1, "replicate"] = 1
    df.loc[1, "sample_rep"] = "S1_1"
    with pytest.raises(ValueError, match="must be unique"):
        validate_metadata(df)


def test_nonnumeric_replicate_raises():
    df = good_metadata()
    df["replicate"] = ["one", "two"]
    with pytest.raises(ValueError, match="'replicate' must be numeric"):
        validate_metadata(df)


def test_duplicate_protein_abundance_name_raises():
    df = good_metadata()
    df.loc[1, "protein_abundance_name"] = "col1"
    with pytest.raises(ValueError, match="must not contain duplicate"):
        validate_metadata(df)

# ————————————————————————————————————————————
#  test validate_proteindata()
# ————————————————————————————————————————————
def test_empty_data_raises():
    df = pd.DataFrame()
    meta = pd.DataFrame({"sample_rep": []})
    with pytest.raises(ValueError, match="data is empty"):
        validate_proteindata(df, meta)


def test_duplicate_protein_ids_case_insensitive():
    df = pd.DataFrame(
        {"s1": [1], "s2": [2]}, index=["GeneA", "genea"]  # duplicate ignoring case
    )
    meta = pd.DataFrame(
        {"sample_rep": ["s1", "s2"], "protein_abundance_name": ["s1", "s2"]}
    )
    with pytest.raises(ValueError, match="Protein identifiers.*unique"):
        validate_proteindata(df, meta)

def test_column_metadata_mismatch_raises():
    df = pd.DataFrame({"S1": [1], "S2": [2]}, index=["Gene1"])
    meta = pd.DataFrame({"sample_rep": ["S1", "S3"]})
    with pytest.raises(ValueError, match="Mismatch between metadata samples"):
        validate_proteindata(df, meta)

def test_non_numeric_values_raises():
    df = pd.DataFrame({"S1": [1], "S2": ["bad"]}, index=["Gene1"])
    meta = pd.DataFrame({"sample_rep": ["S1", "S2"]})
    with pytest.raises(ValueError, match="non-numeric"):
        validate_proteindata(df, meta)

def test_missing_values_raises():
    df = pd.DataFrame({"S1": [1], "S2": [float("nan")]}, index=["Gene1"])
    meta = pd.DataFrame({"sample_rep": ["S1", "S2"]})
    with pytest.warns(UserWarning, match="contains missing"):
        validate_proteindata(df, meta)


# ————————————————————————————————————————————
#  test combine_plots()
# ————————————————————————————————————————————
@pytest.fixture
def image_dir_with_pngs():
    temp_dir = tempfile.mkdtemp()
    subdir1 = os.path.join(temp_dir, "group1")
    subdir2 = os.path.join(temp_dir, "group2")
    os.makedirs(subdir1)
    os.makedirs(subdir2)
    #
    for idx, subdir in enumerate([subdir1, subdir2]):
        img = Image.new("RGB", (100, 100), color=(255, 0, 0))
        img_path = os.path.join(subdir, f"volcano_plot.png")
        img.save(img_path)
    #
    yield temp_dir
    shutil.rmtree(temp_dir)

def test_combine_plots_creates_output(image_dir_with_pngs):
    output_dir = tempfile.mkdtemp()

    combine_plots(
        search_term="volcano_plot.png",
        search_path=image_dir_with_pngs,
        output_dir=output_dir,
        img_size=(100, 100),
        max_cols=2
    )

    # Reconstruct the expected path based on function logic
    combined_path = os.path.join(
        output_dir, "plots", "combined_volcano_plot.png"
    )

    # Assertions
    assert os.path.exists(combined_path)
    assert combined_path.endswith("combined_volcano_plot.png")

    with Image.open(combined_path) as img:
        assert img.size[0] >= 100
        assert img.size[1] >= 100

    shutil.rmtree(output_dir)



# ————————————————————————————————————————————
#  test combine_csv_files()
# ————————————————————————————————————————————
@pytest.fixture
def mock_csv_structure():
    temp_dir = tempfile.mkdtemp()
    os.makedirs(os.path.join(temp_dir, "DrugA"))
    os.makedirs(os.path.join(temp_dir, "DrugB"))

    df1 = pd.DataFrame({
        "gene": ["G1", "G2", "G3"],
        "logFC": [1.5, -0.8, 2.1],
        "pval": [0.01, 0.05, 0.03]
    })
    df2 = pd.DataFrame({
        "gene": ["G4", "G5", "G6"],
        "logFC": [-2.0, 0.7, 1.1],
        "pval": [0.02, 0.06, 0.01]
    })

    df1.to_csv(os.path.join(temp_dir, "DrugA", "top_hits.csv"), index=False)
    df2.to_csv(os.path.join(temp_dir, "DrugB", "top_hits.csv"), index=False)

    yield temp_dir
    shutil.rmtree(temp_dir)

def test_combine_csv_files_outputs_combined(mock_csv_structure):
    output_dir = mock_csv_structure
    combine_csv_files(
        filename="top_hits.csv",
        output_dir=output_dir,
        top_n=2,
        new_column="treatment_pair"
    )

    # Manually construct expected output path
    combined_path = os.path.join(output_dir, "data", "combined_top_hits.csv")

    # Assert file was written
    assert os.path.exists(combined_path)

    # Read the file and verify contents
    df = pd.read_csv(combined_path)
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 4  # 2 rows from each of 2 CSVs
    assert "treatment_pair" in df.columns
    assert set(df["treatment_pair"]) == {"DrugA", "DrugB"}













@pytest.fixture
def create_dummy_csvs(tmp_path):
    def _create_csvs(filenames, subdirs, content_rows=20):
        for subdir in subdirs:
            dir_path = tmp_path / subdir
            dir_path.mkdir(parents=True, exist_ok=True)
            for fname in filenames:
                df = pd.DataFrame({
                    "gene": [f"gene_{i}" for i in range(content_rows)],
                    "logFC": [i for i in range(content_rows)],
                })
                df.to_csv(dir_path / fname, index=False)
        return tmp_path
    return _create_csvs

def test_combine_csv_files_top_abs_logfc(tmp_path):
    # Arrange
    filename = "top_genes.csv"
    subdirs = ["treatA_vs_ctrl", "treatB_vs_ctrl"]
    top_n = 5

    # Create dummy CSVs with mixed positive/negative logFC
    for i, subdir in enumerate(subdirs):
        dir_path = tmp_path / subdir
        dir_path.mkdir(parents=True, exist_ok=True)

        # Shuffle logFC values: include both strong + and - values
        logfc_values = [-8, 2, 0, 9, -3, 4, -10, 1, 5, -7]
        df = pd.DataFrame({
            "gene": [f"{subdir}_gene_{j}" for j in range(len(logfc_values))],
            "logFC": logfc_values
        })
        df.to_csv(dir_path / filename, index=False)

    output_dir = tmp_path
    expected_path = tmp_path / "results" / "data" / f"combined_{filename}"
    expected_path.parent.mkdir(parents=True, exist_ok=True)

    combine_csv_files(
        filename=filename,
        output_dir=str(output_dir),
        output_filename=str(expected_path),
        top_n=top_n,
        new_column="treatment_pair",
        sort_by_logfc=True
    )

    # Assert
    assert expected_path.exists()
    df = pd.read_csv(expected_path)

    # One check already:
    assert df.shape[0] == top_n * len(subdirs)
    assert "treatment_pair" in df.columns
    assert set(df["treatment_pair"]) == set(subdirs)

    # NEW: check for largest absolute logFC values per treatment
    expected_top = sorted([-8, 2, 0, 9, -3, 4, -10, 1, 5, -7], key=lambda x: abs(x), reverse=True)[:top_n]
    for group in subdirs:
        subset = df[df["treatment_pair"] == group]
        assert subset.shape[0] == top_n
        actual_logfc = list(subset["logFC"])
        # Should match top N by absolute value (ignoring order)
        assert sorted(actual_logfc, key=lambda x: abs(x), reverse=True) == expected_top
