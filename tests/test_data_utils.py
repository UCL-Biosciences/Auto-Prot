import sys
from pathlib import Path

import fnmatch
import pandas as pd
import pytest

from PIL import Image


from src.utils.data_utils import get_subset, validate_metadata, validate_proteindata, combine_plots, combine_csv_files

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))

# ————————————————————————————————————————————
# test get_subset()
# ————————————————————————————————————————————


def test_get_subset_happy_path():
    df = pd.DataFrame({"value": [1, 2]}, index=["sample_A", "sample_B"])
    result = get_subset(df, "A")
    assert len(result) == 1
    assert "sample_A" in result.index


def test_get_subset_no_matches():
    df = pd.DataFrame({"value": [1, 2]}, index=["sample_A", "sample_B"])
    with pytest.raises(ValueError, match="No matches found for subset"):
        get_subset(df, "C")


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
def test_duplicate_sample_ids_case_insensitive():
    df = pd.DataFrame({"S1": [1], "s1": [2]}, index=["Gene1", "Gene2"])
    meta = pd.DataFrame(
        {"sample_rep": ["S1", "s1"], "protein_abundance_name": ["S1", "s1"]}
    )
    with pytest.raises(ValueError, match="Sample identifiers.*unique"):
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


# ————————————————————————————————————————————
#  test combine_plots()
# ————————————————————————————————————————————

@pytest.fixture
def create_dummy_images(tmp_path):
    def _create_images(filenames, subdirs=None):
        subdirs = subdirs or [""]
        for subdir in subdirs:
            dir_path = tmp_path / subdir
            dir_path.mkdir(parents=True, exist_ok=True)
            for fname in filenames:
                img_path = dir_path / fname
                img = Image.new("RGB", (100, 100), (255, 0, 0))
                img.save(img_path)
        return tmp_path
    return _create_images


def test_combine_plots_creates_output(tmp_path, create_dummy_images):
    # Arrange
    filenames = ["volcano_plot.png", "volcano_plot.png"]
    search_term = "volcano_plot.png"
    create_dummy_images(filenames, subdirs=["groupA", "groupB"])
    output_dir = tmp_path / "output"
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True)  # <-- ensure plots subfolder exists

    expected_output = plots_dir / f"combined_{search_term.replace('.png', '')}.png"

    # Act
    combine_plots(
        search_term=search_term,
        search_path=str(tmp_path),
        output_dir=str(output_dir)
    )

    # Assert
    assert expected_output.exists()
    img = Image.open(expected_output)
    assert img.size[0] > 0 and img.size[1] > 0

# ————————————————————————————————————————————
#  test combine_csv_files()
# ————————————————————————————————————————————

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
