import sys
from pathlib import Path

import pandas as pd
import pytest

from src.utils.data_utils import get_subset, validate_metadata, validate_proteindata

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))

# ————————————————————————————————————————————
#
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
