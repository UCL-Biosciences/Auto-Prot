import sys
from pathlib import Path

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))

import json
import os
import yaml

import numpy as np
import pandas as pd
import pytest
from pathlib import Path


import src.processing.data_preprocess as dpp
from src.processing.data_processing import (
    clean_data,
    clean_meta,
    clean_prot,
    process_data,
    prot_summary,
)

from src.utils.data_utils import normalise_column_names

# ————————————————————————————————————————————
# Test clean_meta()
# ————————————————————————————————————————————
def test_clean_meta(tmp_path):
    """Test the clean_meta function from data_processing.py."""

    # 1. Build a tiny metadata DataFrame
    df = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S1"],
            "replicate": [1, 1, 2],
            "protein_abundance_name": ["Prot A", "Prot B", "Prot C"],
            "treatment": ["Ctl", "Drug", "Ctl"],
        }
    )
    json_out = os.path.join( tmp_path,  "meta_out.json" )

    # 2. Call clean_meta
    #    - It should return a DataFrame with sample_rep and colours columns
    #    - It should also write a JSON file with summary keys
    #    - The JSON file should contain NUM_SAMPLES, NUM_TREATMENTS, and TREATMENTS keys
    #    - The TREATMENTS string should mention both Ctl and Drug
    result = clean_meta(df.copy(), str(json_out))

    # 3. Check new columns
    assert "protein_abundance_name" in result.columns
    assert list(result["protein_abundance_name"].unique()) == [
        "prot_a",
        "prot_b",
        "prot_c",
    ]
    assert "sample_rep" in result.columns
    # sample_rep should be sample_id_replicate
    assert set(result["sample_rep"]) == {"S1_r1", "S2_r1", "S1_r2"}
    assert "colours" in result.columns
    # There should be exactly as many unique colours as treatments
    assert result["colours"].nunique() == result["treatment"].nunique()

    # 4. JSON file should exist and contain summary keys
    with open(json_out) as f:
        data = json.load(f)
    assert data["NUM_SAMPLES"] == 3
    assert data["NUM_TREATMENTS"] == 2
    # TREATMENTS string should mention both Ctl and Drug
    assert "Ctl:" in data["TREATMENTS"] and "Drug:" in data["TREATMENTS"]


# ————————————————————————————————————————————
# Test clean_prot()
# ————————————————————————————————————————————
def test_clean_prot():
    """Test the clean_prot function from data_processing.py."""
    # 1. Create a prot DataFrame with two realistic protein‐abundance columns plus one extra
    df = pd.DataFrame(
        {
            "sample_1_protein_measured_etc": [10, "x", 30],
            "sample_2_protein_measured_etc": [40, 50, 60],
            "unrelated_column": [7, 8, 9],
        },
        index=["protA", "protB", "protC"],
    )

    # 2. Build corresponding metadata:
    #    - protein_abundance_name matches the long column names
    #    - sample_rep gives the short names we want in the end
    metadata = pd.DataFrame(
        {
            "protein_abundance_name": [
                "sample_1_protein_measured_etc",
                "sample_2_protein_measured_etc",
            ],
            "sample_rep": ["sample_1_1", "sample_2_1"],
        }
    )

    # 3. Call clean_prot
    cleaned_df, _ = clean_prot(df.copy(), metadata)

    # 4. It should drop the unrelated column
    assert "unrelated_column" not in cleaned_df.columns

    # 5. Non‐numeric in the first column becomes NaN
    assert np.isnan(cleaned_df.loc["protB", "sample_1_1"])

    # 6. Columns should be renamed to the short sample_rep names
    assert set(cleaned_df.columns) == {"sample_1_1", "sample_2_1"}

    # 7. Protein IDs (index) remain unchanged
    assert list(cleaned_df.index) == ["protA", "protB", "protC"]


# ————————————————————————————————————————————
# Test prot_summary()
# ————————————————————————————————————————————
def test_prot_summary(tmp_path):
    """Test the prot_summary function from data_processing.py."""
    # 1. Prepare a small DataFrame of processed protein abundances:
    #    Two proteins (rows), two samples (columns)
    df = pd.DataFrame(
        {"sample_1_1": [1.0, 2.0], "sample_2_1": [3.0, 4.0]}, index=["p1", "p2"]
    )
    nrow_original = 3  # pretend one protein was removed earlier

    # 2. Create a starter JSON file with one existing key
    json_out = tmp_path / "prot_summary.json"
    initial = {"FOO": 123}
    json_out.write_text(json.dumps(initial))

    # 3. Call prot_summary
    prot_summary(df, nrow_original, str(json_out))

    # 4. Read the JSON back in
    with open(json_out) as f:
        data = json.load(f)

    # --- Check unchanged keys ---
    assert data["FOO"] == 123

    # --- Check the new summary stats ---
    # Original count
    assert data["NUM_PROTS_OG"] == 3
    # Removed count = original - remaining (3 - 2 = 1)
    assert data["NUM_PROTS_REMOVED"] == 1

    # NUM_PROTS formatted as string with no decimals
    assert data["NUM_PROTS"] == "2"

    # Minimum mean → 2.0, formatted with `:,.0f` → "2"
    assert data["MIN_AVERAGE_ABUNDANCE"] == "2"

    # Maximum mean → 3.0 → "3"
    assert data["MAX_AVERAGE_ABUNDANCE"] == "3"

    # Median of [2.0, 3.0] is 2.5, `:,.0f` rounds to "2"
    assert data["MEDIAN_AVERAGE_ABUNDANCE"] == "2"


# ────────────────────────────────────────────────────
# Test clean_data without metadata
# ────────────────────────────────────────────────────
def test_clean_data_without_metadata():
    """Proteindata branch must error if metadata is None."""
    df = pd.DataFrame({"foo": [1]})
    with pytest.raises(ValueError, match="Metadata is required"):
        clean_data(
            df,
            file_path="my/proteindata.csv",
            metadata=None,
            outPath=".",
            config={},
            json_out="out.json",
        )


# ────────────────────────────────────────────────────
# Test clean_data proteindata happy path
# ────────────────────────────────────────────────────
def test_clean_data(tmp_path, monkeypatch):
    """Test the clean_data function from data_processing.py with proteindata branch."""
    # ────────────────────────────────────────────────────
    # 1) Skip plotting so tests stay headless
    # ────────────────────────────────────────────────────
    monkeypatch.setattr(dpp, "view_prot_distributions", lambda *a, **k: None)

    # ────────────────────────────────────────────────────
    # 2) Tiny metadata DataFrame matching long col names
    # ────────────────────────────────────────────────────
    long_cols = [
        "sample_1_protein_abundance_intensity",
        "sample_2_protein_abundance_intensity",
    ]
    metadata = pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "treatment": ["Ctl", "Drug"],
            "replicate": [1, 1],
            "protein_abundance_name": long_cols,
        }
    )
    # sample_rep will be ["S1_1","S2_1"] when clean_meta runs
    # **Manually create the sample_rep column** (just like clean_meta would)
    metadata["sample_rep"] = (
        metadata["sample_id"] + "_" + metadata["replicate"].astype(str)
    )

    # ────────────────────────────────────────────────────
    # 3) Starter JSON for prot_summary
    # ────────────────────────────────────────────────────
    json_out = Path(tmp_path) / "summary.json"
    with json_out.open("w") as f:
        # Write an empty JSON object to start
        # This will be updated by prot_summary later
        f.write("{}")

    # ────────────────────────────────────────────────────
    # 4) Build proteindata DataFrame with zeros & non-numeric
    # ────────────────────────────────────────────────────
    df = pd.DataFrame(
        {
            long_cols[0]: [1, 0, "x", 6, 6],
            long_cols[1]: [4, 5, 6, 6, 6],
        },
        index=[np.nan, "p2", "p3", "p4", "p4"],
    )

    # ────────────────────────────────────────────────────
    # 5) Pick the imputed DataFrame after process_prot_data, then normalise column names
    # ────────────────────────────────────────────────────
    config = {"df_to_use": "df_imp",
              "data_type": "prot",
              "missing_threshold": 0.5,
              "normalise_method": "sample-median",
              "imputation_method": "hist_grad_boost"}

    df = normalise_column_names(df, file_path="path/to/proteindata.csv", config = config)

    # ────────────────────────────────────────────────────
    # 6) Run the pipeline
    # ────────────────────────────────────────────────────
    cleaned_df = clean_data(
        df.copy(),
        file_path="path/to/proteindata.csv",
        metadata=metadata,
        outPath=str(tmp_path),
        config=config,
        json_out=str(json_out),
    )

    # ────────────────────────────────────────────────────
    # 7) Assertions
    # ────────────────────────────────────────────────────

    # b) Columns renamed to sample_rep: ["S1_1", "S2_1"]
    assert set(cleaned_df.columns) == {"S1_1", "S2_1"}

    # c) Zero→NaN and "x"→NaN, then rows with any NaN are dropped. Duplicate also drop
    #    nan converted to unknown-gene-1. that row remains along with p4
    assert set(cleaned_df.index.tolist()) == {"p4", "Unknown-gene-1"}

    # d) duplicate rows should be removed
    assert cleaned_df.index.is_unique


def write_csv(path, df):
    df.to_csv(path, index=False)


# ────────────────────────────────────────────────────
# Test full clean function with metadata
# ────────────────────────────────────────────────────
def test_process_data_metadata_branch( tmp_path ):
    """Test the process_data function from data_processing.py with metadata branch."""

    # 1) Create a tiny metadata CSV
    long_cols = [
        "sample_1_protein_abundance_intensity",
        "sample_2_protein_abundance_intensity",
    ]
    metadata = pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "treatment": ["Ctl", "Drug"],
            "replicate": [1, 1],
            "protein_abundance_name": long_cols,
        }
    )
    # sample_rep will be ["S1_1","S2_1"] when clean_meta runs
    # **Manually create the sample_rep column** (just like clean_meta would)
    metadata["sample_rep"] = (
        metadata["sample_id"] + "_" + metadata["replicate"].astype(str)
    )
    meta_path = Path(tmp_path) / "my_metadata.csv"
    metadata.to_csv(meta_path)

    # 2) Create a JSON output placeholder
    json_out = Path(tmp_path) / "meta_out.json"
    json_out.write_text("{}")

    # 3) Call process_data on metadata
    result = process_data(str(meta_path), json_out=str(json_out), config={})

    # 4) Should return a DataFrame with the same columns + your clean_meta additions
    assert isinstance(result, pd.DataFrame)
    assert "sample_rep" in result.columns
    # And the JSON file should have been touched by clean_meta
    with json_out.open() as f:
        summary = json.load(f)
    assert "NUM_SAMPLES" in summary
    assert summary["NUM_SAMPLES"] == 2
    assert summary["NUM_TREATMENTS"] == 2
    assert summary["TREATMENTS"] == "Ctl: 1, Drug: 1"


# ────────────────────────────────────────────────────
# Test process protein data
# ────────────────────────────────────────────────────
def test_process_data_proteindata_branch(monkeypatch, tmp_path):
    """Test the process_data function from data_processing.py with proteindata branch."""

    monkeypatch.setattr(dpp, "view_prot_distributions", lambda *a, **k: None)

    # 1) Create a tiny proteindata CSV (two columns)
    long_cols = [
        "sample_1_protein_abundance_intensity",
        "sample_2_protein_abundance_intensity",
    ]
    prot_df = pd.DataFrame(
        {
            long_cols[0]: [1, 0, "x", 6, 6],
            long_cols[1]: [4, 5, 6, 6, 6],
        },
        index=[np.nan, "p2", "p3", "p4", "p4"],
    )

    prot_path = Path(tmp_path) / "my_proteindata.csv"
    prot_df.to_csv(prot_path)

    # 2) Build matching metadata (so clean_data won’t error)
    metadata = pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "treatment": ["Ctl", "Drug"],
            "replicate": [1, 1],
            "protein_abundance_name": long_cols,
        }
    )
    # sample_rep will be ["S1_1","S2_1"] when clean_meta runs
    # **Manually create the sample_rep column** (just like clean_meta would)
    metadata["sample_rep"] = (
        metadata["sample_id"] + "_" + metadata["replicate"].astype(str)
    )

    # 2. Create a starter JSON file with one existing key
    json_out = Path(tmp_path) / "prot_summary.json"
    initial = {"FOO": 123}
    json_out.write_text(json.dumps(initial))

    config = {"df_to_use": "df_imp",
              "data_type": "prot",
              "missing_threshold": 0.5,
              "normalise_method": "sample-median",
              "imputation_method": "hist_grad_boost"}

    # 3) Call process_data on proteindata
    result = process_data(
        str(prot_path),
        metadata=metadata,
        json_out=json_out,
        outPath=str(tmp_path),
        config=config,
    )

    # 4) Should return a DataFrame (the cleaned & renamed data)
    assert isinstance(result, pd.DataFrame)
    # Columns must now be your sample_rep names
    assert set(result.columns) == {"S1_1", "S2_1"}
    ### two prots should be left
    assert ("Unknown-gene-1" in result.index) and ("p4" in result.index)
