import sys
from pathlib import Path

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))

import json
import pandas as pd
import numpy as np
import pytest

import utils.data_processing as dp

from utils.data_processing import clean_meta, clean_prot, prot_summary, clean_data, process_data

# ————————————————————————————————————————————
# Test clean_meta()
# ————————————————————————————————————————————

def test_clean_meta(tmp_path):
    # 1. Build a tiny metadata DataFrame
    df = pd.DataFrame({
        "sample_id": ["S1", "S2", "S1"],
        "replicate": [1, 1, 2],
        "protein_abundance_name": ["Prot A", "Prot B", "Prot C"],
        "treatment": ["Ctl", "Drug", "Ctl"],
    })
    json_out = tmp_path / "meta_out.json"

    # 2. Call clean_meta
    result = clean_meta(df.copy(), str(json_out))

    # 3. Check new columns
    assert "protein_abundance_name" in result.columns
    assert list(result["protein_abundance_name"].unique()) == ["prot_a", "prot_b", "prot_c"]
    assert "sample_rep" in result.columns
    # sample_rep should be sample_id_replicate
    assert set(result["sample_rep"]) == {"S1_1", "S2_1", "S1_2"}
    assert "colours" in result.columns
    # There should be exactly as many unique colours as treatments
    assert result["colours"].nunique() == result["treatment"].nunique()

    # 4. JSON file should exist and contain summary keys
    data = json.loads(json_out.read_text())
    assert data["NUM_SAMPLES"] == 3
    assert data["NUM_TREATMENTS"] == 2
    # TREATMENTS string should mention both Ctl and Drug
    assert "Ctl:" in data["TREATMENTS"] and "Drug:" in data["TREATMENTS"]


# ————————————————————————————————————————————
# Test clean_prot()
# ————————————————————————————————————————————

def test_clean_prot():
    # 1. Create a prot DataFrame with two realistic protein‐abundance columns plus one extra
    df = pd.DataFrame({
        "sample_1_protein_measured_etc": [10, "x", 30],
        "sample_2_protein_measured_etc": [40, 50, 60],
        "unrelated_column": [7, 8, 9]
    }, index=["protA", "protB", "protC"])

    # 2. Build corresponding metadata:
    #    - protein_abundance_name matches the long column names
    #    - sample_rep gives the short names we want in the end
    metadata = pd.DataFrame({
        "protein_abundance_name": [
            "sample_1_protein_measured_etc",
            "sample_2_protein_measured_etc"
        ],
        "sample_rep": ["sample_1_1", "sample_2_1"]
    })

    # 3. Call clean_prot
    cleaned_df, nrow_original = clean_prot(df.copy(), metadata)

    # 4. It should drop the unrelated column
    assert "unrelated_column" not in cleaned_df.columns

    # 5. Non‐numeric in the first column becomes NaN
    assert np.isnan(cleaned_df.loc["protB", "sample_1_1"])

    # 6. Columns should be renamed to the short sample_rep names
    assert set(cleaned_df.columns) == {"sample_1_1", "sample_2_1"}

    # 7. nrow_original should equal the original number of proteins
    assert nrow_original == 3

    # 8. Protein IDs (index) remain unchanged
    assert list(cleaned_df.index) == ["protA", "protB", "protC"]


# ————————————————————————————————————————————
# Test prot_summary()
# ————————————————————————————————————————————

def test_prot_summary(tmp_path):
    # 1. Prepare a small DataFrame of processed protein abundances:
    #    Two proteins (rows), two samples (columns)
    df = pd.DataFrame({
        "sample_1_1": [1.0, 2.0],
        "sample_2_1": [3.0, 4.0]
    }, index=["p1", "p2"])
    nrow_original = 3  # pretend one protein was removed earlier

    # 2. Create a starter JSON file with one existing key
    json_out = tmp_path / "prot_summary.json"
    initial = {"FOO": 123}
    json_out.write_text(json.dumps(initial))

    # 3. Call prot_summary
    prot_summary(df, nrow_original, str(json_out))

    # 4. Read the JSON back in
    data = json.loads(json_out.read_text())

    # --- Check unchanged keys ---
    assert data["FOO"] == 123

    # --- Check the new summary stats ---
    # Original count
    assert data["NUM_PROTS_OG"] == 3
    # Removed count = original - remaining (3 - 2 = 1)
    assert data["NUM_PROTS_REMOVED"] == 1

    # NUM_PROTS formatted as string with no decimals
    assert data["NUM_PROTS"] == "2"

    # Compute expected mean abundances:
    # For prot 1: (1 + 3)/2 = 2.0
    # For prot 2: (2 + 4)/2 = 3.0
    expected_means = [2.0, 3.0]

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
    # ────────────────────────────────────────────────────
    # 1) Skip plotting so tests stay headless
    # ────────────────────────────────────────────────────
    import utils.data_processing as dp
    monkeypatch.setattr(dp.dpp, "view_prot_distributions", lambda *a, **k: None)

    # ────────────────────────────────────────────────────
    # 2) Tiny metadata DataFrame matching long col names
    # ────────────────────────────────────────────────────
    long_cols = [
        "sample_1_protein_abundance_intensity",
        "sample_2_protein_abundance_intensity"
    ]
    metadata = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "treatment": ["Ctl", "Drug"],
        "replicate": [1,    1],
        "protein_abundance_name": long_cols,
    })
    # sample_rep will be ["S1_1","S2_1"] when clean_meta runs
    # **Manually create the sample_rep column** (just like clean_meta would)
    metadata["sample_rep"] = (
        metadata["sample_id"] + "_" + metadata["replicate"].astype(str)
    )

    # ────────────────────────────────────────────────────
    # 3) Starter JSON for prot_summary
    # ────────────────────────────────────────────────────
    json_out = tmp_path / "summary.json"
    json_out.write_text(json.dumps({}))

    # ────────────────────────────────────────────────────
    # 4) Build proteindata DataFrame with zeros & non-numeric
    # ────────────────────────────────────────────────────
    df = pd.DataFrame({
        long_cols[0]: [1,  0, "x", 6, 6],
        long_cols[1]: [4,  5, 6, 6, 6],
    }, index=[np.nan, "p2", "p3", "p4", "p4"])

    # ────────────────────────────────────────────────────
    # 5) Pick the “raw” DataFrame after process_prot_data
    # ────────────────────────────────────────────────────
    config = {"df_to_use": "df_imp"}

    # ────────────────────────────────────────────────────
    # 6) Run the pipeline
    # ────────────────────────────────────────────────────
    cleaned_df, nrow_original = clean_data(
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

    # a) Original count is 3 rows
    assert nrow_original == 5

    # b) Columns renamed to sample_rep: ["S1_1", "S2_1"]
    assert set(cleaned_df.columns) == {"S1_1", "S2_1"}

    # c) Zero→NaN and "x"→NaN, then rows with any NaN are dropped. Duplicate also drop
    #    nan converted to unknown-gene-1. that row remains along with p4
    assert cleaned_df.index.tolist() == ["Unknown-gene-1", "p4"]

    # d) duplicate rows should be removed
    assert cleaned_df.index.is_unique



def write_csv(path, df):
    df.to_csv(path, index=False)


# ────────────────────────────────────────────────────
# Test full clean function with metadata
# ────────────────────────────────────────────────────

def test_process_data_metadata_branch(monkeypatch, tmp_path):
    
    # 1) Create a tiny metadata CSV
    long_cols = [
        "sample_1_protein_abundance_intensity",
        "sample_2_protein_abundance_intensity"
    ]
    metadata = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "treatment": ["Ctl", "Drug"],
        "replicate": [1,    1],
        "protein_abundance_name": long_cols,
    })
    # sample_rep will be ["S1_1","S2_1"] when clean_meta runs
    # **Manually create the sample_rep column** (just like clean_meta would)
    metadata["sample_rep"] = (
        metadata["sample_id"] + "_" + metadata["replicate"].astype(str)
    )
    meta_path = tmp_path / "my_metadata.csv"
    write_csv(meta_path, metadata)

    # 2) Create a JSON output placeholder
    json_out = tmp_path / "meta_out.json"
    json_out.write_text("{}")

    # 3) Call process_data on metadata
    result = process_data(str(meta_path), json_out=str(json_out), config={})

    # 4) Should return a DataFrame with the same columns + your clean_meta additions
    assert isinstance(result, pd.DataFrame)
    assert "sample_rep" in result.columns
    # And the JSON file should have been touched by clean_meta
    summary = json.loads(json_out.read_text())
    assert "NUM_SAMPLES" in summary

# ────────────────────────────────────────────────────
# Test process protein data
# ────────────────────────────────────────────────────

def test_process_data_proteindata_branch(monkeypatch, tmp_path):

    monkeypatch.setattr(dp.dpp, "view_prot_distributions", lambda *a, **k: None)

    # 1) Create a tiny proteindata CSV (two columns)
    long_cols = [
        "sample_1_protein_abundance_intensity",
        "sample_2_protein_abundance_intensity"
    ]
    prot_df = pd.DataFrame({
        long_cols[0]: [1,  0, "x", 6, 6],
        long_cols[1]: [4,  5, 6, 6, 6],
    }, index=[np.nan, "p2", "p3", "p4", "p4"])

    prot_path = tmp_path / "my_proteindata.csv"
    write_csv(prot_path, prot_df)

    # 2) Build matching metadata (so clean_data won’t error)
    metadata = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "treatment": ["Ctl", "Drug"],
        "replicate": [1,    1],
        "protein_abundance_name": long_cols,
    })
    # sample_rep will be ["S1_1","S2_1"] when clean_meta runs
    # **Manually create the sample_rep column** (just like clean_meta would)
    metadata["sample_rep"] = (
        metadata["sample_id"] + "_" + metadata["replicate"].astype(str)
    )
    
    # 2. Create a starter JSON file with one existing key
    json_out = tmp_path / "prot_summary.json"
    initial = {"FOO": 123}
    json_out.write_text(json.dumps(initial))

    # 3) Call process_data on proteindata
    result = process_data(
        str(prot_path),
        metadata=metadata,
        json_out=json_out,
        outPath=str(tmp_path),
        config={"df_to_use": "df"}
    )

    # 4) Should return a DataFrame (the cleaned & renamed data)
    assert isinstance(result, pd.DataFrame)
    # Columns must now be your sample_rep names
    assert set(result.columns) == {"S1_1", "S2_1"}