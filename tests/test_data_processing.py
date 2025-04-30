import sys
from pathlib import Path

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))

import json
import pandas as pd
import numpy as np
import pytest

import utils.data_processing as dp

from utils.data_processing import clean_meta, clean_prot, prot_summary, clean_data

# ————————————————————————————————————————————
# Test clean_meta()
# ————————————————————————————————————————————

def test_clean_meta_creates_expected_columns_and_json(tmp_path):
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

def test_clean_prot_filters_and_renames_realistic():
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

def test_prot_summary_appends_stats_to_json(tmp_path):
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
def test_clean_data_raises_without_metadata():
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
def test_clean_data_selects_and_fixes_index(monkeypatch):
    """
    For proteindata branch, it should:
      - call clean_prot (we stub it to echo back df and count)
      - call process_prot_data (we stub it to return a dict)
      - select dfs[df_to_use]
      - replace NaN in index with 'Unknown-gene-1', drop duplicate labels
    """
    # 1. Input: 3-row DF; actual values don’t matter
    df_in = pd.DataFrame({"anything": [1, 2, 3]})

    # 2. Stub clean_meta so metadata branch never runs
    monkeypatch.setattr(
        "utils.data_processing.clean_meta",
        lambda df, json_out: df,
    )

    # 3. Stub clean_prot to return (df_in, original_count)
    monkeypatch.setattr(
        "utils.data_processing.clean_prot",
        lambda df, md: (df, len(df)),
    )

    # 4. Stub dpp.process_prot_data to return a dict with one DataFrame
    import utils.data_processing as dp
    fake_df = pd.DataFrame({"x": [9, 8, 7]}, index=[np.nan, "G1", "G1"])
    monkeypatch.setattr(
        dp.dpp,
        "process_prot_data",
        lambda df, md, cfg: {"keep_me": fake_df},
    )

    # 5. Stub out plotting and summary to no-ops
    monkeypatch.setattr(dp.dpp, "view_prot_distributions", lambda *a, **k: None)
    monkeypatch.setattr(dp, "prot_summary", lambda *a, **k: None)

    # 6. Call clean_data asking for "keep_me"
    result_df, nrow_original = clean_data(
        df_in,
        file_path="some/proteindata/file.csv",
        metadata=pd.DataFrame({"protein_abundance_name": ["sample_1_1", "sample_2_1"]}),        # not used by our stub
        outPath=".",
        config={"df_to_use": "keep_me"},
        json_out="out.json",
    )

    # 7. Check we got back the right count and DataFrame
    assert nrow_original == 3
    # Index should be ["Unknown-gene-1", "G1"]
    assert result_df.index.tolist() == ["Unknown-gene-1", "G1"]
    
    
# ────────────────────────────────────────────────────
# Test clean_data proteindata error paths
# ────────────────────────────────────────────────────
def test_clean_data_proteindata_errors():
    df = pd.DataFrame({"x":[1]})
    with pytest.raises(ValueError, match="Metadata is required"):
        clean_data(df, file_path="proteindata.csv", metadata=None, outPath=None, config={}, json_out=None)
    bad_md = pd.DataFrame({"foo":[1]}) ## missing protein_abundance_name in the metadata
    with pytest.raises(ValueError, match="missing in the metadata"):
        clean_data(df, file_path="proteindata.csv", metadata=bad_md, outPath=None, config={}, json_out=None)