import pandas as pd
import tempfile
import os
import json
import yaml
import pytest
from pathlib import Path
from src.analysis.analysis import run_analysis

@pytest.fixture
def minimal_input():
    df = pd.DataFrame({
        'sample1': [1.0, 2.0, 3.0],
        'sample2': [2.0, 1.0, 3.0],
        'sample3': [3.0, 3.0, 1.0],
        'sample4': [4.0, 4.0, 4.0],
    }, index=['P1', 'P2', 'P3'])

    metadata = pd.DataFrame({
        'sample_rep': ['sample1', 'sample2', 'sample3', 'sample4'],
        'treatment': ['A', 'A', 'B', 'B'],
        'colours' : ['red', 'red', 'blue', 'blue']
    })

    config = {
        "LFC_threshold": 1.0,
        "FDR_threshold": 0.05,
        "LFC_plot_p_or_FDRp": "Log10_FDR_P_Value"
    }

    return df, metadata, config

def test_run_analysis_end_to_end(minimal_input):
    df, metadata, config = minimal_input

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        json_out = tmpdir_path / "metadata.json"

        os.makedirs(tmpdir_path / "plots" )
        os.makedirs(tmpdir_path / "data" )

        # Start with empty metadata file
        with open(json_out, "w") as f:
            json.dump({}, f)

        results = run_analysis(
            df=df,
            metadata=metadata,
            output_dir=str(tmpdir_path),
            config=config,
            json_out=str(json_out),
        )

        # Basic results structure
        assert "pca" in results
        assert "mds" in results
        assert "heatmap" in results
        assert any(k.startswith("df_lm_") for k in results)

        # Output CSVs
        assert (tmpdir_path / "data/combined_topLFC.csv").exists()

        # Per-pair volcano plots
        volcano_files = list(tmpdir_path.rglob("*volcano_plot.png"))
        assert len(volcano_files) > 0

        # Combined plots
        assert (tmpdir_path / "plots/combined_volcano_plot.png").exists()

        # Version metadata
        with open(json_out) as f:
            meta = yaml.safe_load(f)
        assert "ANALYSIS_VERSION" in meta
