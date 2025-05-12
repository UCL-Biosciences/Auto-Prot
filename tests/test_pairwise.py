import pandas as pd
import numpy as np
import os
from pathlib import Path
from src.analysis.analysis import make_volcano, enrichment_analysis

# ————————————————————————————————————————————
#   test make_volcano()
# ————————————————————————————————————————————

def test_make_volcano(tmp_path):
    # Arrange
    pair_name = "real_test_pair"
    output_dir = tmp_path

    os.makedirs(os.path.join(output_dir, "data", pair_name), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "plots", pair_name), exist_ok=True)

    df_pair = pd.DataFrame({
        "ctrl_1_1": [10, 30],
        "ctrl_2_1": [11, 28],
        "treat_1_1": [50, 10],
        "treat_2_1": [52, 12],
    }, index=["ProteinA", "ProteinB"])

    metadata = pd.DataFrame({
        "sample_rep": ["ctrl_1_1", "ctrl_2_1", "treat_1_1", "treat_2_1"],
        "treatment": ["control", "control", "treatment", "treatment"]
    })

    config = {
        "LFC_threshold": 1.0,
        "FDR_threshold": 0.05,
        "LFC_plot_p_or_FDRp": "Log10_FDR_P_Value"
    }

    # Act
    diff_expr_df = make_volcano(
        df_pair=df_pair,
        output_dir=str(output_dir),
        pair_name=pair_name,
        config=config,
        metadata_pair=metadata,
    )

    # Assert
    assert isinstance(diff_expr_df, pd.DataFrame)
    assert "logFC" in diff_expr_df.columns
    assert "P.Value" in diff_expr_df.columns
    assert "adj.P.Val" in diff_expr_df.columns
    assert "Colour" in diff_expr_df.columns

    plot_path = output_dir / "plots" / pair_name / "volcano_plot.png"
    top_table_path = output_dir / "data" / pair_name / "top_20_by_LFC.csv"
    assert plot_path.exists()
    assert top_table_path.exists()


# ————————————————————————————————————————————
#   test make_volcano() identifies a protein that is DE
# ————————————————————————————————————————————

def test_make_volcano_identifies_clear_de_protein(tmp_path):
    # Arrange
    pair_name = "clear_de_case"
    df_pair = pd.DataFrame({
        "ctrl_1": [10, 20],
        "ctrl_2": [10, 20],
        "treat_1": [100, 20],
        "treat_2": [110, 20],
    }, index=["Protein_DE", "Protein_NDE"])  # only Protein_DE is changing

    metadata = pd.DataFrame({
        "sample_rep": ["ctrl_1", "ctrl_2", "treat_1", "treat_2"],
        "treatment": ["control", "control", "treatment", "treatment"]
    })

    config = {
        "LFC_threshold": 1.0,
        "FDR_threshold": 0.05,
        "LFC_plot_p_or_FDRp": "Log10_FDR_P_Value"
    }

    # Create the expected output path beforehand (as make_volcano will write to it)
    (tmp_path / "data" / pair_name).mkdir(parents=True, exist_ok=True)
    (tmp_path / "plots" / pair_name).mkdir(parents=True, exist_ok=True)

    # Run real analysis
    result_df = make_volcano(
        df_pair=df_pair,
        output_dir=str(tmp_path),
        pair_name=pair_name,
        config=config,
        metadata_pair=metadata
    )

    # Assert
    assert "Colour" in result_df.columns
    assert "Protein_DE" in result_df.index
    assert result_df.loc["Protein_DE", "Colour"] == "blue"
    assert result_df.loc["Protein_NDE", "Colour"] == "gray"

# ————————————————————————————————————————————
#   test enrichment_analysis()
# ————————————————————————————————————————————

def test_enrichment_analysis_runs_and_saves_results(tmp_path):
    
    pair_name = "gprofiler_test"
    output_dir = tmp_path

    # Make sure output folders exist
    (tmp_path / "data" / pair_name).mkdir(parents=True, exist_ok=True)
    (tmp_path / "plots" / pair_name).mkdir(parents=True, exist_ok=True)

    # Example DE results with one real human gene
    df = pd.DataFrame({
        "logFC": [2.5, 0.5, -3.1],
        "P.Value": [0.001, 0.2, 0.03],
        "adj.P.Val": [0.01, 0.4, 0.025]
    }, index=["TP53", "XYZ1", "MTOR__pS2448"])  # Include phospho-form to test cleanup

    config = {
        "LFC_threshold": 1.0,
        "FDR_threshold": 0.05,
        "species": "hsapiens"  # required by g:Profiler
    }

    # Act
    result_df = enrichment_analysis(df, pair_name, config, str(tmp_path))

    # Assert
    # Returns DataFrame or empty DataFrame if nothing enriched
    assert isinstance(result_df, pd.DataFrame)
    if not result_df.empty:
        assert "p_value" in result_df.columns
        assert result_df["p_value"].max() <= 0.05

        plot_path = tmp_path / "plots" / pair_name / f"{pair_name}_pathway_enrichment_plot.png"
        csv_path = tmp_path / "data" / pair_name / f"{pair_name}_pathway_enrichment.csv"
        assert plot_path.exists()
        assert csv_path.exists()
