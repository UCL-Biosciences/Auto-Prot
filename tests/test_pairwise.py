import os
import tempfile
import pandas as pd
import numpy as np
from src.analysis.pairwise import make_volcano, enrichment_analysis

def test_make_volcano_integration_with_spikes():
    np.random.seed(42)

    # Setup: 30 proteins, 6 samples (3 per group)
    proteins = [f"Prot{i}" for i in range(30)]
    samples = [f"S{i}" for i in range(6)]

    # Create baseline expression (proteins = rows, samples = cols)
    expr_matrix = pd.DataFrame(
        np.random.normal(5, 0.5, (30, 6)),
        index=proteins,
        columns=samples
    )

    # Spike 3 proteins with large differences
    spike_proteins = ["Prot0", "Prot1", "Prot2"]
    expr_matrix.loc[spike_proteins, samples[:3]] += 2.5  # Up in group A
    expr_matrix.loc[spike_proteins, samples[3:]] -= 2.5  # Down in group B

    # DON'T Transpose for input (samples as rows, proteins as columns)
    df_pair = expr_matrix

    # Metadata: 3 A, 3 B
    metadata_pair = pd.DataFrame({
        "sample_id": samples,
        "treatment": ["A"] * 3 + ["B"] * 3
    })

    config = {
        "LFC_threshold": 1.0,
        "FDR_threshold": 0.05,
        "LFC_plot_p_or_FDRp": "Log10_FDR_P_Value"
    }

    formula = "~ treatment"
    pair_name = "A_vs_B"

    with tempfile.TemporaryDirectory() as tmpdir:
          # Ensure output directories exist
        os.makedirs(os.path.dirname(tmpdir), exist_ok=True)
        os.makedirs(os.path.join(tmpdir, "data", pair_name), exist_ok=True)
        os.makedirs(os.path.join(tmpdir, "plots", pair_name), exist_ok=True)
        result_df = make_volcano(
            df_pair=df_pair,
            output_dir=tmpdir,
            pair_name=pair_name,
            config=config,
            metadata_pair=metadata_pair,
            formula=formula
        )

        # Output paths
        volcano_path = os.path.join(tmpdir, "plots", pair_name, "volcano_plot.png")
        top20_path = os.path.join(tmpdir, "data", pair_name, "top_20_by_LFC.csv")
        limma_output = os.path.join(tmpdir, "data", pair_name, "limma_output.csv")

        # Check all outputs exist
        assert os.path.exists(volcano_path)
        assert os.path.exists(top20_path)
        assert os.path.exists(limma_output)

        # Check required columns
        assert all(col in result_df.columns for col in ["logFC", "adj.P.Val", "Colour"])

        # ✅ Confirm at least one spiked protein is called "blue"
        spiked = result_df.loc[spike_proteins]
        assert (spiked["Colour"] == "blue").all(), "Spiked proteins not detected as DE"
        not_spiked = result_df.loc[~result_df.index.isin(spike_proteins)]
        assert (not_spiked["Colour"] == "gray").all(), "Spiked proteins not detected as DE"

        # ✅ Confirm top-20 CSV is sorted by abs(logFC)
        top20_df = pd.read_csv(top20_path, index_col=0)
        abs_lfc = top20_df["logFC"].abs().values
        assert (abs_lfc == sorted(abs_lfc, reverse=True)).all(), "Top 20 not sorted by abs(logFC)"


def test_enrichment_analysis_creates_outputs_with_real_genes():
    np.random.seed(42)

    # ✅ Use real human genes known to be in KEGG pathways
    real_genes = [
        "TP53", "EGFR", "MTOR", "PIK3CA", "AKT1", "MAPK1", "KRAS",
        "PTEN", "RB1", "CDK4", "MYC", "ERBB2", "CCND1", "CDKN2A",
        "FOXO3", "BRAF", "MDM2", "CTNNB1", "SMAD4", "STAT3"
    ]

    # Simulate significant values
    df = pd.DataFrame(index=real_genes)
    df["logFC"] = np.random.uniform(1.5, 3.0, size=len(real_genes))  # large FC
    df["adj.P.Val"] = np.random.uniform(0.001, 0.01, size=len(real_genes))  # sig
    df["P.Value"] = df["adj.P.Val"] / 2

    # Add 10 noise genes that are non-significant
    noise_genes = [f"FAKEGENE{i}" for i in range(10)]
    noise_df = pd.DataFrame(index=noise_genes)
    noise_df["logFC"] = np.random.uniform(0, 0.5, size=10)
    noise_df["adj.P.Val"] = 0.9
    noise_df["P.Value"] = 0.8

    # Combine
    full_df = pd.concat([df, noise_df])

    config = {
        "LFC_threshold": 1.0,
        "FDR_threshold": 0.05,
        "species": "hsapiens"
    }

    pair_name = "A_vs_B"

    with tempfile.TemporaryDirectory() as tmpdir:
        os.makedirs(os.path.join(tmpdir, "data", pair_name), exist_ok=True)
        os.makedirs(os.path.join(tmpdir, "plots", pair_name), exist_ok=True)

        # Run real enrichment
        result = enrichment_analysis(full_df, pair_name, config, tmpdir)

        # ✅ Assertions
        assert isinstance(result, pd.DataFrame)
        assert not result.empty
        assert "name" in result.columns
        assert "p_value" in result.columns

        csv_path = os.path.join(tmpdir, "data", pair_name, f"{pair_name}_pathway_enrichment.csv")
        plot_path = os.path.join(tmpdir, "plots", pair_name, f"{pair_name}_pathway_enrichment_plot.png")

        assert os.path.exists(csv_path), "Expected enrichment CSV not created"
        assert os.path.exists(plot_path), "Expected enrichment plot not created"

        
# # ————————————————————————————————————————————
# #   test enrichment_analysis()
# # ————————————————————————————————————————————

# def test_enrichment_analysis_runs_and_saves_results(tmp_path):
    
#     pair_name = "gprofiler_test"
#     output_dir = tmp_path

#     # Make sure output folders exist
#     (tmp_path / "data" / pair_name).mkdir(parents=True, exist_ok=True)
#     (tmp_path / "plots" / pair_name).mkdir(parents=True, exist_ok=True)

#     # Example DE results with one real human gene
#     df = pd.DataFrame({
#         "logFC": [2.5, 0.5, -3.1],
#         "P.Value": [0.001, 0.2, 0.03],
#         "adj.P.Val": [0.01, 0.4, 0.025]
#     }, index=["TP53", "XYZ1", "MTOR__pS2448"])  # Include phospho-form to test cleanup

#     config = {
#         "LFC_threshold": 1.0,
#         "FDR_threshold": 0.05,
#         "species": "hsapiens"  # required by g:Profiler
#     }

#     # Act
#     result_df = enrichment_analysis(df, pair_name, config, str(tmp_path))

#     # Assert
#     # Returns DataFrame or empty DataFrame if nothing enriched
#     assert isinstance(result_df, pd.DataFrame)
#     if not result_df.empty:
#         assert "p_value" in result_df.columns
#         assert result_df["p_value"].max() <= 0.05

#         plot_path = tmp_path / "plots" / pair_name / f"{pair_name}_pathway_enrichment_plot.png"
#         csv_path = tmp_path / "data" / pair_name / f"{pair_name}_pathway_enrichment.csv"
#         assert plot_path.exists()
#         assert csv_path.exists()
