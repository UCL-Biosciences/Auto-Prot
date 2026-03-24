import json
import os
import tempfile

import pandas as pd
import pytest
from PIL import Image

from autoprot.reporting.generate_report import generate_report_html
import autoprot.utils.check_env as ce


@pytest.fixture
def mock_report_setup():
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create directory structure matching what generate_report_html expects
        data_dir = os.path.join(temp_dir, "full_dataset/data")
        plots_dir = os.path.join(temp_dir, "full_dataset/plots")
        report_dir = os.path.join(temp_dir, "report")
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(report_dir, exist_ok=True)

        # Config matching what generate_report_html reads from
        config = {
            "outPath": temp_dir,
            "json_outPath": os.path.join(temp_dir, "metadata.json"),
            "missing_threshold": 0.5,
            "df_to_use": "df_imp",
            "normalise_method": "sample-median",
            "imputation_method": "hist_grad_boost",
            "FDR_threshold": 0.05,
        }

        # Report template at the path get_repo_root() + ./report/report-template.md
        # monkeypatch get_repo_root to return temp_dir so paths resolve correctly
        template_md = os.path.join(report_dir, "report-template.md")
        with open(template_md, "w") as f:
            f.write(
                "# Report\n\n"
                "$PROJECT_NAME\n\n"
                "{{TOP_LFC_PROTS}}\n\n"
                "{{ENRICHMENT_DF}}\n\n"
                "{{ENRICHMENT_PLOT}}\n"
            )

        # Top LFC CSV at expected path
        pd.DataFrame({"gene": ["A", "B"], "logFC": [1.2, -0.8]}).to_csv(
            os.path.join(data_dir, "combined_topLFC.csv"), index=False
        )

        # Enrichment CSV at expected path
        pd.DataFrame({
            "source": ["GO"],
            "native": ["GO:0001"],
            "name": ["test pathway"],
            "p_value": [0.01],
            "term_size": [10],
            "query_size": [5],
            "intersection_size": [2],
            "precision": [0.4],
            "recall": [0.2],
            "treatment_pair": ["A_vs_B"],
        }).to_csv(os.path.join(data_dir, "combined_top_pathway_enrichment.csv"), index=False)

        # Dummy enrichment plot
        img = Image.new("RGB", (100, 100), color=(255, 0, 0))
        img.save(os.path.join(plots_dir, "combined_pathway_enrichment_plot.png"))

        # Starter JSON
        with open(config["json_outPath"], "w") as f:
            json.dump({"PROJECT_NAME": "Test Project"}, f)

        yield config


def test_generate_report_html_generates_file(mock_report_setup, monkeypatch):
    config = mock_report_setup
    
    # make get_repo_root return temp_dir so internal path joins resolve correctly
    monkeypatch.setattr(ce, "get_repo_root", lambda: config["outPath"])
    
    generate_report_html(config=config)

    report_html = os.path.join(config["outPath"], "report-out.html")
    assert os.path.exists(report_html)
    
    with open(report_html) as f:
        html = f.read()
        assert "Protein Abundance Exploratory Analysis" in html
        assert "<img" in html or "No pathway enrichment plot" in html

    with open(config["json_outPath"]) as f:
        meta = json.load(f)
        assert meta["PROJECT_NAME"] == "Test Project"
        assert 'IMP_METHOD' in meta
        assert 'NORM_METHOD' in meta
