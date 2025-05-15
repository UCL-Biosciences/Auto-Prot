import sys
import types
import json
import pandas as pd
from pathlib import Path
from unittest import mock


def test_generate_report_html_patch(tmp_path):
     # Fake the 'markdown' module before anything imports it
    sys.modules["markdown"] = types.SimpleNamespace(markdown=lambda x: f"<p>{x}</p>")
    
    # Now you can import safely
    from src.reporting.generate_report import generate_report_html
    
    # Set up dummy paths
    report_md_path = tmp_path / "template.md"
    report_html_path = tmp_path / "report.html"
    top_lfc_path = tmp_path / "toplfc.csv"
    enrichment_path = tmp_path / "enrich.csv"
    enrichment_plot_path = tmp_path / "enrich.png"
    json_out_path = tmp_path / "data.json"

    # Write template
    report_md_path.write_text("""
    <html>
    <body>
    <h1>Report</h1>
    $GENERATE_REPORT_VERSION
    {{TOP_LFC_PROTS}}
    {{ENRICHMENT_DF}}
    {{ENRICHMENT_PLOT}}
    </body>
    </html>
    """)

    # Write data
    pd.DataFrame({"Protein": ["P1", "P2"], "logFC": [2.1, -1.3]}).to_csv(top_lfc_path, index=False)
    pd.DataFrame({
        "source": ["GO"], "native": ["GO:00001"], "name": ["term1"], "p_value": [0.01],
        "term_size": [50], "query_size": [10], "intersection_size": [5],
        "precision": [0.5], "recall": [0.6], "treatment_pair": ["A_vs_B"]
    }).to_csv(enrichment_path, index=False)
    enrichment_plot_path.write_bytes(b"fake image")
    json_out_path.write_text("{}")

    # Patch markdown and git
    with mock.patch("subprocess.check_output", return_value=b"dummyhash"), \
         mock.patch("markdown.markdown", return_value="<p>Rendered HTML</p>"):
        
        generate_report_html(
            str(report_md_path),
            str(report_html_path),
            str(top_lfc_path),
            str(enrichment_path),
            str(enrichment_plot_path),
            str(json_out_path)
        )

    # Check HTML output
    assert report_html_path.exists()
    html = report_html_path.read_text()
    assert "Rendered HTML" in html
    assert "dummyhash" not in html  # replaced before markdown conversion
    assert "<p>" in html
