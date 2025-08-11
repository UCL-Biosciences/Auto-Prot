import json
import os
import tempfile

import pandas as pd
import pytest

from autoprot.reporting.generate_report import generate_report_html


@pytest.fixture
def mock_report_setup():
    with tempfile.TemporaryDirectory() as temp_dir:
        # Simulate minimal config
        config = {
            "outPath": os.path.join(temp_dir, "output"),
            "json_outPath": os.path.join(temp_dir, "metadata.json"),
        }
        os.makedirs(config["outPath"], exist_ok=True)

        # Create dummy report template with placeholders
        template_md = os.path.join(temp_dir, "template.md")
        with open(template_md, "w") as f:
            f.write(
                "# Report\n\n"
                "$PROJECT_NAME\n\n"
                "{{TOP_LFC_PROTS}}\n\n"
                "{{ENRICHMENT_DF}}\n\n"
                "{{ENRICHMENT_PLOT}}\n"
            )

        # Create dummy top LFC CSV
        top_lfc_csv = os.path.join(temp_dir, "top_lfc.csv")
        pd.DataFrame({"gene": ["A", "B"], "logFC": [1.2, -0.8]}).to_csv(
            top_lfc_csv, index=False
        )

        # Create dummy enrichment CSV
        enrich_csv = os.path.join(temp_dir, "enrich.csv")
        pd.DataFrame(
            {
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
            }
        ).to_csv(enrich_csv, index=False)

        # Create dummy image
        enrich_img = os.path.join(temp_dir, "plot.png")
        with open(enrich_img, "wb") as f:
            f.write(
                b"\x89PNG\r\n\x1a\n" + b"0" * 100
            )  # minimal PNG header + dummy content

        # Create JSON metadata with a PROJECT_NAME key
        with open(config["json_outPath"], "w") as f:
            json.dump({"PROJECT_NAME": "Test Project"}, f)

        # HTML output path
        report_html = os.path.join(config["outPath"], "report.html")

        yield {
            "report_MD": template_md,
            "report_html": report_html,
            "top_LFC_prots_path": top_lfc_csv,
            "enrichment_path": enrich_csv,
            "enrichment_plot_path": enrich_img,
            "json_out": config["json_outPath"],
            "config": config,
        }


def test_generate_report_html_generates_file(mock_report_setup):
    args = mock_report_setup
    generate_report_html(**args)

    # Check HTML file exists
    assert os.path.exists(args["report_html"])
    with open(args["report_html"]) as f:
        html = f.read()
        assert "Test Project" in html
        assert "test pathway" in html
        assert "<img" in html or "No pathway enrichment plot" in html

    # Check that metadata JSON has git hashes recorded
    with open(args["json_out"]) as f:
        meta = json.load(f)
        assert "GENERATE_REPORT_VERSION" in meta
        assert "TEMPLATE_VERSION" in meta


# def test_generate_report_html_patch(tmp_path):
#      # Fake the 'markdown' module before anything imports it
#     sys.modules["markdown"] = types.SimpleNamespace(markdown=lambda x: f"<p>{x}</p>")

#     # Now you can import safely
#     from autoprot.reporting.generate_report import generate_report_html

#     # Set up dummy paths
#     report_md_path = tmp_path / "template.md"
#     report_html_path = tmp_path / "report.html"
#     top_lfc_path = tmp_path / "toplfc.csv"
#     enrichment_path = tmp_path / "enrich.csv"
#     enrichment_plot_path = tmp_path / "enrich.png"
#     json_out_path = tmp_path / "data.json"

#     # Write template
#     report_md_path.write_text("""
#     <html>
#     <body>
#     <h1>Report</h1>
#     $GENERATE_REPORT_VERSION
#     {{TOP_LFC_PROTS}}
#     {{ENRICHMENT_DF}}
#     {{ENRICHMENT_PLOT}}
#     </body>
#     </html>
#     """)

#     # Write data
#     pd.DataFrame({"Protein": ["P1", "P2"], "logFC": [2.1, -1.3]}).to_csv(top_lfc_path, index=False)
#     pd.DataFrame({
#         "source": ["GO"], "native": ["GO:00001"], "name": ["term1"], "p_value": [0.01],
#         "term_size": [50], "query_size": [10], "intersection_size": [5],
#         "precision": [0.5], "recall": [0.6], "treatment_pair": ["A_vs_B"]
#     }).to_csv(enrichment_path, index=False)
#     enrichment_plot_path.write_bytes(b"fake image")
#     json_out_path.write_text("{}")

#     # Patch markdown and git
#     with mock.patch("subprocess.check_output", return_value=b"dummyhash"), \
#          mock.patch("markdown.markdown", return_value="<p>Rendered HTML</p>"):

#         generate_report_html(
#             str(report_md_path),
#             str(report_html_path),
#             str(top_lfc_path),
#             str(enrichment_path),
#             str(enrichment_plot_path),
#             str(json_out_path)
#         )

#     # Check HTML output
#     assert report_html_path.exists()
#     html = report_html_path.read_text()
#     assert "Rendered HTML" in html
#     assert "dummyhash" not in html  # replaced before markdown conversion
#     assert "<p>" in html
