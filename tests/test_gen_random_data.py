import os
import pandas as pd
import tempfile
import shutil
from src.utils.gen_random_data import generate_protein_data  # update path if needed

def test_generate_protein_data_creates_files():
    # Create temporary directory to simulate repo_root
    with tempfile.TemporaryDirectory() as temp_root:
        # Prepare dummy gene list
        gene_list_path = os.path.join(temp_root, "input", "data")
        os.makedirs(gene_list_path, exist_ok=True)
        gene_df = pd.DataFrame({"Symbol": [f"Gene{i}" for i in range(200)]})
        gene_df.to_csv(os.path.join(gene_list_path, "human_genes.txt"), index=False)

        # Run generator
        generate_protein_data(
            repo_root=temp_root,
            human_genes_file="input/data/human_genes.txt",
            n_genes=50,
            n_treatments=2,
            output_prefix="input/data/",
            random_seed=42,
        )

        # Check output files
        protein_file = os.path.join(temp_root, "input", "data", "proteindata.csv")
        metadata_file = os.path.join(temp_root, "input", "data", "metadata.csv")

        assert os.path.exists(protein_file)
        assert os.path.exists(metadata_file)

        # Read files and check basic structure
        proteindf = pd.read_csv(protein_file, index_col=0)
        metadf = pd.read_csv(metadata_file)

        assert "PTM.ModificationTitle" in proteindf.columns
        assert "PTM.SiteLocation" in proteindf.columns
        assert "PTM.SiteAA" in proteindf.columns

        assert "sample_id" in metadf.columns
        assert "replicate" in metadf.columns
        assert "protein_abundance_name" in metadf.columns
        assert len(proteindf) == 50