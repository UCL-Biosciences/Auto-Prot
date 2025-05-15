import tempfile
import os
import pandas as pd
from src.utils.gen_random_data import generate_protein_data  # replace with actual module name

def test_generate_dummy_protein_data_creates_valid_files():
    # Create a temporary directory to simulate the repo root
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a dummy gene list file
        gene_file_path = os.path.join(tmpdir, "input/data/human_genes.txt")
        os.makedirs(os.path.dirname(gene_file_path), exist_ok=True)
        pd.DataFrame({"Symbol": [f"GENE{i}" for i in range(200)]}).to_csv(gene_file_path, index=False)

        n_genes = 100
        samples = list('ABCDEF')
        n_replicates = 4

        # Run the data generator
        generate_protein_data(repo_root=tmpdir,
                              n_genes = n_genes,
                              samples = samples,
                              n_replicates = n_replicates)

        # Check the output files exist
        data_file = os.path.join(tmpdir, "input/data/proteindata.csv")
        metadata_file = os.path.join(tmpdir, "input/data/metadata.csv")
        assert os.path.exists(data_file)
        assert os.path.exists(metadata_file)

        # Load the files and test basic expectations
        df = pd.read_csv(data_file, index_col=0)
        meta = pd.read_csv(metadata_file)

        assert df.shape[0] == n_genes  # 100 genes
        assert "prot_abundance_A_1" in df.columns  # test for a column
        assert meta.shape[0] == len(samples) * n_replicates  # 6 samples × 5 replicates
        assert meta["protein_abundance_name"].nunique() == df.shape[1]
