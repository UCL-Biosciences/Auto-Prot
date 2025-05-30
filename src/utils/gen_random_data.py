import sys
from pathlib import Path

# Add parent of `src/` to sys.path to import code.utils.check_env
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import os
import pandas as pd
import numpy as np
import random

from src.utils.check_env import get_repo_root

def generate_protein_data(
    repo_root,
    human_genes_file="input/data/human_genes.txt",
    n_treatments=2,
    n_genes=1000,
    samples=None,
    n_replicates=5,
    spike_multiplier=4,
    na_fraction=0.05,
    random_seed=0,
    output_prefix="input/data/"
):
    
    random.seed(random_seed)
    np.random.seed(random_seed)

    if samples is None:
        samples = list("ABCDEF")

    # Prepare sample replicate columns
    sample_replicates = [
        f"{sample}_{replicate}"
        for sample in samples
        for replicate in range(1, n_replicates + 1)
    ]
    col_names = [f"prot_abundance_{sr}" for sr in sample_replicates]
    n_total_samp_rep = len(samples) * n_replicates
    n_samples_per_treatment = int(n_total_samp_rep / n_treatments)

    treatments = ["positive"] * n_samples_per_treatment + ["negative"] * n_samples_per_treatment

    # Read and sample genes
    genes = pd.read_csv(os.path.join(repo_root, human_genes_file))["Symbol"].dropna().unique().tolist()
    sampled_genes = random.sample(genes, n_genes)

    # Generate data
    gene_means = np.random.uniform(3, 5, size=n_genes)
    gene_sigmas = np.random.uniform(0.3, 1.5, size=n_genes)
    data = np.array([
        np.random.lognormal(mean=mu, sigma=sigma, size=n_total_samp_rep)
        for mu, sigma in zip(gene_means, gene_sigmas)
    ])
    df = pd.DataFrame(data, index=sampled_genes, columns=col_names)
    df.index.name = "PG.genes"
    # Spiking
    n_spiked = int(0.05 * n_genes)
    spike_set_1 = np.random.choice(sampled_genes, size=n_spiked, replace=False)
    remaining_genes = list(set(sampled_genes) - set(spike_set_1))
    spike_set_2 = np.random.choice(remaining_genes, size=n_spiked, replace=False)

    df.loc[spike_set_1, col_names[:n_samples_per_treatment]] *= spike_multiplier
    df.loc[spike_set_2, col_names[n_samples_per_treatment:]] *= spike_multiplier

    # Add missing values
    mask = np.random.rand(*df.shape) < na_fraction
    df = df.mask(mask)

    # Add PTM columns
    n_rows = df.shape[0]
    
    # PTM.ModificationTitle: 90% "Phospho (STY)", 10% "Carbamidomethyl (C)"
    mod_titles = np.random.choice(
        ["Phospho (STY)", "Carbamidomethyl (C)"],
        size=n_rows,
        p=[0.9, 0.1]
    )

    # Add to DataFrame
    df["PTM.ModificationTitle"] = mod_titles
    df["PTM.SiteAA"] = ["S"] * n_rows
    df["PTM.SiteLocation"] =  np.random.randint(1, 1001, size=n_rows)



    # Save protein data
    os.makedirs(os.path.join(repo_root, output_prefix), exist_ok=True)
    df.to_csv(os.path.join(repo_root, output_prefix, "proteindata.csv"))

    # Save metadata
    metadata_df = pd.DataFrame({
        "sample_id": [s for s in samples for _ in range(n_replicates)],
        "replicate": list(range(1, n_replicates + 1)) * len(samples),
        "treatment": treatments,
        "protein_abundance_name": col_names,
    })
    metadata_df.to_csv(os.path.join(repo_root, output_prefix, "metadata.csv"), index=False)

    print("Dummy protein data and metadata saved.")

if __name__ == "__main__":
    from pathlib import Path

    # Example: get the repo root as the parent of this script
    REPO_ROOT = get_repo_root()  # adjust as needed
    generate_protein_data(repo_root=REPO_ROOT)
