###### Generate Random Protein Abundance Dataset ######

### for a given number of samples and treatments
### gene names downloaded from NCBI and are randomly assigned

##### load libraries
import random
import pandas as pd
import os
import numpy as np
from check_env import get_repo_root 

#### set number parameters ####
n_treatments = 2
n_genes = 1000
samples = list('ABCDEFGHIJ')
n_replicates = 2

## good solution for combining. Print sample_replicate and loop through samples and 1:max number of replicates (exclusive)
sample_replicates = [f"{sample}_{replicate}" for sample in samples for replicate in range(1, n_replicates + 1)]
col_names = list(map(lambda x:'prot_abundance_'+str(x),sample_replicates)) # col names have a string before as we expect real data will

### total number of samples = samples * n_replicates
n_total_samp_rep = len(samples * n_replicates)
n_samples_per_treatment = int(n_total_samp_rep/n_treatments)

# Create the treatment list
treatments = ['positive'] * n_samples_per_treatment + ['negative'] * n_samples_per_treatment

##### read and randomly sample gene list
REPO_ROOT=get_repo_root()

human_genes_file = os.path.join(REPO_ROOT, 'input/data/human_genes.txt')
genes = pd.read_csv(human_genes_file)['Symbol'].dropna().unique().tolist()
random.seed(0) # for repeatability
sampled_genes = random.sample(genes, n_genes)

##### generate df with random numbers (set seed for repeatability)
np.random.seed(0) ; df =  pd.DataFrame(np.random.randint(10,10000, 
                                                         size=(len(sampled_genes),
                                                               len(sample_replicates))),
                                                               index=sampled_genes,
                                                               columns=col_names )

### we are going to randomly change some numbers to create some differences to test pathway enrichment
# Define the number of spiked proteins (5% of total genes)
n_spiked = int(0.05 * len(sampled_genes))

# Randomly select two mutually exclusive sets of genes for spiking
spike_set_1 = np.random.choice(sampled_genes, size=n_spiked, replace=False)
remaining_genes = list(set(sampled_genes) - set(spike_set_1))
spike_set_2 = np.random.choice(remaining_genes, size=n_spiked, replace=False)

#### we will spike gene set 1 in treatment group 1 and gene set 2 in treatment group 2
cols_for_set_1 = list(map(lambda x:'prot_abundance_'+str(x),sample_replicates[:n_samples_per_treatment]))
cols_for_set_2 = list(map(lambda x:'prot_abundance_'+str(x),sample_replicates[n_samples_per_treatment:]))

# Apply spiking
df.loc[spike_set_1, cols_for_set_1] *= 5  # Spike first set in first treatment group
df.loc[spike_set_2, cols_for_set_2] *= 5  # Spike second set in second treatment group

# Convert to integers
df = df.astype(int)

# save to file
df.to_csv(os.path.join(REPO_ROOT, 'input/data/proteindata.csv'), index=False)

###### create metadata
metadata_df = pd.DataFrame({
    "sample_id" : [item for item in samples for _ in range(n_replicates)], # sample ID is each sample (for item in samples) as many times as there are replicates ( for _ in range(n_replicates) )
    "replicate": list(range(1, n_replicates + 1 )) * len(samples), # for each sample, replicate will be 1:n_replicates
    "treatment": treatments,
    "protein_abundance_name": col_names
})

metadata_df.to_csv(os.path.join(REPO_ROOT, 'input/data/metadata.csv'), index=False)

print(f"Generated protein data and metadata and saved to file. All done.")