# Setup and Input Requirements

This page describes everything you need to get started: environments, input file formats, and where to put your data.

## Environments

### Main pipeline (Python/conda)
Users do not need to install separately the many packages used by the pipeline. The packages are stored in a set of files in the `/configs` directory in config files e.g. `configs/auto-prot-env-windowsOS.yml`. Recreate with `conda env create --name --file=configs/auto-prot-env-windowsOS.yml` and activate with `conda activate auto-proteomics`. A corresponding file for mac OS is also in configs. You will need to create environments for the general pipeline (configs/auto-prot-env-windowsOS.yml) and R functions (configs/auto-prot-env-limma-windowsOS.yml).

For more information on conda environments, see links in the [conda homepage](https://docs.conda.io/en/latest).
 
This project includes an automatic Conda environment check before running the report generator. It ensures the active environment matches the environments/auto-prot-env.yml file. The check function is saved in `autoprot/utils/check_env.py`.

If the environments differ, you will receive a warning before continuing. Note this doesn't happen in notebooks in VSC.

### R comparison script (renv, no conda)
The script `autoprot/r_scripts/DE-with-MSstats-Limma.R` is independent of conda and manages its own R packages via `renv`. See [compare-analysis-methods.md](compare-analysis-methods.md) for full setup instructions. In short: install R 4.5 or above from https://cran.r-project.org, add it to your PATH, then run the script directly with `Rscript`. All packages install automatically into the repo on first run — nothing touches your system R libraries.

#### Conda Issues
I ran into a few common issues with the environment worth noting. They mostly come down to making sure that your computer knows where to look for the packages installed by conda. E.g.:
* make sure the conda executable is in $PATH.
* make sure the location where the environments are saved is in $PATH. It can be different to the location of conda installation
* make sure the python version running in any script is the correct one i.e. it matches the version in the conda env.

## Input File Format
**Please please please never mess with raw data**. Everything in this pipeline runs automatically. You can always re-run it to recover outputs. If you have a single copy of your raw data and overwrite it, it might not be possible to retrieve the raw data! It would be good to have a copy of the original data safely stored (and backed up) in a location that will never be touched by this pipeline (or any other). E.g. for UCL people, store the raw data on the [Research Data Storage Service](https://www.ucl.ac.uk/isd/research-data-storage-service) and make a local copy as input for the pipeline.

This tool expects two input files in CSV format:
- `proteindata.csv`: protein abundance table
- `metadata.csv`: sample metadata table

You can place these in the `/data/` directory or specify a different path in the config file.

---

### 1. Protein Data (`proteindata.csv`)

| Protein_ID | Sample1 | Sample2 | Sample3 | ... |
|------------|---------|---------|---------|-----|
| P01234     | 1283    | 1092    | 1189    | ... |
| P05678     | 2049    | 1983    | 1877    | ... |

- **Rows** = proteins or genes (IDs can be UniProt, HGNC symbols, etc.)
- **Columns** = sample names
- All abundance values must be **numeric**
- Missing values should be represented as `NA` (not blank cells)
- Any column with 'gene' in the column name will be picked to name the rows

Needs only the raw protein intensity data for each sample, with column names matching the values in the `protein_abundance_name` in the metadata. If the genes associated with the proteins are available, include the gene names in a column with "genes" in the name, and this will be used throughout to label the proteins.

The name of the protein data by default is proteindata.csv in /data/input. You can rename the file, as long as you update the protPath field in the config (see config section). The protein data must have the word "protein" (lower case) in the file name as this determines the data cleaning steps, which are different for protein data and metadata.

Some protein data need filtering for the target species, quality etc. This isn't handled by the pipeline and should be done before running the pipeline.

It is useful to generate an index of protein or gene names so we can describe patterns for individual proteins. In order to be flexible, the pipeline can proceed without a specified column describing each protein/gene. By default it will use any column with "gene" in the name because this is what works with the data we have tried so far. If you don't have that, you could rename a column in the protein data to include "gene". Missing values in this column are replaced with "unknown-gene-X", with X being a unique number for each of the rows with missing gene name. If the names in the "gene" column are duplicated, the pipeline will break. We suggest you rename the genes (e.g. gene, gene to gene-1, gene-2) before running the pipeline.

#### Phosphoproteomic data
The pipeline can be used with phosphoproteomic data. There are a few things to be aware of before running the pipeline.

1. In phospho data, each protein can be included multiple times (in multiple rows), with each row including data for a different PTM. In this case, the gene name won't be unique and the pipeline will crash. The `data_type` field in the config should be set to "phospho" and `phoshpho_row_id` should include the names of columns used to create unique row IDs. The phospho_row_id follows the format below. "fields" tells us which columns to use to generate unique row IDs, "missing_value" is what to do if some rows have missing data in the "fields" columns.

```
"phospho_row_id" : {
        "fields": ["Amino acid", "Position", "Multiplicity"],
        "missing_value": "NA"   
    },
```
2. Phoshoproteomic data can include different types of post-translational modifications and we may not want to include them all in the analysis. It is up to the user to remove the unwanted modifications before running the pipeline (i.e. no filtering takes place within the pipeline).
---

### 2. Metadata (`metadata.csv`)

| sample_id | treatment | replicate | timepoint | protein_abundance_name |
|-----------|-----------|-----------|-----------|------------------------|
| Sample1   | control   | 1         | 1         | Sample1                |
| Sample2   | treated   | 1         | 1         | Sample2                |

- `Sample_ID`: used to identify samples
- `treatment`: treatment group for each sample. E.g. positive, negative, healthy etc. required for group comparison and plotting
- `replicate`: integer; replicates for sample IDs. If there are no replicates, enter 1 for all rows.
- `timepoint`: optional, can be included in model for differential expression or used to subset samples
- `protein_abundance_name`:  Very important column linking metadata to protein abundance. This must contain the exact name of the column containing protein data for each sample, replicate etc

The pipeline will run with just these four columns in metadata (timepoint is optional). However, sample_id + replicate combinations must be unique. If there is a time point variable to be included, it should be under in a column called `timepoint`. This will then be added to sample_id and replicate to generate a unique identifier.

Additional columns can be included in the metadata. They will mostly be ignored unless you want to include them in the differential expression model or use them to subset the data (see section on config file for how to do this).

---

## Example Files

To generate example input data, `python autoprot/utils/gen_random_data.py` will populate /input/data with valid `proteindata.csv` and `metadata.csv` files.

An example of the report generated is available in `docs/example-report-out.html`.

## Config file
We want this tool to be accessible for people with limited coding experience. Therefore, key parameters can be set without the need to change any code. These parameters are controlled from a config file in json format. By default: `configs/auto-prot-config.json`. The main functions in the code find this file and extract the parameters as required.

To make the code work, you must enter the correct parameters and combination of parameters. For example, if you change the input file but not the output file, it will overwrite any previous output in the output file.

### Run Settings

**process_prot_data.** Whether to complete the full protein data processing. If you have already done the data processing and have a `data/proteinAbundance.csv` file to use, set this to false to save a few seconds. If true (default), process protein data from scratch. If false, load from `output/data/proteinAbundance.csv`. 


### Input Data

**data_type.** Type of data. Options:
- "prot" → Proteomics (default)
- "phospho" → Phosphoproteomics

**phospho_row_id.fields.** Columns used to create unique row IDs if using phosphoproteomics. Ignore if `data_type = "prot"`. Example given in the config.yaml file.

**phospho_row_id.missing_value**. Value to use if any of the above fields are missing. Default: "NA"

**missing_threshold.** Proportion threshold for missingness per group (0–1). Example: 0.75 means proteins must appear in at least 75% of samples per group. Default: 0.75

### Normalisation and imputation

**normalise_method.** Normalisation method. Options: "vsn" (variance stabilisation, recommended), "sample-median"

**imputation_method.**: Imputation method. Options: "pimms_collabfilter", "hist_grad_boost"

**df_to_use.** Which dataset to analyse. Options: "df" → Raw (not recommended), "df_log2" → Log2 transformed, "df_norm" → Normalised, "df_imp" → Normalised + imputed. Default: df_imp

**z_score_for_clustering.** Whether to z-score data for PCA, MDS and heatmap clustering (true/false). Uses scikit-learn's [StandardScaler](https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html#sklearn.preprocessing.StandardScaler). 

### Annotation

**species.** Species for pathway enrichment. Examples: "hsapiens", "mmusculus". Default: "hsapiens"

#FILE PATHS

**protPath**: Path to protein data file. Example: input/data/proteindata.csv

**metaPath**: Path to metadata file. Example: input/data/metadata.csv

**outPath**: Directory for all outputs. Example: output

**json_outPath**: JSON file for storing report-related info. Example: output/data/data_for_report.json

#### Subsetting

**analyse_full_dataset.** Analyse full dataset i.e. all samples. true/false. Default: true

**analyse_subsets.** Analyse subsets based on a metadata variable. true/false. Default: false

**subset_variable.** Metadata column to subset on. Example: "timepoint"

**subsets.** Specific subset values to run. Empty = all values of subset_variable. Example: "timepoint 1" or ["treatment A", "treatment B"]

#### Differential Expression Calculation

**DE_full_formula.** Formula for DE on full dataset (R-style model formula). Example: "~ treatment"

**DE_subset_formula.** Formula for DE on subsets. Example: "~ treatment". Ignore if not running the analysis on any subsets.

#### Figures

**LFC_plot_p_or_FDRp.** Which p-value metric to use to define differential expression, including as y-axis on volcano plot. Options: "Log10_unadjusted_p_Value", "Log10_FDR_P_Value" (default). If there are no proteins meeting thresholds for LFC and FDR p-value, you might want to use an unadjsuted (i.e. less conservative) p-value. Be cautious.

**LFC_threshold.** Log2 fold-change threshold for highlighting points. Default: 1

**FDR_threshold.** FDR cut-off for colouring points. Default: 0.05

#### Pathway Enrichment

**enrichment_pathways.** Pathway database for enrichment. Options: KEGG, Reactome, GO. Default: KEGG

## Need help?
If your data isn’t working or you'd like to check formatting, feel free to get in touch. Contact details in main README.
