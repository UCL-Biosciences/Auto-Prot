# Setup and Input Requirements

This page describes everything you need to get started: conda environments, input file formats, and where to put your data.

## Environments
Users do not need to install separately the many packages used by the pipeline. The packages are stored in a set of files in the `/configs` directory in config files e.g. `configs/auto-prot-env-windowsOS.yml`. Recreate with `conda env create --name --file=configs/auto-prot-env-windowsOS.yml` and activate with `conda activate auto-proteomics`. A corresponding file for mac OS is also in configs. You will need to create environments for the general pipeline (configs/auto-prot-env-windowsOS.yml), R functions (configs/auto-prot-env-limma-windowsOS.yml) and markdown functions (configs/auto-prot-env-markdown-windowsOS.yml).

For more information on conda environments, see links in the [conda homepage](https://docs.conda.io/en/latest).
 
This project includes an automatic Conda environment check before running the report generator. It ensures the active environment matches the environments/auto-prot-env.yml file. The check function is saved in `src/utils/check_env.py`.

If the environments differ, you will receive a warning before continuing. Note this doesn't happen in notebooks in VSC.

#### Conda Issues
I ran into a few common issues with the environment worth noting. They mostly come down to making sure that your computer knows where to look for the packages installed by conda. E.g.:
* make sure the conda executable is in $PATH.
* make sure the location where the environments are saved is in $PATH. It can be different to the location of conda installation
* make sure the python version running in any script is the correct one i.e. it matches the version in the conda env.

## Input File Format

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

---

### 2. Metadata (`metadata.csv`)

| sample_id | treatment | replicate | timepoint | abundance_col |
|-----------|-----------|-----------|-----------|----------------|
| Sample1   | control   | 1         | 1         | Sample1        |
| Sample2   | treated   | 1         | 1         | Sample2        |

- `Sample_ID`: matches a column in `proteindata.csv`
- `abundance_col`: must exactly match the corresponding column name in `proteindata.csv`
- `treatment`: required for group comparison and plotting
- `replicate`: integer; used in visualisation
- `timepoint`: optional, can be included in model for differential expression or used to subset samples

Additional columns can be included in the metadata. They will mostly be ignored unless you want to include them in the differential expression model or use them to subset the data (see section on config file for how to do this).

---

## Example Files

To generate example input data, `python src/utils/gen_random_data.py` will populate /input/data with valid `proteindata.csv` and `metadata.csv` files.

## Example Directory Layout
/data/
  ├── proteindata.csv
  └── metadata.csv

/configs/
  ├── auto-prot-env-windowsOS.yml
  ├── auto-prot-env-limma-windowsOS.yml
  └── auto-prot-env-markdown-windowsOS.yml

the same files are available for macOS users.

## Need help?
If your data isn’t working or you'd like to check formatting, feel free to get in touch: james.d.gilbert@ucl.ac.uk.
