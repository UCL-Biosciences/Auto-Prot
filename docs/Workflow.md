# Workflow Overview
Here we discuss the workflow in more detail, describing the scripts, their roles and some important details.

## Input
Please please please never mess with raw data. Everything in this pipeline runs automatically. You can always re-run it to recover outputs. If you have a single copy of your raw data and overwrite it, it might not be possible to retrieve the raw data! It would be good to have a copy of the original data safely stored (and backed up) in a location that will never be touched by this pipeline (or any other). E.g. for UCL people, store the raw data on the [Research Data Storage Service](https://www.ucl.ac.uk/isd/research-data-storage-service) and make a local copy as input for the pipeline.

### Metadata
Requires four columns:
- sample_id. unique ID for each sample
- replicate. replicates for sample IDs. If there are no replicates, these values will all be 1.
- treatment. treatment group for each sample. E.g. positive, negative, healthy etc.
- protein_abundance_name. Very important column linking metadata to protein abundance. This must contain the exact name of the column containing protein data for each sample, replicate etc

The pipeline will run with just these four columns in metadata. However, sample_id + replicate combinations must be unique. If there is a time point variable to be included, it should be under in a column called `timepoint`. This will then be added to sample_id and replicate to generate a unique identifier.

### Protein data
Needs only the raw protein intensity data for each sample, with column names matching the values in the `protein_abundance_name` in the metadata. If the genes associated with the proteins are available, include the gene names in a column with "genes" in the name, and this will be used throughout to label the proteins.

The name of the protein data by default is proteindata.csv in /data/input. You can rename the file, as long as you update the protPath field in the config (see below). The protein data must have the word "protein" (lower case) in the file name as this determines the data cleaning steps, which are different for protein data and metadata.

Some protein data need filtering for the target species, quality etc. This isn't handled by the pipeline and should be done before running the pipeline.

It is useful to generate an index of protein or gene names so we can describe patterns for individual proteins. In order to be flexible, the pipeline can proceed without a specified column describing each protein/gene. By default it will use any column with "gene" in the name because this is what works with the data we have tried so far. If you don't have that, you could rename a column in the protein data to include "gene". Missing values in this column are replaced with "unknown-gene-X", with X being a unique number each of the rows with missing gene name.

Proteins need unique gene names. If multiple proteins have the same associated gene name (in a column with "gene" in the title), they will be renamed to Gene-1 and Gene-2 automatically.

#### Phosphoproteomic data
The pipeline can be used with phosphoproteomic data. There are a few things to be aware of before running the pipeline.

1. In phospho data, each protein can be included multiple times (in multiple rows), with each row including data for a different PTM. In this case, the gene name won't be unique and the pipeline will crash. The `data_type` field in the config should be set to "phospho" and `phoshpho_row_id` should include the names of columns used to create unique row IDs. The phospho_row_id follows the following format. "fields" tells us which columns to use to generate unique row IDs, "missing_value" is what to do if some rows have missing data in the "fields" columns.

```
"phospho_row_id" : {
        "fields": ["Amino acid", "Position", "Multiplicity"],
        "missing_value": "NA"   
    },
```


### Config file
We want this tool to be accessible for people with limited coding experience. Therefore, key parameters can be set without the need to change any code. These parameters are controlled from a config file in json format. By default: `configs/auto-prot-config.json`. The main functions find this file and extract the parameters as required.

To make the code work, you must enter the correct parameters and combination of parameters. For example, if you change the input file but not the output file, it will overwrite any previous output in the output file.

- data_type. "prot" = proteomics, this is the default use. "phospho" = phosphoproteomic data. This invokes some additional processing steps e.g. filtering to keep only phosphorylation of STY amino acids, other post-translational modifications are removed.
- phospho_row_id. If phosphorylation dataset, protein and gene names will be duplicated in multiple rows, with rows having different PTM info. This parameter determines which other columns will be used to create a unique ID for each row.
- missing_threshold. threshold for missingness (0-1). Protein must have been observed in more than this threshold per group. If 0.8, protein must be in 80% of samples per group or will be removed.
- normalise_method. "vsn" = variance stabilisation normalisation from the vsn R package, as recommended by Välikangas et al (2018). "sample-median" will normalise by subtracting each sample's median value from all of that sample's protein intensities.
- df_to_use. Which dataset to use in the pipeline. Options are raw data ("df_to_use" : "df" NOT REOMMENDED), log2 transformed data ("df_to_use" : "df_log2"), normalised data ("df_to_use" : "df_norm"), or normalised and imputed data ("df_to_use" : "df_imp").
- species. which species to use when looking up the annotation information for pathway enrichment analysis.
- protPath. path to protein data
- metaPath. path to metadata
- outPath. path to general output directory
- json_outPath. Path to a file (json format) to store some information that is referenced in the html report
- analyse_full_dataset. whether to analyse the full dataset, which means without filtering any proteins.
- analyse_subsets. whether to run the analysis on subsets of samples.
- subset_variable. which variable to subset the data on when running on subsets.
- subsets. if empty, will run on all unique values of the subset_variable. If you want to run on specific subsets, enter them here.

## Data Processing

### Subsets
You might wish to re-run the analysis on specific subsets of samples. For example, you might wish to select on samples from a given timepoint. You can choose whether to analyse any subsets by setting the config field `analyse_subsets` to `true`. You can define which subsets to analyse using the `subset_variable` (which variable to apply the subset too) and `subsets` fields. `subsets` needs a list of values to analyse e.g. ["1", "2"] if you want to analyse separately time points 1 and 2. By default, this will run the whole analysis for each subset. There is no default option to filter based on >1 variable.

## Analysis


## Creating the summary report
The report is generated in two main of steps:
1. `python main.py` creates the output
2. `src/reporting/generate_report.py` creates an html output file using the following:
     * a template for the report `./reporting/report-template.md`
     * outputs (tables, plots) from `main.py`
     * values in `output/data/data_for_report.json` that go directly into the markdown text
     * `src/reporting/generate_report.py` which converts the markdown template into an html file.


## Refs
Välikangas, Tommi, Tomi Suomi, and Laura L. Elo. "A systematic evaluation of normalization methods in quantitative label-free proteomics." Briefings in bioinformatics 19.1 (2018): 1-11.


