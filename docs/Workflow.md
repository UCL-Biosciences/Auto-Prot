# Workflow Overview
Here we discuss the workflow in more detail, describing the scripts, their roles and some important details.

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


