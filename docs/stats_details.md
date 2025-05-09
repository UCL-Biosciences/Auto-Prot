# Stats and Calculations
A bit more detail on important calculations.

## Pre Processing
### Stats

Include details of packages used and others available

This tool is designed to produce a reasonable first exploration of your data. However, no single approach works optimally for all datasets. You should check the data look sensible before drawing strong conclusions. Are the distributions as you would expect after logging, normalising and imputing? Do any technical replicates cluster tightly? Is there variation within or among treatments? Do the different clustering methods agree? Answers to these questions will tell you how stable the results are. There are four datasets you can use in the analysis: raw values, log2 values, sample-median normalised, or normalised and imputed. The default is normalised and imputed. The datast can be changed by using the "df_to_use" field in the `configs/auto-prot-config.json` file: "df", "df_log2", "df_median_norm", "df_imp" 

You are of course welcome (and encouraged) to adjust the code to work best for your dataset. [AlphaPepStats](https://github.com/MannLabs/alphapeptstats) has a more comprehensive set of options for pre-processing data, which can be added to this tool by editing the `data_processing.py` script.

As always, we appreciate feedback and suggestions. Please create issues or get in touch - details on main README file.

#### Normalisation
We saw considerable variability in mean and total abundance among samples and treatments. This can look like lots of proteins are overexpressed in a given treatment group, but reflects overall sample abundance (i.e. a sampling artefact) rather than genuine differences per protein. To address this, we normalise by dividing all values by the sample median. In most cases, we are looking for relative differences among treatments, so this strategy is appropriate. However, if most proteins really do go up/down in one condition, the effect will be reduced or lost. Similarly, you won't get a clear idea of the total abundance per sample/treatment.

#### Imputation
Previous studies have found Random Forest imputation methods perform best (see below for refs). The pipeline can use python functions `HistGradientBoostingRegressor` and `IterativeImputer` but with more than 1,000 proteins, it becomes very slow with default parameters. There are some we have tweaked to reduce imputation time:
- HistGradientBoostingRegressor
  -   `max_iter`
  -   `min_samples_leaf`
  -   `max_depth`
-   IterativeImputer
    - `max_iter`
    - `n_nearest_features`

##### Imputation Details 
To impute missing values in protein abundance data, we use iterative multivariate imputation via the `IterativeImputer` from scikit-learn. This method estimates each missing value by modelling it as a function of the most correlated features, cycling through all features over multiple iterations to refine the estimates.

To estimate each value, we use the `HistGradientBoostingRegressor` (HGBR) — a fast, regularised "gradient boosting" method that trains shallow decision trees sequentially to improve predictions. HGBR accelerates training using histogram-based feature binning and handles missing values natively.

While previous studies have found Random Forests particularly effective for imputation due to their robustness and ability to capture nonlinear feature interactions, we follow the approach used in AlphaPepStats by applying HGBR, which is similar to Random Forest but can offer higher predictive power.

## Refs
Kokla, Marietta, et al. "Random forest-based imputation outperforms other methods for imputing LC-MS metabolomics data: a comparative study." BMC bioinformatics 20 (2019): 1-11.

Krismer, Elena, et al. "AlphaPeptStats: an open-source Python package for automated and scalable statistical analysis of mass spectrometry-based proteomics." Bioinformatics 39.8 (2023): btad461.

Jin, Liang, et al. "A comparative study of evaluating missing value imputation methods in label-free proteomics." Scientific reports 11.1 (2021): 1760.



