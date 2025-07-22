# Stats and Calculations
This tool is designed to produce a reasonable first exploration of your data. However, no single approach works optimally for all datasets. You should check the data look sensible before drawing strong conclusions. Are the distributions as you would expect after logging, normalising and imputing? Do any technical replicates cluster tightly? Is there variation within or among treatments? Do the different clustering methods agree? Answers to these questions will tell you how stable the results are. There are four datasets you can use in the analysis: raw values, log2 values, sample-median normalised, or normalised and imputed. The default is normalised and imputed. The datast can be changed by using the "df_to_use" field in the `configs/auto-prot-config.json` file: "df", "df_log2", "df_median_norm", "df_imp" (see setup doc for more info).

You are of course welcome (and encouraged) to adjust the code to work best for your dataset. [AlphaPepStats](https://github.com/MannLabs/alphapeptstats) has a more comprehensive set of options for pre-processing data, which can be added to this tool by editing the `data_processing.py` script.

As always, we appreciate feedback and suggestions. Please create issues or get in touch - details on main README file.

## Pre Processing
### Log Transformation
If sample-median normalisation is used, values are log2-transformed. vsn normalisation requires raw values.

### Normalisation
Normalisation is controlled in the config field "normalise_method". "vsn" will use variance stablisation normalisation, "sample-median" will subtract all values from the respective sample's median. Explained below.

The default normalisation is variance stablisiation normalisation, as recommended by Välikangas et al (2018): "We found that variance stabilization normalization (Vsn) reduced variation the most between technical replicates in all examined data sets. Vsn also performed consistently well in the differential expression analysis." Note, "the Vsn normalization performs a transformation similar to the log transformation and requires the input data to be untransformed." More details in Välikangas et al.

"sample-median" will normalise simply by subtracting all values from the median value for each protein. We saw considerable variability in mean and total abundance among samples and treatments. This can look like lots of proteins are overexpressed in a given treatment group, but reflects overall sample abundance (i.e. a sampling artefact) rather than genuine differences per protein. Subtracting the median can help with this but is not as effective as other methods.

### Imputation
#### Collaborative Filter
When you pick the “pimms_collabfiltering” option, the pipeline treats your protein-by-sample matrix like a recommender system: it first learns a compact set of latent factors for proteins and samples via alternating‐least‐squares (ALS) matrix factorisation, then uses those factors in a simple iterative regressor to fill in missing values one protein at a time. You can speed it up by asking it to find fewer latent factors or to make fewer ALS passes over the data, and by limiting how many other proteins it considers when regressing each one, or you can improve accuracy (at the cost of longer runtimes) by increasing those same settings. To speed up the collaborative‐filtering imputer, you can:
- Reduce the number of latent factors (factors) so the ALS solver learns a smaller, simpler representation.
- Lower the ALS iteration count (iterations) so it makes fewer passes over the data.
- Decrease the neighbor count (n_nearest_features) so each regression only looks at a handful of other proteins.

To improve accuracy (at the expense of runtime), you do the opposite:

#### Gradient Boosting
Previous studies have found Random Forest imputation methods perform best (see below for refs). If the "imputation_method" config field is set to "hist_grad_boost", the pipeline will use python functions `HistGradientBoostingRegressor` and `IterativeImputer`. This is a robust machine learning imputation method similar to random forest, but with more than 1,000 proteins, it becomes very slow with default parameters. There are some we have tweaked to reduce imputation time:
- HistGradientBoostingRegressor
  -   `max_iter`
  -   `min_samples_leaf`
  -   `max_depth`
-   IterativeImputer
    - `max_iter`
    - `n_nearest_features`

#### Imputation Details 
To impute missing values in protein abundance data, we use iterative multivariate imputation via the `IterativeImputer` from scikit-learn. This method estimates each missing value by modelling it as a function of the most correlated features, cycling through all features over multiple iterations to refine the estimates.

To estimate each value, we use the `HistGradientBoostingRegressor` (HGBR) — a fast, regularised "gradient boosting" method that trains shallow decision trees sequentially to improve predictions. HGBR accelerates training using histogram-based feature binning and handles missing values natively.

While previous studies have found Random Forests particularly effective for imputation due to their robustness and ability to capture nonlinear feature interactions, we follow the approach used in AlphaPepStats by applying HGBR, which is similar to Random Forest but can offer higher predictive power.

## Exploratory Analysis
In proteomics data analysis, before diving into pairwise treatment comparisons or complex modelling, it's critical to understand the** overall structure** of the dataset. Exploratory visualisations such as PCA, MDS, and heatmaps offer intuitive ways to:
- Assess **variation between biological conditions or treatments**.
- Evaluate **replicate consistency** (do replicates cluster together?).
- Detect** outliers or unexpected batch effects**.
- Identify **global patterns** in protein expression across all samples.

These methods are unsupervised and data-driven—they do not rely on predefined hypotheses, making them powerful tools for data QC and hypothesis generation.

**What does clustering tell us?**
Clustering groups samples or proteins that behave similarly across the dataset. For example:
- If all replicates from a treatment cluster together, it suggests consistent treatment effects.
- If samples from different treatments intermingle, it may indicate weak or overlapping biological signals—or a need to inspect metadata or preprocessing.

Visualising these groupings in 2D (e.g. PCA, MDS) or as a clustered heatmap helps you visually **validate assumptions** and **guide downstream analysis**.

### PCA (Principal Component Analysis)
**Purpose:** Display thousands of proteins using plots with two axes. Reduce dimensionality while retaining as much variance as possible. Highlights major axes of variation in the data.

**Calculation:**
- `sklearn.decomposition.PCA` is run twice
  - PCA is run to determine how many components are needed to explain 95% of the total variance.
  - PCA is then re-run using that optimal number of components.
- A 2D scatterplot of the first two components (PC1 and PC2) is generated.
- Group centroids and 95% confidence ellipses are added (via custom function using covariance - ellipses are scaled to 4x standard deviation ≈95% CI).
- Log-transformed, (VSN-)normalised data are not re-scaled for PCA as VSN stablises variance. If you use a different normalisation, you might need to scale the data.

### MDS (Metric Multidimensional Scaling)
**Purpose:** Visualise similarity or dissimilarity between samples based on distance metrics (here, Euclidean distance on protein abundance).

**Calculation:** 
- Distance matrix is computed via `scipy.spatial.distance.pdist`.
- MDS is performed using `sklearn.manifold.MDS` with `dissimilarity="precomputed"` and `metric=True`.
- Again, plots show ellipses and centroids for treatments.

### Heatmap
**Purpose:** Visualise the full protein abundance matrix, clustering both proteins and samples to reveal global patterns.

**Calculation:**
- Uses `seaborn.clustermap` to perform hierarchical clustering.
- Dendrograms reflect similarity between samples/proteins.
- Column colours are dynamically re-ordered to match clustering.
- Heatmap uses the raw input orientation (samples = rows, proteins = columns), unlike PCA and MDS, which require transposed data.

## Differential Expression
We use the Limma R Package (Ritchie et al. 2015) to quantify differences between treatment groups for all proteins (following e.g. Lou et al. 2022). Differential expression is calculated for each pair of treatments. 

By default, it will fit a model comparing all samples in each treatment group (i.e. `~ treatment`). You can specify additional variables for the model using the config file. You should be able to specify any [valid R formula](https://www.datacamp.com/tutorial/r-formula-tutorial) contained in brackets. The DE_full_formula field specifies the formula for the full dataset, and can include any additional variables in the metadata. For example, you could fit a model with  `~ treatment + timepoint` to account for differences among timepoints.

A separate formula is supplied for any subsets that you want to analyse DE for. This is because you might want to create subsets for a variable that is contained in the full dataset's model - e.g. timepoint - which might break the analysis.

Be careful interpreting the direction of the differential expression calculations. By default, the pairs of treatment groups compared are all unique combinations of the treatment column in metadata. 

## Pathway Enrichment
This analysis step performs pathway enrichment using the g:Profiler tool, focusing on proteins with large fold changes and low adjusted p-values, based on thresholds defined in the config. It uses over-representation analysis (ORA) to identify enriched biological processes or pathways. By default, we use KEGG ontology. This can be changed but currently needs to be changed in `src/analysis/pairwise.py`, `source = ["source"]` where source can be "KEGG", "GO", or "REAC". Results are saved as a CSV file and visualised as a bubble plot showing the top 20 pathways ranked by p-value. Note, there are many cases where relatively few proteins are differentially expressed and no enrichment terms are found. In this case you can investigate whether the code is running as expected and can reduce the thresholds that define differential expression through config fields `LFC_threshold` (default = 1) and `FDR_threshold` (default = 0.05).

## Refs
Kokla, Marietta, et al. "Random forest-based imputation outperforms other methods for imputing LC-MS metabolomics data: a comparative study." BMC bioinformatics 20 (2019): 1-11.

Krismer, Elena, et al. "AlphaPeptStats: an open-source Python package for automated and scalable statistical analysis of mass spectrometry-based proteomics." Bioinformatics 39.8 (2023): btad461.

Lou, Ronghui, et al. "Benchmarking commonly used software suites and analysis workflows for DIA proteomics and phosphoproteomics." Nature communications 14.1 (2023): 94.

Jin, Liang, et al. "A comparative study of evaluating missing value imputation methods in label-free proteomics." Scientific reports 11.1 (2021): 1760.

Ritchie, Matthew E., et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies." Nucleic acids research 43.7 (2015): e47-e47.

Välikangas, Tommi, Tomi Suomi, and Laura L. Elo. "A systematic evaluation of normalization methods in quantitative label-free proteomics." Briefings in bioinformatics 19.1 (2018): 1-11.



