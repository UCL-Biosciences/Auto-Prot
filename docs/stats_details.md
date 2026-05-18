# Stats and Calculations
This tool is designed to produce a reasonable first exploration of your data. However, no single approach works optimally for all datasets. You should check the data look sensible before drawing strong conclusions. Are the distributions as you would expect after logging, normalising and imputing? Do any technical replicates cluster tightly? Is there variation within or among treatments? Do the different clustering methods agree? Answers to these questions will tell you how stable the results are. There are four datasets you can use in the analysis: raw values, log2 values, sample-median normalised, or normalised and imputed. The default is normalised and imputed. The datast can be changed by using the "df_to_use" field in the `configs/auto-prot-config.yaml` file: "df", "df_log2", "df_norm", "df_imp" (see setup doc for more info).

You are of course welcome (and encouraged) to adjust the code to work best for your dataset. [AlphaPepStats](https://github.com/MannLabs/alphapeptstats) has a more comprehensive set of options for pre-processing data, which can be added to this tool by editing the `data_processing.py` script.

As always, we appreciate feedback and suggestions. Please create issues or get in touch - details on main README file.

## Pre Processing
### Zero Values and Missing Data
Important note: mass spectrometry instruments record the intensities of proteins that are detected above a certain threshold. Below that threshold, it is not clear whether a protein is completely absent or whether it just wasn't detected in the sample. In this pipeline, we assume all zeros are missing values and convert them to NA before further processing.

Low-variance proteins are removed prior to analysis. The interquartile range of normalised intensities is calculated for all proteins, and the bottom XX% are removed. The threshold is determined by the "IQR_threshold" config parameter. The default value is 0.1 (10%) and can be set to 0 if you don't want to remove any low-variance proteins.

### Log Transformation
If sample-median normalisation is used, values are log2-transformed. vsn normalisation requires raw values.

### Normalisation
Normalisation is controlled in the config field "normalise_method". "vsn" will use variance stablisation normalisation, "sample-median" will subtract all values from the respective sample's median. Explained below.

The default normalisation is variance stablisiation normalisation, as recommended by Välikangas et al (2018): "We found that variance stabilization normalization (Vsn) reduced variation the most between technical replicates in all examined data sets. Vsn also performed consistently well in the differential expression analysis." Note, "the Vsn normalization performs a transformation similar to the log transformation and requires the input data to be untransformed." More details in Välikangas et al.

"sample-median" will normalise simply by subtracting all values from the median value for each protein. We saw considerable variability in mean and total abundance among samples and treatments. This can look like lots of proteins are overexpressed in a given treatment group, but reflects overall sample abundance (i.e. a sampling artefact) rather than genuine differences per protein. Subtracting the median can help with this but is not as effective as other methods.

There is also an option to scale values before doing clustering analyses (PCA, MDS, and heatmap) using scikit-learn's [`StandardScaler`](https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html#sklearn.preprocessing.StandardScaler). By default, this is set to true.

### Imputation
In MS proteomics, it is common to have a lot of missing values. Some analyses (e.g. PCA, some machine learning techniques) can't handle missing values. That means you must either remove all proteins with any missing values or impute the missing values. Whether to do that is an important decision and requires consideration from researchers.

Different imputation methods vary in speed, scalability, and how well they preserve biological structure. Gradient boosting offers high accuracy but can be slow on large datasets, while PIMMS is faster and scales better, especially with thousands of proteins, though it may be slightly less precise in capturing fine-grained patterns. PIMMS' collaborative filtering model is tested on a dataset with ~ 50 samples and may not perform well with smaller datasets. Imputation options are controlled in the `config` file. See `docs/setup.md` for details.

Note, gradient boosting normalisation requires normalised values, while PIMMS works with the raw or log-transformed values. I tried normalising after PIMMS imputing and the LFCs were all tiny - it removed too much natural variation. I removed it and now the PIMMS pipeline includes no normalisation. It isn't clear from the PIMMS paper whether they normalise afterwards or not.

#### Gradient Boosting Imputation
Previous studies have found Random Forest imputation methods perform best (see below for refs). The pipeline can use python functions `HistGradientBoostingRegressor` and `IterativeImputer` but with more than 1,000 proteins, it becomes slow with default parameters. You can see the parameters we used to improve performance in `autoprot/processing/data_preprocess.py` (see `impute_prot_data_histgradboost` function).

To fill in missing values in our protein abundance data, we use a method called _iterative multivariate imputation_ with scikit-learn’s `IterativeImputer`. This approach predicts each missing value based on the values of correlated proteins. It goes through the data several times, gradually improving its guesses by learning patterns between proteins.

For each prediction, we use a model called `HistGradientBoostingRegressor` (HGBR). This is a fast and efficient machine learning algorithm that builds a series of simple decision trees, where each new tree helps correct the errors of the previous ones. HGBR is especially useful here because it can handle missing values directly and speeds up training by grouping similar values into bins.

While previous studies have found Random Forests particularly effective for imputation due to their robustness and ability to capture nonlinear feature interactions, we follow the approach used in AlphaPepStats (Krismer et al., 2023) by applying HGBR, which is similar to Random Forest but can offer higher predictive power.

#### PIMMS Imputation
PIMMS (Proteomics Imputation Modeling Mass Spectrometry; Webel et al. 2024) is a deep learning-based framework designed to impute missing values in MS-based proteomics data. Unlike traditional methods that rely on distributional assumptions or local similarity (e.g. random forest or k-nearest neighbours), PIMMS uses self-supervised models trained directly on intensity data to learn complex patterns between proteins and samples.

PIMMS includes three model types: collaborative filtering (CF), denoising autoencoders (DAE), and variational autoencoders (VAE). All are trained to reconstruct observed intensities and predict missing ones, without requiring normalisation. These models scale well to high-dimensional proteomics data and can outperform or match traditional imputation methods on multiple datasets.

In our pipeline, we use the collaborative filtering model implemented in the `impute_pimms_cf` function (see autoprot/processing/data_preprocess.py). This method learns latent representations of samples and proteins and predicts missing values using their interactions. It runs efficiently on medium to large datasets and supports GPU acceleration. Compared to other methods, PIMMS-CF offers a good balance of scalability, performance, and robustness across a wide range of intensity values and missingness patterns.

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
This analysis step performs pathway enrichment using the g:Profiler tool, focusing on proteins with large fold changes and low adjusted p-values, based on thresholds defined in the config. It uses over-representation analysis (ORA) to identify enriched biological processes or pathways. By default, we use KEGG ontology. This can be changed in the config file to "KEGG", "GO", or "REAC". Results are saved as a CSV file and visualised as a bubble plot showing the top 20 pathways ranked by p-value. Note, there are many cases where relatively few proteins are differentially expressed and no enrichment terms are found. In this case you can investigate whether the code is running as expected and can reduce the thresholds that define differential expression through config fields `LFC_threshold` (default = 1) and `FDR_threshold` (default = 0.05).

## Refs
Kokla, Marietta, et al. "Random forest-based imputation outperforms other methods for imputing LC-MS metabolomics data: a comparative study." BMC bioinformatics 20 (2019): 1-11.

Krismer, Elena, et al. "AlphaPeptStats: an open-source Python package for automated and scalable statistical analysis of mass spectrometry-based proteomics." Bioinformatics 39.8 (2023): btad461.

Lou, Ronghui, et al. "Benchmarking commonly used software suites and analysis workflows for DIA proteomics and phosphoproteomics." Nature communications 14.1 (2023): 94.

Jin, Liang, et al. "A comparative study of evaluating missing value imputation methods in label-free proteomics." Scientific reports 11.1 (2021): 1760.

Ritchie, Matthew E., et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies." Nucleic acids research 43.7 (2015): e47-e47.

Välikangas, Tommi, Tomi Suomi, and Laura L. Elo. "A systematic evaluation of normalization methods in quantitative label-free proteomics." Briefings in bioinformatics 19.1 (2018): 1-11.

Webel, Henry, et al. "Imputation of label-free quantitative mass spectrometry-based proteomics data using self-supervised deep learning." Nature Communications 15.1 (2024): 5405.

