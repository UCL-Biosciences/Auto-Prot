# Methods boilerplate text
If you publish results generated with this tool, please cite according to instructions on the main README. Here, we provide a template for describing what auto-prot does.

## Short
Note this would come after the description of how protein intensities/abundances are generated. Sections can be chosen and edited as required. Version can be checked locally or in the latest conda environment files in `.../configs`. R packages may reflect local system package setup - please check what you used.

General summary
_We used auto-prot version X (REF to auto-prot) to process the protein abundances and estimate differences between groups. Full details are given at https://github.com/UCL-Biosciences/Auto-Prot._ 

Processing
_After converting zero values to missing, we removed proteins not observed in at least XX% in all groups, leaving XX proteins. We then log2-transformed, normalised the data using the R package vsn (version XX) and imputed with the CollaborativeFilteringTransformer function from python package pimms (version XX)/HistGradientBoostingRegressor from the python package scikit-learn (version XX)._

Clustering
_To understand relationships among samples, we generated PCA and MDS plots using scikit-learn (version XX) and scipy (vXX). We made a heatmap showing relationships among samples and proteins using the clustermap function from sseaborn (vXX). Proteins with missing data were not included in making these plots._

Differential Expression
_Differental expression calculations were done using the R package limma (version XX). For each pair of treatments, we ran limma, including Bayesian variance shrinkage, and generated volcano plots with seaborn (version XX). We used Benajmini-Hochberg False Discovery Rate adjustment to control for multiple testing. Proteins were considered differentially expressed if they had an FDR-adjusted p-value <XX and an absolute LFC > XX._

Pathway Enrichment
_We performed pathway enrichment using a python implementation of Gprofiler (version XX). We queried against the Gprofiler database for species XX and report enrichment of GO/KEGG terms. We applied Gprofiler's g_SCS correction for multiple testing and considered terms enriched if they had an adjusted p-value < 0.05. We plot pathway enrichment data with seaborn (version XX)._




