# Methods boilerplate text
If you publish results generated with this tool, please cite according to instructions on the main README. Here, we provide a template for describing what auto-prot does.

## Short
Note this would come after the description of how protein intensities/abundances are generated. Sections can be chosen and edited as required. Version can be checked locally or in the latest conda environment files in `.../configs`. R packages may reflect local system package setup - please check what you used.

General summary
_We used auto-prot version X (REF to auto-prot) to process the protein abundances and estimate differences between groups. Full details are given at https://github.com/UCL-Biosciences/Auto-Prot. Analyses were carried out using python vXX and vXX._ 

Processing
_After converting zero values to missing, we removed proteins not observed in at least XX% of samples in all groups, leaving XX proteins. We then log2-transformed, normalised the data using the R package vsn (version XX) and imputed with the CollaborativeFilteringTransformer function from python package pimms (version XX)/HistGradientBoostingRegressor from the python package scikit-learn (version XX)._

Clustering
_To understand relationships among samples, we generated PCA and MDS plots using scikit-learn (version XX) and scipy (vXX). We made a heatmap showing relationships among samples and proteins using the clustermap function from sseaborn (vXX). Proteins with missing data were not included in making these plots._

Differential Expression
_Differental expression calculations were done using the R package limma (version XX). For each pair of treatments, we ran limma, including Bayesian variance shrinkage, and generated volcano plots with seaborn (version XX). We used Benajmini-Hochberg False Discovery Rate adjustment to control for multiple testing. Proteins were considered differentially expressed if they had an FDR-adjusted p-value <XX and an absolute log fold change of at least XX._

Pathway Enrichment
_We performed pathway enrichment using a python implementation of Gprofiler (version XX). We queried against the Gprofiler database for species XX and report enrichment of GO/KEGG terms. We applied Gprofiler's g_SCS correction for multiple testing and considered terms enriched if they had an adjusted p-value < 0.05. We plot pathway enrichment data with seaborn (version XX)._

## Package references (delete as appropriate)
vsn: Huber, W. et al. (2002). “Variance Stabilization Applied to Microarray Data Calibration and to the Quantification of Differential Expression.” Bioinformatics, 18 Suppl. 1, S96-S104.
pimms-learn: Webel, H. et al. Imputation of label-free quantitative mass spectrometry-based proteomics data using self-supervised deep learning. Nat Commun 15, 5405 (2024). https://doi.org/10.1038/s41467-024-48711-5
scikit-learn: Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.
scipy: Virtanen, P., et al. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272. DOI: 10.1038/s41592-019-0686-2.
seaborn: Waskom, M. L., (2021). seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021, https://doi.org/10.21105/joss.03021.
limma: Ritchie, M. E. et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47–e47.
gprofiler: Kolberg, L. et al. "g: Profiler—interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update)." Nucleic acids research 51.W1 (2023): W207-W212.








