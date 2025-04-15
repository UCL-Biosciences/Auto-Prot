# Stats and Calculations
A bit more detail on important calculations.

## Pre Processing
### Stats

Include details of packages used and others available

This tool is designed to produce a reasonable first exploration of your data. However, no single approach works optimally for all datasets. You should check the data look sensible before drawing strong conclusions. Are the distributions as you would expect after logging, normalising and imputing? Do any technical replicates cluster tightly? Is there variation within or among treatments? Do the different clustering methods agree? Answers to these questions will tell you how stable the results are.

You are of course welcome (and encouraged) to adjust the code to work best for your dataset. [AlphaPepStats](https://github.com/MannLabs/alphapeptstats) has a more comprehensive set of options for pre-processing data, which can be added to this tool by editing the `data_processing.py` script.

As always, we appreciate feedback and suggestions. Please create issues or get in touch - details on main README file.

#### Normalisation
We saw considerable variability in mean and total abundance among samples and treatments. This can look like lots of proteins are overexpressed in a given treatment group, but reflects overall sample abundance (i.e. a sampling artefact) rather than genuine differences per protein. To address this, we normalise by dividing all values by the sample median. In most cases, we are looking for relative differences among treatments, so this strategy is appropriate. However, if most proteins really do go up/down in one condition, the effect will be reduced or lost. Similarly, you won't get a clear idea of the total abundance per sample/treatment.

#### Imputation
Alphastats recommends using a Random Forest approach. The pipeline can use python functions `HistGradientBoostingRegressor` and `IterativeImputer` but with more than 1,000 proteins, it becomes very slow with default parameters. There are some we have tweaked to reduce imputation time:
- HistGradientBoostingRegressor
  -   `max_iter` reduced to 30 (default 100)
  -   `min_samples_leaf` to 5 (default 31)
  -   `max_depth` to 4 (default is unconstrained)
-   IterativeImputer
    - `max_iter` reduced to 3 (default 10)
    - `n_nearest_features` to 30 (default all features). How many other proteins to use when imputing. This takes the 50 most correlated proteins, with the assumption that they will be the most informative. Improves performance but also removes influence of noisy/unrelated proteins.

## Fold change
We add +1 to each group when calculating fold change. Avoids /0 problems. Adding the minimum count to avoid influencing FC.
