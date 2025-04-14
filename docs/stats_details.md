# Stats and Calculations
A bit more detail on important calculations.

## Pre Processing
### Stats

Include details of packages used and others available

#### Imputation
Alphastats recommends using a Random Forest approach. The pipeline can use python functions `HistGradientBoostingRegressor` and `IterativeImputer` but with more than 1,000 proteins, it becomes very slow with default parameters. There are some we have tweaked to reduce imputation time:
- HistGradientBoostingRegressor
  -   `max_iter` reduced to 50 (default 100)
  -   `min_samples_leaf` to 3 (default 31)
  -   `max_depth` to 5 (default is unconstrained)
-   IterativeImputer
    - `max_iter` reduced to 5 (default 10)
    - `n_nearest_features` to 50 (default all features). How many other proteins to use when imputing. This takes the 50 most correlated proteins, with the assumption that they will be the most informative. Improves performance but also removes influence of noisy/unrelated proteins.


## Fold change

We add +1 to each group when calculating fold change. Avoids /0 problems. Adding the minimum count to avoid influencing FC.
