# Compare analysis methods
Given the range of options available for analysing proteomics data, we thought it would be useful to do a quick demo of the similarities and differences between methods and packages for calculating differential expression. This aims to show two things:
1. Auto-Prot isn't crazy. Auto-Prot uses the R package Limma for DE calculation - it should produce similar results to results produced by a simple Limma pipeline in R.
2. Choice of method can have big impact on results.

It is not a comprehensive benchmarking exercise. It is supposed to compare light-touch analyses by Auto-Prot, Limma if exectued in minimal R scripts, and MS stats (another R package for proteomics data).

There are four main steps:

1. Set up the R env for running the R script with Limma and MS stats functions
2. Run the R script - this creates DE output for Limma (once on log-transformed data, once after normalising and imputing) and MS stats
3. Run auto-prot twice. Once analysing the log2-transformed data, once analysing the normalised and imputed data
4. Work through the comparison notebooks to see the log fold changes and p-values produced by the different methods.

## Usage
Note, this should be straight forward if you are confident with running programmes from the command line. Experience with Auto-Prot also helps. It might be tricky for people with less experience. Have a go and let us know if you have problems via contact details on main repo.

### Prerequisites: install R
The R script is independent of conda and manages its own packages via `renv`. You need system R installed and on your PATH — do **not** run this inside a conda environment.

1. Download and install R 4.5 or above from https://cran.r-project.org
2. Add R to your PATH:
   - **Windows (Git Bash)** — add to `~/.bashrc` (replace `R-4.x.x` with your version):
     ```bash
     echo 'export PATH="/c/Program Files/R/R-4.x.x/bin:$PATH"' >> ~/.bashrc
     source ~/.bashrc
     ```
   - **macOS** — add to `~/.zshrc`:
     ```bash
     echo 'export PATH="/Library/Frameworks/R.framework/Resources/bin:$PATH"' >> ~/.zshrc
     source ~/.zshrc
     ```
3. Install system build tools needed to compile some R packages:
   - **Windows**: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
   - **macOS**: run `xcode-select --install` in Terminal
4. Verify: `Rscript --version`

### Run R Script
Deactivate any conda environment first (`conda deactivate`), then run:

```bash
Rscript autoprot/r_scripts/DE-with-MSstats-Limma.R
```

On first run, all required R packages are installed automatically into `renv/library/` inside the repo — nothing is written to your system R libraries. Packages are cached globally so a second clone won't re-download them. This first run may take several minutes.

The script:
- downloads data from an online tutorial
- processes and filters
- calculates DE with MS stats
- runs Limma on the log2-transformed data
- runs Limma after normalising (vsn) and imputing with QRILC, a quantile regression approach for the imputation of left-censored missing data in quantitative proteomics.

All of the outputs are saved to file.

### Run Auto-Prot twice
We want to compare Limma-in-R and Auto-Prot twice; once with just log-transformed data where we know the results should be the same, and once after applying normalisation (both VSN) and imputation (different in Limma-in-R and Auto-Prot). This allows us to see where differences arise. Therefore, we run Auto-Prot twice, once to calculate DE with log-transformed data and once with imputed data.

1. DE with log2-transformed data. To run this, just change three fields in the config file and run with `python main.py`:
- `df_to_use: df_log2`
- `outPath: output/log2`
- `json_outPath: output/log2/data/data_for_report.json`

2. DE with imputed data. Make the following changes in the config file and run `python main.py`:
- `df_to_use: df_imp` 
- `outPath: output/imp_GB`
- `json_outPath: output/imp_GB/data/data_for_report.json`

make sure you have `imputation_method: hist_grad_boost`.

Note these must be correct as they specify the outpath paths used by the jupyter notebook that compares the results.

### Compare the results
Open the jupyter notebook `autoprot/reporting/compare-autoprot-msstats.ipynb` (you will need something that runs jupyter notebooks). Use the `auto-proteomics` kernel. This is where the fun starts, if it works :'D

First the **log2-transformed** data. At this point, both analyses have received the same data and applied log2-transformation, and both use Limma to calculate DE. The results should be the same.

Now we compare DE results from **imputed data**. This contains two types of data. Proteins that had no missing data and have therefore been normalised but have no imputed values. These should be the same in Auto-Prot and Limma-in-R because they have the same transformation and normalisation. Look at the blue points - are they correlated? There are some proteins with missing values; these have undergone different imputation methods and the values will be different. How different are they? Look at the orange points.

What does this tell you about imputation methods?

Finally we compare Auto-Prot with MS Stats, which uses a different protocol completely. The LFC __should__ be similar but the overall p-values and proteins that are DE may well be different.

## Conclusions
Pretty clear - choices matter! Imputation introduces noise and can significantly impact results. Choose carefully, get to know your data well, and, if possible, validate results experimentally!
