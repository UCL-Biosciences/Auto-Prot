##### Calculate differential protein expression using limma #####
## We process some proteomics data in python
## But limma has a great differential expression/abundance calculation
## This script loads the protein abundances for a treatment pair,
## calculates the DE, and saves to a file
## which can be read again by the main pipeline.

# Define project-local library
find_repo_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, ".git"))) return(cur)
    parent <- dirname(cur)
    if (parent == cur) stop("No git repository found above ", start)
    cur <- parent
  }
}

repo_root <- find_repo_root()
proj_lib <- file.path(repo_root, "output/r_libs")
if (!dir.exists(proj_lib)) dir.create(proj_lib, recursive = TRUE)

# Prepend to library search path
.libPaths(c(proj_lib, .libPaths()))

# Ensure vsn is available in that local path
if (!requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma",
    lib = proj_lib, ask = FALSE, update = FALSE)
}

### load package
library(limma)

# run_limma.R

#### Load inputs
#### note all inputs come via the python subprocess command
#### where normally we would load from file

## The args are those written at the end of the subprocess command
args <- commandArgs(trailingOnly = TRUE)
## args 1 has the intensities for the pair to be analysed
expr_file <- args[1]
# Check input file exists
if (!file.exists(expr_file)) {
  stop("ERROR: Abundance file does not exist.")
}
# metedata (including sample_id, treatment and any other vars to be included in the model)
meta_file <- args[2]
output_file <- args[3]
### the formula to be used in the model
formula_str <- args[4]

# Load data
df <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
mode(df) <- "numeric"

# load metadata
meta <- read.csv(meta_file)

# design is taken directly from the input string
design <- model.matrix(as.formula(formula_str), data = meta)

# Fit model
print( paste0("running model with formula: ", formula_str))
fit <- lmFit(df, design)
print("running ebayes")
fit <- eBayes(fit)

# Output result
res <- topTable(fit, coef = 2, number = Inf)
write.csv(res, file = output_file)

cat("limma ran successfully. DE table written to:\n")
print(output_file)