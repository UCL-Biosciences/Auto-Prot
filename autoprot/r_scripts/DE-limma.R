##### Calculate differential protein expression using limma #####
## We process some proteomics data in python
## But limma has a great differential expression/abundance calculation
## This script loads the protein abundances for a treatment pair,
## calculates the DE, and saves to a file
## which can be read again by the main pipeline.

# Set seed for reproducibility
set.seed(123)

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
proj_lib <- file.path(repo_root, "r_libs")
if (!dir.exists(proj_lib)) dir.create(proj_lib, recursive = TRUE)

# Prepend to library search path
.libPaths(c(proj_lib, .libPaths()))

# Ensure vsn is available in that local path
if (!requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma",
    lib = proj_lib, ask = FALSE, update = FALSE)
}

### load package
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

# Define project-local library
proj_lib <- file.path(repo_root, "r_libs")
.libPaths(c(proj_lib, .libPaths()))

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", lib = proj_lib, version = "3.20")
}
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
ref_level <- args[5]  # base level to use in limma calcs - so it is consistent with figure labelling

# Load data
df <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
mode(df) <- "numeric"

# load metadata
meta <- read.csv(meta_file)
meta$treatment <- relevel(factor(meta$treatment), ref = ref_level)

# design is taken directly from the input string
design <- model.matrix(as.formula(formula_str), data = meta)

# Fit model
print( paste0("running model with formula: ", formula_str))
fit <- lmFit(df, design)
print("running ebayes")
fit <- eBayes(fit)

## which coefficient to use for the DE calculation?
# two steps - first find which metadata column is the treatment to contrast
# then set the coefficient to the treatment level of interest (which is the second level of the factor)

## find the grouping column
# Identify the grouping variable as the first term in the formula.
# By convention the config formula must list the group variable first
# (e.g. "~ treatment" or "~ treatment + batch"), which avoids passing the
# column name in as a separate argument.
group_col <- trimws(strsplit(sub("^\\s*~\\s*", "", formula_str), "[+*]")[[1]][1])

if (!group_col %in% colnames(meta)) {
  stop("Grouping variable '", group_col, "' from formula '", formula_str,
       "' is not a column in the metadata. Metadata columns are: ",
       paste(colnames(meta), collapse = ", "))
}

meta[[group_col]] <- relevel(factor(meta[[group_col]]), ref = ref_level)

## Then define the contrast level for the group column
# ref_level is passed as arg so level to contrast is the other level in the pair
other_level <- setdiff(levels(meta[[group_col]]), ref_level)
stopifnot(length(other_level) == 1)  # Python passes one pair at a time

# Now we can find the coefficient name in the design matrix
# The coefficient name is typically in the form "group_colOtherLevel", e.g., "treatmentDrugA" if group_col is "treatment" and other_level is "DrugA".
coef_name <- paste0(group_col, other_level)
if (!coef_name %in% colnames(design)) {
  stop("Expected coefficient '", coef_name, "' not found in design matrix. ",
       "Design columns are: ", paste(colnames(design), collapse = ", "))
}

res <- topTable(fit, coef = coef_name, number = Inf)

# Output result
res <- topTable(fit, coef = 2, number = Inf)
write.csv(res, file = output_file)

cat("limma ran successfully. DE table written to:\n")
print(output_file)