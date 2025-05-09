# cat("R version: ", R.version.string, "\n")
# cat("limma version: ", as.character(packageVersion("limma")), "\n")
# cat("Working directory: ", getwd(), "\n")
# cat("R_HOME: ", Sys.getenv("R_HOME"), "\n")


##### Calculate differential protein expression using limma #####
## We process some proteomics data in python
## But limma has a great differential expression/abundance calculation
## This script loads the protein abundances for a treatment pair,
## calculates the DE, and saves to a file
## which can be read again by the main pipeline.

# run_limma.R
args <- commandArgs(trailingOnly = TRUE)
expr_file <- args[1]
# Check input file exists
if (!file.exists(expr_file)) {
  stop("ERROR: Abundance file does not exist.")
}
meta_file <- args[2]
output_file <- args[3]

# cat("Expression file: ", expr_file, "\n")
# cat("Metadata file: ", meta_file, "\n")
# cat("Output file: ", output_file, "\n")


### load package
library(limma)

# Load data
df <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
mode(df) <- "numeric"
# cat("Expression data loaded. Dimensions: ", dim(df)[1], "proteins x", dim(df)[2], "samples\n")
# cat("Missing values in data:", sum(is.na(df)), "\n")

# load metadata
meta <- read.csv(meta_file)

# Define groups# Load metadata into a named vector: names = sample IDs, values = treatment
treatment_map <- setNames(as.character(meta$treatment), meta$sample_rep)

# Use column names of the expression matrix to look up treatments
treatment <- factor(treatment_map[colnames(df)])

# Check alignment with number of samples
if (length(treatment) != ncol(df)) {
  stop(paste("ERROR: Number of treatment labels (", length(treatment), ") does not match number of samples (", ncol(df), ")"))
}

# cat("Treatment levels:\n")
# print(levels(treatment))
# cat("Treatment vector:\n")
# print(treatment)

# Build design matrix and run limma
design <- model.matrix(~ treatment)

# Fit model
fit <- lmFit(df, design)
fit <- eBayes(fit)

# Output result
res <- topTable(fit, coef = 2, number = Inf)
write.csv(res, file = output_file)

cat("limma ran successfully. DE table written to:\n")
print(output_file)