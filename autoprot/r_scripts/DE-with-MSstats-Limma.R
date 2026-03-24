#!/usr/bin/env Rscript
#
# Proteomics analysis to compare with Auto-Prot
#
# This script downloads and processes proteomic data,
# then performs differential expression analysis using
# MSstats and limma.
#
# It uses `renv` to ensure reproducible environments:
# - project-local library only (no system packages)
# - packages recorded in renv.lock
# - reproducible installation from CRAN/Bioconductor/GitHub
#

## -------------------------------
## Locate repo root
## -------------------------------
### if this doesn't work, set your wd (setwd()) to the Auto-Prot dir
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

## -------------------------------
## Bootstrap renv (fully isolated)
## -------------------------------
# All packages install into renv/library/ inside the repo.
# Nothing is written to your system or user R libraries.
#
# First run: installs all packages automatically (takes a few minutes).
# Subsequent runs: instant (packages already cached).
#
# Requires: R installed from https://cran.r-project.org and on your PATH.
# Do NOT run this script with conda's R — activate no conda env beforehand.

if (grepl("conda|envs", R.home(), ignore.case = TRUE) ||
    nzchar(Sys.getenv("CONDA_DEFAULT_ENV"))) {
  stop(
    "This script must be run with system R (https://cran.r-project.org), ",
    "not conda's R.\n",
    "Deactivate any conda environment first: conda deactivate"
  )
}

source(file.path(repo_root, "renv", "activate.R"))  # bootstraps renv if needed, isolates library

if (getRversion() < "4.5.0") stop("R >= 4.5.0 is required. Download from https://cran.r-project.org")

message("Restoring R packages from renv.lock (first run may take several minutes)...")
renv::restore(project = repo_root, prompt = FALSE)

## -------------------------------
## Load packages
## -------------------------------
library(dplyr)
library(tidyr)
library(vsn)
library(MSstats)
library(limma)
library(MSnbase)
library(imputeLCMD)

## -------------------------------
## Set dirs
## -------------------------------
out_dir = file.path(repo_root, "output/compare-msstats"); dir.create(out_dir, recursive = T, showWarnings = F)

## -------------------------------
## Analysis pipeline
## -------------------------------

### Download and prepare data ###

#
# Example dataset from Galaxy training material:
# https://training.galaxyproject.org/training-material/topics/proteomics/tutorials/maxquant-msstats-dda-lfq/tutorial.html
#

## -------------------------------
## Download data
## -------------------------------
print("downloading data")
evidence_url <- "https://zenodo.org/record/4896554/files/MaxQuant_Evidence.tabular"
pg_url       <- "https://zenodo.org/record/4896554/files/MaxQuant_proteingroups.tabular"

options(timeout = 300)  # set to 5 minutes
evi_full <- read.delim(evidence_url, stringsAsFactors = FALSE)
pg_full  <- read.delim(pg_url, stringsAsFactors = FALSE)

## -------------------------------
## Basic filtering
## -------------------------------
print("processing data")
pg_full <- pg_full %>%
  filter(
    Potential.contaminant != "+",
    Reverse != "+",
    Only.identified.by.site != "+"
  )

# Replace zeros with NA (MaxQuant LFQ often encodes missing as 0)
pg_full[pg_full == 0] <- NA

## -------------------------------
## Metadata: sample groups
## -------------------------------
lfq_cols <- grep("^LFQ\\.intensity", colnames(pg_full), value = TRUE)

samples <- data.frame(
  sample = lfq_cols
) %>%
  mutate(
    group = gsub("^LFQ\\.intensity\\.", "", sample) %>%
      sub("_.*$", "", .) %>%
      factor()
  )

## -------------------------------
## Filtering by missing values
## -------------------------------
# Threshold: require ≥75% non-missing values in *every* group
missing_threshold <- 0.75

# Indices of LFQ columns per group
group_idx <- split(seq_along(lfq_cols), samples$group)

# Row-wise filter function
keep <- apply(pg_full[lfq_cols], 1, function(x) {
  all(sapply(group_idx, function(idx) mean(!is.na(x[idx])) >= missing_threshold))
})

pg_full <- pg_full[keep, ]

### write protein groups to file so they can be used with auto-prot
write.csv(pg_full %>% dplyr::rename(Genes = Protein.IDs),
          file.path(repo_root, "input/data/proteindata.csv"),
          row.names = F)


## -------------------------------
## Output check
## -------------------------------
head(pg_full[lfq_cols])

## -------------------------------
## MS stats analysis
## -------------------------------
print("run MS Stats")
# Build annotation: Raw.file must match evi$Raw.file exactly
# Columns required: Raw.file, Condition, BioReplicate, Run
msstats_meta <- evi_full %>%
  distinct(Raw.file) %>% # one row per file
  arrange(Raw.file) %>% # sort by name
  mutate(
    Condition = sub("_.*$", "", Raw.file), # e.g. "metast" from "metast_1.raw"
    BioReplicate = ave(seq_along(Raw.file), Condition, FUN = seq_along), # 1, 2, ... per condition
    Run = 1 # dummy column (not used in this design)
  )

# Convert to MSstats format
quant <- MaxQtoMSstatsFormat(
  evidence = evi_full, 
  annotation = msstats_meta,
  proteinGroups = pg_full,
  useUniquePeptide = TRUE, # only unique peptides 
  removeFewMeasure = TRUE, # filter proteins with few measurements
  removeProtein_with1Feature = TRUE, # filter proteins with only 1 peptide
  summaryforMultipleRows = max) # summarise multiple intensity values per peptide

# Process to protein level
# dataProcess normalises, imputes, and performs QC
# see ?dataProcess for options
processed <- dataProcess(quant)

### DE analysis ###
levels(processed$ProteinLevelData$GROUP)
comparison = matrix( c(-1, 1), nrow = 1 )
row.names(comparison) <- "metast-RDEB%"
colnames(comparison) = c("metast", "RDEB" ) 

# Tests for differentially abundant proteins with models:
testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=processed)

# save results
write.csv( file = file.path(out_dir, "MSstats-out.csv"),
           testResultOneComparison$ComparisonResult )

## -------------------------------
## Get missingness info
## -------------------------------

### print number missing per protein
# Select LFQ columns
lfq_cols <- grep("^LFQ", names(pg_full), value = TRUE)

# Calculate proportion missing per row
pg_full$proportion_missing <- rowMeans(is.na(pg_full[lfq_cols]))

# Build output table
prop_missing <- pg_full[, c("Protein.IDs", "proportion_missing")]

write.csv(file = file.path(out_dir, "prop_missing.csv"),
          prop_missing, row.names = F)

metadata = data.frame(
  protein_abundance_name = pg_full %>% select(starts_with("LFQ")) %>% colnames()) %>%
  mutate( 
    sample_id = gsub("Intensity.", "", protein_abundance_name ) %>% sub("\\.raw\\..*$", "", .),
    replicate = 1,
    timepoint = 1,
    treatment  = gsub("LFQ.intensity.", "", protein_abundance_name ) %>% sub("_.*$", "", .)
  )

write.csv(metadata,
          file.path(repo_root, "input/data/metadata.csv"),
          row.names = F)

## -------------------------------
## Limma analysis
## -------------------------------
print("running limma")
# Select LFQ columns
lfq <- pg_full[, lfq_cols]
rownames(lfq) <- pg_full$Protein.IDs

# log2 transform
lfq_log <- log2(lfq)

# design matrix
# requires treatment groups and samples to be in same order as lfq columns
design <- model.matrix(~0 + samples$group)
colnames(design) <- levels(samples$group)

### fit the model without normalisation or imputation
## we will first compare Auto-Prot and this "raw" limma resulton log2 data
fit <- lmFit(lfq_log, design) # fit the model
cont.matrix <- makeContrasts(metast - RDEB, levels=design) # define contrast matrix
fit2 <- contrasts.fit(fit, cont.matrix) # 
fit2 <- eBayes(fit2) # apply ebayes shrinkage

res <- topTable(fit2, adjust="BH", number=Inf) # get all results including BH FDR adjustment
head(res)
write.csv(res, file.path(out_dir, "Rstudio-limma-out-log2.csv"))

### Now with normalisation and imputation
# normalise
# Convert to matrix
## note vsn requires raw positive values, i.e. no log transformation
expr <- as.matrix(lfq)

# apply VSN
lfq_vsn <- vsn2(expr)

#--- Simple imputation ---
lfq_ms <- MSnSet(exprs = as.matrix(lfq_vsn))

# imputation
lfq_imp <- MSnbase::impute(lfq_ms, method = "QRILC")
exprs_imp = exprs(lfq_imp)

### refit the model as above but with imputed data
fit <- lmFit(lfq_imp, design)
cont.matrix <- makeContrasts(metast - RDEB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

res2 <- topTable(fit2, adjust="BH", number=Inf)
res2[order(rownames(res2)), ] %>% head
write.csv(res2, file.path(out_dir, "Rstudio-limma-out-imputed.csv"))

print("finished creating limma and MS Stats output")