### MSstats ###

library(renv)

# Initialise once per project if not already done
#renv::init()

# Install MSstats (dependencies auto-installed into project library)
# renv::install("bioc::MSstats", "dplyr", "tidyr")

library(MSstats)
library(dplyr)
library(tidyr)

## tutorial https://www.bioconductor.org/packages/release/bioc/vignettes/MSstats/inst/doc/MSstats.html

DIR = "C:/Projects/inProgress/OneDrive_Exile/Auto-prot-validation/MSstats/Auto-Prot"

out_dir = file.path(DIR, "output/compare-msstats"); dir.create(out_dir, recursive = T)

# Load the example
# MSstats DIA data are too small - only 2 samples:
# data("DIARawData", package = "MSstats")

# so we use an example from galaxy training: https://training.galaxyproject.org/training-material/topics/proteomics/tutorials/maxquant-msstats-dda-lfq/tutorial.html?utm_source=chatgpt.com
# Download
evidence_url <- "https://zenodo.org/record/4896554/files/MaxQuant_Evidence.tabular"
pg_url       <- "https://zenodo.org/record/4896554/files/MaxQuant_proteingroups.tabular"

# Read into R
evi_full <- read.delim(evidence_url, stringsAsFactors = FALSE)
pg_full  <- read.delim(pg_url, stringsAsFactors = FALSE) %>%
  filter( Potential.contaminant != "+",
          Reverse != "+",
          Only.identified.by.site != "+")

# convert 0 to NA
pg_full[pg_full == 0] <- NA

#--- Metadata ---
samples <- data.frame(
  sample = colnames(pg_full %>% select(starts_with("LFQ")) )
) %>% 
  mutate(
    group = gsub("LFQ\\.intensity\\.", "", sample) %>% sub("_.*$", "", .) %>% as.factor()
  )

#--- Filter: keep proteins with >=75% non-missing values per group ---
# indices per group
group_idx <- split(seq_len(ncol(pg_full %>% select(starts_with("LFQ")))), samples$group)

# keep if every group has ≥XX% non-missing
missing_threshold = 0.75
# to the columns of each group, apply a function
keep <- apply(pg_full %>% select(starts_with("LFQ")), 1, function(x) {
  
  # for each group in group_idx
  all(sapply(group_idx, function(idx) {
    
    mean(!is.na(x[idx])) >= missing_threshold
    
  }))
  
})


pg_full <- pg_full[keep, ]

pg_full %>% select(starts_with("LFQ"))

# evi = read.delim("https://raw.githubusercontent.com/MannLabs/alphapeptstats/main/testfiles/maxquant_go/evidence.txt", 
#                  stringsAsFactors = F)
# 
# pg = read.delim("https://raw.githubusercontent.com/MannLabs/alphapeptstats/main/testfiles/maxquant_go/proteinGroups.txt", 
#                  stringsAsFactors = F)


# 2) Build annotation: Raw.file must match evi$Raw.file exactly
# Columns required: Raw.file, Condition, BioReplicate, Run
msstats_meta <- evi_full %>%
  distinct(Raw.file) %>%
  arrange(Raw.file) %>%
  mutate(
    Condition = sub("_.*$", "", Raw.file),
    BioReplicate = ave(seq_along(Raw.file), Condition, FUN = seq_along),
    Run = 1
  )

# 3) Convert to MSstats format
quant <- MaxQtoMSstatsFormat(
  evidence = evi_full,
  annotation = msstats_meta,
  proteinGroups = pg_full,
  useUniquePeptide = TRUE,
  removeFewMeasure = TRUE,
  removeProtein_with1Feature = TRUE,
  summaryforMultipleRows = max)

# 4) Process to protein level
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


### write protein groups to file so they can be used with auto-prot
write.csv(pg_full,
          file.path(DIR, "input/data/proteindata.csv"),
          row.names = F)

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
          file.path(DIR, "input/data/metadata.csv"),
          row.names = F)


###########
## Limma ##
###########
# following : https://www.bioconductor.org/packages/devel/bioc/vignettes/proDA/inst/doc/data-import.html

# renv::install("bioc::limma", "bioc::MSnbase", "bioc::imputeLCMD")
library(limma)
library(MSnbase )
library(vsn)

lfq_cols <- grep("^LFQ\\.intensity\\.", colnames(pg_full))
lfq <- pg_full[, lfq_cols]
rownames(lfq) <- pg_full$Protein.IDs
lfq[lfq == 0] <- NA

#--- log2 transform ---
lfq_log <- log2(lfq)

#--- limma model ---
design <- model.matrix(~0 + samples$group)
colnames(design) <- levels(samples$group)

fit <- lmFit(lfq_log, design)
cont.matrix <- makeContrasts(metast - RDEB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, adjust="BH", number=Inf)
head(res)
write.csv(res, file.path(out_dir, "Rstudio-limma-out-log2.csv"))

### Now with normalisation and imputation

#--- normalise ---
# Convert to matrix
expr <- as.matrix(lfq_log)

# VSN
lfq_vsn <- vsn2(expr)

#--- Simple imputation ---
lfq_ms <- MSnSet(exprs = as.matrix(lfq_vsn))

# imputation
lfq_imp <- MSnbase::impute(lfq_ms, method = "bpca")

### refit the model
fit <- lmFit(lfq_imp, design)
cont.matrix <- makeContrasts(metast - RDEB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

res2 <- topTable(fit2, adjust="BH", number=Inf)
head(res2)
write.csv(res2, file.path(out_dir, "Rstudio-limma-out-imputed.csv"))
