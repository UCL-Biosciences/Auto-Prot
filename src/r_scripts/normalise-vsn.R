
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
meanSdPlot_path = args[3]

print(input_file)
file.exists(input_file)
print(output_file)
file.exists(output_file)
print(meanSdPlot_path)
file.exists(meanSdPlot_path)

library(vsn)
library(tibble)

# Read input
data <- read.csv(input_file, check.names = FALSE)
rownames(data) <- data[[1]]
data <- data[, -1]

print("data loaded OK")

# Convert to matrix
expr <- as.matrix(data)

# VSN
print( "fitting vsn" )
vsn_fit <- vsn2(expr)
# Save mean-SD plot to PNG
### ### mean and sd shouldn't be correlated! ### ###
png(meanSdPlot_path, width = 800, height = 600)
meanSdPlot(vsn_fit)
dev.off()
vsn_data <- predict(vsn_fit, expr)

# Convert back to data.frame with protein IDs
vsn_df <- as_tibble(vsn_data, rownames = "genes")

print("saving output")
# Save output
write.csv(vsn_df, output_file, row.names = FALSE)
