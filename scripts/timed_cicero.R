################## 2nd part
################## 2nd part
################## 2nd part
################## 2nd part

run_cicero <- function(cds,
                       genomic_coords,
                       window = 500000,
                       silent = FALSE,
                       sample_num = 100) {
  # Check input
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }

  if (!silent) print("Starting Cicero")
  if (!silent) print("Calculating distance_parameter value")
  distance_parameters <- cicero::estimate_distance_parameter(cds, window=window,
                                  maxit = 100, sample_num = sample_num,
                                   distance_constraint = 250000,
                                   distance_parameter_convergence = 1e-22,
                                   genomic_coords = genomic_coords)


  # Calculate mean distance parameter
  mean_distance_parameter <- mean(unlist(distance_parameters))

  print(paste("Mean distance parameter: ", mean_distance_parameter))
  #write.csv(as.matrix(mean_distance_parameter), snakemake@params$distance_parameter)

  # Run models
  if (!silent) print("Running models")
  cicero_out <-
    cicero::generate_cicero_models(cds,
                           distance_parameter = mean_distance_parameter,
                           window = window,
                           genomic_coords = genomic_coords)

  # Assemble connections
  if (!silent) print("Assembling connections")
  all_cons <- cicero::assemble_connections(cicero_out, silent=silent)

  if (!silent) print("Done")
  all_cons
  }

# functions that need to be renamed
int_elementMetadata <- SingleCellExperiment::int_elementMetadata
counts <- SingleCellExperiment::counts

# obtain chromosome sizes
genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
chromosome_sizes <- data.frame(V1 = genome@seqinfo@seqnames,
                                V2 = genome@seqinfo@seqlengths)

# functions that need to be renamed
int_elementMetadata <- SingleCellExperiment::int_elementMetadata
counts <- SingleCellExperiment::counts

print('read cds rds')
cicero_cds <- readRDS(snakemake@input$cicero_cds_rds)
print('run cicero')
cicero <- run_cicero(
  cds = cicero_cds, # Infer peak-links
  genomic_coords = chromosome_sizes,
  window = snakemake@params$distance_threshold,             # Default = 5e+05
  silent = FALSE,             # Default = FALSE
  sample_num = snakemake@params$sample_num) # Default = 100

threshold <- -1
# Remove NAs, double edges, and edges with coaccess score <=0
# Check for coaccess = NA
if (length(which(is.na(cicero$coaccess))) > threshold) {
  cicero <- cicero[which(!is.na(cicero$coaccess)), ]  # Remove NAs
}
cicero$temp <- NA  # Helper column to check and remove double edges
my_cols <- which(as.character(cicero$Peak1) <= as.character(cicero$Peak2))
cicero$temp[my_cols] <- paste(cicero$Peak1[my_cols],
                              cicero$Peak2[my_cols],
                              sep = ";")

my_cols <- which(as.character(cicero$Peak1) > as.character(cicero$Peak2))
cicero$temp[my_cols] <- paste(cicero$Peak2[my_cols],
                              cicero$Peak1[my_cols],
                              sep = ";")

# Sort table according to temp-column (each entry appears double)
cicero <- cicero[with(cicero, order(temp, decreasing = TRUE)), ]
rownames(cicero) <- c(1:dim(cicero)[1])
A <- as.character(cicero$Peak1[seq(1, dim(cicero)[1], 2)])
Anum <- round(cicero$coaccess[seq(1, dim(cicero)[1], 2)], 10)
B <- as.character(cicero$Peak2[seq(2, dim(cicero)[1], 2)])
Bnum <- round(cicero$coaccess[seq(2, dim(cicero)[1], 2)], 10)
# length(which(A==B & Anum==Bnum)) 
# Each edge appears twice with same coaccess score (rounded to 10 digits)
cicero <- cicero[seq(1, dim(cicero)[1], 2), ] # Remove double edges
cicero$temp <- NULL # Remove helper column
cicero <- cicero[with(cicero, order(cicero$coaccess,
                                    decreasing = TRUE)), ]  # Sort
rownames(cicero) <- c(1:dim(cicero)[1])
cicero$Peak1 <- gsub("_", "-", cicero$Peak1)
# Peak names 2x"-" to match bipartites
cicero$Peak2 <- gsub("_", "-", cicero$Peak2)
# Peak names 2x"-" to match bipartites ? 2x"-" or 2x"_"

peak_network <- cicero[which(cicero$coaccess > threshold), ]
# Remove edges with coaccess score <= threshold

cat("\n", dim(peak_network)[1], "peak edges with a coaccess score >",
      threshold, "were found.\n")

write.table(peak_network, snakemake@output$network)
