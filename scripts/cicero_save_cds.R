

reduction_method <- "UMAP"
seed <- 2
number_cells_per_clusters <- snakemake@params$number_cells_per_clusters

############# SCRIPT #############
atac_path <- snakemake@input$atac_path

atac <- t(read.table(atac_path, header = TRUE, row.names = 1))
atac[1:5, 1:5]


# Get scATAC-seq data
scATAC <- as.matrix(atac)
rownames(scATAC) <- gsub("_", "-", rownames(scATAC))

# Matrix to edgelist
acc <- reshape2::melt(scATAC)
colnames(acc) <- c("V1", "V2", "V3")

# Prepare cicero input
input_cds <- cicero::make_atac_cds(acc, binarize = TRUE) # Create CDS object
set.seed(seed)
# It is required that there is no empty cell
if (length(which(colSums(as.matrix(monocle3::exprs(input_cds))) == 0)) == 0
) {
# Calculating size factors using default method = mean-geometric-mean-total
  input_cds <- monocle3::estimate_size_factors(input_cds)
  # Preprocessing using LSI
  input_cds <- monocle3::preprocess_cds(input_cds, method = "LSI")
  # Dimensionality reduction using UMAP
  input_cds <- monocle3::reduce_dimension(
                                input_cds,
                                reduction_method = reduction_method,
                                preprocess_method = "LSI")
} else {
  print("Error: there is at least one cell with no signal.")
}
# Get reduced (UMAP) coordinates
umap_coords <- SingleCellExperiment::reducedDims(input_cds)$UMAP
# Compute pseudocells
cicero_cds <- cicero::make_cicero_cds(
  input_cds,  # Create a Cicero CDS object
  reduced_coordinates = umap_coords,
  k = number_cells_per_clusters,  #number neighbors/ Default = 50
  summary_stats = NULL,         # Default
  size_factor_normalize = TRUE, # Default
  silent = FALSE)               # Default


write.table(
  as.matrix(cicero_cds@assays@data$counts),
  snakemake@output$cicero_cds_tsv)

saveRDS(cicero_cds, snakemake@output$cicero_cds_rds)