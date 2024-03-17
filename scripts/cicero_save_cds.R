#library(SeuratDisk)
library(HuMMuS)


########### FUNCTION WRAPPER ###########
run_cicero_wrapper <- function(
    hummus,
    atac_assay = "peaks",
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    window = 500000,
    number_cells_per_clusters = 50,
    sample_num = 100,
    seed = 2025,
    verbose = 1,
    threshold = 0,
    reduction_method = "UMAP"
    ) {

    # functions that need to be renamed
    int_elementMetadata <- SingleCellExperiment::int_elementMetadata
    counts <- SingleCellExperiment::counts

    # obtain chromosome sizes
    chromosome_sizes <- data.frame(V1 = genome@seqinfo@seqnames,
                                   V2 = genome@seqinfo@seqlengths)
#    chromosome_sizes <- chromosome_sizes[1:10, 1:2]
#    chromosome_sizes[1, 2] <- 4100000
#    chromosome_sizes[2, 2] <- 4100000
#    chromosome_sizes[3, 2] <- 4100000
#    chromosome_sizes[4, 2] <- 4100000
#    chromosome_sizes[5, 2] <- 4100000
#    chromosome_sizes[6, 2] <- 4100000
#    chromosome_sizes[7, 2] <- 4100000
#    chromosome_sizes[8, 2] <- 4100000
#    chromosome_sizes[9, 2] <- 4100000
#    chromosome_sizes[10,2] <- 4100000

    # Get scATAC-seq data
    scATAC <- as.matrix(hummus@assays[[atac_assay]]@counts)
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
      snakemake@output$cicero_cds)

    cicero <- cicero::run_cicero(
      cds = cicero_cds, # Infer peak-links
      genomic_coords = chromosome_sizes,
      window = window,             # Default = 5e+05
      silent = FALSE,             # Default = FALSE
      sample_num = sample_num) # Default = 100

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

    if (verbose > 0) {
      cat("\n", dim(peak_network)[1], "peak edges with a coaccess score >",
          threshold, "were found.\n")
    }



    # Return peak network including edges with positive coaccess score
    return(peak_network)
}

############# SCRIPT #############
atac_path <- snakemake@input$atac_path

atac <- t(read.table(atac_path, header = TRUE, row.names = 1))
atac[1:5, 1:5]
# Transform seurat in signac object
print("Add ATAC to Seurat object")
rna_X <- atac[1:2, ]
seurat_object <- SeuratObject::CreateSeuratObject(rna_X)
print(seurat_object)


seurat_object@assays$peaks <- Signac::CreateChromatinAssay(
    counts = atac,
    sep = c("_", "_")
    )
print(seurat_object)

# create the HuMMuS object
print("Create HuMMuS object")
hummus <- as(seurat_object, "hummus_object")

genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10

## Compute peak network
print("Create peak network")
cicero_results <- run_cicero_wrapper(
            hummus,
            window = snakemake@params$distance_threshold,
            atac_assay = "peaks",
            number_cells_per_clusters = snakemake@params$number_cells_per_clusters,
            genome = genome,
            threshold = -1,
            verbose = 1,
            )

write.table(cicero_results, snakemake@output$network)
#write.table(hummus@multilayer@multiplex$peaks@networks$peak_network_cicero, snakemake@output$network)