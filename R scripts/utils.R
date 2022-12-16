###############
## Functions ##
###############

message("load data preprocessing functions...")

#################################
## Mapping metadata to objects ##
#################################
mapping_vec_function <- function(reik, metadata) {
  n <- length(reik$barcode)
  mapping_vec <- rep(NA, n)

  for(i in 1:n){
    if(i %% floor(n/10) == 0) cat('*')
    tmp <- intersect(which(metadata$barcode == reik$barcode[i]),
                     which(metadata$sample == reik$sample[i]))

    if(length(tmp) > 1) stop()
    if(length(tmp) > 0) mapping_vec[i] <- tmp
  }

  return(mapping_vec)
}

###############################################
## Converting row names in scATAC-seq object ##
###############################################
# Adjust rownames of dataset: Convert chr-xxxx-xxxx to chr:xxxx-xxxx
seuat_adjust_row_names <- function(seurat_obj) {
  # Set the default assay as ATAC
  SetAssayData(object = seurat_obj, assay = "ATAC")

  # Convert the row names of the Seurat object to a data frame
  DF <- as.data.frame(row.names(seurat_obj))

  # Define a function to adjust the row names by replacing hyphens with colons
  adjust_row_names <- function(x) {
    sub("-", ":", x)
  }

  # Apply the adjustment function to the row names
  adjusted <- sapply(DF, adjust_row_names)

  # Set the adjusted row names as the row names for the ATAC assay counts
  row.names(seurat_obj@assays$ATAC@counts) <- adjusted

  # Return the modified Seurat object
  rm(DF)
  rm(adjusted)
  return(seurat_obj)
}


##########################################################################
## Binarise the scATAC-seq count matrix and save the matrix for SCENIC+ ##
##########################################################################
# Function to set all non-zero values in the ATAC assay counts to 1, convert the counts
# to a data frame, add a column containing the row names, and write the data frame to a feather file
process_atac_assay_counts <- function(seurat_obj, output_file,metadata_output_file) {

  # Check if the output file starts with a /
  if (!startsWith(output_file, "/") || !startsWith(metadata_output_file,"/")) {
    stop("The output file must start with a /")
  }

  # Check if the output file is a string
  if (!is.character(output_file) || !is.character(metadata_output_file)) {
    stop("The output file must be a string")
  }

  # Set all non-zero values in the ATAC assay counts to 1
  seurat_obj@assays$ATAC@counts[seurat_obj@assays$ATAC@counts > 0] <- 1

  # Convert the ATAC assay counts to a data frame
  ct <- as.data.frame(seurat_obj@assays$ATAC@counts)

  # Add a column containing the row names
  ct$NAMES <- row.names(ct)

  # Write the data frame to a feather file
  write_feather(ct, paste0=directory$objects,output_file)

  # remove the count matrix object
  rm(ct)

  # collect metadata and write out as a csv
  seurat_obj<-seurat_obj@meta.data
  write.csv(seurat_obj,file=paste0(directory$objects,metadata_output_file),row.names=FALSE)

  # Return the modified Seurat object
  print("Created the count matrix feather object and metadata.csv files for the processed data")
  return(seurat_obj)
}
