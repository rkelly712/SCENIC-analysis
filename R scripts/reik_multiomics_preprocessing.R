###########################
## SCENIC+ preprocessing ##
###########################

###################
## Load datasets ##
###################

source("/data/kellyrc/SCENIC+/scripts/r_script/settings.r")
rm(list=ls())

reik <- readRDS(paste0(directory$rna,"/seurat.rds"))
metadata <- read.csv(paste0(directory$rna"/sample_metadata.txt"),
                     sep = "\t")

#############################
## Generate mapping vector ##
#############################

n <- ncol(reik)
mapping_vec <- rep(NA, n)
for(i in 1:n){
  if(i %% floor(n/10) == 0) cat('*')
  tmp <- intersect(which(metadata$barcode == reik$barcode[i]),
                   which(metadata$sample == reik$sample[i]))

  if(length(tmp) > 1) stop()
  if(length(tmp) > 0) mapping_vec[i] <- tmp
}

reik$stage <- metadata$stage[mapping_vec]
reik$genotype <- metadata$genotype[mapping_vec]
reik$pass_rnaQC <- metadata$pass_rnaQC[mapping_vec]
reik$doublet_score <- metadata$doublet_score[mapping_vec]
reik$doublet_call <- metadata$doublet_call[mapping_vec]
reik$celltype <- metadata$celltype[mapping_vec]
reik$celltype.score <- metadata$celltype.score[mapping_vec]
reik$TSSEnrichment_atac <- metadata$TSSEnrichment_atac[mapping_vec]
reik$ReadsInTSS_atac <- metadata$ReadsInTSS_atac[mapping_vec]
reik$PromoterRatio_atac <- metadata$PromoterRatio_atac[mapping_vec]
reik$NucleosomeRatio_atac <- metadata$NucleosomeRatio_atac[mapping_vec]
reik$BlacklistRatio_atac <- metadata$BlacklistRatio_atac[mapping_vec]
reik$pass_atacQC <- metadata$pass_atacQC[mapping_vec]
reik$celltype.predicted <- metadata$celltype.predicted[mapping_vec]
Seurat::DefaultAssay(reik) <- "RNA"
reik

############################
## Load scATAC-seq object ##
############################

reik_atac <- ArchR::loadArchRProject(paste0(directory$atac, "/archR/"))
ArchR::getAvailableMatrices(reik_atac)
reik_atac <- ArchR::addPeakMatrix(reik_atac, binarize = FALSE, force = TRUE)# this takes ages
# This is a SummarizedExperiment object

se_obj <- ArchR::getMatrixFromProject(ArchRProj = reik_atac,
                                      useMatrix = "PeakMatrix") # takes around 25 minutes

peak_granges <- se_obj@colData
n <- ncol(reik)
mapping_vec <- rep(NA, n)
for(i in 1:n){
  if(i %% floor(n/10) == 0) cat('*')
  tmp <- intersect(which(peak_granges$barcode == reik$barcode[i]),
                   which(peak_granges$sample == reik$sample[i]))

  if(length(tmp) > 1) stop()
  if(length(tmp) > 0) mapping_vec[i] <- tmp
}

keep_vec <- rep(1, ncol(reik))
keep_vec[which(is.na(mapping_vec))] <- 0
reik$keep <- keep_vec
reik <- subset(reik, keep == 1)

mapping_vec2 <- mapping_vec[which(!is.na(mapping_vec))]
mat <- se_obj@assays@data$PeakMatrix[,mapping_vec2]
chromosome_vec <- as.character(se_obj@rowRanges@seqnames)
range_vec <- as.character(se_obj@rowRanges@ranges)
p <- length(chromosome_vec)
row_vec <- sapply(1:p, function(j){
  paste0(chromosome_vec[j], "_", range_vec[j])
})
rownames(mat) <- row_vec
reik[["ATAC"]] <- Seurat::CreateAssayObject(counts = mat)

##########
## Save ##
##########

# main object
SaveH5Seurat(reik, filename = paste0(directory$objects, "/reik_multiomics.h5Seurat"), overwrite = TRUE)

# Subsetted objects for WT and KO
reik_WT<- subset(x = reik, subset = sample == "E8.5_CRISPR_T_WT")
reik_KO<- subset(x = reik, subset = sample == "E8.5_CRISPR_T_KO")
SaveH5Seurat(reik_WT, filename = paste0(directory$objects, "/reik_WT.h5Seurat"), overwrite = TRUE)
SaveH5Seurat(reik_KO, filename = paste0(directory$objects,"/reik_KO.h5Seurat"), overwrite = TRUE)

###########################################
## Convert from R to python for analysis ##
###########################################


message("Converting rownames for WT...")
DF<-as.data.frame(row.names(reik_WT))
colnames(DF)<-"Name"
adjusted<-sapply(DF,function(x) sub('-', ':', DF$Name))
row.names(reik_WT@assays$ATAC@counts)<-adjusted

print("Converting rownames for KO...")
DF<-as.data.frame(row.names(reik_KO))
colnames(DF)<-"Name"
adjusted<-sapply(DF,function(x) sub('-', ':', DF$Name))
row.names(reik_KO@assays$ATAC@counts)<-adjusted

rm(DF)
rm(adjusted)

# Binarise the matrix
reik_WT@assays$ATAC@counts[reik_WT@assays$ATAC@counts>0]<-1
reik_KO@assays$ATAC@counts[reik_KO@assays$ATAC@counts>0]<-1

# Save count matrices
ct <- reik_WT@assays$ATAC@counts
ct <- as.data.frame(ct)
ct$NAMES <- row.names(ct)
write_feather(ct, paste0(directory$object,"/WT_count_matrix.feather"))
rm(ct)

ct <- reik_KO@assays$ATAC@counts
ct <- as.data.frame(ct)
ct$NAMES <- row.names(ct)
write_feather(ct, paste0(directory$object,"/KO_count_matrix.feather"))

rm(ct)

###################
## Save metadata ##
###################

WT_meta<-reik_WT@meta.data
KO_meta<-reik_KO@meta.data
WT_meta$NAMES <- row.names(WT_meta)
KO_meta$NAMES <- row.names(KO_meta)

write.csv(WT_meta,file=paste0(directory$object,"/WT_celldata.csv"),row.names=FALSE)
write.csv(WT_meta,file=paste0(directory$object,"/KO_celldata.csv"),row.names=FALSE)


