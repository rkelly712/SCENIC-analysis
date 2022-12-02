# Load packages, set seed and directory
if (!require("pacman")) install.packages("pacman")


message("Loading Packages")

pacman::p_load(tibble,arrow,SeuratDisk,SeuratData,cisTopic,Seurat,Signac,ArchR,tiltedCCA, install=TRUE)

print("Setting directories")

directory <- list()
if (grepl("FNM-S105936",Sys.info()['nodename'])) {
  directory$basedir <- "/data/kellyrc/SCENIC+"
  directory$rna <- "/data/kellyrc/SCENIC+/rna"
  directory$atac <-"/data/kellyrc/SCENIC+/atac"
  directory$output<-"/data/kellyrc/SCENIC+/output"
  directory$r_scripts<-"/data/kellyrc/SCENIC+/scripts/r_scripts"
  directory$bash_scripts<-"/data/kellyrc/SCENIC+/scripts/bash_scripts"
  directory$objects<-"/data/kellyrc/SCENIC+/objects"
  print("Home Computer")
  setwd(directory$basedir)
  paste0("Working Directory:",getwd())
} else if (grepl("*cn*", Sys.info()['nodename'])) {
  if (grepl("kellyrc", Sys.info()['effective_user'])) {
  directory$basedir <- "/data/kellyrc/SCENIC+"
  directory$rna <- "/data/kellyrc/SCENIC+/rna"
  directory$atac <-"/data/kellyrc/SCENIC+/atac"
  directory$output<-"/data/kellyrc/SCENIC+/output"
  directory$r_scripts<-"/data/kellyrc/SCENIC+/scripts/r_scripts"
  directory$bash_scripts<-"/data/kellyrc/SCENIC+/scripts/bash_scripts"
  directory$objects<-"/data/kellyrc/SCENIC+/objects"
  setwd(directory$basedir)
  print("You are working in the Biowulf Cluster")
  paste0("Working Directory:", getwd())
  }

} else {
  stop("Computer/node not recognised. Please check the naming convention in the settings.r file")
}
