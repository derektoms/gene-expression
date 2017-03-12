# 2017-02-24
# Trying out getting and analyzing the date from GEO

if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()

db = src_sqlite("GEOmetadb.sqlite")
src_tbls(db)
gse = tbl(db, 'gse')
gse_gpl = tbl(db, 'gse_gpl')
gpl = tbl(db, 'gpl') 
gsm = tbl(db, 'gsm')

# Equivalent human platform is GPL570 with 127 514 samples
# HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
biocLite("hgu133plus2.db")
gse_gsm = tbl(db, 'gse_gsm')

# diabetes
search_string = "\\bpancreas\\b|\\bislets?\\b|\\bbeta-cells\\b|\\bdiabetes\\b"
# search_string = "\\bpancreas?(tic)?\\b&(\\bislets?\\b|\\bbeta-cells?\\b|\\bdiabetes\\b)"
pancreatic_gpl570 = gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  filter(gpl == "GPL570") %>% 
  collect() %>% 
  filter(str_detect(title, search_string))
write_csv(pancreatic_gpl570, 'pancreatic_gpl570.csv')

biocLite(c('MergeMaid','GEOquery','inSilicoMerging','affy','sva','Rtsne','metaArray','testthat'))
library(MergeMaid)
library(GEOquery)
library(readr)
library(testthat)
library(dplyr)
library(metaArray)
library(inSilicoMerging) ## problem here
library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(sva)
library(limma)
library(affy)
library(stringr)

# gpl1261_metadata = read_csv("retinal_gpl1261_with_tissue.csv")

if(!file.exists("gpl570_datasets.rda")) { 
  gse = lapply(pancreatic_gpl570$gse, getGEO)
  gse = lapply(gse, "[[", 1)
  names(gse) = pancreatic_gpl570$gse
  save(gse, file = "gpl570_datasets.rda")
} else {
  load("gpl570_datasets.rda")
}

pancreatic_gpl570_metadata = read_csv("pancreatic_gpl570_metadata.csv")

# UI connects to a,b,c
a<-'islets'
b<-'pancreas'
c<-'beta-cells'

keep_types <- c(a,b,c)
includes = pancreatic_gpl570_metadata$gse[pancreatic_gpl570_metadata$include == "yes"] %>% na.omit()
tissue_keeps = pancreatic_gpl570_metadata$gse[pancreatic_gpl570_metadata$tissue %in% keep_types]
gse_to_keep = intersect(includes, tissue_keeps)

#get_files = TRUE
get_files = FALSE
if (get_files) {
  # get raw CEL files
  rawFilePaths = lapply(gse_to_keep, function(x) {
    if (!dir.exists(x)) {
      getGEOSuppFiles(x)
    }
  })
  
  tar_dirs = list.files(pattern = "GSE")
  tar_files = lapply(tar_dirs, list.files, pattern = ".tar", full.names = TRUE)
  
  # these are downloaded as tar archives so need to extract
  lapply(tar_files, function(x) {
    untar(x, exdir = dirname(x)) 
  })
  
  # create a vector of file names for each study
  affy_files = lapply(tar_files, function(x) {
    unlist(list.files(dirname(x), pattern = "[Cc][Ee][Ll].gz", full.names = TRUE))
  })
  
  # load up all the array data
  affy_batches = lapply(affy_files, function(x) ReadAffy(filenames = x))
  
  # convert to expression values using the RMA method from the affy package
  # this will correct the background and normalize
  esets = lapply(affy_batches, rma)
  names(esets) = tar_dirs
  
  save(esets, file = "esets_2016-02-12.rda")
} else {
  load("esets_2016-02-12.rda")
}