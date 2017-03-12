# 2017-02-23
# working through the required data collection and processing; took care of "get_geo_info.Rmd"
# will look into processing the data tomorrow


# Setting up the workspace:
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("limma","annotate", "mouse4302.db","genefilter","ComplexHeatmap","pheatmap","cowplot","GEOmetadb","hgu133plus2.db")
library(limma) #via biocLite
library(annotate) # via biocLite
library(mouse4302.db) #via biocLite
library(genefilter) #via biocLite
library(ComplexHeatmap) #via biocLite
library(pheatmap) #via biocLite
library(cowplot) #via biocLite
library(GEOmetadb) #via biocLite

#install.packages(c('dplyr','tidyr','ggplot2','RColorBrewer','readr','stringr','shiny','shinythemes','shinyjs'))

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)

library(shiny)
library(shinythemes)
library(shinyjs)


# Here is a copy of the code used to search the GEO database to find the relevant datasets.  Note that resulting GEO datasets were manually filtered and labelled according to tissue type prior to analysis.

# get the latest GEO sqlite file (this is gigantic, ~350 MB)
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()

db = src_sqlite("GEOmetadb.sqlite")
src_tbls(db)
gse = tbl(db, 'gse')
gse_gpl = tbl(db, 'gse_gpl')
gpl = tbl(db, 'gpl') 
gsm = tbl(db, 'gsm')

# Search for specific platform: GPL1261 with 50 437 samples
# [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array
search_string = "\\beye\\b|\\bretinas?l?\\b|\\bcones?\\b|\\brods?\\b|\\bphotoreceptors?\\b"

retinal_gpl1261 = gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  filter(gpl == "GPL1261") %>% 
  collect() %>% 
  filter(str_detect(title, search_string))


all_retinal = gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  left_join(gpl %>% select(gpl, organism), copy = TRUE) %>% collect() %>% 
  filter(organism == "Mus musculus" ) %>% 
  filter(str_detect(title, search_string))

write_csv(retinal_gpl1261, "retinal_gpl1261.csv")
write_csv(all_retinal, "all_eye.csv")
write_csv(gse_info, "gse_info.csv")

### New stuff:
# Equivalent human platform is GPL570 with 127 514 samples
# HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
biocLite("hgu133plus2.db")
gse_gsm = tbl(db, 'gse_gsm')

# diabetes
search_string = "\\bpancreas\\b|\\bislets?\\b|\\bbeta\\b|\\bdiabetes\\b"

pancreatic_gpl1261 = gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  filter(gpl == "GPL1261") %>% 
  collect() %>% 
  filter(str_detect(title, search_string))
write_csv(pancreatic_gpl1261, "pancreatic_gpl1261_03.csv")

# diabetes_2
search_string = "\\bpancreas?(tic)?\\b|\\bislets?\\b|\\bbeta-cells?\\b|\\bdiabetes\\b"
pancreatic_gpl570 = gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  filter(gpl == "GPL570") %>% 
  collect() %>% 
  filter(str_detect(title, search_string))
write_csv(pancreatic_gpl570, "pancreatic_gpl570.csv")

# human retina
retinal_gpl570 = gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  filter(gpl == "GPL570") %>% 
  collect() %>% 
  filter(str_detect(title, search_string))
write_csv(retinal_gpl570, "retinal_gpl570.csv")