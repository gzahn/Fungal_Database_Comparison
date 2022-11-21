# Build phyloseq objects ####
library(tidyverse)
library(phyloseq)
library(readxl)
library(dada2)


# parse filenames, make sure they match ####
path <- "./output/output_ps_data"
seqtables <- list.files(path, pattern = "seqtab_nochim.RDS",full.names = TRUE)
metatables <- list.files(path = "./Data/study_metadata",pattern = "SraRunTable",full.names = TRUE)
seqtable_ids <- basename(seqtables) %>% str_remove_all("_seqtab_nochim.RDS")
metatable_ids <- basename(metatables) %>% str_remove_all("_SraRunTable.txt.txt")
metatables <- metatables[metatable_ids %in% seqtable_ids]
metatable_ids <- metatable_ids[metatable_ids %in% seqtable_ids]
UNITE_tax_tables <- list.files(path = "./output/output_ps_data",pattern = "UNITE.RDS",full.names = TRUE)
UNITE_Euk_tax_Tables <- list.files(path = "./output/output_ps_data",pattern = "UNITE_Euk.RDS",full.names = TRUE)

# read in metadata ####
meta <- list()
for(i in seq_along(seqtables)){
  meta[[i]] <- read_delim(metatables[i]) %>% as.data.frame()
}
names(meta) <- metatable_ids

# change "Run" to "SampleID" and rename rows by this value in each metadata df
for(i in seqtable_ids){
  names(meta[[i]]) <- names(meta[[i]]) %>% str_replace("Run","SampleID")
  row.names(meta[[i]]) <- meta[[i]]$SampleID
}



# read in sequence tables ####
seqtabs <- list()
for(i in seq_along(metatable_ids)){
  seqtabs[[i]] <- readRDS(seqtables[i])  
}
names(seqtabs) <- metatable_ids


# read in taxonomy ####
unite_tax <- list()
for(i in seq_along(metatable_ids)){
  unite_tax[[i]] <- readRDS(UNITE_tax_tables[i])
}
names(unite_tax) <- metatable_ids

unite_euk_tax <- list()
for(i in seq_along(metatable_ids)){
  unite_euk_tax[[i]] <- readRDS(UNITE_Euk_tax_Tables[i])
}
names(unite_euk_tax) <- metatable_ids


# combine for standard UNITE phyloseq objects and write to disk

for(i in metatable_ids){
  x <- phyloseq(otu_table(seqtabs[[i]],taxa_are_rows = FALSE),
           sample_data(meta[[i]]),
           tax_table(unite_tax[[i]]))  
  print(i)
  print(x)
  ps_name <- file.path(path,paste0(i,"_UNITE_ps_object.RDS"))
  saveRDS(x,ps_name)
}

for(i in metatable_ids){
  x <- phyloseq(otu_table(seqtabs[[i]],taxa_are_rows = FALSE),
                sample_data(meta[[i]]),
                tax_table(unite_euk_tax[[i]]))  
  print(i)
  print(x)
  ps_name <- file.path(path,paste0(i,"_UNITE_Euk_ps_object.RDS"))
  saveRDS(x,ps_name)
}


