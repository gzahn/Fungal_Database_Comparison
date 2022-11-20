# ASSIGN TAXONOMY ####

# load packages
library(dada2)
library(ShortRead)
library(Biostrings)
library(tidyverse)

# paths for databases
unite.ref <- "/uufs/chpc.utah.edu/common/home/u6033249/Data/databases/fungi/sh_general_release_dynamic_s_10.05.2021.fasta.gz"
unite.euk.ref <- "/uufs/chpc.utah.edu/common/home/u6033249/Data/databases/fungi/sh_general_release_dynamic_s_all_10.05.2021.fasta.gz"
unite.insd <- "/uufs/chpc.utah.edu/common/home/u6033249/Data/databases/fungi/UNITE_INSD.fasta.gz"

# find sequence table filepaths ####


directories <- readLines("/uufs/chpc.utah.edu/common/home/u6033249/Clayton_SRA/dir_paths.txt")
directory <- directories[1]
for(directory in directories){
  
  output.dir <- file.path(directory,"output")
  study.name <- directory %>% str_split("/") %>% map_chr(7) %>% unique()
  seqtab.path <- file.path(output.dir,paste0(study.name,"_seqtab_nochim.RDS"))

# load sequence table    
seqtab.nochim <- readRDS(seqtab.path)

# for each one, assign taxonomy in various ways

# UNITE Fungal
taxa_UNITE <- assignTaxonomy(seqtab.nochim, unite.ref, multithread=20, tryRC = TRUE)
# UNITE + Euk
taxa_UNITE_Euk <- assignTaxonomy(seqtab.nochim, unite.euk.ref, multithread=20, tryRC = TRUE)
# NCBI nr
taxa_UNITE_INSD <- assignTaxonomy(seqtab.nochim, unite.insd, multithread=20, tryRC = TRUE)

# save each
saveRDS(taxa_UNITE,file.path(output.dir,paste0(study.name,"_taxa_UNITE.RDS")))
saveRDS(taxa_UNITE_Euk,file.path(output.dir,paste0(study.name,"_taxa_UNITE_Euk.RDS")))
saveRDS(taxa_UNITE_INSD,file.path(output.dir,paste0(study.name,"_taxa_UNITE_INSD.RDS")))

}




# next script:
# import metadata and build a new version of study ps object for each taxonomy assignment