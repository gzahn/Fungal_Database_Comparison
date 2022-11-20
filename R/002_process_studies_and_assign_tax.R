# SETUP ####
library(tidyverse)
library(dada2)
library(phyloseq)

unite.ref <- "/uufs/chpc.utah.edu/common/home/u6033249/Data/databases/fungi/sh_general_release_dynamic_16.10.2022.fasta.gz"
unite.euk.ref <- "/uufs/chpc.utah.edu/common/home/u6033249/Data/databases/fungi/sh_general_release_dynamic_all_16.10.2022.fasta.gz"


ps_list_notax <- readRDS("/scratch/general/lustre/Zahn/Clayton_SRA_Data/ps_list_notax.RDS")

# PREPARE DATA ####

# remove any faulty studies
ps_list_notax <- ps_list_notax[which(!is.na(ps_list_notax))]

# get accessions and label list elements
names(ps_list_notax) <- 
ps_list_notax %>% map(sample_data) %>% map("bio_project") %>% map(unique) %>% unlist()

# combine all studies into one ps object
full <- 
  ps_list_notax %>% 
  reduce(merge_phyloseq)

taxa <- c(
"CTTGGTCATTTAGAGGAAGTAAAAGTCATAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTAGAGAACGAGCTTTGCTCGATTATCTATATAACACATGTGAGAAGTGGCTTCGGCCTTATTACCTTTATGTTTGAATGTATCAGCATAACAAAAGAAAAACTTTCAACAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCCTGTCTCTTATACACATCTCCGAGCCCACGAGACCG",                                                                                                                                                                                                  
"CTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTCTCCGTAGGTGAACCTGCGGAGGGATCATTACCGAGCGAGGGCGGACGCCGCCCGACCTCCCACCCCGTGTCTGTCGGACTACCACGTTGCCTCGGGGGGCCCCCCGGCCCCCCGGCCGGCCACCGTTGACCCACCCTGCGTCCCCGCGCGTCGGAGTAAATGGGCTTGGACCCCAAATTGATGGAAAAACTTTCAACAACGGATCTCTTGGT"
)
test_set <- 
  full %>% 
  subset_samples(bio_project == "PRJNA316729")

# pull out OTU table
otu <- otu_table(full)

# ASSIGN TAXONOMY ####

# testing process first...
taxa_UNITE <- assignTaxonomy(test_set@otu_table@.Data[,1:10], unite.ref, multithread=16, tryRC = TRUE)
test_set@otu_table@.Data <- test_set@otu_table@.Data[,1:10]
phyloseq(test_set,tax_table(taxa_UNITE))

# UNITE Fungal
taxa_UNITE <- assignTaxonomy(otu, unite.ref, multithread=16, tryRC = TRUE)
ps_UNITE <- phyloseq(full,tax_table(taxa_UNITE))

# UNITE + Euk
taxa_UNITE_Euk <- assignTaxonomy(otu, unite.euk.ref, multithread=16, tryRC = TRUE)
ps_UNITE_Euk <- phyloseq(full,tax_table(taxa_UNITE_Euk))


saveRDS(ps_UNITE,"/scratch/general/lustre/Zahn/Clayton_SRA_Data/ps_UNITE.RDS")
saveRDS(ps_UNITE_Euk,"/scratch/general/lustre/Zahn/Clayton_SRA_Data/ps_UNITE_Euk.RDS")