# Process ITS reads #

# load packages ####
library(dada2)
library(ShortRead)
library(Biostrings)
library(tidyverse)


# for-loop to build seqtab for each SRA accession ####
directories <- readLines("/uufs/chpc.utah.edu/common/home/u6033249/Clayton_SRA/dir_paths.txt")

for(directory in directories){
  
  output.dir <- file.path(directory,"output")
  dir.create(output.dir) # make output directory
  
  # find files ####
  path <- directory
  fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
  dir.create(output.dir) # make output directory
  
  # create filtered file paths
  filtFs <- fnFs %>% str_replace_all("_1.fastq.gz","_1_trimmed.fastq.gz")
  
  # get.sample.name <- function(fname){strsplit(basename(fname), "_")[[1]][1]}
  study.name <- fnFs %>% str_split("/") %>% map_chr(7) %>% unique()  
  sample.names <- fnFs %>% str_split("/") %>% map_chr(8) %>% str_split("_") %>% map_chr(1)
 
  
  print(path)
  print(output.dir)
  print(study.name)
  print(sample.names)
  print(file.path(output.dir,paste0(study.name,"_seqtab_nochim.RDS")))
}