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

        
    # filter and trim ####
    out <- filterAndTrim(fnFs, filtFs, maxN = 0, maxEE = c(2, 2), 
                         truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
      
    # learn errors ####
    errF <- learnErrors(filtFs, multithread = TRUE)
    errorplot <- plotErrors(errF, nominalQ = TRUE)
    ggsave(file.path(output.dir,paste0(study.name,"_dada2_errorplot.png")),
           plot = errorplot,
           device = "png")
      
    # de-replicate
    derepFs <- derepFastq(filtFs, verbose = TRUE)
      
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
      
    # sample inferrence ####
    dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
      
    # save intermediate output
    saveRDS(dadaFs,file.path(output.dir,paste0(study.name,"_dada.RDS")))
      
    # make sequence table ####
    seqtab <- makeSequenceTable(dadaFs) # FWD seqs only
      
    # remove chimeras ####
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
      
    # track reads through pipeline ####
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoised", "nonchim")
    rownames(track) <- sample.names
      
    # save sequence table and tracked reads
    saveRDS(seqtab.nochim,file.path(output.dir,paste0(study.name,"_seqtab_nochim.RDS")))
    write_csv(track,file.path(output.dir,paste0(study.name,"_trackedreads.csv")))
  
}



