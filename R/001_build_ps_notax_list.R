
# Run dada2 and build seq table for each study
library(tidyverse)
library(dada2)
library(phyloseq)
df_list <- readRDS("/uufs/chpc.utah.edu/common/home/u6033249/Clayton_SRA/study_metadata/df_list.RDS")

# function to build phyloseq object for each study (without taxonomy)
build_seqtab <- function(x){
  samplenames <- x$fwd_filename %>% str_remove(".fastq.gz")
  filts <- x$fwd_filepath %>% str_replace(".fastq.gz","_filt.fastq.gz")

# remove missing files
  fwds <- x$fwd_filepath[file.exists(x$fwd_filepath)]
  filts <- filts[file.exists(x$fwd_filepath)]
  x <- x[file.exists(x$fwd_filepath), ]

  out <- filterAndTrim(x$fwd_filepath, filts, # fnRs, filtRs,
                       maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=FALSE)
  filts <- sort(list.files(dirname(filts) %>% unique(), full.names = TRUE,pattern="_filt.fastq.gz"))
  errF <- learnErrors(filts, multithread=TRUE, MAX_CONSIST = 20)
  derep <- derepFastq(filts, verbose=FALSE)
  if(length(derep) != length(samplenames)){
    samplenames <- unlist(map(strsplit(basename(filts), "_filt"), 1))
  }
  names(derep) <- samplenames
  dadaFs <- dada(derep, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")
  seqtab <- makeSequenceTable(dadaFs)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  x <- x[x$run %in% row.names(seqtab.nochim),] # remove metadata for missing files
  meta <- data.frame(bio_project=x$bio_project,
                     run=x$run,
                     host=x$host,
                     habitat=x$habitat,
                     filename=x$fwd_filename)
  
  met <- sample_data(meta)
  sample_names(met) <- meta$run
  ps <- phyloseq(otu_table(seqtab.nochim,taxa_are_rows = FALSE),
                 met)
  return(ps)  
}


# apply function to list of studies

ps_list_notax <- map(df_list,build_seqtab)

# save output (list of ps objects without taxonomy)
saveRDS(ps_list_notax,"/scratch/general/lustre/Zahn/Clayton_SRA_Data/ps_list_notax.RDS")

