library(tidyverse)
library(janitor)

# functions ####
read_clean <- function(x){read_csv(x) %>% clean_names()}

# load metadata files
meta_files <- list.files("./Data/study_metadata/SRA_Metadata", 
                         pattern = ".txt$",
                         full.names = TRUE)

meta <- map(meta_files,read_clean)

# find and remove empty metadata files
keepers <- map(meta,1) %>% map(is.null) %>% map(isFALSE) %>% unlist()
meta <- meta[keepers]
meta_files <- meta_files[keepers]


# make directories
dir_list <- meta %>% map("bio_project") %>% map(unique) %>% map(1) %>% unlist()
# for(i in dir_list){
#   dir.create(file.path("./Data/fastq",i))
# }

write_lines(dir_list,"./Data/dir_list")
# use sra-toolkit to download straight to the CHPC


# combine all metadata tables
to_char <- function(x){x %>% mutate_all(as.character)}
df <- map(meta,to_char) %>% 
  reduce(full_join)
glimpse(df)

# find prokaryote samples and remove them
find_prokaryote_samples <- function(x){grep("prokary|16S|bacter|mock",ignore.case = TRUE,x)}
bact_samples <- apply(df,2,find_prokaryote_samples) %>% unlist() %>% unique()
df <- df[-bact_samples,]
df <- df %>% arrange(bio_project)

# add metadata about each remaining sample from publications
meta <- readxl::read_xlsx("./Data/study_metadata/Articles for Dr. Zahn Publication.xlsx")
meta <- meta[meta$Accession_Number %>% str_squish() %>% grep("PRJNA",.),]
valid_accessions <- meta$Accession_Number %>% str_squish() %>% str_split(", ") %>% unlist()

# remove df rows that don't match metadata
df <- df[df$bio_project %in% valid_accessions,]

# add metadata from publication (curated caterories)
dict_accession <- meta$Accession_Number
dict_habitat <- meta$Habitat
dict_host <- meta$Host

df$habitat <- str_replace_all(df$bio_project, setNames(dict_habitat,dict_accession))
df$host <- str_replace_all(df$bio_project, setNames(dict_host,dict_accession))

# get file names
df$fwd_filename <- paste0(df$run,".fastq.gz")
df$fwd_filepath <- paste0("/scratch/general/lustre/Zahn/Clayton_SRA_Data/",df$bio_project,"/",df$fwd_filename)

# remove samples that have no file associated in the HPC
df <- df[file.exists(df$fwd_filepath),]


df_list <- df %>% group_by(bio_project) %>% group_split()
saveRDS(df_list,"./Data/study_metadata/df_list.RDS")

