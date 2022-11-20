# packages and setup ####
library(tidyverse)
library(phyloseq)
library(patchwork)
library(microbiome)
library(skimr)
source("./R/plot_bar2.R")
source("./R/palettes.R")

theme_set(theme_minimal())

# collect file paths 
path <- "./output/output_ps_data"
UNITE_list <- list.files(path,pattern = "UNITE_ps_object.RDS",full.names = TRUE)
EUK_list <- list.files(path,pattern = "UNITE_Euk_ps_object.RDS",full.names = TRUE)
StudyIDs <- UNITE_list %>% str_remove_all("./output/output_ps_data/") %>% str_remove_all("_UNITE_ps_object.RDS")


agreement_dfs_list <- list()
# make the rest of this into a for-loop to cycle through all studies
# save diversity values (discrepancies) in a list for each study

# set up lists for filling....
total_metadata_skim <- list()
total_metadata <- list()
alpha_comparisons <- list()
agreement_plots <- list()
kingdom_barplots <- list()
phylum_barplots <- list()
class_barplots <- list()

# clean up log file between debugging runs
file.remove("./output/06_compare_alpha_div.log")

for(i in seq_along(StudyIDs)){

# import ps objects
UNITE_ps <- readRDS(UNITE_list[i])
EUK_ps <- readRDS(EUK_list[i])
studyname <- UNITE_list[i] %>% 
  str_remove_all("./output/output_ps_data/") %>% 
  str_split("_") %>% 
  map_chr(1)

# print info to log file for debugging
sink("./output/06_compare_alpha_div.log",append = TRUE)
print(studyname)
print(UNITE_ps)
print(rank_names(UNITE_ps))
print(EUK_ps)
print(rank_names(EUK_ps))
sink(NULL)
# Look at available metadata
total_metadata_skim[[studyname]] <- UNITE_ps@sam_data %>% meta() %>% skim()
total_metadata[[studyname]] <- UNITE_ps@sam_data %>% meta()
# Remove non-fungi ####
nonfungikingdoms <- EUK_ps@tax_table[,1] %>% unique() != "k__Fungi"
nonfungikingdoms <- unique(EUK_ps@tax_table[,1])[which(nonfungikingdoms != FALSE | is.na(nonfungikingdoms))]
row.names(nonfungikingdoms) <- NULL
nonfungikingdoms <- nonfungikingdoms %>% as.data.frame() %>% as.matrix()
nonfungikingdoms <- nonfungikingdoms[,1]
EUK_ps_nf <- EUK_ps %>% 
  subset_taxa(!Kingdom %in% nonfungikingdoms) %>% 
  subset_taxa(!is.na(Kingdom))

# check to see if there are any tax assignments at all...
EUK_ps %>% tax_table()
taxtab_debug <- skim(tax_table(UNITE_ps)) %>% as.data.frame()
completed <- !taxtab_debug$complete_rate == 0
proceed <- as.numeric(completed) %>% sum()

if(proceed == 1){

# Stacked barplots ####

# kingdom-level
bpa <- UNITE_ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Kingdom") + ggtitle("UNITE") + theme(axis.text.x = element_blank()) +
  labs(y="Relative abundance") +
scale_fill_manual(values=pal.discrete)

bpb <- EUK_ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Kingdom") + ggtitle("UNITE+Euk") +
  labs(y="Relative abundance") +
  scale_fill_manual(values=pal.discrete) +
  theme(legend.title = element_blank(),
        axis.text.x = element_blank())

bpc <- EUK_ps_nf %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Kingdom") + ggtitle("UNITE+Euk (non-fungi removed)") +
  labs(y="Relative abundance") +
  scale_fill_manual(values=pal.discrete) +
  theme(legend.title = element_blank())


kingdom_barplot <- bpa / bpb / bpc
kingdom_barplots[[studyname]] <- kingdom_barplot + plot_annotation(title = studyname)
# Taxonomic differences ####

# pull tax tables out
UNITE_Tax <- tax_table(UNITE_ps)
EUK_Tax <- tax_table(EUK_ps)
# clean up names for easier reading
row.names(UNITE_Tax) <- NULL
row.names(EUK_Tax) <- NULL

# what's the first point of disagreement?
mat <- UNITE_Tax == EUK_Tax %>% 
  as.data.frame() %>% as.matrix()
mat[is.na(mat)] <- FALSE  
agreement_dfs <- as.data.frame(mat) %>% 
  mutate(ASV=paste0("ASV_",1:nrow(EUK_Tax))) %>% 
  pivot_longer(-ASV,names_to="Rank",values_to="Agreement") %>%
  mutate(Rank = factor(Rank,levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))

agreement_dfs_list[i] <- agreement_dfs



agreement_plots[[studyname]] <- as.data.frame(mat) %>% 
  mutate(ASV=paste0("ASV_",1:nrow(EUK_Tax))) %>% 
  pivot_longer(-ASV,names_to="Rank",values_to="Agreement") %>%
  mutate(Rank = factor(Rank,levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>% 
  ggplot(aes(x=Rank,fill=Agreement)) +
  geom_bar() +
  scale_fill_manual(values=pal.discrete[c(1,10)]) +
  labs(y="Proportional agreement between databases",x="Rank",fill="Agreement",
       title=studyname)

} else if(proceed > 1){
  # kingdom-level
  bpa <- UNITE_ps %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    plot_bar2(fill="Kingdom") + ggtitle("UNITE") + theme(axis.text.x = element_blank()) +
    labs(y="Relative abundance") +
    scale_fill_manual(values=pal.discrete)
  
  bpb <- EUK_ps %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    plot_bar2(fill="Kingdom") + ggtitle("UNITE+Euk") +
    labs(y="Relative abundance") +
    scale_fill_manual(values=pal.discrete) +
    theme(legend.title = element_blank(),axis.text.x = element_blank())
  
  bpc <- EUK_ps_nf %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    plot_bar2(fill="Kingdom") + ggtitle("UNITE+Euk (non-fungi removed)") +
    labs(y="Relative abundance") +
    scale_fill_manual(values=pal.discrete) +
    theme(legend.title = element_blank())
  
  
  kingdom_barplot <- bpa / bpb / bpc
  kingdom_barplots[[studyname]] <- kingdom_barplot + plot_annotation(title = studyname)
  # Taxonomic differences ####
  
  # pull tax tables out
  UNITE_Tax <- tax_table(UNITE_ps)
  EUK_Tax <- tax_table(EUK_ps)
  # clean up names for easier reading
  row.names(UNITE_Tax) <- NULL
  row.names(EUK_Tax) <- NULL
  
  # what's the first point of disagreement?
  mat <- UNITE_Tax == EUK_Tax %>% 
    as.data.frame() %>% as.matrix()
  mat[is.na(mat)] <- FALSE  
  agreement_dfs <- as.data.frame(mat) %>% 
    mutate(ASV=paste0("ASV_",1:nrow(EUK_Tax))) %>% 
    pivot_longer(-ASV,names_to="Rank",values_to="Agreement") %>%
    mutate(Rank = factor(Rank,levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
  
  agreement_plots[[studyname]] <- as.data.frame(mat) %>% 
    mutate(ASV=paste0("ASV_",1:nrow(EUK_Tax))) %>% 
    pivot_longer(-ASV,names_to="Rank",values_to="Agreement") %>%
    mutate(Rank = factor(Rank,levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>% 
    ggplot(aes(x=Rank,fill=Agreement)) +
    geom_bar() +
    scale_fill_manual(values=pal.discrete[c(1,10)]) +
    labs(y="Proportional agreement between databases",x="Rank",fill="Agreement",
         title=studyname)

    # phylum-level
  bp1 <- UNITE_ps %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Phylum") + ggtitle("UNITE") + theme(axis.text.x = element_blank()) +
    labs(y="Relative abundance") +
    scale_fill_manual(values=pal.discrete) +
    theme(axis.title.x = element_blank())
  
  bp2 <- EUK_ps %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    plot_bar2(fill="Phylum") + ggtitle("UNITE+Euk") + theme(axis.text.x = element_blank()) +
    labs(y="Relative abundance") +
    scale_fill_manual(values=pal.discrete) +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank())
  
  bp3 <- EUK_ps_nf %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    plot_bar2(fill="Phylum") + ggtitle("UNITE+Euk (non-fungi removed)") +
    labs(y="Relative abundance") +
    scale_fill_manual(values=pal.discrete) +
    theme(legend.title = element_blank())
  
  phylum_barplot <- bp1 / bp2 / bp3
  phylum_barplots[[studyname]] <- phylum_barplot + plot_annotation(title = studyname)
  phylum_barplot + plot_annotation(title = studyname)
  ggsave(paste0("./output/figs/",studyname,"_phylum_comparison_barplot.png"), dpi = 300, height = 8, width = 12)
  
  UNITE_div <- UNITE_ps %>% 
    tax_glom("Phylum") %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>%
    estimate_richness(measures = "Shannon") %>% 
    rename(UNITE=Shannon)
  EUK_div <- EUK_ps %>% 
    tax_glom("Phylum") %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>%
    estimate_richness(measures = "Shannon") %>% 
    rename(UNITE_Euk=Shannon)
  EUK_div_nf <- EUK_ps_nf %>% 
    tax_glom("Phylum") %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>%
    estimate_richness(measures = "Shannon") %>% 
    rename(UNITE_Euk_Fungi=Shannon)
  cbind(UNITE_div,EUK_div,EUK_div_nf) %>% 
    pivot_longer(dplyr::contains("UNITE"),
                 names_to="Database",
                 values_to="Shannon_div") %>% 
    ggplot(aes(x=Database,y=Shannon_div,fill=Database)) +
    geom_boxplot() +
    geom_jitter(alpha=.25,width = .1) +
    scale_fill_manual(values = pal.discrete) +
    labs(y="Shannon diversity",title=paste0(studyname,": Phylum-level alpha diversity"))
  
  phylum_shannon <- cbind(UNITE_div,EUK_div,EUK_div_nf) %>% 
    pivot_longer(dplyr::contains("UNITE"),
                 names_to="Database",
                 values_to="Shannon_div") %>% 
    mutate(Tax_level = "Phylum")
  alpha_comparisons[[studyname]] <- phylum_shannon
  
} else {
  total_metadata_skim[[studyname]] <- NA
  alpha_comparisons[[studyname]] <- NA
  agreement_plots[[studyname]] <- NA
  kingdom_barplots[[studyname]] <- NA
  phylum_barplots[[studyname]] <- NA
  class_barplots[[studyname]] <- NA
}
}


# clean up strays
rm(list=ls(pattern = "^bp"))
rm(list=ls(pattern = "^EUK"))
rm(list=ls(pattern = "^UNITE"))
rm(list=ls(pattern = "shannon"))
rm(mat)

# Look at taxonomic agreement plots ####
map(agreement_plots,plot)

for(i in names(agreement_plots)){
  agreement_plots[[i]]
  ggsave(file.path("./output/figs",paste0(i,"_agreement_plot.png")),
         dpi=300,height = 6,width = 10)
}

# p1 <- agreement_plots[[1]]$data

agreement_data_full <- list()
for(i in 1:length(agreement_plots)){
  
  agreement_data_full[[i]] <- 
    
  agreement_plots[[i]]$data %>% 
    group_by(Rank,Agreement) %>% 
    summarize(N=n()) %>% 
    mutate(TOTAL = sum(N)) %>%
    ungroup() %>% 
    filter(Agreement == FALSE) %>% 
    mutate(Percent_Agreement = (TOTAL - N) / TOTAL) %>% 
    mutate(Study = agreement_plots[[i]]$labels$title)
  
  
}


df <- agreement_data_full %>% 
  map(selectem) %>% 
  reduce(full_join)


df %>% 
  ggplot(aes(x=Rank,y=Percent_Agreement)) +
  geom_boxplot(fill=pal.discrete[1]) +
  labs(y="Percent agreement\n",
       x="\nRank")
ggsave("./output/figs/all_studies_agreement_plot.png",dpi=300,width = 8,height = 8)


agreement_dfs_list



# Look at alpha diversity for all studies ####
for(i in names(alpha_comparisons)){
  alpha_comparisons[[i]] <- alpha_comparisons[[i]] %>% 
    mutate(Accession=i)
}

full_comparisons <- reduce(alpha_comparisons,full_join)

full_comparisons %>% 
  ggplot(aes(x=Tax_level,y=Shannon_div,fill=Database)) +
  geom_boxplot() +
  facet_wrap(~Accession,nrow = 2,scales = 'free') +
  scale_fill_manual(values = pal.discrete) +
  labs(y="Shannon diversity",y="")

ggsave("./output/figs/all_studies_alpha_comparisons.png",
       height = 6,width = 12,dpi=300)

total_metadata_skim[[1]]
