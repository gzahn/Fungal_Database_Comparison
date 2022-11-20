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

UNITES <- map(UNITE_list,readRDS)
EUKS <- map(EUK_list,readRDS)

full_UNITE_ps <- reduce(UNITES,merge_phyloseq)
full_EUK_ps <- reduce(EUKS,merge_phyloseq)



# alpha diversity ####
# remove non-fungi and estimate alpha div
UNITE_alpha <- 
  full_UNITE_ps %>% 
  subset_taxa(Kingdom == "k__Fungi") %>% 
  estimate_richness(measures = c("Observed","Shannon"))

EUK_alpha <- 
  full_EUK_ps %>% 
  subset_taxa(Kingdom == "k__Fungi") %>% 
  estimate_richness(measures = c("Observed","Shannon"))

UNITE_alpha$Database <- "UNITE"
UNITE_alpha$Study <- full_UNITE_ps@sam_data$BioProject
EUK_alpha$Database <- "UNITE+EUK"
EUK_alpha$Study <- full_EUK_ps@sam_data$BioProject

full_alpha <- 
UNITE_alpha %>% 
  full_join(EUK_alpha) %>% 
  pivot_longer(c("Observed","Shannon"),names_to = "Measure",values_to = "Value")
head(full_alpha)


# plot alpha div
full_alpha %>% 
  filter(Measure == "Observed") %>% 
  ggplot(aes(x=Database,y=Value,fill=Database)) +
  geom_boxplot(position = "dodge") +
  facet_wrap(~Study,scales = 'free_y') +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(face='bold',size=14),
        strip.text = element_text(face='bold',size=14),
        legend.title = element_text(face='bold',size=14)) +
  labs(y="Richness") +
  scale_fill_manual(values = c("Orange","DarkGreen"))
ggsave("./output/manuscript_figs/richness_difference_by_study.png",dpi=300,height = 8,width = 12)

full_alpha %>% 
  filter(Measure == "Shannon") %>% 
  ggplot(aes(x=Database,y=Value,fill=Database)) +
  geom_boxplot(position = "dodge") +
  facet_wrap(~Study,scales = 'free_y') +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(face='bold',size=14),
        strip.text = element_text(face='bold',size=14),
        legend.title = element_text(face='bold',size=14)) +
  labs(y="Shannon diversity") +
  scale_fill_manual(values = c("Orange","DarkGreen"))
ggsave("./output/manuscript_figs/shannon_difference_by_study.png",dpi=300,height = 8,width = 12)


blast_tax <- 
  read_delim("./taxonomy/All_Seqs_Top_BLAST_Hit_UNITE_EUK.txt",col_names = FALSE) %>% 
  select(X1,X2) %>% 
  rename("ASV"="X1","BLAST"="X2")

blast_tax$BLAST

seqs <- euk_tax %>% row.names()
names(seqs) <- paste0("ASV_",seq_along(seqs))


df <- data.frame(ASV = names(seqs),
           Sequence = seqs) %>% 
  full_join(blast_tax)

df$BLAST %>% is.na() %>% sum()
df$ASV

df <- 
df %>% 
  unique.data.frame() %>% 
  mutate(
    UNITE_taxonomy = unite_tax %>% 
      mutate(UNITE_taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus,Species,sep=";")) %>% 
      pluck("UNITE_taxonomy"),
    UNITE_EUK_taxonomy = euk_tax %>% 
      mutate(EUK_taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus,Species,sep=";")) %>% 
      pluck("EUK_taxonomy")
  )

df[6954,]

nf <- 
df$BLAST %>% grep(pattern = "k__Fungi",invert = TRUE)


df[nf,"UNITE_taxonomy"] %>% str_split(";") %>% map_chr(1) %>% table()
bl <- df[nf,"BLAST"]
bl[is.na(bl)] <- "NA|NA|NA|NA|NA;NA;NA;NA;NA;NA;NA"

bl %>% str_split("\\|") %>% map_chr(5) %>% str_split(";") %>% map_chr(1) %>% str_remove("k__") %>% table

getmeta <- function(x){sample_data(x) %>% meta()}
full_meta <- map(UNITES,getmeta)
map(full_meta,glimpse)


# ordinate
nmds <- full_EUK_ps %>% 
  ordinate(method = "NMDS")

ord_by_study <- plot_ordination(physeq = full_EUK_ps, ordination = nmds,color="BioProject") +
  coord_cartesian(xlim = c(-75,-25),ylim = c(-75,-20)) +
  scale_color_viridis_d(end=.8)

ord_by_feature <- plot_ordination(physeq = full_EUK_ps, ordination = nmds,color="Organism") +
  coord_cartesian(xlim = c(-75,-25),ylim = c(-75,-20))

ord_by_study / ord_by_feature


# Add "DatabaseMethod" and merge everything

full_UNITE_ps@sam_data$Database <- "UNITE"
full_EUK_ps@sam_data$Database <- "UNITE+Euk"
sample_names(full_EUK_ps) <- paste0(sample_names(full_EUK_ps),"_UNITE_Euk")
everything <- merge_phyloseq(full_UNITE_ps,full_EUK_ps)



# Ordinate
full_nmds <- everything %>% 
  # tax_glom("Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS")

full_ord_plot <- plot_ordination(everything, full_nmds, color="BioProject")
full_ord_plot
full_ord_plot +
  coord_cartesian(xlim = c(-12,0),ylim = c(-8,2)) +
  facet_wrap(~Database)
