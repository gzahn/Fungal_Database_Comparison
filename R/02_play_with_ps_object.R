#### Setup ####
library(tidyverse)
library(phyloseq)
library(readxl)
library(patchwork)
library(modelr)

plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
                       title = NULL, facet_grid = NULL) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


#### Load Data ####
# load phyloseq object
ps <- readRDS("./output/ps_object_not_cleaned.RDS")

# load metadata (just to check)
meta <- read_xlsx("Data/PRJNA664779.xlsx")


#### Explore ps object ####
sample.names(ps)
rank_names(ps)

otu_table(ps) %>% rowSums() %>% plot()
sam_data(ps)$disease_status
tax_table(ps)[,2] %>% table()


#### Remove non-fungi ####
fung <- subset_taxa(ps, Kingdom == "k__Fungi")
tax_table(fung)[,1] %>% table()

#### Remove empty ASVs ####
colSums(otu_table(fung)) %>% summary()
# no emptyies found

rowSums(otu_table(fung)) %>% summary()
# no empty samples found

data.frame(SamplingEffort = rowSums(otu_table(fung)),
           Disease_Status = fung@sam_data$disease_status) %>% 
  ggplot(aes(x=factor(Disease_Status),y=SamplingEffort)) +
  geom_boxplot() + geom_jitter()


#### remove NA values from disease_status ####
fung <- subset_samples(fung, !is.na(disease_status))




# transform to relative abundance
fung_ra <- transform_sample_counts(fung, function(x){x/sum(x)})

# bar plot of phyla
plot_bar(fung_ra,fill = "Phylum")

fung@sam_data$disease_status

# bar plot for phylum
fung %>% 
  merge_samples("disease_status") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Phylum") +
  theme_minimal() +
  labs(y="Relative abundance")

# bar plot for class
fung %>% 
  merge_samples("disease_status") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Class") +
  theme_minimal() +
  labs(y="Relative abundance") +
  scale_fill_viridis_d()

# alpha diversity
fung %>% 
  plot_richness(x="disease_status",measures = c("Observed","Shannon")) +
  geom_boxplot()


# beta diversity
DCA <- fung %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "DPCoA")

plot_ord2 <- function(i){
  p <- plot_ordination(fung %>% transform_sample_counts(function(x){x/sum(x)}),
                  DCA,
                  color=i) +
    theme_minimal()
  return(p)
}

p1 <- plot_ord2("BMI")
p2 <- plot_ord2("sex")
p3 <- plot_ord2("disease_status")

p1 + p2 + p3



# Add new column for disease (binary)
fung@sam_data$disease_status[fung@sam_data$disease_status == "PD"] <- 1
fung@sam_data$disease_status[fung@sam_data$disease_status == "Control"] <- 0
fung@sam_data$disease_status <- as.numeric(fung@sam_data$disease_status)


fung@sam_data$Shannon <- estimate_richness(fung,measures = "Shannon")$Shannon

# model disease presence as function of shannon diversity
mod <- glm(data=fung@sam_data %>% as("data.frame"),
           formula = disease_status ~ Shannon,
           family = "binomial")

add_predictions(data=fung@sam_data %>% as("data.frame"),
                model = mod,
                type = "response") %>% 
  ggplot(aes(x=Shannon,y=pred,color=sex)) +
  geom_point()
summary(mod)




ps@sam_data$Bristol