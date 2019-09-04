---
title: "Triclosan_Tetracycline_Analysis"
author: "Stephen Techtmann"
date: "9/4/2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Load Appropriate packages
```{r}
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(csv)
library(vegan)
library(DESeq2)
```

##Set plotting them
```{r}
theme_set(theme_bw())
Best_colors <-c("darkgray", "blue","goldenrod", "green", "magenta", "orange", 'cyan',"red", 'grey0')
```

## Set working directory whre the `.rds` files from dada2 output are
```{r}
setwd("~/Documents/Manuscripts/AZIZ_TCS/Analysis/Microbial_Response_Triclosan_tetracycline/")
```

## Import sample metadata
```{r}
sdata=as.csv("Combined_metadata.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
samdata=sample_data(sdata)
```

## Import chimera checked ASV table from dada2
```{r}
seqtab.nochim=readRDS("seqtabnochim.rds")
seqtab=otu_table(seqtab.nochim,taxa_are_rows = FALSE)
```

## Import taxa assignments for ASVs from dada2
```{r}
taxasp=readRDS("taxasp.rds")
taxa=tax_table(taxasp)
```

## Merge imported files into a phyloseq object

```{r}
TC_phylo=phyloseq(otu_table(seqtab),tax_table(taxa),sample_data(samdata))
TC_phylo
```

## Rarify ASV table

```{r}
samplesover1000_all <- subset_samples(TC_phylo, sample_sums(TC_phylo) > 1000)
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)
set.seed(81)
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))
TCS_rare=rarefy_samplesover1000_all
min(sample_sums(TCS_rare))
```

## PCoA of all samples colored by treatment and shape by location

```{r}
MD.ord <- ordinate(TCS_rare,"PCoA", "bray")
plotr1=plot_ordination(TCS_rare, MD.ord)#, type = "site", color= "SITE", shape= "OIL_TYPE")
plotr1+
  geom_point(aes(fill = location, shape= treatment_factor_2),size=4, color="black", stroke=1) + 
  scale_shape_manual(values=c(21,22,23,24,25,20))+
  #scale_shape_discrete(breaks=c("CONTROL","BAKKEN","BITUMEN")) +
  scale_fill_manual(values=Best_colors)+
  theme_set(theme_bw())+
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+
  labs(title = "PCoA All Samples")+
  guides(fill=guide_legend(override.aes=list(shape=21)))
```

## Subset phyloseq object into different regions
```{r}
GB_rare=subset_samples(TCS_rare,region=="Green_Bay")
Houghton_rare=subset_samples(TCS_rare,region=="Houghton")
HM_rare=subset_samples(TCS_rare,region=="Huron_Mountains")
```

## Plot PCoA of Green Bay
```{r}
MD.ord <- ordinate(GB_rare,"PCoA", "bray")
plotr1= plot_ordination(GB_rare, MD.ord)#, type = "site", color= "location", shape= "treatment_factor_2")
plotr1+
  geom_point(aes(fill = time, shape= treatment_factor_2),size=4, color="black", stroke=1) + 
  scale_shape_manual(values=c(21,22,23,24,25,20))+
  #scale_shape_discrete(breaks=c("CONTROL","BAKKEN","BITUMEN")) +
  scale_fill_manual(values=Best_colors)+
  theme_set(theme_bw())+
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+
  labs(title = "PCoA All Samples")+
  guides(fill=guide_legend(override.aes=list(shape=21)))
```


## Plot PCoA of Houghton
```{r}
MD.ord <- ordinate(Houghton_rare,"PCoA", "bray")
plotr1= plot_ordination(Houghton_rare, MD.ord)#, type = "site", color= "location", shape= "treatment_factor_2")
plotr1+
  geom_point(aes(fill = time, shape= treatment_factor_2),size=4, color="black", stroke=1) + 
  scale_shape_manual(values=c(21,22,23,24,25,20))+
  #scale_shape_discrete(breaks=c("CONTROL","BAKKEN","BITUMEN")) +
  scale_fill_manual(values=Best_colors)+
  theme_set(theme_bw())+
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+
  labs(title = "PCoA All Samples")+
  guides(fill=guide_legend(override.aes=list(shape=21)))
```


## Plot PCoA of Huron Mountains
```{r}
MD.ord <- ordinate(HM_rare,"PCoA", "bray")
plotr1= plot_ordination(HM_rare, MD.ord)#, type = "site", color= "location", shape= "treatment_factor_2")
plotr1+
  geom_point(aes(fill = time, shape= treatment_factor_2),size=4, color="black", stroke=1) + 
  scale_shape_manual(values=c(21,22,23,24,25,20))+
  #scale_shape_discrete(breaks=c("CONTROL","BAKKEN","BITUMEN")) +
  scale_fill_manual(values=Best_colors)+
  theme_set(theme_bw())+
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+
  labs(title = "PCoA All Samples")+
  guides(fill=guide_legend(override.aes=list(shape=21)))
```

## PERMANOVA Analysis for each location

```{r}
##Stats for all regions
all_brayr <- phyloseq::distance(TCS_rare, method = "bray")
sampledf <- data.frame(sample_data(TCS_rare))
adonis(all_brayr~treatment_factor_2, data=sampledf ,permutations = 999, method = "bray")
adonis(all_brayr~treatment_factor_2*time, data=sampledf ,permutations = 999, method = "bray")
adonis(all_brayr~treatment_factor_2*time*location, data=sampledf ,permutations = 999, method = "bray")
```

## Stats for Green_Bay
```{r}
GB_brayr <- phyloseq::distance(GB_rare, method = "bray")
sampledf_GB <- data.frame(sample_data(GB_rare))
adonis(GB_brayr~treatment_factor_2, data=sampledf_GB ,permutations = 999, method = "bray")
adonis(GB_brayr~treatment_factor_2*time, data=sampledf_GB ,permutations = 999, method = "bray")
```

## Stats for Houghton
```{r}
Houghton_brayr <- phyloseq::distance(Houghton_rare, method = "bray")
sampledf_Houghton <- data.frame(sample_data(Houghton_rare))
adonis(Houghton_brayr~treatment_factor_2, data=sampledf_Houghton ,permutations = 999, method = "bray")
adonis(Houghton_brayr~treatment_factor_2*time, data=sampledf_Houghton ,permutations = 999, method = "bray")
```

## Stats for Huron_Mountains
```{r}
HM_brayr <- phyloseq::distance(HM_rare, method = "bray")
sampledf_HM <- data.frame(sample_data(HM_rare))
adonis(HM_brayr~treatment_factor_2, data=sampledf_HM ,permutations = 999, method = "bray")
adonis(HM_brayr~treatment_factor_2*time, data=sampledf_HM ,permutations = 999, method = "bray")
```

## Pairwise PERMANOVAs

```{r}
calc_pairwise_permanovas(Houghton_brayr,sampledf_Houghton,"location")
```
## Determine Alpha diversity from rarifying table 100 times and calculating mean of shannon, inverse simpsons, and richness(observed ASVs)

```{r}
min_lib <- min(sample_sums(prune_samplesover1000_all))

nsamp = nsamples(prune_samplesover1000_all)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(prune_samplesover1000_all)

shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(shannon) <- sample_names(prune_samplesover1000_all)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(prune_samplesover1000_all)


# It is always important to set a seed when you subsample so your result is replicable 
set.seed(81)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(prune_samplesover1000_all, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
  
  # Calculate shannon
  shan <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
  shannon[ ,i] <- shan
}


#---------------mean and std -----------------#
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(shannon)
mean <- apply(shannon, 1, mean)
sd <- apply(shannon, 1, sd)
measure <- rep("Shannon", nsamp)
Shannon_stats <- data.frame(SampleID, mean, sd, measure)

##estimate richness###
alpha <- rbind(rich_stats, even_stats,Shannon_stats)

s <- data.frame(sample_data(prune_samplesover1000_all))
names(s)[names(s) == "X.SampleID"] <- "SampleID"
alphadiv <- cbind(alpha,s)

rarefy_alpha_rich=subset(alphadiv,measure=="Richness")
rarefy_alpha_even=subset(alphadiv,measure=="Inverse Simpson")
rarefy_alpha_shannon=subset(alphadiv,measure=="Shannon")
```

## Plot Richness Diversity

```{r}
ggplot(rarefy_alpha_rich, aes(x = time, y = mean)) +
  geom_boxplot()+
  geom_point(aes(fill = treatment_factor_2, shape = treatment_factor_2),color="black", size = 5) + 
  geom_line(size = 0.8) +
  facet_grid(treatment_factor_2~location)+
  #facet_wrap(~measure, ncol = 1, scales = "free") +
  scale_shape_manual(values=c(21,22,23,24,25,20))+
  scale_fill_manual(values = Best_colors) +
  scale_x_discrete(
    breaks = c("0H", "24H", "7D", "14D","21D","28D", "29D","35D","36D", "42D","43D","49D","50D", "56D"),
    #labels = c("Jul", "Aug", "Sep", "Oct"), 
    drop = FALSE
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  guides(fill = guide_legend(override.aes=list(shape=21)))
```

## Statistics for alpha diversity
```{r}
anova_shannon_treatment_all <- aov(rarefy_alpha_shannon$mean ~ treatment_factor_2, rarefy_alpha_shannon)
summary(anova_shannon_treatment_all)
```
## Subset samples by region
```{r}
shannon_GB=subset(rarefy_alpha_shannon,region=="Green_Bay")
shannon_Houghton=subset(rarefy_alpha_shannon,region=="Houghton")
shannon_HM=subset(rarefy_alpha_shannon,region=="Huron_Mountains")
```

## ANOVA for GB

```{r}
anova_shannon_GB <- aov(shannon_GB$mean ~ treatment_factor_2, shannon_GB)
summary(anova_shannon_GB)
```
```{r}
anova_shannon_GB_treat_time <- aov(shannon_GB$mean ~ treatment_factor_2*time, shannon_GB)
summary(anova_shannon_GB_treat_time)
```
## ANOVA for Houghton
```{r}
anova_shannon_Houghton <- aov(shannon_Houghton$mean ~ treatment_factor_2, shannon_Houghton)
summary(anova_shannon_Houghton)
```
```{r}
anova_shannon_Houghton_treat_time <- aov(shannon_Houghton$mean ~ treatment_factor_2*time, shannon_Houghton)
summary(anova_shannon_Houghton_treat_time)
```


## ANOVA for Huron Mountains
```{r}
anova_shannon_HM <- aov(shannon_HM$mean ~ treatment_factor_2, shannon_HM)
summary(anova_shannon_HM)
```
```{r}
anova_shannon_HM_treat_time <- aov(shannon_HM$mean ~ treatment_factor_2*time, shannon_HM)
summary(anova_shannon_HM_treat_time)
```

## Taxa Plots

```{r}
# melt a ps object with absolute counts (we don't want log10 here)
TCS_rare_melt <- TCS_rare %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# taxa plotting - define plotting parameters
number.taxa<-length(tax_table(TCS_rare))

microbiomeseq_cols <- function(){
  colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00",
               "#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000",
               "#FFFF00",grey.colors(1000));
  return(colours)
}

colours <- microbiomeseq_cols()


theme_set(theme_bw())
pal = "Set3"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

Tims_Colors <- c("darkgray", "blue", "aquamarine", "deepskyblue", "chocolate", "goldenrod","green","hotpink3","magenta","orange","olivedrab3","navy","seagreen","royalblue","red","purple","peachpuff2","tan1","orangered")
```

```{r}
# taxaplot
p.pre<-ggplot(data=TCS_rare_melt,aes(x=TCS_rare_melt$time,fill=TCS_rare_melt$Phylum, y=TCS_rare_melt$Abundance))+
  geom_bar(stat="identity", position="fill")+#,color="black") +
  facet_grid(TCS_rare_melt$treatment_factor_2~location, drop=TRUE,scale="free",space="free_x")

p.pre<- p.pre+ guides(fill=guide_legend(ncol=1))+
  #scale_color_manual(values="black")+
  scale_fill_manual(values=colours[1:(number.taxa+1)])+
  theme_bw()+
  xlab("")+
  ylab("Relative Abundance (% all reads)") +
  ggtitle("Bacterial Community Comparison (Sample Type vs Cruise)")
```

```{r}
p.pre + 
  labs(fill = "Bacterial Class")+
  theme(axis.title = 
          element_text(size = rel(2)) 
  )+
  theme(axis.text.x = 
          element_text(size = rel(1.5)) 
  )+
  theme(axis.text.y = 
          element_text(size = rel(1.5)) 
  )+ 
  theme(legend.text=element_text(size=20))
```

## DESeq


### Subset into two factor tables.
```{r}
Control_TCS2ppm=subset_samples(prune_samplesover1000_all,treatment_factor_2=="Control"|treatment_factor_2=="Triclosan_2ppm")
Control_TCS6ppm=subset_samples(prune_samplesover1000_all,treatment_factor_2=="Control"|treatment_factor_2=="Triclosan_6ppm")
Control_Tet=subset_samples(prune_samplesover1000_all,treatment_factor_2=="Control"|treatment_factor_2=="Tetracycline")
```
### Convert phyloseq object to DESeq object non-rarified data
```{r}
sample_sums(Control_TCS2ppm)
All_treat_DESeq <- subset_samples(Control_TCS2ppm, treatment_factor_2 ="Control")

diagddsALL = phyloseq_to_deseq2(All_treat_DESeq, ~ treatment_factor_2)
diagddsALL$treatment_factor_2 <- relevel( diagddsALL$treatment_factor_2, "Control" )
```

### Performe DESeq analysis
```{r}
as.data.frame(colData(diagddsALL))
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagddsALL), 1, gm_mean)
diagddsALL = estimateSizeFactors(diagddsALL, geoMeans = geoMeans)
diagddsALL = DESeq(diagddsALL, fitType="local")
head(diagddsALL)

resultsALL = results(diagddsALL, cooksCutoff = FALSE)
alpha = 0.05
sigtabdsALL = resultsALL[which(resultsALL$padj < alpha), ]
sigtabdsALL = cbind(as(sigtabdsALL, "data.frame"), as(tax_table(ALLDESeq)[rownames(sigtabdsALL), ], "matrix"))

head(sigtabdsALL)
class(sigtabdsALL)
dim(sigtabdsALL)

normalized_countsALL <- counts(diagddsALL, normalized = TRUE)
dim(normalized_countsALL)
normalized_countsALL <- t(normalized_countsALL)
dim(normalized_countsALL)

write.csv(sigtabdsALL,"Location_sigtabdsALL.csv")
sigtabdsALL


########Volcano plot#########
## my data is called GB_res.  This is not hte sig table, it is the output from DESeq.
## You should have a command similar to this in your DESeq code
GB_res = cbind(as(resultsALL, "data.frame"), as(tax_table(prune_samplesover1000_all)[rownames(resultsALL), ], "matrix"))

##this command will add a column with labels for the OTUs that are significantly enriched.
GB_res$threshold = as.factor(abs(GB_res$log2FoldChange) > 2 & GB_res$padj < 0.05)
#This command will rename those columns to indicate which OTUs were enriched or not enriched
GB_res$threshold=as.numeric(GB_res$threshold)
GB_res$threshold <- as.character(GB_res$threshold)
GB_res$threshold[GB_res$threshold == "1"] <- "Not Enriched"
GB_res$threshold[GB_res$threshold == "2"] <- "Enriched"

#This is a color palette that I like
cbbPalette <- c("orange","blue","darkgray", "black","#009E73","#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#This makes a volcano plot
g = ggplot(data=GB_res, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(colour=threshold, fill=threshold), size=3) +
  #opts(legend.position = "none") +
  #xlim(c(-10, 10)) + ylim(c(0, 15)) +
  scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  #scale_colour_gradient2(mid="light grey",high="dark blue")+
  #scale_fill_gradient2(mid="light grey",high="dark blue")+
  scale_shape_manual(values=c(21,22,23,24,25,20))+
  labs(title = "Oil vs Control",
       x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed") 
g
```
