library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(purrr)
library(scales)


common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white", size = 1),
  strip.text = element_text(size = 6),
  legend.key=element_blank(),
  legend.direction="vertical",
  legend.background = element_rect(colour ="white", size = .3),
  legend.text.align = 0,
  legend.title = element_text(size=10, face="bold"),
  legend.title.align = 0.5,
  legend.margin = margin(c(2,2,2,2)),
  legend.key.height= unit(.3,"cm"),
  legend.key.width = unit(.3,"cm"),
  legend.text = element_text(size = 8),
  axis.line = element_line(colour = "black", size = 0.5),
  axis.text = element_text(size = 6, colour = "black"),
  axis.title = element_text(size = 7,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))


# Load utility functions
source("code/utility.R")
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# Load the metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", sep =",", header = T)
rownames(metadata.df) <- metadata.df$Index

# AND load the processed metadata that is unfiltered. This is for the qPCR results.
metadata_unfiltered.df <- read.csv("Result_tables/processed_unfiltered_metadata.csv", sep =",", header = T)
metadata_unfiltered.df <- metadata_unfiltered.df[metadata_unfiltered.df$Cohort == "organ transplant recipient",]
rownames(metadata_unfiltered.df) <- metadata_unfiltered.df$Index
metadata_unfiltered.df <- metadata_unfiltered.df[is.na(metadata_unfiltered.df$qPCR_exclude),]
metadata_unfiltered.df %>% group_by(Cohort, Sample_type) %>% tally()
metadata_unfiltered.df <- metadata_unfiltered.df[,c("Index","Sample_type","Swab_ID", "Cohort",
                                                    "Length_of_immunosuppression_group_1","Length_of_immunosuppression_group_2",
                                                    "Sample_type_colour","Sample_type_shape","S_aureus_qPCR","Staph_spp_qPCR", "qPCR_16S")]

# Remove negative samples (if present)
metadata.df <- metadata.df[!metadata.df$Sample_type == "negative",]
metadata_unfiltered.df <- metadata_unfiltered.df[! metadata_unfiltered.df$Sample_type == "negative",]

# Load abundance data
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)
genus.m <- df2m(read.csv("Result_tables/count_tables/Genus_counts.csv", header = T))


# Separate taxonomy string and create taxonomy label
genus_data.df <- separate(genus_rel.df, "taxonomy_genus", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus"), remove =F, sep = ";")
genus_data.df$taxonomy_label <- with(genus_data.df, paste0(Domain,";", Class,";",Family, ";", Genus))
genus_data.df$taxonomy_label <- gsub("[a-z]__", "", genus_data.df$taxonomy_label)

# Melt and combine with metadata
genus_data.df <- melt(genus_data.df, variable.name = "Index", value.name = "Relative_abundance")
genus_data.df <- left_join(genus_data.df, metadata.df, by = "Index")

# Set levels for sample type and suppression groups (probably need to do again later)
genus_data.df$Sample_type <- factor(genus_data.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))
genus_data.df$Length_of_immunosuppression_group_1 <- 
  with(genus_data.df, factor(Length_of_immunosuppression_group_1, levels = rev(c("2 to 6", "7 to 15", "16 and higher"))))
genus_data.df$Length_of_immunosuppression_group_2 <- 
  with(genus_data.df, factor(Length_of_immunosuppression_group_2, levels = rev(c("2 to 8", "9 to 20", "21 and higher"))))

# Filter abundance data to IS cohort
genus_data.df <- genus_data.df[genus_data.df$Cohort == "organ transplant recipient",]

genus_summary_length_IS_group_1.df <- 
  genus_data.df %>% 
  dplyr::group_by(Length_of_immunosuppression_group_1, Sample_type) %>%
  dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
  dplyr::group_by(taxonomy_genus,taxonomy_label,Length_of_immunosuppression_group_1, Sample_type) %>%
  # dplyr::group_by(taxonomy_label, Sample_type) %>%
  dplyr::mutate(Samples_present = sum(Relative_abundance > 0, na.rm = TRUE)) %>%
  dplyr::summarise(
    Samples_present = max(Samples_present),
    Samples_total = max(Samples_total),
    # Mean_relative_abundance2 = round(mean(sum(Relative_abundance)/max(Samples_total)), 5),
    Min_relative_abundance = round(min(Relative_abundance),5),
    Max_relative_abundance = round(max(Relative_abundance),5),
    Median_relative_abundance = round(median(Relative_abundance), 5),
    Mean_relative_abundance = round(mean(Relative_abundance), 5),
    
  ) %>%
  arrange(dplyr::desc(Mean_relative_abundance))

genus_summary_length_IS_group_2.df <- 
  genus_data.df %>% 
  dplyr::group_by(Length_of_immunosuppression_group_2, Sample_type) %>%
  dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
  dplyr::group_by(taxonomy_genus,taxonomy_label,Length_of_immunosuppression_group_2, Sample_type) %>%
  # dplyr::group_by(taxonomy_label, Sample_type) %>%
  dplyr::mutate(Samples_present = sum(Relative_abundance > 0, na.rm = TRUE)) %>%
  dplyr::summarise(
    Samples_present = max(Samples_present),
    Samples_total = max(Samples_total),
    # Mean_relative_abundance2 = round(mean(sum(Relative_abundance)/max(Samples_total)), 5),
    Min_relative_abundance = round(min(Relative_abundance),5),
    Max_relative_abundance = round(max(Relative_abundance),5),
    Median_relative_abundance = round(median(Relative_abundance), 5),
    Mean_relative_abundance = round(mean(Relative_abundance), 5),
    
  ) %>%
  arrange(dplyr::desc(Mean_relative_abundance))

# Identify the top genus for each lesion type
number_of_top_taxa <- 5

top_genus_summary_length_IS_group_1.df <-
  genus_summary_length_IS_group_1.df %>% dplyr::group_by(Length_of_immunosuppression_group_1,Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)
top_genus_summary_length_IS_group_2.df <-
  genus_summary_length_IS_group_2.df %>% dplyr::group_by(Length_of_immunosuppression_group_2,Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)

# Create a unique list of the top taxa
# top_genus_length_of_IS_group_1.v <- unique(c(top_genus_summary_length_IS_group_1.df$taxonomy_label))
# top_genus_length_of_IS_group_2.v <- unique(c(top_genus_summary_length_IS_group_2.df$taxonomy_label))
top_genus_length_of_IS_group_1.v <- unique(c(top_genus_summary_length_IS_group_1.df$taxonomy_genus))
top_genus_length_of_IS_group_2.v <- unique(c(top_genus_summary_length_IS_group_2.df$taxonomy_genus))

# sort(top_genus_length_of_IS_group_2.v) == sort(top_genus_length_of_IS_group_1.v)
top_genus <- unique(c(top_genus_length_of_IS_group_1.v, top_genus_length_of_IS_group_2.v))

# ------
# Create palette based on unique set
# publication_palette_10 <- c("#d35238","#6ab74d","#8562cc","#c3ab41","#688bcd","#c07b44","#4bb193","#c361aa","#6d8038","#c9566e")
# 
# # top_genus[!top_genus %in% both_cohorts_lesions_top_genus]
# genus_palette <- setNames(publication_palette_10[1:length(top_genus)], top_genus)
# 
# # Set Other colour to grey
# genus_palette["Other"] <- "grey"
# 

# Load pre-existing palette
taxa_palette.df <- read.csv("Result_tables/genus_palette.csv", header = T)

# This should be all TRUE
summary(top_genus %in% taxa_palette.df$taxonomy_genus)
taxa_palette.df <- taxa_palette.df[taxa_palette.df$taxonomy_genus %in% top_genus,]

# Finally, change palette taxonomy names to match labels
taxa_palette.df <-separate(taxa_palette.df, "taxonomy_genus", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus"), remove =F, sep = ";")
taxa_palette.df$taxonomy_label <- with(taxa_palette.df, paste0(Domain,";", Class,";",Family,";", Genus))
taxa_palette.df$taxonomy_label <- gsub("[a-z]__", "", taxa_palette.df$taxonomy_label)

genus_palette <- with(taxa_palette.df, setNames(Colour, taxonomy_label))

# Set Other colour to grey
genus_palette["Other"] <- "grey"


# ----------------------------------------

# -----

# Re-label any taxa not in the top taxa (across ALL cohorts/groups) to "Other" for each summary table
# genus_summary_length_IS_group_1.df[!genus_summary_length_IS_group_1.df$taxonomy_label %in% top_genus,]$taxonomy_label <- "Other"
# genus_summary_length_IS_group_2.df[!genus_summary_length_IS_group_2.df$taxonomy_label %in% top_genus,]$taxonomy_label <- "Other"
genus_summary_length_IS_group_1.df[!genus_summary_length_IS_group_1.df$taxonomy_genus %in% top_genus,]$taxonomy_label <- "Other"
genus_summary_length_IS_group_2.df[!genus_summary_length_IS_group_2.df$taxonomy_genus %in% top_genus,]$taxonomy_label <- "Other"
genus_summary_length_IS_group_1.df[!genus_summary_length_IS_group_1.df$taxonomy_genus %in% top_genus,]$taxonomy_genus <- "Other"
genus_summary_length_IS_group_2.df[!genus_summary_length_IS_group_2.df$taxonomy_genus %in% top_genus,]$taxonomy_genus <- "Other"

# Normalise the Mean_relative_abundance values within each sample type
sum(genus_summary_length_IS_group_1.df$Mean_relative_abundance)
sum(genus_summary_length_IS_group_1.df$Normalised_mean_relative_abundance)
genus_summary_length_IS_group_1.df <- 
  genus_summary_length_IS_group_1.df %>% 
  dplyr::group_by(Length_of_immunosuppression_group_1, Sample_type) %>% 
  dplyr::mutate(Normalised_mean_relative_abundance = Mean_relative_abundance/sum(Mean_relative_abundance)) %>% 
  as.data.frame()
sum(genus_summary_length_IS_group_1.df$Mean_relative_abundance)
sum(genus_summary_length_IS_group_1.df$Normalised_mean_relative_abundance)

sum(genus_summary_length_IS_group_2.df$Mean_relative_abundance)
sum(genus_summary_length_IS_group_2.df$Normalised_mean_relative_abundance)
genus_summary_length_IS_group_2.df <- 
  genus_summary_length_IS_group_2.df %>% 
  dplyr::group_by(Length_of_immunosuppression_group_2, Sample_type) %>% 
  dplyr::mutate(Normalised_mean_relative_abundance = Mean_relative_abundance/sum(Mean_relative_abundance)) %>% 
  as.data.frame()
sum(genus_summary_length_IS_group_2.df$Mean_relative_abundance)
sum(genus_summary_length_IS_group_2.df$Normalised_mean_relative_abundance)


genus_summary_length_IS_group_1.df[genus_summary_length_IS_group_1.df$taxonomy_label == "Other","Mean_relative_abundance"] <- NA
genus_summary_length_IS_group_2.df[genus_summary_length_IS_group_2.df$taxonomy_label == "Other","Mean_relative_abundance"] <- NA

# Need a single entry for the Other group, so collapse the normalised abundance values for Other by summing
genus_summary_length_IS_group_1.df <-
  genus_summary_length_IS_group_1.df %>% 
  # group_by(Length_of_immunosuppression_group_1, Sample_type, taxonomy_label) %>% 
  group_by(Length_of_immunosuppression_group_1, Sample_type, taxonomy_genus, taxonomy_label) %>% 
  dplyr::summarise(Mean_relative_abundance = max(Mean_relative_abundance), 
                   Normalised_mean_relative_abundance = sum(Normalised_mean_relative_abundance)) %>% 
  as.data.frame()

genus_summary_length_IS_group_2.df <-
  genus_summary_length_IS_group_2.df %>% 
  # group_by(Length_of_immunosuppression_group_2, Sample_type, taxonomy_label) %>%
  group_by(Length_of_immunosuppression_group_2, Sample_type, taxonomy_genus, taxonomy_label) %>%
  dplyr::summarise(Mean_relative_abundance = max(Mean_relative_abundance), 
                   Normalised_mean_relative_abundance = sum(Normalised_mean_relative_abundance)) %>% 
  as.data.frame()
# sum(genus_summary_length_IS_group_2.df$Mean_relative_abundance)
# sum(genus_summary_length_IS_group_2.df$Normalised_mean_relative_abundance)

genus_summary_length_IS_group_1.df[genus_summary_length_IS_group_1.df$taxonomy_label == "Other","Mean_relative_abundance"] <- NA
genus_summary_length_IS_group_2.df[genus_summary_length_IS_group_2.df$taxonomy_label == "Other","Mean_relative_abundance"] <- NA

Mean_qpcr_totals_length_IS_group_1.df <-
  metadata.df %>% 
  filter(Cohort == "organ transplant recipient") %>%
  group_by(Cohort, Length_of_immunosuppression_group_1, Sample_type) %>% 
  dplyr::summarise(N_samples = n_distinct(Index), #!!!
                   Mean_qPCR_16S = mean(qPCR_16S, na.rm = T)
  ) %>% 
  as.data.frame()

Mean_qpcr_totals_length_IS_group_2.df <-
  metadata.df %>% 
  filter(Cohort == "organ transplant recipient") %>%
  group_by(Cohort, Length_of_immunosuppression_group_2, Sample_type) %>% 
  dplyr::summarise(N_samples = n_distinct(Index), #!!!
                   Mean_qPCR_16S = mean(qPCR_16S, na.rm = T)
  ) %>% 
  as.data.frame()

genus_summary_length_IS_group_1.df$Cohort <- "organ transplant recipient"
genus_summary_length_IS_group_1.df$Mean_relative_abundance_qPCR_16S_proportional <- NA
genus_summary_length_IS_group_1.df <- left_join(genus_summary_length_IS_group_1.df, Mean_qpcr_totals_length_IS_group_1.df, by = c("Cohort" = "Cohort",
                                                                                                                               "Length_of_immunosuppression_group_1" = "Length_of_immunosuppression_group_1",
                                                                                                                               "Sample_type" = "Sample_type"))

genus_summary_length_IS_group_2.df$Cohort <- "organ transplant recipient"
genus_summary_length_IS_group_2.df$Mean_relative_abundance_qPCR_16S_proportional <- NA
genus_summary_length_IS_group_2.df <- left_join(genus_summary_length_IS_group_2.df, Mean_qpcr_totals_length_IS_group_2.df, by = c("Cohort" = "Cohort",
                                                                                                                               "Length_of_immunosuppression_group_2" = "Length_of_immunosuppression_group_2",
                                                                                                                               "Sample_type" = "Sample_type"))

# Calculate Normalised_mean_relative_abundance value proportional to qPCR totals
genus_summary_length_IS_group_1.df$Mean_relative_abundance_qPCR_16S_proportional <- 
  with(genus_summary_length_IS_group_1.df, Normalised_mean_relative_abundance * Mean_qPCR_16S)

genus_summary_length_IS_group_2.df$Mean_relative_abundance_qPCR_16S_proportional <- 
  with(genus_summary_length_IS_group_2.df, Normalised_mean_relative_abundance * Mean_qPCR_16S)

# Order taxonomy by the abundance. This is only approximate.
# The Other group should be last.
genus_summary_length_IS_group_1.df <- genus_summary_length_IS_group_1.df %>% group_by(Sample_type) %>% arrange(Normalised_mean_relative_abundance) %>% as.data.frame()
my_levels <- c(unique(genus_summary_length_IS_group_1.df$taxonomy_label)[unique(genus_summary_length_IS_group_1.df$taxonomy_label) != "Other"], "Other")
genus_summary_length_IS_group_1.df$taxonomy_label <- factor(genus_summary_length_IS_group_1.df$taxonomy_label, levels = my_levels)
genus_summary_length_IS_group_1.df$value_label <- as.character(lapply(genus_summary_length_IS_group_1.df$Normalised_mean_relative_abundance, function(x) ifelse(x >= 0.05, paste0(round(x*100), "%"), "")))

genus_summary_length_IS_group_2.df <- genus_summary_length_IS_group_2.df %>% group_by(Sample_type) %>% arrange(Normalised_mean_relative_abundance) %>% as.data.frame()
my_levels <- c(unique(genus_summary_length_IS_group_2.df$taxonomy_label)[unique(genus_summary_length_IS_group_2.df$taxonomy_label) != "Other"], "Other")
genus_summary_length_IS_group_2.df$taxonomy_label <- factor(genus_summary_length_IS_group_2.df$taxonomy_label, levels = my_levels)
genus_summary_length_IS_group_2.df$value_label <- as.character(lapply(genus_summary_length_IS_group_2.df$Normalised_mean_relative_abundance, function(x) ifelse(x >= 0.05, paste0(round(x*100), "%"), "")))

# Factorise taxonomy_label
# genus_summary_length_IS_group_1.df$taxonomy_label <-
#   with(genus_summary_length_IS_group_1.df, factor(taxonomy_label, levels = unique(taxonomy_label)))
# 
# genus_summary_length_IS_group_2.df$taxonomy_label <-
#   with(genus_summary_length_IS_group_2.df, factor(taxonomy_label, levels = unique(taxonomy_label)))


# Set levels for sample type and suppression groups (probably need to do again later)
genus_summary_length_IS_group_1.df$Sample_type <- factor(genus_summary_length_IS_group_1.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))
genus_summary_length_IS_group_1.df$Length_of_immunosuppression_group_1 <- 
  with(genus_summary_length_IS_group_1.df, factor(Length_of_immunosuppression_group_1, levels = rev(c("2 to 6", "7 to 15", "16 and higher"))))

genus_summary_length_IS_group_2.df$Sample_type <- factor(genus_summary_length_IS_group_2.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))
genus_summary_length_IS_group_2.df$Length_of_immunosuppression_group_2 <- 
  with(genus_summary_length_IS_group_2.df, factor(Length_of_immunosuppression_group_2,levels = rev(c("2 to 8", "9 to 20", "21 and higher"))))


# Create grouping label
genus_summary_length_IS_group_1.df$Group_label <- 
  with(genus_summary_length_IS_group_1.df, paste0(Length_of_immunosuppression_group_1, " (n = ", N_samples, ")"))
genus_summary_length_IS_group_2.df$Group_label <- 
  with(genus_summary_length_IS_group_2.df, paste0(Length_of_immunosuppression_group_2, " (n = ", N_samples, ")"))

# Factorise grouping label to be correct order
genus_summary_length_IS_group_1.df$Group_label <- factor(genus_summary_length_IS_group_1.df$Group_label, as.character(unique(genus_summary_length_IS_group_1.df$Group_label[order(genus_summary_length_IS_group_1.df$Length_of_immunosuppression_group_1)])))
genus_summary_length_IS_group_2.df$Group_label <- factor(genus_summary_length_IS_group_2.df$Group_label, as.character(unique(genus_summary_length_IS_group_2.df$Group_label[order(genus_summary_length_IS_group_2.df$Length_of_immunosuppression_group_2)])))

# -----------------------------------------------------------------------------
# Format metadata for 16S qpcr panels

# metadata_unfiltered.df
metadata_unfiltered.df <- 
  metadata_unfiltered.df %>% 
  dplyr::group_by(Length_of_immunosuppression_group_1,Cohort, Sample_type) %>%
  dplyr::mutate(N_samples_LoIS_group_1 = n_distinct(Index),
                LoIS_group_1_label = paste0(Length_of_immunosuppression_group_1, " (n = ", n_distinct(Index), ")")) %>%
  dplyr::group_by(Length_of_immunosuppression_group_2,Cohort, Sample_type) %>%
  dplyr::mutate(N_samples_LoIS_group_2 = n_distinct(Index),
                LoIS_group_2_label = paste0(Length_of_immunosuppression_group_2, " (n = ", n_distinct(Index), ")")) %>%
  as.data.frame()


metadata_unfiltered.df$Sample_type <- factor(metadata_unfiltered.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))
metadata_unfiltered.df$Length_of_immunosuppression_group_1 <- 
  with(metadata_unfiltered.df, factor(Length_of_immunosuppression_group_1, levels = rev(c("2 to 6", "7 to 15", "16 and higher"))))
metadata_unfiltered.df$Length_of_immunosuppression_group_2 <- 
  with(metadata_unfiltered.df, factor(Length_of_immunosuppression_group_2,levels = rev(c("2 to 8", "9 to 20", "21 and higher"))))
metadata_unfiltered.df$LoIS_group_1_label <-
  with(metadata_unfiltered.df, 
       factor(LoIS_group_1_label, as.character(unique(LoIS_group_1_label[order(Length_of_immunosuppression_group_1)]))))
metadata_unfiltered.df$LoIS_group_2_label <-
  with(metadata_unfiltered.df, 
       factor(LoIS_group_2_label, as.character(unique(LoIS_group_2_label[order(Length_of_immunosuppression_group_2)]))))

LoIS_group_1_16S_mean.df <- 
  metadata_unfiltered.df %>%
  dplyr::group_by(Length_of_immunosuppression_group_1,Sample_type,Sample_type_colour,Sample_type_shape, LoIS_group_1_label) %>%
  dplyr::summarise(mean_qPCR_16S = mean(qPCR_16S),
                   sd_qPCR_16S = sd(qPCR_16S),
                   ) %>%
  as.data.frame()
LoIS_group_2_16S_mean.df <- metadata_unfiltered.df %>%
  dplyr::group_by(Length_of_immunosuppression_group_2,Sample_type,Sample_type_colour,Sample_type_shape, LoIS_group_2_label) %>%
  dplyr::summarise(mean_qPCR_16S = mean(qPCR_16S),
                   sd_qPCR_16S = sd(qPCR_16S),
  ) %>%
  as.data.frame()


LoIS_group_1_16S_mean.df$LoIS_group_1_label <- factor(LoIS_group_1_16S_mean.df$LoIS_group_1_label, 
                                                   as.character(unique(LoIS_group_1_16S_mean.df$LoIS_group_1_label[order(LoIS_group_1_16S_mean.df$Length_of_immunosuppression_group_1)])))

LoIS_group_2_16S_mean.df$LoIS_group_2_label <- factor(LoIS_group_2_16S_mean.df$LoIS_group_2_label, 
                                                   as.character(unique(LoIS_group_2_16S_mean.df$LoIS_group_2_label[order(LoIS_group_2_16S_mean.df$Length_of_immunosuppression_group_2)])))

# ------------------------------------------------------------------------------------------
# Calculate the whether qPCR values are significantly different between sample types


LoIS_group_1_kruskal.df <- metadata_unfiltered.df %>% 
  group_by(Sample_type) %>%
  nest() %>%
  mutate(p.value = map(data, ~kruskal.test(.x$qPCR_16S~.x$LoIS_group_1_label)$p.value)) %>%
  unnest(c(p.value,data) ) %>%
  # select(c(Sample_type, Length_of_immunosuppression_group_1, LoIS_group_1_label, p.value)) %>%
  select(c(Sample_type,p.value)) %>%
  mutate(Length_of_immunosuppression_scheme = "Length_of_immunosuppression_group_1") %>%
  unique() %>%
  as.data.frame()

LoIS_group_2_kruskal.df <- metadata_unfiltered.df %>% 
  group_by(Sample_type) %>%
  nest() %>%
  mutate(p.value = map(data, ~kruskal.test(.x$qPCR_16S~.x$LoIS_group_2_label)$p.value)) %>%
  unnest(c(p.value,data) ) %>%
  # select(c(Sample_type, Length_of_immunosuppression_group_2, LoIS_group_2_label, p.value)) %>%
  select(c(Sample_type,p.value)) %>%
  mutate(Length_of_immunosuppression_scheme = "Length_of_immunosuppression_group_2") %>%
  unique() %>%
  as.data.frame()

LoIS_kruskal.df <- rbind(LoIS_group_1_kruskal.df, LoIS_group_2_kruskal.df)
LoIS_kruskal.df$P_value_label <- as.character(lapply(LoIS_kruskal.df[,"p.value"], function(x) ifelse(x <= 0.0001, "****",
                                                                                             ifelse(x <= 0.001, "***", 
                                                                                                    ifelse(x <= 0.01, "**", 
                                                                                                           ifelse(x <= 0.05, "*", "ns"))))))
# LoIS_kruskal.df <- LoIS_kruskal.df[LoIS_kruskal.df$P_value_label != "ns",]
write.csv(x = LoIS_kruskal.df, file = "Result_tables/dunn_tests/qPCR_length_of_immunosuppression_sampletype_kruskal.csv", quote = F, row.names = F)


# temp <- metadata_unfiltered.df %>% 
#   group_by(Sample_type) %>%
#   nest() %>%
#   mutate(model = map(data, ~kruskal.test(.x$qPCR_16S~.x$LoIS_group_1_label)$p.value))%>%
#   unnest(model) %>%
#   select(-data) %>%
#   as.data.frame()

# metadata_unfiltered.df %>%
#   split(., list(.$Sample_type)) %>%
#   map(~dunnTest(.x$qPCR_16S~.x$LoIS_group_1_label, method = "bh",alpha = 0.05)$res) %>%
#   map_dfr(~ as.data.frame(.), .id = "Sample_type") %>%
#   select(-Z) %>%
#   mutate(Comparison = gsub("[()]", "", Comparison)) %>%
#   separate(Comparison, into = c("Group_1", "Group_2"), sep = " - ") %>%
#   as.data.frame()

LoIS_group_1_dunn.df <- 
  metadata_unfiltered.df %>% 
  group_by(Sample_type) %>%
  nest() %>%
  mutate(model = map(data, ~dunnTest(.x$qPCR_16S~.x$LoIS_group_1_label,method = "bh", alpha = 0.05)$res)) %>%
  unnest(model) %>%
  select(-data,-Z) %>%
  mutate(Comparison = gsub("[()]", "", Comparison)) %>%
  separate(Comparison, into = c("Group_1", "Group_2"), sep = " - ") %>%
  as.data.frame()

LoIS_group_2_dunn.df <- 
  metadata_unfiltered.df %>% 
  group_by(Sample_type) %>%
  nest() %>%
  mutate(model = map(data, ~dunnTest(.x$qPCR_16S~.x$LoIS_group_2_label,method = "bh", alpha = 0.05)$res)) %>%
  unnest(model) %>%
  select(-data,-Z) %>%
  mutate(Comparison = gsub("[()]", "", Comparison)) %>%
  separate(Comparison, into = c("Group_1", "Group_2"), sep = " - ") %>%
  as.data.frame()

LoIS_group_1_dunn.df$Length_of_immunosuppression_scheme <- "Length_of_immunosuppression_group_1"
LoIS_group_2_dunn.df$Length_of_immunosuppression_scheme <- "Length_of_immunosuppression_group_2"

LoIS_dunn.df <- rbind(LoIS_group_1_dunn.df,LoIS_group_2_dunn.df)
LoIS_dunn.df <- LoIS_dunn.df[,c("Length_of_immunosuppression_scheme", grep("Length_of_immunosuppression_scheme", names(LoIS_dunn.df), value =T, invert = T))]

LoIS_dunn.df$P_value_label <- as.character(lapply(LoIS_dunn.df[,"P.adj"], function(x) ifelse(x <= 0.0001, "****",
                                                                             ifelse(x <= 0.001, "***", 
                                                                                    ifelse(x <= 0.01, "**", 
                                                                                           ifelse(x <= 0.05, "*", "ns"))))))
# LoIS_dunn.df <- LoIS_dunn.df[LoIS_dunn.df$P_value_label != "ns",]
write.csv(x = LoIS_dunn.df, file = "Result_tables/dunn_tests/qPCR_length_of_immunosuppression_sampletype_dunn.csv", quote = F, row.names = F)

# ----------------------------
# Plot abundance bar graphs

LOI_group_1_just_legend_plot <- 
  ggplot(genus_summary_length_IS_group_1.df, aes(x = Length_of_immunosuppression_group_1, y = Normalised_mean_relative_abundance, fill = taxonomy_label)) +
  geom_bar(stat = "identity", colour = "black", lwd = .2) +
  coord_flip() +
  scale_fill_manual(values = genus_palette, name = "Domain;Class;Family;Genus", guide = guide_legend(title.position = "top",nrow= 4)) +
  common_theme

LOI_group_2_just_legend_plot <- 
  ggplot(genus_summary_length_IS_group_2.df, aes(x = Length_of_immunosuppression_group_2, y = Normalised_mean_relative_abundance, fill = taxonomy_label)) +
  geom_bar(stat = "identity", colour = "black", lwd = .2) +
  coord_flip() +
  scale_fill_manual(values = genus_palette, name = "Domain;Class;Family;Genus", guide = guide_legend(title.position = "top",nrow= 4)) +
  common_theme

LOI_group_1_abundance_plot <- 
  ggplot(genus_summary_length_IS_group_1.df, aes(Group_label, Normalised_mean_relative_abundance*100, fill =taxonomy_label)) +
  geom_bar(stat = "identity", colour = "black", lwd = .2) +
  # geom_text(aes(label = value_label), position = position_stack(vjust = 0.5), size = 1.5,color = "grey10") +
  coord_flip() +
  scale_fill_manual(values = genus_palette, guide = F) +
  scale_y_continuous(breaks = seq(0,100, by = 10), limits = c(0,101),expand = c(0, 0)) +
  xlab("Years of immunosuppression (# samples)") +
  ylab("Normalised mean relative abundance (%)") +
  common_theme + 
  facet_wrap(~Sample_type, ncol = 1, scales = "free_y")
LOI_group_1_abundance_plot

LOI_group_2_abundance_plot <- 
  ggplot(genus_summary_length_IS_group_2.df, aes(Group_label, Normalised_mean_relative_abundance*100, fill =taxonomy_label)) +
  geom_bar(stat = "identity", colour = "black", lwd = .2) +
  # geom_text(aes(label = value_label), position = position_stack(vjust = 0.5), size = 1.5,color = "grey10") +
  coord_flip() +
  scale_fill_manual(values = genus_palette, guide = F) +
  scale_y_continuous(breaks = seq(0,100, by = 10), limits = c(0,101),expand = c(0, 0)) +
  xlab("Years of immunosuppression (# samples)") +
  ylab("Normalised mean relative abundance (%)") +
  common_theme + 
  facet_wrap(~Sample_type, ncol = 1, scales = "free_y")

# Qpcr plots
# LOI_group_1_bar_plot <-
#   ggplot(LoIS_group_1_16S_mean.df,
#          aes(x = LoIS_group_1_label, 
#              y = mean_qPCR_16S, 
#              fill = Sample_type_colour,
#              shape = Sample_type_shape
#          )
#   ) +
#   geom_bar(stat = "identity",colour = "black") +
#   geom_errorbar(aes(ymin=mean_qPCR_16S-sd_qPCR_16S, ymax=mean_qPCR_16S+sd_qPCR_16S), position = "dodge")+
#   scale_fill_identity() + 
#   scale_shape_identity() + 
#   facet_wrap(~Sample_type, ncol = 1, scales = "free_y")+
#   coord_flip() +
#   scale_y_continuous(expand = c(0,0)) +
#   common_theme 

# LoIS_group_1_16S_mean.df %>% filter(Sample_type == "SCC")
LOI_group_1_qpcr_plot <- 
  ggplot(metadata_unfiltered.df, 
         aes(x = LoIS_group_1_label,
             y = qPCR_16S,
             # colour = Sample_type_colour,
             fill = Sample_type_colour,
             shape = Sample_type_shape)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(show.legend = F,size = 1,stroke = .3, alpha = 1,position = position_jitter(width = .15)) +
  # geom_errorbar(data = LoIS_group_1_16S_mean.df[,c("LoIS_group_1_label","Sample_type", "mean_qPCR_16S")],
  #               aes(ymax = mean_qPCR_16S, ymin = mean_qPCR_16S,x = LoIS_group_1_label),inherit.aes = F,
  #               width = .4, lwd =.6, linetype = "solid", colour = "black") +
  # stat_summary(data = LoIS_group_1_16S_mean.df[,c("LoIS_group_1_label","Sample_type", "Sample_type_colour","mean_qPCR_16S")],
  #              aes(y = mean_qPCR_16S), colour = "grey2", geom = "point",
  #              shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  stat_summary(fun = "mean", fill = "black", geom = "point",
               shape = 21,size = 1, position = position_dodge(width = .75),show.legend = F) +
  scale_y_log10(limits=c(10^-1, 10^5),breaks=10^(-1:5),labels = trans_format('log10',math_format(10^.x))) +
  # scale_y_continuous(limits = c(0,5010), breaks = seq(0,5000,500))+
  scale_colour_identity() + 
  scale_fill_identity() + 
  scale_shape_identity() +
  # scale_shape_manual(values = variable_shapes) +
  # labs(title = "Organ transplant recipient") +
  xlab("Length of immunosuppression (years)") +
  ylab("SSU rRNA equivalents per sampling area") +
  facet_wrap(~Sample_type, ncol = 1, scales =  "free_y") +
  coord_flip() +
  common_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6))
LOI_group_1_qpcr_plot


LOI_group_2_qpcr_plot <- 
  ggplot(metadata_unfiltered.df, 
         aes(x = LoIS_group_2_label,
             y = qPCR_16S,
             # colour = Sample_type_colour,
             fill = Sample_type_colour,
             shape = Sample_type_shape)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(show.legend = F,size = 1,stroke = .3, alpha = 1,position = position_jitter(width = .15)) +
  # geom_errorbar(data = LoIS_group_2_16S_mean.df[,c("LoIS_group_2_label","Sample_type", "mean_qPCR_16S")],
  #               aes(ymax = mean_qPCR_16S, ymin = mean_qPCR_16S,x = LoIS_group_2_label),inherit.aes = F,
  #               width = .4, lwd =.6, linetype = "solid", colour = "black") +
  # stat_summary(data = LoIS_group_2_16S_mean.df[,c("LoIS_group_2_label","Sample_type", "Sample_type_colour","mean_qPCR_16S")],
  #              aes(y = mean_qPCR_16S), colour = "grey2", geom = "point",
  #              shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  stat_summary(fun = "mean", fill = "black", geom = "point",
               shape = 21,size = 1, position = position_dodge(width = .75),show.legend = F) +
  scale_y_log10(limits=c(10^-1, 10^5),breaks=10^(-1:5),labels = trans_format('log10',math_format(10^.x))) +
  # scale_y_continuous(limits = c(0,5010), breaks = seq(0,5000,500))+
  scale_colour_identity() + 
  scale_fill_identity() + 
  scale_shape_identity() +
  # scale_shape_manual(values = variable_shapes) +
  # labs(title = "Organ transplant recipient") +
  xlab("Length of immunosuppression (years)") +
  ylab("SSU rRNA equivalents per sampling area") +
  facet_wrap(~Sample_type, ncol = 1, scales =  "free_y") +
  coord_flip() +
  common_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6))
LOI_group_2_qpcr_plot

# Extract the legend
# Will be the same for both length of suppression groupings
my_legend_taxa <- cowplot::get_legend(LOI_group_1_just_legend_plot + 
                                        theme(
                                          legend.position = "right",
                                          legend.text = element_text(size = 6),
                                          legend.title = element_text(size=7, face="bold"),
                                          legend.justification = "center",
                                          legend.direction = "horizontal",
                                          legend.box.just = "bottom",
                                          plot.margin = unit(c(0, 0, 0, 0), "cm")
                                        )
)

# Grouping 1
LOI_group_1_grid_plot <- plot_grid(plotlist = list(LOI_group_1_abundance_plot+
                                         theme(plot.margin = unit(c(0,25,0,0),"pt"),
                                               axis.text.x = element_text(size = 7),
                                               axis.text.y = element_text(size = 7),
                                               axis.title = element_text(size = 9),
                                               plot.title = element_text(size = 12),
                                               strip.text = element_text(size =7, face = "bold")
                                               ),
                                       NULL,
                                       LOI_group_1_qpcr_plot +
                                         xlab("") +
                                         theme(plot.margin = unit(c(0,0,0,0),"pt"),
                                               axis.text.x = element_text(size = 7),
                                               axis.text.y = element_text(size = 7),
                                               axis.title = element_text(size = 9),
                                               strip.text = element_text(size =7, face = "bold")
                                               )
                                       ),
                                   ncol = 3, nrow=1,  rel_widths = c(1,-0.05,1),align = "hv")

LOI_group_1_grid_plot <- plot_grid(LOI_group_1_grid_plot, my_legend_taxa, rel_heights = c(1,0.4), ncol = 1, nrow=2)

# Grouping 2
LOI_group_2_grid_plot <- plot_grid(plotlist = list(LOI_group_2_abundance_plot+
                                                     theme(plot.margin = unit(c(0,25,0,0),"pt"),
                                                           axis.text.x = element_text(size = 7),
                                                           axis.text.y = element_text(size = 7),
                                                           axis.title = element_text(size = 9),
                                                           plot.title = element_text(size = 12),
                                                           strip.text = element_text(size =7, face = "bold")
                                                           ),
                                                   NULL,
                                                   LOI_group_2_qpcr_plot +
                                                     xlab("") +
                                                     theme(plot.margin = unit(c(0,0,0,0),"pt"),
                                                           axis.text.x = element_text(size = 7),
                                                           axis.text.y = element_text(size = 7),
                                                           axis.title = element_text(size = 9),
                                                           strip.text = element_text(size =7, face = "bold")
                                                           )
                                                   ),
                                   ncol = 3, nrow=1,  rel_widths = c(1,-0.05,1),align = "hv")

LOI_group_2_grid_plot <- plot_grid(LOI_group_2_grid_plot, my_legend_taxa, rel_heights = c(1,0.4), ncol = 1, nrow=2)


# Add (a) (b) figure labels
text_par <- grid::gpar(col = "black", fontsize = 13, 
                       fontface = "bold", lineheight = 0.9, alpha = 1)
text.grob_a <- grid::textGrob("(a)", x = grid::unit(0.5, "npc"), 
                              y = grid::unit(0.5, "npc"), hjust = 0.5, vjust = 0.5, 
                              rot = 0, gp = text_par)
text.grob_b <- grid::textGrob("(b)", x = grid::unit(0.5, "npc"), 
                              y = grid::unit(0.5, "npc"), hjust = 0.5, vjust = 0.5, 
                              rot = 0, gp = text_par)

LOI_group_1_grid_plot <- LOI_group_1_grid_plot +
  annotation_custom(text.grob_a, xmin = 0.02, xmax = 0.02, ymin = .98, ymax=.98) +
  annotation_custom(text.grob_b, xmin = 0.52, xmax = 0.52, ymin = .98, ymax=.98)

LOI_group_2_grid_plot <- LOI_group_2_grid_plot +
  annotation_custom(text.grob_a, xmin = 0.02, xmax = 0.02, ymin = .98, ymax=.98) +
  annotation_custom(text.grob_b, xmin = 0.52, xmax = 0.52, ymin = .98, ymax=.98)


ggsave(filename = "Result_figures/abundance_analysis_plots/length_of_immunosuppression_grouping_1.pdf", 
       plot = LOI_group_1_grid_plot, width = 25, height = 15, units = "cm")
ggsave(filename = "Result_figures/abundance_analysis_plots/length_of_immunosuppression_grouping_2.pdf", 
       plot = LOI_group_2_grid_plot, width = 25, height = 15, units = "cm")

ggsave(filename = "Result_figures/abundance_analysis_plots/length_of_immunosuppression_grouping_1.svg", 
       plot = LOI_group_1_grid_plot, width = 25, height = 15, units = "cm",device = "svg")
ggsave(filename = "Result_figures/abundance_analysis_plots/length_of_immunosuppression_grouping_2.svg", 
       plot = LOI_group_2_grid_plot, width = 25, height = 15, units = "cm",device = "svg")


