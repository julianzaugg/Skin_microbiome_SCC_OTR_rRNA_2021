library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
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


source("code/utility.R")

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


# Load the processed metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", sep =",", header = T)
rownames(metadata.df) <- metadata.df$Index

# AND load the processed metadata that is unfiltered. This is for the qPCR results.
metadata_unfiltered.df <- read.csv("Result_tables/processed_unfiltered_metadata.csv", sep =",", header = T)
rownames(metadata_unfiltered.df) <- metadata_unfiltered.df$Index
metadata_unfiltered.df <- metadata_unfiltered.df[is.na(metadata_unfiltered.df$qPCR_exclude),]
metadata_unfiltered.df %>% group_by(Cohort, Sample_type) %>% tally()
metadata_unfiltered.df <- metadata_unfiltered.df[,c("Index","Sample_type","Swab_ID", "Cohort","Sample_type_colour","Sample_type_shape","S_aureus_qPCR","Staph_spp_qPCR", "qPCR_16S")]

# Remove negative samples (if present)
metadata.df <- metadata.df[!metadata.df$Sample_type == "negative",]
metadata_unfiltered.df <- metadata_unfiltered.df[! metadata_unfiltered.df$Sample_type == "negative",]

# List of forearm swabs
forearm_swab_ids_IC = c("522","523","564","565","678","679","740","741", "1172", "1200","1201","1322","1323",
                        "1358","1359","1492","1493")

forearm_swab_ids_IS <- c("1382","1383","1384", "1385","1470","1471","1561","1562",
                         "1599","1600","1649", "1650")


# Load abundance data
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)
genus.m <- df2m(read.csv("Result_tables/count_tables/Genus_counts.csv", header = T))

# Separate taxonomy string and create taxonomy label
genus_data.df <- separate(genus_rel.df, "taxonomy_genus", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus"), remove =F, sep = ";")
genus_data.df$taxonomy_label <- with(genus_data.df, paste0(Domain,";", Class,";",Family,";", Genus))
genus_data.df$taxonomy_label <- gsub("[a-z]__", "", genus_data.df$taxonomy_label)

# Melt and combine with metadata
genus_data.df <- melt(genus_data.df, variable.name = "Index", value.name = "Relative_abundance")
genus_data.df <- left_join(genus_data.df, metadata.df, by = "Index")

# Set levels for sample type
genus_data.df$Sample_type <- factor(genus_data.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))


# Create cohort specific datasets
immunosuppressed_genus_data.df <- subset(genus_data.df, Cohort == "organ transplant recipient")
immunocompetent_genus_data.df <- subset(genus_data.df, Cohort == "immunocompetent")
immunosuppressed_forearm_genus_data.df <- subset(genus_data.df, Cohort == "organ transplant recipient" & Swab_ID %in% forearm_swab_ids_IS)
immunocompetent_forearm_genus_data.df <- subset(genus_data.df, Cohort == "immunocompetent" & Swab_ID %in% forearm_swab_ids_IC)

# forearm_swab_ids_IS[!forearm_swab_ids_IS %in% immunosuppressed_forearm_genus_data.df$Swab_ID]
# forearm_swab_ids_IC[!forearm_swab_ids_IC %in% immunocompetent_forearm_genus_data.df$Swab_ID]

# Set factor levels for sample type for specific cohort data frames
immunosuppressed_genus_data.df$Sample_type <- factor(immunosuppressed_genus_data.df$Sample_type, levels = rev(c("NS","PDS", "AK", "SCC_PL", "SCC")))
immunocompetent_genus_data.df$Sample_type <- factor(immunocompetent_genus_data.df$Sample_type, levels = rev(c("PDS", "AK", "SCC_PL", "SCC")))
immunosuppressed_forearm_genus_data.df$Sample_type <- factor(immunosuppressed_forearm_genus_data.df$Sample_type, levels = rev(c("SCC_PL", "SCC")))
immunocompetent_forearm_genus_data.df$Sample_type <- factor(immunocompetent_forearm_genus_data.df$Sample_type, levels = rev(c("SCC_PL", "SCC")))


# Generate full genus summary for each sample type
immunosuppressed_genus_summary.df <- 
  immunosuppressed_genus_data.df %>% 
  dplyr::group_by(Sample_type) %>%
  dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
  dplyr::group_by(taxonomy_genus, taxonomy_label, Sample_type) %>%
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

immunocompetent_genus_summary.df <- 
  immunocompetent_genus_data.df %>% 
  dplyr::group_by(Sample_type) %>%
  dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
  dplyr::group_by(taxonomy_genus,taxonomy_label, Sample_type) %>%
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

immunosuppressed_forearm_genus_summary.df <- 
  immunosuppressed_forearm_genus_data.df %>% 
  dplyr::group_by(Sample_type) %>%
  dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
  dplyr::group_by(taxonomy_genus,taxonomy_label, Sample_type) %>%
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

immunocompetent_forearm_genus_summary.df <- 
  immunocompetent_forearm_genus_data.df %>% 
  dplyr::group_by(Sample_type) %>%
  dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
  dplyr::group_by(taxonomy_genus,taxonomy_label, Sample_type) %>%
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

immunosuppressed_top_genus_summary.df <-
  immunosuppressed_genus_summary.df %>% dplyr::group_by(Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)

immunocompetent_top_genus_summary.df <-
  immunocompetent_genus_summary.df %>% dplyr::group_by(Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)

immunosuppressed_forearm_top_genus_summary.df <-
  immunosuppressed_forearm_genus_summary.df %>% dplyr::group_by(Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)

immunocompetent_forearm_top_genus_summary.df <-
  immunocompetent_forearm_genus_summary.df %>% dplyr::group_by(Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)

# Create a unique list of the top taxa across all the cohorts / groups
# both_cohorts_lesions_top_genus <- unique(c(immunosuppressed_top_genus_summary.df$taxonomy_label,
#                                            immunocompetent_top_genus_summary.df$taxonomy_label,
#                                            immunosuppressed_forearm_top_genus_summary.df$taxonomy_label,
#                                            immunocompetent_forearm_top_genus_summary.df$taxonomy_label))
# table(immunosuppressed_genus_summary.df$taxonomy_label)[table(immunosuppressed_genus_summary.df$taxonomy_label) > 5]

both_cohorts_lesions_top_genus <- unique(c(immunosuppressed_top_genus_summary.df$taxonomy_genus,
                                           immunocompetent_top_genus_summary.df$taxonomy_genus,
                                           immunosuppressed_forearm_top_genus_summary.df$taxonomy_genus,
                                           immunocompetent_forearm_top_genus_summary.df$taxonomy_genus))

# ----------------------------------------
# Create palette based on unique set
# publication_palette_10 <- c("#d35238","#6ab74d","#8562cc","#c3ab41","#688bcd","#c07b44","#4bb193","#c361aa","#6d8038","#c9566e")
# publication_palette_12 <- c("#61b64e","#d64080","#b7ba3e","#6975c9","#797b34","#d09348","#55b2d4","#57aa7a","#cf5235","#ba5f62","#c67bb7","#aa55c6")

# both_cohorts_genus_palette <- setNames(publication_palette_10[1:length(both_cohorts_lesions_top_genus)], both_cohorts_lesions_top_genus)
# both_cohorts_genus_palette <- setNames(publication_palette_12[1:length(both_cohorts_lesions_top_genus)], both_cohorts_lesions_top_genus)


# Load pre-existing palette
taxa_palette.df <- read.csv("Result_tables/genus_palette.csv", header = T)

# This should be all TRUE
summary(both_cohorts_lesions_top_genus %in% taxa_palette.df$taxonomy_genus)
taxa_palette.df <- taxa_palette.df[taxa_palette.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]

# Finally, change palette taxonomy names to match labels
taxa_palette.df <-separate(taxa_palette.df, "taxonomy_genus", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus"), remove =F, sep = ";")
taxa_palette.df$taxonomy_label <- with(taxa_palette.df, paste0(Domain,";", Class,";",Family,";", Genus))
taxa_palette.df$taxonomy_label <- gsub("[a-z]__", "", taxa_palette.df$taxonomy_label)

both_cohorts_genus_palette <- with(taxa_palette.df, setNames(Colour, taxonomy_label))

# Set Other colour to grey
both_cohorts_genus_palette["Other"] <- "grey"


# ----------------------------------------

# Re-label any taxa not in the top taxa (across ALL cohorts/groups) to "Other" for each summary table
# immunosuppressed_genus_summary.df[!immunosuppressed_genus_summary.df$taxonomy_label %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"
# immunocompetent_genus_summary.df[!immunocompetent_genus_summary.df$taxonomy_label %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"
# immunosuppressed_forearm_genus_summary.df[!immunosuppressed_forearm_genus_summary.df$taxonomy_label %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"
# immunocompetent_forearm_genus_summary.df[!immunocompetent_forearm_genus_summary.df$taxonomy_label %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"

immunosuppressed_genus_summary.df[!immunosuppressed_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"
immunocompetent_genus_summary.df[!immunocompetent_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"
immunosuppressed_forearm_genus_summary.df[!immunosuppressed_forearm_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"
immunocompetent_forearm_genus_summary.df[!immunocompetent_forearm_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_label <- "Other"
immunosuppressed_genus_summary.df[!immunosuppressed_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_genus <- "Other"
immunocompetent_genus_summary.df[!immunocompetent_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_genus <- "Other"
immunosuppressed_forearm_genus_summary.df[!immunosuppressed_forearm_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_genus <- "Other"
immunocompetent_forearm_genus_summary.df[!immunocompetent_forearm_genus_summary.df$taxonomy_genus %in% both_cohorts_lesions_top_genus,]$taxonomy_genus <- "Other"

# Normalise the Mean_relative_abundance values within each sample type
sum(immunosuppressed_genus_summary.df$Mean_relative_abundance)
sum(immunosuppressed_genus_summary.df$Normalised_mean_relative_abundance)
immunosuppressed_genus_summary.df <- 
  immunosuppressed_genus_summary.df %>% 
  dplyr::group_by(Sample_type) %>% 
  dplyr::mutate(Normalised_mean_relative_abundance = Mean_relative_abundance/sum(Mean_relative_abundance)) %>% 
  as.data.frame()

immunocompetent_genus_summary.df <- 
  immunocompetent_genus_summary.df %>% 
  dplyr::group_by(Sample_type) %>% 
  dplyr::mutate(Normalised_mean_relative_abundance = Mean_relative_abundance/sum(Mean_relative_abundance)) %>% 
  as.data.frame()

immunosuppressed_forearm_genus_summary.df <-
  immunosuppressed_forearm_genus_summary.df %>%
  dplyr::group_by(Sample_type) %>%
  dplyr::mutate(Normalised_mean_relative_abundance = Mean_relative_abundance/sum(Mean_relative_abundance)) %>%
  as.data.frame()

immunocompetent_forearm_genus_summary.df <-
  immunocompetent_forearm_genus_summary.df %>%
  dplyr::group_by(Sample_type) %>%
  dplyr::mutate(Normalised_mean_relative_abundance = Mean_relative_abundance/sum(Mean_relative_abundance)) %>%
  as.data.frame()
sum(immunosuppressed_genus_summary.df$Mean_relative_abundance)
sum(immunosuppressed_genus_summary.df$Normalised_mean_relative_abundance)

immunosuppressed_genus_summary.df[immunosuppressed_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA
immunocompetent_genus_summary.df[immunocompetent_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA
immunosuppressed_forearm_genus_summary.df[immunosuppressed_forearm_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA
immunocompetent_forearm_genus_summary.df[immunocompetent_forearm_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA

# Need a single entry for the Other group, so collapse the normalised abundance values for Other by summing
immunosuppressed_genus_summary.df <-
  immunosuppressed_genus_summary.df %>% 
  group_by(Sample_type, taxonomy_genus, taxonomy_label) %>%
  # group_by(Sample_type, taxonomy_genus) %>%
  dplyr::summarise(Mean_relative_abundance = max(Mean_relative_abundance), 
                   Normalised_mean_relative_abundance = sum(Normalised_mean_relative_abundance)) %>% 
  as.data.frame()

immunocompetent_genus_summary.df <- 
  immunocompetent_genus_summary.df %>% 
  group_by(Sample_type, taxonomy_genus, taxonomy_label) %>%
  # group_by(Sample_type, taxonomy_genus) %>% 
  dplyr::summarise(Mean_relative_abundance = max(Mean_relative_abundance), 
                   Normalised_mean_relative_abundance = sum(Normalised_mean_relative_abundance)) %>% 
  as.data.frame()

immunosuppressed_forearm_genus_summary.df <-
  immunosuppressed_forearm_genus_summary.df %>%
  group_by(Sample_type, taxonomy_genus, taxonomy_label) %>%
  # group_by(Sample_type, taxonomy_genus) %>% 
  dplyr::summarise(Mean_relative_abundance = max(Mean_relative_abundance),
                   Normalised_mean_relative_abundance = sum(Normalised_mean_relative_abundance)) %>%
  as.data.frame()

immunocompetent_forearm_genus_summary.df <-
  immunocompetent_forearm_genus_summary.df %>%
  group_by(Sample_type, taxonomy_genus, taxonomy_label) %>%
  # group_by(Sample_type, taxonomy_genus) %>% 
  dplyr::summarise(Mean_relative_abundance = max(Mean_relative_abundance),
                   Normalised_mean_relative_abundance = sum(Normalised_mean_relative_abundance)) %>%
  as.data.frame()


immunosuppressed_genus_summary.df[immunosuppressed_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA
immunocompetent_genus_summary.df[immunocompetent_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA
immunosuppressed_forearm_genus_summary.df[immunosuppressed_forearm_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA
immunocompetent_forearm_genus_summary.df[immunocompetent_forearm_genus_summary.df$taxonomy_genus == "Other","Mean_relative_abundance"] <- NA


# Calculate the mean qPCR totals for each cohort + lesion type.
# These are used for the mean bars in the scatter plots
# THESE ARE ONLY QPCR VALUES FOR SAMPLES PASSING QC
Mean_qpcr_totals.df <-
  metadata.df %>% 
  group_by(Cohort, Sample_type) %>% 
  dplyr::summarise(N_samples = n_distinct(Index), #!!!
                   Mean_qPCR_16S = mean(qPCR_16S, na.rm = T)
  ) %>% 
  as.data.frame()

Mean_qpcr_totals_forearm.df <-
  metadata.df %>% 
  filter(Swab_ID %in% c(forearm_swab_ids_IS,forearm_swab_ids_IC)) %>%
  group_by(Cohort, Sample_type) %>% 
  dplyr::summarise(N_samples = n_distinct(Index), #!!!
                   Mean_qPCR_16S = mean(qPCR_16S, na.rm = T)
  ) %>% 
  as.data.frame()

immunosuppressed_genus_summary.df$Cohort <- "organ transplant recipient"
immunocompetent_genus_summary.df$Cohort <- "immunocompetent"
immunosuppressed_forearm_genus_summary.df$Cohort <- "organ transplant recipient"
immunocompetent_forearm_genus_summary.df$Cohort <- "immunocompetent"


immunosuppressed_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- NA
immunocompetent_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- NA
immunosuppressed_forearm_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- NA
immunocompetent_forearm_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- NA

immunosuppressed_genus_summary.df <- left_join(immunosuppressed_genus_summary.df, Mean_qpcr_totals.df, by = c("Cohort" = "Cohort", "Sample_type" = "Sample_type"))
immunocompetent_genus_summary.df <- left_join(immunocompetent_genus_summary.df, Mean_qpcr_totals.df, by = c("Cohort" = "Cohort", "Sample_type" = "Sample_type"))
immunosuppressed_forearm_genus_summary.df <- left_join(immunosuppressed_forearm_genus_summary.df, Mean_qpcr_totals_forearm.df, by = c("Cohort" = "Cohort", "Sample_type" = "Sample_type"))
immunocompetent_forearm_genus_summary.df <- left_join(immunocompetent_forearm_genus_summary.df, Mean_qpcr_totals_forearm.df, by = c("Cohort" = "Cohort", "Sample_type" = "Sample_type"))

# Calculate Normalised_mean_relative_abundance value proportional to qPCR totals
# We don't use this for the final plot
immunosuppressed_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- 
  with(immunosuppressed_genus_summary.df, Normalised_mean_relative_abundance * Mean_qPCR_16S)
immunocompetent_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- 
  with(immunocompetent_genus_summary.df, Normalised_mean_relative_abundance * Mean_qPCR_16S)

immunosuppressed_forearm_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- 
  with(immunosuppressed_forearm_genus_summary.df, Normalised_mean_relative_abundance * Mean_qPCR_16S)
immunocompetent_forearm_genus_summary.df$Mean_relative_abundance_qpcr_16S_proportional <- 
  with(immunocompetent_forearm_genus_summary.df, Normalised_mean_relative_abundance * Mean_qPCR_16S)


# Order taxonomy by the abundance. This is only approximate.
# The Other group should be last.
immunosuppressed_genus_summary.df <- immunosuppressed_genus_summary.df %>% group_by(Sample_type) %>% arrange(Normalised_mean_relative_abundance) %>% as.data.frame()
my_levels <- c(unique(immunosuppressed_genus_summary.df$taxonomy_label)[unique(immunosuppressed_genus_summary.df$taxonomy_label) != "Other"], "Other")
immunosuppressed_genus_summary.df$taxonomy_label <- factor(immunosuppressed_genus_summary.df$taxonomy_label, levels = my_levels)
immunosuppressed_genus_summary.df$value_label <- as.character(lapply(immunosuppressed_genus_summary.df$Normalised_mean_relative_abundance, function(x) ifelse(x >= 0.05, paste0(round(x*100), "%"), "")))

immunocompetent_genus_summary.df <- immunocompetent_genus_summary.df %>% group_by(Sample_type) %>% arrange(Normalised_mean_relative_abundance) %>% as.data.frame()
my_levels <- c(unique(immunocompetent_genus_summary.df$taxonomy_label)[unique(immunocompetent_genus_summary.df$taxonomy_label) != "Other"], "Other")
immunocompetent_genus_summary.df$taxonomy_label <- factor(immunocompetent_genus_summary.df$taxonomy_label, levels = my_levels)
immunocompetent_genus_summary.df$value_label <- as.character(lapply(immunocompetent_genus_summary.df$Normalised_mean_relative_abundance, function(x) ifelse(x >= 0.05, paste0(round(x*100), "%"), "")))

immunosuppressed_forearm_genus_summary.df <- immunosuppressed_forearm_genus_summary.df %>% group_by(Sample_type) %>% arrange(Normalised_mean_relative_abundance) %>% as.data.frame()
my_levels <- c(unique(immunosuppressed_forearm_genus_summary.df$taxonomy_label)[unique(immunosuppressed_forearm_genus_summary.df$taxonomy_label) != "Other"], "Other")
immunosuppressed_forearm_genus_summary.df$taxonomy_label <- factor(immunosuppressed_forearm_genus_summary.df$taxonomy_label, levels = my_levels)
immunosuppressed_forearm_genus_summary.df$value_label <- as.character(lapply(immunosuppressed_forearm_genus_summary.df$Normalised_mean_relative_abundance, function(x) ifelse(x >= 0.05, paste0(round(x*100), "%"), "")))

immunocompetent_forearm_genus_summary.df <- immunocompetent_forearm_genus_summary.df %>% group_by(Sample_type) %>% arrange(Normalised_mean_relative_abundance) %>% as.data.frame()
my_levels <- c(unique(immunocompetent_forearm_genus_summary.df$taxonomy_label)[unique(immunocompetent_forearm_genus_summary.df$taxonomy_label) != "Other"], "Other")
immunocompetent_forearm_genus_summary.df$taxonomy_label <- factor(immunocompetent_forearm_genus_summary.df$taxonomy_label, levels = my_levels)
immunocompetent_forearm_genus_summary.df$value_label <- as.character(lapply(immunocompetent_forearm_genus_summary.df$Normalised_mean_relative_abundance, function(x) ifelse(x >= 0.05, paste0(round(x*100), "%"), "")))


# Set levels for sample type
immunosuppressed_genus_summary.df$Sample_type <- factor(immunosuppressed_genus_summary.df$Sample_type, levels = rev(c("NS", "PDS", "AK", "SCC_PL","SCC")))
immunocompetent_genus_summary.df$Sample_type <- factor(immunocompetent_genus_summary.df$Sample_type, levels = rev(c("PDS", "AK", "SCC_PL","SCC")))
immunosuppressed_forearm_genus_summary.df$Sample_type <- factor(immunosuppressed_forearm_genus_summary.df$Sample_type, levels = rev(c("SCC_PL","SCC")))
immunocompetent_forearm_genus_summary.df$Sample_type <- factor(immunocompetent_forearm_genus_summary.df$Sample_type, levels = rev(c("SCC_PL","SCC")))

# Combine the forearm with the full dataset and set forearm status
immunosuppressed_forearm_genus_summary.df$Location <- "Forearm"
immunosuppressed_genus_summary.df$Location <- "Forearm"
immunosuppressed_genus_summary.df[immunosuppressed_genus_summary.df$Sample_type %in% c("SCC_PL", "SCC"),]$Location <- "All body sites"
immunosuppressed_genus_summary_combined_forearm.df <- rbind(immunosuppressed_genus_summary.df, immunosuppressed_forearm_genus_summary.df)

immunocompetent_forearm_genus_summary.df$Location <- "Forearm"
immunocompetent_genus_summary.df$Location <- "Forearm"
immunocompetent_genus_summary.df[immunocompetent_genus_summary.df$Sample_type %in% c("SCC_PL", "SCC"),]$Location <- "All body sites"
immunocompetent_genus_summary_combined_forearm.df <- rbind(immunocompetent_genus_summary.df, immunocompetent_forearm_genus_summary.df)

# Create sample_type_label
immunosuppressed_genus_summary_combined_forearm.df$Sample_type_label <- 
  with(immunosuppressed_genus_summary_combined_forearm.df, paste0(Sample_type, "\n", Location, "\nn = ", N_samples))

immunocompetent_genus_summary_combined_forearm.df$Sample_type_label <- 
  with(immunocompetent_genus_summary_combined_forearm.df, paste0(Sample_type, "\n", Location, "\nn = ", N_samples))

immunosuppressed_genus_summary_combined_forearm.df$Sample_type_label <- 
  factor(immunosuppressed_genus_summary_combined_forearm.df$Sample_type_label, 
         levels = rev(unique(immunosuppressed_genus_summary_combined_forearm.df[with(immunosuppressed_genus_summary_combined_forearm.df, order(Location, Sample_type)),]$Sample_type_label)))

immunocompetent_genus_summary_combined_forearm.df$Sample_type_label <- 
  factor(immunocompetent_genus_summary_combined_forearm.df$Sample_type_label, 
         levels = rev(unique(immunocompetent_genus_summary_combined_forearm.df[with(immunocompetent_genus_summary_combined_forearm.df, order(Location, Sample_type)),]$Sample_type_label)))

# -----------------------------------------------------------------------------
# Format metadata for scatter plots
metadata_unfiltered_forearm.df <- metadata_unfiltered.df %>% 
  filter(Swab_ID %in% c(forearm_swab_ids_IS,forearm_swab_ids_IC) | Sample_type %in% c("NS", "PDS","AK"))
metadata_unfiltered_all_sites.df <- metadata_unfiltered.df %>% 
  filter(Sample_type %in% c("SCC_PL", "SCC"))

metadata_unfiltered_forearm.df$Location <- "Forearm"
metadata_unfiltered_all_sites.df$Location <- "All body sites"


metadata_unfiltered_forearm.df <- 
  metadata_unfiltered_forearm.df %>% 
  dplyr::group_by(Cohort, Sample_type) %>%
  dplyr::mutate(N_samples = n_distinct(Index)) %>%
  as.data.frame()

metadata_unfiltered_all_sites.df <- 
  metadata_unfiltered_all_sites.df %>% 
  dplyr::group_by(Cohort, Sample_type) %>%
  dplyr::mutate(N_samples = n_distinct(Index)) %>%
  as.data.frame()
unique(metadata_unfiltered_forearm.df[,c("Cohort","Sample_type", "N_samples")])
unique(metadata_unfiltered_all_sites.df[,c("Cohort","Sample_type", "N_samples")])

# Create sample type labels in the metadata
metadata_unfiltered_forearm.df$Sample_type_label <- 
  with(metadata_unfiltered_forearm.df, paste0(Sample_type, "\n", Location, "\nn = ", N_samples))
metadata_unfiltered_all_sites.df$Sample_type_label <- 
  with(metadata_unfiltered_all_sites.df, paste0(Sample_type, "\n", Location, "\nn = ", N_samples))

# Factorise sample type in the metadata
metadata_unfiltered_forearm.df$Sample_type <- factor(metadata_unfiltered_forearm.df$Sample_type, 
                                                     levels = rev(c("NS", "PDS", "AK", "SCC_PL","SCC")))
metadata_unfiltered_all_sites.df$Sample_type <- factor(metadata_unfiltered_all_sites.df$Sample_type, 
                                                       levels = rev(c("SCC_PL","SCC")))
# Order the labels by location and sample type
metadata_unfiltered_forearm.df$Sample_type_label <- 
  factor(metadata_unfiltered_forearm.df$Sample_type_label, 
         levels = rev(unique(metadata_unfiltered_forearm.df[with(metadata_unfiltered_forearm.df, order(Location, Sample_type)),]$Sample_type_label)))

metadata_unfiltered_all_sites.df$Sample_type_label <- 
  factor(metadata_unfiltered_all_sites.df$Sample_type_label, 
         levels = rev(unique(metadata_unfiltered_all_sites.df[with(metadata_unfiltered_all_sites.df, order(Location, Sample_type)),]$Sample_type_label)))

# Combine metadata for both cohorts
metadata_unfiltered_combined.df <- rbind(metadata_unfiltered_forearm.df, metadata_unfiltered_all_sites.df)
metadata_unfiltered_combined.df$Sample_type <- factor(metadata_unfiltered_combined.df$Sample_type, 
                                                      levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))


write.csv(metadata_unfiltered_combined.df %>% mutate(Sample_type_label = gsub("\n", ";",Sample_type_label)),
          file = "Result_tables/qPCR_data_samples_in_16S_qPCR_plot.csv",quote = F,row.names = F)

# ------------------------------------------------------------------------------------------
# Calculate the whether qPCR values are significantly different between sample types
kw_suppressed <- kruskal.test(qPCR_16S~Sample_type_label, data = subset(metadata_unfiltered_combined.df, Cohort == "organ transplant recipient"))
kw_suppressed_forearm <- kruskal.test(qPCR_16S~Sample_type_label, data = subset(metadata_unfiltered_combined.df, Cohort == "organ transplant recipient" & Location == "Forearm"))
kw_suppressed_allbodysites <- kruskal.test(qPCR_16S~Sample_type_label, data = subset(metadata_unfiltered_combined.df, Cohort == "organ transplant recipient" & Location == "All body sites"))
kw_competent <- kruskal.test(qPCR_16S~Sample_type_label, data = subset(metadata_unfiltered_combined.df, Cohort == "immunocompetent"))
kw_competent_forearm <- kruskal.test(qPCR_16S~Sample_type_label, data = subset(metadata_unfiltered_combined.df, Cohort == "immunocompetent" & Location == "Forearm"))
kw_competent_allbodysites <- kruskal.test(qPCR_16S~Sample_type_label, data = subset(metadata_unfiltered_combined.df, Cohort == "immunocompetent" & Location == "All body sites"))

immunosuppressed_qPCR_dunn <- dunnTest(x = qPCR_16S~Sample_type_label, 
                                       data = subset(metadata_unfiltered_combined.df, Cohort == "organ transplant recipient"), 
                                       method = "bh", alpha = 0.05)$res

immunocompetent_qPCR_dunn <-  dunnTest(x = qPCR_16S~Sample_type_label, 
                                       data = subset(metadata_unfiltered_combined.df, Cohort == "immunocompetent"), 
                                       method = "bh", alpha = 0.05)$res
immunosuppressed_qPCR_dunn$Cohort <- "organ transplant recipient"
immunocompetent_qPCR_dunn$Cohort <- "immunocompetent"
dunn <- rbind(immunosuppressed_qPCR_dunn,immunocompetent_qPCR_dunn)
dunn$Comparison <- gsub("\\n", " ", dunn$Comparison)
dunn <- separate(dunn, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Cohort","Group_1","Group_2","P.unadj","P.adj")]
dunn$P_value_label <- as.character(lapply(dunn[,"P.adj"], function(x) ifelse(x <= 0.0001, "****",
                                                                             ifelse(x <= 0.001, "***", 
                                                                                    ifelse(x <= 0.01, "**", 
                                                                                           ifelse(x <= 0.05, "*", "ns"))))))
dunn <- dunn[dunn$P_value_label != "ns",]
write.csv(x = dunn, file = "Result_tables/dunn_tests/qPCR_sampletype_dunn.csv", quote = F, row.names = F)

metadata_unfiltered_combined.df %>% 
  dplyr::group_by(Cohort, Sample_type,Sample_type_label) %>%
  tally()
metadata_unfiltered_combined.df %>% 
  dplyr::group_by(Cohort, Sample_type,Sample_type_label) %>%
  tally()

mean_qpcr_values.df <- metadata_unfiltered_combined.df %>%
  dplyr::group_by(Cohort, Sample_type,Sample_type_label) %>%
  dplyr::summarise(Mean_qPCR_16S = mean(qPCR_16S,na.rm = T)) %>%
  as.data.frame()
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
#                                       PLOTTING 

# Combine processed data for both cohorts
all_combined.df <- rbind(immunosuppressed_genus_summary_combined_forearm.df, immunocompetent_genus_summary_combined_forearm.df)

# Add colours for sample types
all_combined.df <- unique(left_join(all_combined.df, metadata.df[,c("Sample_type","Sample_type_colour")], by = "Sample_type"))

temp <- all_combined.df
temp$Sample_type_label <- gsub("\\n", ";", temp$Sample_type_label)
write.csv(x = temp %>% select(-Sample_type_colour) %>% as.data.frame(), 
          file = "Result_tables/abundance_analysis/genus_abundances_combined_publication.csv", 
          quote = F,
          row.names = F)

all_combined.df %>% 
  dplyr::group_by(Cohort, Sample_type,Sample_type_label) %>%
  tally()

# Compare SCC and SCC_PL from all body sites to forearm only samples, 
# to determine change in mean relative abundances
temp <- subset(all_combined.df, Sample_type %in% c("SCC", "SCC_PL") & taxonomy_label != "Other")
temp$Cohort_sampletype_taxonomy_label <- with(temp, paste0(Cohort, "__", Sample_type, "__", taxonomy_label))
temp$Cohort_sampletype_taxonomy_genus <- with(temp, paste0(Cohort, "__", Sample_type, "__", taxonomy_genus))

differences <- list("immunocompetent" = c(),
                    "immunosuppressed" = c())
for (tax in unique(temp$Cohort_sampletype_taxonomy_genus)){
  cohort <- strsplit(tax, "__")[[1]][1]
  f_abundance <- round(subset(temp, Cohort_sampletype_taxonomy_genus == tax & Location == "Forearm")$Mean_relative_abundance*100,3)
  abs_abundance <- round(subset(temp, Cohort_sampletype_taxonomy_genus == tax & Location == "All body sites")$Mean_relative_abundance*100,3)
  difference <- f_abundance - abs_abundance
  if (cohort == "immunocompetent"){
    differences$immunocompetent <- c(differences$immunocompetent, difference)
  }else{
    differences$immunosuppressed <- c(differences$immunosuppressed, difference)
  }
  
  print(paste(tax, f_abundance, abs_abundance, f_abundance - abs_abundance))
}
mean(abs(differences$immunocompetent))
sd(abs(differences$immunocompetent))
median(abs(differences$immunocompetent))
mean(abs(differences$immunosuppressed))
sd(abs(differences$immunosuppressed))
median(abs(differences$immunosuppressed))

# ------------------------------------------------------------------------------------------
# Plot using qPCR scatter plots
# undo order-by-abundance and just use all_combined order
immunosuppressed_genus_summary_combined_forearm.df$taxonomy_label <- 
  factor(immunosuppressed_genus_summary_combined_forearm.df$taxonomy_label, levels = levels(all_combined.df$taxonomy_label))
immunocompetent_genus_summary_combined_forearm.df$taxonomy_label <- 
  factor(immunocompetent_genus_summary_combined_forearm.df$taxonomy_label, levels = levels(all_combined.df$taxonomy_label))


just_legend_plot <- 
  ggplot(all_combined.df, 
       aes(x = Sample_type_label, 
           y = Normalised_mean_relative_abundance, 
           fill = taxonomy_label)) +
  geom_bar(stat = "identity", colour = "black", lwd = .1) +
  # coord_flip() +
  scale_fill_manual(values = both_cohorts_genus_palette, name = "Domain;Class;Family;Genus", 
                    guide = guide_legend(title.position = "top",nrow= 4)) +
  xlab("Sample site") +
  ylab("Normalised mean relative abundance (%)") +
  common_theme +
  # theme(plot.title = element_text(hjust = .5, face = "bold"))
  theme(plot.title = element_text(hjust = .5, face = "bold"),
        axis.text.y = element_blank(),
        axis.title.x = element_blank())


immunosuppressed_abundance_plot <- 
  ggplot(immunosuppressed_genus_summary_combined_forearm.df, 
         aes(x = Sample_type_label, 
             y = Normalised_mean_relative_abundance*100, 
             fill = taxonomy_label)) +
  geom_bar(stat = "identity", colour = "black", lwd = .1) +
  # geom_text(aes(label = value_label), position = position_stack(vjust = 0.5), size = 2,color = "grey10") +
  scale_fill_manual(values = both_cohorts_genus_palette, guide = F) +
  scale_y_continuous(breaks = seq(0,100, by = 10), limits = c(0,101), expand = c(0, 0)) +
  # scale_y_continuous(breaks = c(seq(0,100, by = 10),1000),limits =c(0,100.01), expand = c(0, 0)) +
  # coord_cartesian(ylim = c(0,100.01)) +
  xlab("Sample type") +
  ylab("Normalised mean relative abundance (%)") +
  # ylab(expression(bold(atop("Normalised mean relative abundance (%)","",sep ="\n")))) +
  common_theme +
  theme(plot.title = element_text(hjust = .5, face = "bold"),
        axis.text.x = element_blank(),
        # axis.title.y = element_text(lineheight = 2),
        axis.title.x = element_blank())
immunosuppressed_abundance_plot

immunocompetent_abundance_plot <- 
  ggplot(immunocompetent_genus_summary_combined_forearm.df,
         aes(x = Sample_type_label, 
             y = Normalised_mean_relative_abundance*100,
             fill = taxonomy_label)) +
  geom_bar(stat = "identity", colour = "black", lwd = .1) +
  # geom_text(aes(label = value_label), position = position_stack(vjust = 0.5), size = 2,color = "grey10") +
  scale_fill_manual(values = both_cohorts_genus_palette, guide = F) +
  scale_y_continuous(breaks = seq(0,100, by = 10), limits = c(0,101), expand = c(0, 0)) +
  xlab("Sample type") +
  # ylab("Normalised mean relative abundance (%)") +
  common_theme +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"))
immunocompetent_abundance_plot

immunosuppressed_qpcr_plot <- 
  ggplot(subset(metadata_unfiltered_combined.df, Cohort == "organ transplant recipient"), 
         aes(x = Sample_type_label,
             y = qPCR_16S,
             # colour = Sample_type_colour,
             fill = Sample_type_colour,
             shape = Sample_type_shape)) +
  geom_errorbar(data = subset(mean_qpcr_values.df, Cohort == "organ transplant recipient")[,c("Sample_type_label", "Mean_qPCR_16S")],
                aes(ymax = Mean_qPCR_16S, ymin = Mean_qPCR_16S,x = Sample_type_label),inherit.aes = F,
                width = .4, lwd =.6, linetype = "solid", colour = "black") +
  geom_jitter(show.legend = F,size = 1.5,stroke = .3, alpha = 1,position = position_jitter(width = .15)) +
  scale_y_log10(limits=c(10^-1, 10^5),breaks=10^(-1:5),labels = trans_format('log10',math_format(10^.x))) +
  scale_colour_identity() + 
  scale_fill_identity() + 
  scale_shape_identity() +
  # scale_shape_manual(values = variable_shapes) +
  # labs(title = "Organ transplant recipient") +
  xlab("Sample type") +
  ylab("SSU rRNA equivalents per sampling area") +
  common_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6))
immunosuppressed_qpcr_plot

immunocompetent_qpcr_plot <- 
  ggplot(subset(metadata_unfiltered_combined.df, Cohort == "immunocompetent"), 
         aes(x = Sample_type_label,
             y = qPCR_16S,
             # colour = Sample_type_colour,
             fill = Sample_type_colour,
             shape = Sample_type_shape)) +
  geom_errorbar(data = subset(mean_qpcr_values.df, Cohort == "immunocompetent")[,c("Sample_type_label", "Mean_qPCR_16S")],
                aes(ymax = Mean_qPCR_16S, ymin = Mean_qPCR_16S,x = Sample_type_label),inherit.aes = F,
                width = .4, lwd =.6, linetype = "solid", colour = "black") +
  geom_jitter(show.legend = F,size = 1.5,stroke = .3, alpha = 1,position = position_jitter(width = .15)) +
  scale_y_log10(limits=c(10^-1, 10^5),breaks=10^(-1:5),labels = trans_format('log10',math_format(10^.x))) +
  scale_colour_identity() + 
  scale_fill_identity() + 
  scale_shape_identity() + 
  # scale_shape_manual(values = variable_shapes) +
  # labs(title = "Immunocompetent") +
  xlab("Sample type") +
  ylab("SSU rRNA equivalents per sampling area") +
  common_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6))


# Extract the legend
my_legend_taxa <- cowplot::get_legend(just_legend_plot + 
                                        theme(
                                          legend.position = "right",
                                          legend.text = element_text(size = 7),
                                          legend.title = element_text(size = 8, face="bold"),
                                          legend.justification = "center",
                                          legend.direction = "horizontal",
                                          legend.box.just = "bottom",
                                          plot.margin = unit(c(0, 0, 0, 0), "cm")
                                        )
)
# now add the title
title <- ggdraw() + 
  draw_label(
    "",
    fontface = 'bold',
    size = 8
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0)
  )



grid_plot_immunosuppressed <- plot_grid(plotlist = list(immunosuppressed_abundance_plot+
                                                          labs(title = "Organ transplant recipient") +
                                                          theme(plot.margin = unit(c(5.5,5.5,5.5,14),"pt"),
                                                                axis.text.x = element_text(size = 7),
                                                                axis.text.y = element_text(size = 7),
                                                                axis.title.y = element_text(size = 9),
                                                                plot.title = element_text(size = 12)),
                                                        immunosuppressed_qpcr_plot +
                                                          theme(plot.margin = unit(c(5.5,5.5,5.5,12),"pt"),
                                                                axis.text.x = element_text(size = 7),
                                                                axis.text.y = element_text(size = 7),
                                                                axis.title = element_text(size = 9))
),
ncol = 1, rel_widths = c(1,1),scale = c(1,1),axis= "r")
grid_plot_immunosuppressed

grid_plot_immunocompetent <- plot_grid(plotlist = list(immunocompetent_abundance_plot +
                                                         labs(title = "Immunocompetent") +
                                                         theme(plot.title = element_text(size = 12),
                                                               axis.text.x = element_text(size = 7)), 
                                                       immunocompetent_qpcr_plot +
                                                         theme(axis.text.x = element_text(size = 7),
                                                               axis.text.y = element_blank(),
                                                               axis.title.x = element_text(size = 9),
                                                               axis.title.y = element_blank())),
                                       ncol = 1, rel_widths = c(1,1))
grid_plot_immunocompetent
grid_plot <- plot_grid(grid_plot_immunosuppressed, 
                       grid_plot_immunocompetent,
                       rel_heights = c(1,1),rel_widths = c(1,.8),scale = c(1,1), ncol = 2, nrow=1)

grid_plot <- plot_grid(grid_plot, my_legend_taxa, rel_heights = c(1,0.1), ncol = 1, nrow=2)
grid_plot


# Add (a) (b) figure labels
text_par <- grid::gpar(col = "black", fontsize = 16, 
                       fontface = "bold", lineheight = 0.9, alpha = 1)
text.grob_a <- grid::textGrob("(a)", x = grid::unit(0.5, "npc"), 
                              y = grid::unit(0.5, "npc"), hjust = 0.5, vjust = 0.5, 
                              rot = 0, gp = text_par)
text.grob_b <- grid::textGrob("(b)", x = grid::unit(0.5, "npc"), 
                              y = grid::unit(0.5, "npc"), hjust = 0.5, vjust = 0.5, 
                              rot = 0, gp = text_par)

grid_plot <- grid_plot +
  annotation_custom(text.grob_a, xmin = 0.02, xmax = 0.02, ymin = .98, ymax=.98) +
  annotation_custom(text.grob_b, xmin = 0.02, xmax = 0.02, ymin = .54, ymax=.54) 

# annotation_custom(rect_grob, xmin = 0.09, xmax = 0.14, ymin = .11, ymax=.144)

ggsave(filename = "Result_figures/abundance_analysis_plots/IS_IC_sample_type_relative_abundance_and_rRNA_qPCR_scatter.pdf", 
       plot = grid_plot, width = 26, 
       height = 25, units = "cm")
ggsave(filename = "Result_figures/abundance_analysis_plots/IS_IC_sample_type_relative_abundance_and_rRNA_qPCR_scatter.svg", 
       plot = grid_plot, width = 26, height = 25, units = "cm",device = "svg")
ggsave(plot =my_legend_taxa,"Result_figures/abundance_analysis_plots/taxa.svg", 
       width = 27, units = "cm", device = "svg")

