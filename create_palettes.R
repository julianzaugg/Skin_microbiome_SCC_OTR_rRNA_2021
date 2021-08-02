# Create palettes for publication figures
# Specifically, unique and consistent colours for the taxa

source("code/utility.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", sep =",", header = T)
rownames(metadata.df) <- metadata.df$Index

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

# Calculate summary stats of top taxa. All figures using colours for taxonomy are genus based,
# and depict the top N genera per various groupings. Below we identify the top genus
# analysed in companion scripts and assign unique colours to each. These genus should
# be only those depicted in the figures!

# Cohort - sample type summary
cohort_sampletype_genus_summary.df <-
  genus_data.df %>% 
  dplyr::group_by(Cohort,Sample_type) %>%
  dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
  dplyr::group_by(taxonomy_genus,Cohort, Sample_type) %>%
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

loi_sampletype_genus_summary.df <-
  genus_data.df %>% 
    filter(Cohort == "organ transplant recipient") %>% 
    melt(measure.vars = c("Length_of_immunosuppression_group_1", "Length_of_immunosuppression_group_2"), 
         variable.name = "LOI_type", value.name = "LOI_group") %>%
    dplyr::group_by(LOI_type, LOI_group, Sample_type) %>%
    dplyr::mutate(Samples_total = n_distinct(Index)) %>% # number of unique samples/index
    dplyr::group_by(taxonomy_genus,LOI_type, LOI_group, Sample_type) %>%
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

top_cohort_sampletype_genus_summary.df <-
  cohort_sampletype_genus_summary.df %>% 
  dplyr::group_by(Cohort, Sample_type) %>% 
  top_n(10, Mean_relative_abundance)

top_loi_sampletype_genus_summary.df <-
  cohort_sampletype_genus_summary.df %>% 
  dplyr::group_by(Cohort, Sample_type) %>% 
  top_n(10, Mean_relative_abundance)

top_genus <- unique(c(top_cohort_sampletype_genus_summary.df$taxonomy_genus,
                      top_loi_sampletype_genus_summary.df$taxonomy_genus))
top_genus <- sort(top_genus)


colour_palette_25_soft <- c("#d843ae","#85b533","#e27c8e","#df926b","#b88ed9","#914f9c","#507533","#48b5a3","#a9ac66","#9d486f","#ab4a4b","#dc9334","#926b2e","#db3d7f","#63a7de","#cb5326","#b5ab37","#64b871","#de3b48","#bf53d2","#7457d6","#556dad","#4ab943","#db80bc","#6277dc")

genus_palette <- setNames(colour_palette_25_soft[1:length(top_genus)], top_genus)

# Manually specify some genus (a bit trial-and-error until figures look good...)
names(genus_palette)[grepl("Micro", names(genus_palette))] 
genus_palette["d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus"] <- "#ab4a4b"
genus_palette["d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus"] <- "#e38d3d"
genus_palette["d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Micrococcales;f__Micrococcaceae;g__Micrococcus"] <- "#1fc2a7"



# Set Other colour to grey
genus_palette["Other"] <- "grey"

palette.df <-
  data.frame("taxonomy_genus" = names(genus_palette), 
             "Colour" = as.character(genus_palette))

write.csv(palette.df, file = "Result_tables/genus_palette.csv", quote = F, row.names = F)

