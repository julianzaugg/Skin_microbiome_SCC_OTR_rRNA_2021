library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(lemon)

# Calculate the significance values for taxa between multiple groups
calculate_taxa_significances_multiple <- function(mydata, variable_column, value_column, taxonomy_column){
  results.df <- data.frame("Taxonomy" = character(),
                           "Variable" = character(),
                           "Group_1" = character(),
                           "Group_2" = character(),
                           "Dunn_pvalue" = character(),
                           "Dunn_padj" = character(),
                           "KrusW_pvalue" = character()
  )
  for (taxa in unique(mydata[,taxonomy_column])){ # For each taxa in the taxonomy column
    taxa_data <- subset(mydata, get(taxonomy_column) == taxa)
    n_groups = length(as.character(unique(taxa_data[,variable_column])))
    if (any(is.na(taxa_data[,variable_column]))){
      return()
    }
    if (all(taxa_data[,value_column] == 0)){
      next()
    }
    # print(taxa)
    # print(n_groups)
    if (n_groups > 2){
      kw <- kruskal.test(get(value_column)~get(variable_column), data = taxa_data)
      dunn <- dunnTest(x = get(value_column)~get(variable_column), data = taxa_data, method = "bh", alpha = 0.05)
      dunn <- separate(dunn$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
      names(dunn) <- c("Group_1","Group_2","Dunn_pvalue","Dunn_padj")
      dunn$Taxonomy <- taxa
      dunn$KrusW_pvalue <- kw$p.value
      dunn$Variable <- variable_column
      results.df <- rbind(results.df, dunn)
    }
  }
  results.df[,c("Taxonomy", "Variable", "Group_1","Group_2", "Dunn_pvalue", "Dunn_padj", "KrusW_pvalue")]
}

generate_p_labels <- function(sig_table){
  for (sig_column in c("Dunn_padj")){
    metric = strsplit(sig_column, "_")[[1]][1]
    sig_table[,paste0(metric, "_p_label")] <-
      as.character(lapply(sig_table[,sig_column], 
                          function(x) ifelse(x <= 0.001, "***", 
                                             ifelse(x <= 0.01, "**",
                                                    ifelse(x <= 0.05, "*", "ns")))))
  }
  sig_table
}
common_theme <- theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
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
  legend.key.height=unit(.4,"cm"),
  legend.text = element_text(size = 8),
  axis.text = element_text(size = 8, colour = "black"),
  axis.title = element_text(size = 10,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))



source("code/utility.R")


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Load the processed metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", sep =",", header = T)

# Filter to immunocomprised only
metadata.df <- metadata.df %>% filter(Cohort == "organ transplant recipient")

# Create LOIS group labels
metadata.df <-
  metadata.df %>%
  dplyr::group_by(Length_of_immunosuppression_group_1,Cohort, Sample_type) %>%
  dplyr::mutate(N_samples_LoIS_group_1 = n_distinct(Index),
                LoIS_group_1_label = paste0(Length_of_immunosuppression_group_1, " (n = ", n_distinct(Index), ")")) %>%
  dplyr::group_by(Length_of_immunosuppression_group_2,Cohort, Sample_type) %>%
  dplyr::mutate(N_samples_LoIS_group_2 = n_distinct(Index),
                LoIS_group_2_label = paste0(Length_of_immunosuppression_group_2, " (n = ", n_distinct(Index), ")")) %>%
  as.data.frame()

# Load abundance data
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)
genus_rel.df <- genus_rel.df[,c("taxonomy_genus", metadata.df$Index)]

# Calculate the mean relative abundances for each cohort + sample type
genus_mean_abundances.df <-
  genus_rel.df %>%
  pivot_longer(data = .,
               cols = -taxonomy_genus,
               names_to = "Index",
               values_to = "Relative_abundance") %>%
  left_join(metadata.df, by = "Index") %>% 
  pivot_longer(data = .,
               cols = c("Length_of_immunosuppression_group_1","Length_of_immunosuppression_group_2"),
               names_to = c("LOIS_group_scheme"),
               values_to = c("LOIS_group")) %>%
  dplyr::group_by(Cohort, Sample_type,LOIS_group_scheme,LOIS_group, taxonomy_genus,.drop = F) %>% 
  dplyr::summarise(Mean_relative_abundance = round(mean(Relative_abundance)*100, 5),
                   N_samples = n_distinct(Index)) %>%
  arrange(dplyr::desc(Mean_relative_abundance)) %>%
  as.data.frame()

# Get list of unique taxa
unique_top_genus.v <- genus_mean_abundances.df %>% 
  group_by(Cohort, Sample_type,LOIS_group_scheme,LOIS_group) %>% 
  dplyr::slice_head(n = 5) %>% 
  pull(taxonomy_genus) %>% 
  unique()

# Filter abundance matrix to top taxa
genus_rel.df <- genus_rel.df[genus_rel.df$taxonomy_genus %in% unique_top_genus.v,]

# We are going to look at specific taxa (commensals) later...
specific_taxa.v <- c("Staphylococcus","Paracoccus","Cutibacterium","Malassezia","Micrococcus","Pseudomonas")

# ...for now, will look at the top # taxa
# Combine abundances for top taxa with the metadata
genus_rel_processed.df <- 
  genus_rel.df %>%
  pivot_longer(data = .,
               cols = -taxonomy_genus,
               names_to = "Index",
               values_to = "Relative_abundance") %>%
  left_join(metadata.df, by = "Index") %>% 
  pivot_longer(data = .,
               cols = c("Length_of_immunosuppression_group_1","Length_of_immunosuppression_group_2"),
               names_to = c("LOIS_group_scheme"),
               values_to = c("LOIS_group")) %>%
  select(Index, Sample_type,
         Cohort, Relative_abundance,taxonomy_genus,LOIS_group_scheme,LOIS_group) %>%
  as.data.frame()

# Create genus label
genus_rel_processed.df$Genus <- gsub("g__", "", unlist(lapply(as.character(genus_rel_processed.df$taxonomy_genus), function(x) unlist(strsplit(x,";"))[6])))

# And multiply abundance
genus_rel_processed.df$Relative_abundance <- genus_rel_processed.df$Relative_abundance * 100

# Order with Staphylococcus first
genus_rel_processed.df$Genus <- factor(genus_rel_processed.df$Genus, levels = c("Staphylococcus", 
                                                                                sort(unique(genus_rel_processed.df$Genus)[unique(genus_rel_processed.df$Genus) != "Staphylococcus"])))

# Factorise variable columns
genus_rel_processed.df$Sample_type <- factor(genus_rel_processed.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL", "SCC"))
genus_rel_processed.df$LOIS_group <- factor(genus_rel_processed.df$LOIS_group,
                                            levels = c("2 to 6", "2 to 8", "7 to 15", "9 to 20", "16 and higher", "21 and higher"))

# Calculate significances for both cohorts and merge with abundance data
genus_significances.df <-
  # Abundance data for top taxa
  genus_rel_processed.df %>% 
  group_by(Cohort,Sample_type,LOIS_group_scheme, taxonomy_genus) %>%
  nest() %>%
  # ------------------------------
  # Filter out nested groups where the max abundance is 0 (otherwise error thrown)
  mutate(max_abundance = map(data, ~max(.x$Relative_abundance))) %>%
  filter(max_abundance > 0) %>%
  # ------------------------------
  # Calculate p-values within nested groups
  mutate(KrusW_pvalue = map(data, ~kruskal.test(.x$Relative_abundance~.x$LOIS_group)$p.value),
         dunn_model = map(data, ~dunnTest(.x$Relative_abundance~.x$LOIS_group,method = "bh", alpha = 0.05)$res)) %>%
  # ------------------------------
  unnest(c(KrusW_pvalue,dunn_model)) %>%
  select(-data,-Z) %>%
  # ------------------------------
  # Process "Comparison" column from dunnTest
  mutate(Comparison = gsub("[()]", "", Comparison)) %>% 
  separate(Comparison, into = c("Group_1", "Group_2"), sep = " - ") %>% #
  # ------------------------------
  select(c(taxonomy_genus, Cohort, Group_1, Group_2,KrusW_pvalue, P.unadj, P.adj)) %>%
  dplyr::rename(Dunn_pvalue = P.unadj, Dunn_padj ="P.adj") %>% 
  filter(Dunn_padj <= 0.05, KrusW_pvalue <= 0.05) %>%
  # ------------------------------
  # Now add abundance/N samples data for each group
  pivot_longer(data = .,
             cols = c("Group_1","Group_2"),
             names_to = c("Group"),
             values_to = c("LOIS_group")) %>%
  left_join(genus_mean_abundances.df, by = c("Cohort", "Sample_type","LOIS_group_scheme","LOIS_group", "taxonomy_genus")) %>%
  pivot_wider(data = .,
              names_from = c(Group),
              values_from = c(LOIS_group, N_samples, Mean_relative_abundance)
  ) %>%
  
  # Create genus label
  dplyr::mutate(Genus = gsub("g__", "", unlist(lapply(as.character(taxonomy_genus), function(x) unlist(strsplit(x,";"))[6])))) %>%
  
  # Reorder columns
  relocate(taxonomy_genus, Genus, Cohort, Sample_type,LOIS_group_scheme,
           LOIS_group_Group_1, LOIS_group_Group_2,
           N_samples_Group_1, N_samples_Group_2,
           Mean_relative_abundance_Group_1, Mean_relative_abundance_Group_2) %>%
  unique() %>% 
  as.data.frame()


# ------------------------------------------------------------------------
# Testing (alternative method)

# temp_sig <- genus_significances.df
# genus_significances.df <- data.frame()
# for (st in unique(genus_rel_processed.df$Sample_type)){
#   for (lgs in unique(genus_rel_processed.df$LOIS_group_scheme)){
#     data_subset.df <- genus_rel_processed.df %>% filter(Sample_type == st, LOIS_group_scheme == lgs) %>% as.data.frame()
#     temp <- calculate_taxa_significances_multiple(mydata = data_subset.df,
#                                                   variable_column = "LOIS_group",
#                                                   value_column = "Relative_abundance",
#                                                   taxonomy_column = "taxonomy_genus")
#     temp$Sample_type <- st
#     temp$LOIS_group_scheme <- lgs
#     genus_significances.df <- rbind(genus_significances.df, temp)
#   }
# }
# 
# genus_significances.df <- genus_significances.df %>% filter(Dunn_padj <= 0.05, KrusW_pvalue <= 0.05)
# sort(temp_sig$Dunn_padj) == sort(genus_significances.df$Dunn_padj)
# 
# genus_significances.df$Cohort <- "organ transplant recipient"
# names(genus_significances.df)[names(genus_significances.df) == "Taxonomy"] <- "taxonomy_genus"
# genus_significances.df$Variable <- NULL
# 
# genus_significances.df <-
#   genus_significances.df %>%
#   pivot_longer(data = .,
#                cols = c("Group_1","Group_2"),
#                names_to = c("Group"),
#                values_to = c("LOIS_group")) %>%
#   left_join(genus_mean_abundances.df, by = c("Cohort", "Sample_type","LOIS_group_scheme","LOIS_group", "taxonomy_genus")) %>%
#   pivot_wider(data = .,
#               names_from = c(Group),
#               values_from = c(LOIS_group, N_samples, Mean_relative_abundance)) %>%
#   as.data.frame()
# Create genus label
# genus_significances.df$Genus <- gsub("g__", "", unlist(lapply(as.character(genus_significances.df$taxonomy_genus), function(x) unlist(strsplit(x,";"))[6])))
# ------------------------------------------------------------------------

# Generate p-value labels
genus_significances.df <- generate_p_labels(genus_significances.df)

# Write results to file
write.csv(file = "Result_tables/abundance_analysis/sample_type_lois__genus_abundance_significances.csv",
          x = genus_significances.df, row.names = F,quote = F)

# Also, filter to taxa of interest and write to file
write.csv(x = subset(genus_significances.df, Genus %in% specific_taxa.v),
          file = "Result_tables/abundance_analysis/sample_type_lois__genus_specific_abundance_significances.csv", row.names = F, quote = F)


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Boxplots

sample_type_colours <- unique(metadata.df[,c("Sample_type", "Sample_type_colour")])
sample_type_colours <- setNames(as.character(sample_type_colours[,"Sample_type_colour"]), sample_type_colours[,"Sample_type"])
sample_type_shapes <- unique(metadata.df[,c("Sample_type", "Sample_type_shape")])
sample_type_shapes <- setNames(as.numeric(sample_type_shapes[,"Sample_type_shape"]), sample_type_shapes[,"Sample_type"])

sample_type_outline_colours <- unique(metadata.df[,c("Sample_type", "Sample_type_colour")])
sample_type_outline_colours <- setNames(darken(as.character(sample_type_outline_colours[,"Sample_type_colour"])), sample_type_outline_colours[,"Sample_type"])

# Now filter to the specific commensals
genus_rel_specific.df <- genus_rel_processed.df %>% filter(Genus %in% specific_taxa.v)

lois_group_1_plot <- 
  ggplot(genus_rel_specific.df %>% filter(LOIS_group_scheme == "Length_of_immunosuppression_group_1"), 
         aes(x = LOIS_group, 
             y = Relative_abundance, 
             fill = Sample_type, 
             shape = Sample_type)) +
  geom_boxplot(position = position_dodge(width =.75), 
               outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .2,
                                              dodge.width = .75)) +
  scale_shape_manual(values = sample_type_shapes,name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1,
               position = position_dodge(width = .75),show.legend = F) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 20)) + 
  ylab("Relative abundance") +
  lemon::facet_rep_wrap(Genus~Sample_type, scales = "free_x", ncol = 5,repeat.tick.labels = F) +
  common_theme +
  theme(
    # axis.text.x = element_text(face = "italic"),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0, "cm")
  )
lois_group_1_plot 

lois_group_2_plot <- 
  ggplot(genus_rel_specific.df %>% filter(LOIS_group_scheme == "Length_of_immunosuppression_group_2"), 
         aes(x = LOIS_group, 
             y = Relative_abundance, 
             fill = Sample_type, 
             shape = Sample_type)) +
  geom_boxplot(position = position_dodge(width =.75), 
               outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .2,
                                              dodge.width = .75)) +
  scale_shape_manual(values = sample_type_shapes,name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1,
               position = position_dodge(width = .75),show.legend = F) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 20)) + 
  xlab("Years of immunosuppression") +
  ylab("Relative abundance") +
  lemon::facet_rep_wrap(Genus~Sample_type, scales = "free_x", ncol = 5,repeat.tick.labels = F) +
  common_theme +
  theme(
    # axis.text.x = element_text(face = "italic"),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0, "cm")
  )
lois_group_2_plot

ggsave(plot = lois_group_1_plot, 
       filename = "Result_figures/abundance_analysis_plots/lois_scheme1__abundance_boxplot_publication.pdf",
       width = 30, height = 25, units = "cm"
)

ggsave(plot = lois_group_2_plot, 
       filename = "Result_figures/abundance_analysis_plots/lois_scheme2__abundance_boxplot_publication.pdf",
       width = 30, height = 25, units = "cm"
)



