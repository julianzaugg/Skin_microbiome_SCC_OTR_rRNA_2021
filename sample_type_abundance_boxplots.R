library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)

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

# List of forearm swabs
forearm_swab_ids_IC = c("522","523","564","565","678","679","740","741", "1172", "1200","1201","1322","1323",
                        "1358","1359","1492","1493")

forearm_swab_ids_IS <- c("1382","1383","1384", "1385","1470","1471","1561","1562",
                         "1599","1600","1649", "1650")

metadata.df$Location <- "All body sites"
metadata_forearm.df <- metadata.df[metadata.df$Swab_ID %in% c(forearm_swab_ids_IC,forearm_swab_ids_IS),]
metadata_forearm.df$Sample_type <- paste0(metadata_forearm.df$Sample_type, "_forearm")
# metadata.df <- rbind(metadata.df,metadata_forearm.df)
metadata.df[metadata.df$Swab_ID %in% c(forearm_swab_ids_IC,forearm_swab_ids_IS),]$Location <- "Forearm"
 
# Make label for sample type
metadata.df <-
  metadata.df %>%
  dplyr::group_by(Cohort, Sample_type) %>%
  dplyr::mutate(N_samples_Sample_type = n_distinct(Index),
                Sample_type_label = paste0(Sample_type, " (n = ", n_distinct(Index), ")")) %>%
  as.data.frame()


# Load abundance data
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)

# Calculate the mean relative abundances for each cohort + sample type
genus_mean_abundances.df <-
  melt(genus_rel.df, variable.name = "Index", value.name = "Relative_abundance") %>%
  left_join(metadata.df, by =  "Index") %>%
  dplyr::group_by(Cohort, Sample_type,taxonomy_genus,.drop = F) %>% 
  dplyr::summarise(Mean_relative_abundance = round(mean(Relative_abundance)*100, 5),
                   N_samples = n_distinct(Index)) %>%
  arrange(dplyr::desc(Mean_relative_abundance)) %>%
  dplyr::group_by(Cohort, Sample_type) %>% 
  # dplyr::slice_head(n = 30) %>%
  as.data.frame()

# ...and again specifically for SCC/SCC_PL samples, accounting for forearm specific samples
genus_mean_abundances_location.df <-
  melt(genus_rel.df, variable.name = "Index", value.name = "Relative_abundance") %>%
  left_join(metadata.df, by =  "Index") %>%
  filter(Sample_type %in% c("SCC", "SCC_PL")) %>%
  dplyr::group_by(Cohort, Sample_type,taxonomy_genus,.drop = F) %>% 
  dplyr::summarise(Mean_relative_abundance = round(mean(Relative_abundance)*100, 5),
                   Mean_relative_abundance_forearm = round(mean(Relative_abundance[Location == "Forearm"])*100, 5),
                   N_samples = n_distinct(Index),
                   N_samples_forearm = n_distinct(Index[Location == "Forearm"])
                   ) %>%
  arrange(dplyr::desc(Mean_relative_abundance)) %>%
  dplyr::group_by(Cohort, Sample_type) %>% 
  as.data.frame()


# Get list of unique taxa
unique_top_genus.v <- genus_mean_abundances.df %>% group_by(Cohort, Sample_type) %>% dplyr::slice_head(n = 30) %>% pull(taxonomy_genus) %>% unique()

# What is the mean absolute difference in mean abundances for SCC/SCC_PL from all body sites vs only forearm 
# This will be different to the values reported in the abundance_qpcr_stacked_scatter script (different set of taxa)
genus_mean_abundances_location.df %>% filter(taxonomy_genus %in% unique_top_genus.v) %>%
  summarise(Mean_relative_abundance_difference = Mean_relative_abundance - Mean_relative_abundance_forearm) %>%
  pull(Mean_relative_abundance_difference) %>%
  mean()


# Filter abundance matrix to top taxa from both cohorts
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
  # melt(variable.name = "Index", value.name = "Relative_abundance") %>%
  left_join(metadata.df, by = "Index") %>% 
  select(Index, Sample_type, Sample_type_colour,
         Cohort,Location, Relative_abundance,taxonomy_genus) 

# Create genus label
genus_rel_processed.df$Genus <- gsub("g__", "", unlist(lapply(as.character(genus_rel_processed.df$taxonomy_genus), function(x) unlist(strsplit(x,";"))[6])))

# And multiply abundance
genus_rel_processed.df$Relative_abundance <- genus_rel_processed.df$Relative_abundance * 100

# Order with Staphylococcus first
genus_rel_processed.df$Genus <- factor(genus_rel_processed.df$Genus, levels = c("Staphylococcus", 
                                                                                sort(unique(genus_rel_processed.df$Genus)[unique(genus_rel_processed.df$Genus) != "Staphylococcus"])))

# Factorise variable columns
genus_rel_processed.df$Sample_type <- factor(genus_rel_processed.df$Sample_type,
                                             levels = c("NS", "PDS", "AK", "SCC_PL", "SCC"))
# genus_rel_processed.df$Sample_type <- factor(genus_rel_processed.df$Sample_type,
                                             # levels = c("NS", "PDS", "AK","SCC_PL_forearm", "SCC_forearm", "SCC_PL", "SCC"))

# Calculate significances for both cohorts and merge with abundance data
genus_significances.df <-
  # Abundance data for top taxa
  genus_rel_processed.df %>% 
  group_by(Cohort,taxonomy_genus) %>%
  nest() %>%
  # ------------------------------
  # Filter out nested groups where the max abundance is 0 (otherwise error thrown)
  mutate(max_abundance = map(data, ~max(.x$Relative_abundance))) %>%
  filter(max_abundance > 0) %>%
  # ------------------------------
  # Calculate p-values within nested groups
  mutate(KrusW_pvalue = map(data, ~kruskal.test(.x$Relative_abundance~.x$Sample_type)$p.value),
         dunn_model = map(data, ~dunnTest(.x$Relative_abundance~.x$Sample_type,method = "bh", alpha = 0.05)$res)) %>%
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
  unique() %>% 
  # ------------------------------
  # Now add abundance/N samples data for each group
  pivot_longer(data = .,
               cols = c("Group_1","Group_2"),
               names_to = c("Group"),
               values_to = c("Sample_type")) %>%
  left_join(genus_mean_abundances.df, by = c("Cohort", "Sample_type", "taxonomy_genus")) %>%
  pivot_wider(data = .,
              names_from = c(Group),
              values_from = c(Sample_type, N_samples, Mean_relative_abundance)
              ) %>%
  
  # Create genus label
  dplyr::mutate(Genus = gsub("g__", "", unlist(lapply(as.character(taxonomy_genus), function(x) unlist(strsplit(x,";"))[6])))) %>%
  
  # Reorder columns
  relocate(taxonomy_genus, Genus,Cohort, 
           Sample_type_Group_1, Sample_type_Group_2,
           N_samples_Group_1, N_samples_Group_2,
           Mean_relative_abundance_Group_1, Mean_relative_abundance_Group_2) %>%

  as.data.frame()

# ------------------------------------------------------------------------
# Testing (alternative method)
# IS_genus_significances.df <-
#   calculate_taxa_significances_multiple(mydata = genus_rel_processed.df %>% filter(Cohort == "organ transplant recipient"),
#                                       variable_column = "Sample_type",
#                                       value_column = "Relative_abundance",
#                                       taxonomy_column = "taxonomy_genus") %>% 
#   filter(Dunn_padj <= 0.05, KrusW_pvalue <= 0.05)
# 
# IC_genus_significances.df <-
#   calculate_taxa_significances_multiple(mydata = genus_rel_processed.df %>% filter(Cohort == "immunocompetent"),
#                                         variable_column = "Sample_type",
#                                         value_column = "Relative_abundance",
#                                         taxonomy_column = "taxonomy_genus") %>% 
#   filter(Dunn_padj <= 0.05, KrusW_pvalue <= 0.05)

# IS_genus_significances.df$Cohort <- "organ transplant recipient"
# IC_genus_significances.df$Cohort <- "immunocompetent"

# genus_significances.df <- rbind(IS_genus_significances.df,IC_genus_significances.df)
# genus_significances.df <-
#   genus_significances.df %>%
#   pivot_longer(data = .,
#                cols = c("Group_1","Group_2"),
#                names_to = c("Group"),
#                values_to = c("Sample_type")) %>%
#   left_join(genus_mean_abundances.df, by = c("Cohort", "Sample_type", "Taxonomy" = "taxonomy_genus")) %>%
#   pivot_wider(data = .,
#               names_from = c(Group),
#               values_from = c(Sample_type, N_samples, Mean_relative_abundance)
#   ) %>%
#   as.data.frame() 

# Create genus label
# genus_significances.df$Genus <- gsub("g__", "", unlist(lapply(as.character(genus_significances.df$taxonomy_genus), function(x) unlist(strsplit(x,";"))[6])))
# ------------------------------------------------------------------------

# Generate p-value labels
genus_significances.df <- generate_p_labels(genus_significances.df)

# Write results to file
write.csv(file = "Result_tables/abundance_analysis/sample_type__genus_abundance_significances.csv",
          x = genus_significances.df, row.names = F,quote = F)

# Also, filter to taxa of interest and write to file
write.csv(x = subset(genus_significances.df, Genus %in% specific_taxa.v),
          file = "Result_tables/abundance_analysis/sample_type__genus_specific_abundance_significances.csv", row.names = F, quote = F)


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
IS_genus_rel_specific.df <- genus_rel_processed.df %>% filter(Genus %in% specific_taxa.v, Cohort == "organ transplant recipient")
IC_genus_rel_specific.df <- genus_rel_processed.df %>% filter(Genus %in% specific_taxa.v, Cohort == "immunocompetent")


IS_genus_plot <- ggplot(IS_genus_rel_specific.df, aes(x = Genus, 
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
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 10)) + 
  ylab("Relative abundance") +
  common_theme +
  theme(axis.text.x = element_text(face = "italic"))
IS_genus_plot


ggsave(plot = IS_genus_plot, 
       filename = "Result_figures/abundance_analysis_plots/IS_abundance_boxplot_publication.pdf",
       width = 18, height = 12, units = "cm"
)



IC_genus_plot <- ggplot(IC_genus_rel_specific.df, aes(x = Genus, 
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
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, by = 10)) + 
  ylab("Relative abundance") +
  common_theme +
  theme(axis.text.x = element_text(face = "italic"))
IC_genus_plot

ggsave(plot = IC_genus_plot, 
       filename = "Result_figures/abundance_analysis_plots/IC_abundance_boxplot_publication.pdf",
       width = 18, height = 12, units = "cm"
)

# Extract the legend
my_legend <- cowplot::get_legend(IS_genus_plot + 
                                   theme(
                                     legend.position = "right",
                                     legend.text = element_text(size = 7),
                                     legend.title = element_text(size =8, face="bold"),
                                     legend.justification = "center",
                                     legend.direction = "horizontal",
                                     legend.box.just = "bottom",
                                     plot.margin = unit(c(0, 0, 0, 0), "cm")
                                   )
)


grid_plot <- plot_grid(plotlist = list(IS_genus_plot + 
                                         labs(title = "Organ transplant recipient") + 
                                         theme(legend.position = "none",
                                               plot.title = element_text(size = 10, 
                                                                         face = "bold",
                                                                         hjust = 0.5)), 
                                       IC_genus_plot + labs(title = "Immunocompetent") + 
                                         theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               legend.position = "none",
                                               plot.title = element_text(size = 10, 
                                                                         face = "bold",
                                                                         hjust = 0.5
                                               ))),
                       ncol =2,rel_widths = c(1,.85))
grid_plot <- plot_grid(plotlist = list(grid_plot, my_legend), 
                       nrow = 2, rel_heights = c(1,.1))

ggsave(plot = grid_plot, 
       filename = "Result_figures/abundance_analysis_plots/IS_IC_abundance_boxplot_publication.pdf",
       width = 30, height = 12, units = "cm"
)

ggsave(plot = grid_plot, 
       filename = "Result_figures/abundance_analysis_plots/IS_IC_abundance_boxplot_publication.svg",
       width = 30, height = 12, units = "cm",device = "svg"
)


# -----------------------------------------------------------------
# Testing - calculate SCC and SCC_PL vs forearm equivalents separately.

# genus_significances_location_IS.df <- data.frame()
# for (cohort in c("organ transplant recipient")){
#   for (tax in unique(genus_rel_processed.df$taxonomy_genus)){
#     for (st in c("SCC", "SCC_PL")){
#       forearm_abundances.v <- 
#         genus_rel_processed.df %>% 
#         filter(Location == "Forearm",
#                Sample_type == st, 
#                Cohort == cohort,
#                taxonomy_genus == tax) %>% 
#         pull(Relative_abundance)
#       
#       all_body_sites_abundances.v <- 
#         genus_rel_processed.df %>% 
#         filter(Sample_type == st, 
#                Cohort == cohort,
#                taxonomy_genus == tax) %>% 
#         pull(Relative_abundance)
#       p_value <- wilcox.test(forearm_abundances.v,all_body_sites_abundances.v, exact = F)$p.value
#       # print(cohort)
#       # print(st)
#       # print(tax)
#       # print(p_value)
#       genus_significances_location_IS.df <- rbind(genus_significances_location_IS.df, data.frame(cohort,st,tax,p_value))
#     }
#   }
#   genus_significances_location_IS.df <- genus_significances_location_IS.df[!is.na(genus_significances_location_IS.df$p_value),]
#   genus_significances_location_IS.df$padj <- round(p.adjust(genus_significances_location_IS.df$p_value,method = "BH"),6)
# }
# 
# genus_significances_location_IC.df <- data.frame()
# for (cohort in c("organ transplant recipient")){
#   for (tax in unique(genus_rel_processed.df$taxonomy_genus)){
#     for (st in c("SCC", "SCC_PL")){
#       forearm_abundances.v <- 
#         genus_rel_processed.df %>% 
#         filter(Location == "Forearm",
#                Sample_type == st, 
#                Cohort == cohort,
#                taxonomy_genus == tax) %>% 
#         pull(Relative_abundance)
#       
#       all_body_sites_abundances.v <- 
#         genus_rel_processed.df %>% 
#         filter(Sample_type == st, 
#                Cohort == cohort,
#                taxonomy_genus == tax) %>% 
#         pull(Relative_abundance)
#       p_value <- wilcox.test(forearm_abundances.v,all_body_sites_abundances.v, exact = F)$p.value
#       # print(cohort)
#       # print(st)
#       # print(tax)
#       # print(p_value)
#       genus_significances_location_IC.df <- rbind(genus_significances_location_IC.df, data.frame(cohort,st,tax,p_value))
#     }
#   }
#   genus_significances_location_IC.df <- genus_significances_location_IC.df[!is.na(genus_significances_location_IC.df$p_value),]
#   genus_significances_location_IC.df$padj <- round(p.adjust(genus_significances_location_IC.df$p_value,method = "BH"),6)
# }

# Mann Whitney comparison between locations for SCC/SCC_PL
# genus_significances_location.df <-
#   genus_rel_processed.df %>%
#   filter(Sample_type %in% c("SCC", "SCC_PL")) %>%
#   group_by(Cohort,Sample_type,taxonomy_genus) %>%
#   nest() %>%
#   mutate(max_abundance = map(data, ~max(.x$Relative_abundance))) %>%
#   filter(max_abundance > 0) %>%
#   mutate(pvalue = map(data, ~wilcox.test(Relative_abundance ~ Location, exact = F, data = .)$p.value)) %>%
#   unnest(pvalue) %>%
#   select(-data, -max_abundance) %>%
#   mutate(padj = round(p.adjust(pvalue,method = "BH"),6)) %>%
#   # filter(padj <= 0.05) %>%
#   as.data.frame()
# genus_significances_location.df[order(genus_significances_location.df$padj, decreasing = F),] %>% filter(grepl(paste(specific_taxa.v, collapse = "|"), taxonomy_genus))

# -----------------------------------------------------------------

