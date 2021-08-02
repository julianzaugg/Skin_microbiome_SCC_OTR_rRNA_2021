
library(tidyverse)
library(ComplexHeatmap)


genus_relabeller_function <- function(my_labels){
  unlist(lapply(my_labels, 
                function(x) {
                  phylostring <- unlist(strsplit(x, split = ";"))
                  paste0("(",phylostring[2],") ",phylostring[5],";",phylostring[6])
                }))
}


source("code/utility.R")


# Load the processed metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", sep =",", header = T)
rownames(metadata.df) <- metadata.df$Index


# AND load the processed metadata that is unfiltered. This is for the qPCR results.
metadata_unfiltered.df <- read.csv("Result_tables/processed_unfiltered_metadata.csv", sep =",", header = T)
rownames(metadata_unfiltered.df) <- metadata_unfiltered.df$Index
metadata_unfiltered.df <- metadata_unfiltered.df[is.na(metadata_unfiltered.df$qPCR_exclude),]
metadata_unfiltered.df %>% group_by(Cohort, Sample_type) %>% tally()
metadata_unfiltered.df <- metadata_unfiltered.df[,c("Index","Subject" ,"Sample_type","Swab_ID", "Cohort","Sample_type_colour","Sample_type_shape","S_aureus_qPCR","Staph_spp_qPCR", "qPCR_16S")]

# metadata.df %>% filter(Subject == "MS007")
# metadata_unfiltered.df %>% filter(Subject == "MS007")
# Remove negative samples (if present)
metadata.df <- metadata.df[!metadata.df$Sample_type == "negative",]
metadata_unfiltered.df <- metadata_unfiltered.df[! metadata_unfiltered.df$Sample_type == "negative",]

# Calculate Mean_qPCR_16S per subject
mean_qPCR_16S.df <-
  metadata_unfiltered.df %>% 
  dplyr::group_by(Subject, Sample_type) %>%
  dplyr::summarise(Mean_qPCR_16S = mean(qPCR_16S)) %>%
  as.data.frame()

# Add to metadata.df
metadata.df <- left_join(metadata.df, mean_qPCR_16S.df, by  = c("Subject", "Sample_type"))


# Set levels for sample type and Cohort
metadata.df$Sample_type <- factor(metadata.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))
metadata.df$Cohort <- factor(metadata.df$Cohort, levels = c("immunocompetent", "organ transplant recipient"))


# Load abundance data
genus_rel.df <- read.csv("Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", header = T)
genus.m <- df2m(read.csv("Result_tables/count_tables/Genus_counts.csv", header = T))
genus_rel.m <- df2m(genus_rel.df)

# Create taxonomy label
genus_data.df <- separate(genus_rel.df, "taxonomy_genus", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus"), remove =F, sep = ";")

# Melt and combine with metadata
genus_data.df <- melt(genus_data.df, variable.name = "Index", value.name = "Relative_abundance")
genus_data.df <- left_join(genus_data.df, metadata.df, by = "Index")

# Set levels for sample type
genus_data.df$Sample_type <- factor(genus_data.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))

# Create taxonomy label
genus_data.df$taxonomy_label <- with(genus_data.df, paste0(Domain,";", Class,";",Family,";", Genus))
genus_data.df$taxonomy_label <- gsub("[a-z]__", "", genus_data.df$taxonomy_label)
# table(genus_data.df$taxonomy_label)[table(genus_data.df$taxonomy_label) != 388]

# Create cohort specific datasets
immunosuppressed_genus_data.df <- subset(genus_data.df, Cohort == "organ transplant recipient")
immunocompetent_genus_data.df <- subset(genus_data.df, Cohort == "immunocompetent")


# Generate full genus summary for each sample type, for each cohort
immunosuppressed_genus_summary.df <- 
  immunosuppressed_genus_data.df %>% 
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

# Generate genus summary for each subject + sample type
subject_genus_summary.df <-
  genus_data.df %>%
  dplyr::group_by(Subject,Sample_type, Cohort,taxonomy_genus, taxonomy_label) %>%
  dplyr::summarise(
  # Mean_relative_abundance2 = round(mean(sum(Relative_abundance)/max(Samples_total)), 5),
  Min_relative_abundance = round(min(Relative_abundance),5),
  Max_relative_abundance = round(max(Relative_abundance),5),
  Median_relative_abundance = round(median(Relative_abundance), 5),
  Mean_relative_abundance = round(mean(Relative_abundance), 5),
  ) %>%
  arrange(dplyr::desc(Mean_relative_abundance))

# Identify the top genus for each lesion type for both cohorts
number_of_top_taxa <- 10

immunosuppressed_top_genus_summary.df <-
  immunosuppressed_genus_summary.df %>% dplyr::group_by(Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)

immunocompetent_top_genus_summary.df <-
  immunocompetent_genus_summary.df %>% dplyr::group_by(Sample_type) %>% top_n(number_of_top_taxa, Mean_relative_abundance)

# Create a unique list of the top taxa across all the cohorts / groups
both_cohorts_lesions_top_genus <- unique(c(immunosuppressed_top_genus_summary.df$taxonomy_genus,
                                           immunocompetent_top_genus_summary.df$taxonomy_genus))

# ------------------------------------------------------------
# Generate filtered abundance matrices
genus_top_genus_rel.m <- genus_rel.m[both_cohorts_lesions_top_genus,]
immunocompetent_genus_top_genus_rel.m <- genus_rel.m[both_cohorts_lesions_top_genus, metadata.df %>% filter(Cohort == "immunocompetent") %>% pull(Index) ]
immunosuppressed_genus_top_genus_rel.m <- genus_rel.m[both_cohorts_lesions_top_genus, metadata.df %>% filter(Cohort == "organ transplant recipient") %>% pull(Index) ]

create_subject_sample_type_matrix <- function(sample_type){
  subject_genus_summary.df %>% 
    dplyr::filter(Sample_type == sample_type, taxonomy_genus %in% both_cohorts_lesions_top_genus) %>%
    dplyr::ungroup() %>%
    dplyr::select(Subject, taxonomy_genus, Mean_relative_abundance) %>%
    spread(Subject, Mean_relative_abundance,fill = 0) %>% as.data.frame() %>% df2m
}

create_subject_sample_type_metadata <- function(sample_type){
  subject_metadata.df <- unique(metadata.df[,c("Subject", "Subject_colour", "Cohort","Cohort_colour", "Sample_type", "Sample_type_colour","Mean_qPCR_16S"), drop = F])
  subject_metadata.df <- subject_metadata.df %>% filter(Sample_type == sample_type)
  rownames(subject_metadata.df) <- subject_metadata.df$Subject
  subject_metadata.df <- subject_metadata.df[order(subject_metadata.df$Subject),]
  subject_metadata.df
}

subject_scc.m <- create_subject_sample_type_matrix("SCC")
subject_scc_pl.m <- create_subject_sample_type_matrix("SCC_PL")
subject_ak.m <- create_subject_sample_type_matrix("AK")
subject_pds.m <- create_subject_sample_type_matrix("PDS")
subject_ns.m <- create_subject_sample_type_matrix("NS")
# dim(subject_scc.m)
# dim(subject_scc_pl.m)
# dim(subject_ak.m)
# dim(subject_scc.m)
# dim(subject_pds.m)

# Create metadata sheets for each subset
subject_scc_metadata.df <- create_subject_sample_type_metadata("SCC")
subject_scc_pl_metadata.df <- create_subject_sample_type_metadata("SCC_PL")
subject_ak_metadata.df <- create_subject_sample_type_metadata("AK")
subject_pds_metadata.df <- create_subject_sample_type_metadata("PDS")
subject_ns_metadata.df <- create_subject_sample_type_metadata("NS")

# ------------------------------------------------------------
# genus_row_labels.df <- data.frame(rownames(genus_rel.m), genus_relabeller_function(rownames(genus_rel.m)))
genus_row_labels.df <- unique(genus_data.df[, c("taxonomy_genus", "taxonomy_label")])
source("code/utility.R")
temp <- make_heatmap(log(subject_scc.m*100+0.0000001,2),
                     row_labels.df =genus_row_labels.df,
                     # col_labels = genus_col_labels.df,
                     metadata.df = subject_scc_metadata.df,
                     annotation_variables = c("Cohort"),
                     annotation_palette = NULL,
                     show_top_annotation = T,
                     annotation_bar_name_size = 6,
                     
                     column_split = sort(subject_scc_metadata.df$Cohort),
                     column_gap = unit(.5, "cm"),
                     plot_title = "SCC",
                     
                     grid_colour = "grey20",
                     # grid_colour = "grey10",
                     grid_thickness = .1,
                     
                     
                     row_title = "Domain;Class;Family;Genus",
                     column_title = "Subject",
                     plot_height = 6,
                     plot_width = 10,
                     my_padding = unit(c(0,0,0,0),"cm"),
                     cluster_columns = F,
                     cluster_rows = T,
                     show_column_dend = F,
                     show_row_dend = F,
                     
                     legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, ">= 60"),
                     my_breaks = unlist(lapply(c(0.000001, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,log,2)),
                     discrete_legend = T,
                     
                     legend_title = expression("Log"[2]~"(Mean relative abundance %)"),
                     default_palette_choice = 'blue',
                     # heatmap_palette = colorRampPalette(c("white","#17468a","#ffdd47","#99113a"))(10),
                     # default_palette_choice = 'bluered',
                     # heatmap_palette = c("cyan", "blue", "violet","red"),
                     row_dend_width = unit(4, "cm"),
                     column_dend_height= unit(1, "cm"),
                     simple_anno_size = unit(.25, "cm"),
                     column_names_rot = 45,
                     col_name_size = 6,
                     show_cell_values = F,
                     
                     filename = "Result_figures/heatmaps/subject_scc_genus_relative_abundance.pdf",
)

make_heatmap(log(subject_scc_pl.m*100+0.0000001,2),
             row_labels.df = genus_row_labels.df,
             # col_labels = genus_col_labels.df,
             metadata.df = subject_scc_pl_metadata.df,
             show_top_annotation = T,
             annotation_variables = c("Cohort"),
             annotation_palette = NULL,
             annotation_bar_name_size = 6,
             
             column_split = sort(subject_scc_pl_metadata.df$Cohort),
             column_gap = unit(.5, "cm"),
             plot_title = "SCC_PL",
             
             grid_colour = "grey20",
             grid_thickness = .1,
             
             row_title = "Domain;Class;Family;Genus",
             column_title = "Subject",
             plot_height = 6,
             plot_width = 10,
             my_padding = unit(c(0,0,0,0),"cm"),
             cluster_columns = F,
             cluster_rows = T,
             show_column_dend = T,
             show_row_dend = F,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, ">= 60"),
             my_breaks = unlist(lapply(c(0.000001, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,log,2)),
             discrete_legend = T,
             
             legend_title = expression("Log"[2]~"(Mean relative abundance %)"),
             default_palette_choice = 'blue',
             #heatmap_palette = colorRampPalette(c("white","#17468a","#ffdd47","#99113a"))(10),
             row_dend_width = unit(4, "cm"),
             column_dend_height= unit(1, "cm"),
             simple_anno_size = unit(.25, "cm"),
             column_names_rot = 45,
             col_name_size = 6,
             
             filename = "Result_figures/heatmaps/subject_scc_pl_genus_relative_abundance.pdf",
)


temp <- make_heatmap(log(subject_ak.m*100+0.0000001,2),
             row_labels.df = genus_row_labels.df,
             # col_labels = genus_col_labels.df,
             metadata.df = subject_ak_metadata.df,
             show_top_annotation = T,
             annotation_variables = c("Cohort"),
             annotation_palette = NULL,
             annotation_bar_name_size = 6,
             
             column_split = sort(subject_ak_metadata.df$Cohort),
             column_gap = unit(.5, "cm"),
             plot_title = "AK",
             
             grid_colour = "grey20",
             grid_thickness = .1,
             
             row_title = "Domain;Class;Family;Genus",
             column_title = "Subject",
             plot_height = 6,
             plot_width = 10,
             my_padding = unit(c(0,0,0,0),"cm"),
             cluster_columns = F,
             cluster_rows = T,
             show_column_dend = T,
             show_row_dend = F,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, ">= 60"),
             my_breaks = unlist(lapply(c(0.000001, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,log,2)),
             discrete_legend = T,
             
             legend_title = expression("Log"[2]~"(Mean relative abundance %)"),
             default_palette_choice = 'blue',
             #heatmap_palette = colorRampPalette(c("white","#17468a","#ffdd47","#99113a"))(10),
             row_dend_width = unit(4, "cm"),
             column_dend_height= unit(1, "cm"),
             simple_anno_size = unit(.25, "cm"),
             column_names_rot = 45,
             col_name_size = 6,
             
             filename = "Result_figures/heatmaps/subject_ak_genus_relative_abundance.pdf",
)

make_heatmap(log(subject_pds.m*100+0.0000001,2),
             row_labels.df = genus_row_labels.df,
             # col_labels = genus_col_labels.df,
             metadata.df = subject_pds_metadata.df,
             show_top_annotation = T,
             annotation_variables = c("Cohort"),
             annotation_palette = NULL,
             annotation_bar_name_size = 6,
             
             column_split = sort(subject_pds_metadata.df$Cohort),
             column_gap = unit(.5, "cm"),
             plot_title = "PDS",
             
             grid_colour = "grey20",
             grid_thickness = .1,
             
             row_title = "Domain;Class;Family;Genus",
             column_title = "Subject",
             plot_height = 6,
             plot_width = 11,
             my_padding = unit(c(0,0,0,0),"cm"),
             cluster_columns = F,
             cluster_rows = T,
             show_column_dend = T,
             show_row_dend = F,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, ">= 60"),
             my_breaks = unlist(lapply(c(0.000001, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,log,2)),
             discrete_legend = T,
             
             legend_title = expression("Log"[2]~"(Mean relative abundance %)"),
             default_palette_choice = 'blue',
             #heatmap_palette = colorRampPalette(c("white","#17468a","#ffdd47","#99113a"))(10),
             row_dend_width = unit(4, "cm"),
             column_dend_height= unit(1, "cm"),
             simple_anno_size = unit(.25, "cm"),
             column_names_rot = 45,
             col_name_size = 6,
             
             filename = "Result_figures/heatmaps/subject_pds_genus_relative_abundance.pdf",
)

make_heatmap(log(subject_ns.m*100+0.0000001,2),
             row_labels.df = genus_row_labels.df,
             # col_labels = genus_col_labels.df,
             metadata.df = subject_ns_metadata.df,
             show_top_annotation = T,
             annotation_variables = c("Cohort"),
             annotation_palette = NULL,
             annotation_bar_name_size = 6,
             
             column_split = sort(subject_ns_metadata.df$Cohort),
             column_gap = unit(.5, "cm"),
             plot_title = "NS",
             
             grid_colour = "grey20",
             grid_thickness = .1,
             
             row_title = "Domain;Class;Family;Genus",
             column_title = "Subject",
             plot_height = 5,
             plot_width = 5.5,
             my_padding = unit(c(0,0,0,0),"cm"),
             cluster_columns = F,
             cluster_rows = T,
             show_column_dend = T,
             show_row_dend = F,
             
             legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.5,.1))*100, ">= 60"),
             my_breaks = unlist(lapply(c(0.000001, 0.001, 0.005,0.05, seq(.1,.6,.1))*100,log,2)),
             discrete_legend = T,
             
             legend_title = expression("Log"[2]~"(Mean relative abundance %)"),
             default_palette_choice = 'blue',
             #heatmap_palette = colorRampPalette(c("white","#17468a","#ffdd47","#99113a"))(10),
             row_dend_width = unit(4, "cm"),
             column_dend_height= unit(1, "cm"),
             simple_anno_size = unit(.25, "cm"),
             column_names_rot = 45,
             col_name_size = 6,
             
             filename = "Result_figures/heatmaps/subject_ns_genus_relative_abundance.pdf",
)

# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------
# TESTING
# hm <- Heatmap(subject_ak.m,cluster_columns = F)
# ha <- HeatmapAnnotation(Cohort = as.character(subject_ak_metadata.df[colnames(subject_ak.m),c("Cohort")]),
#                   Mean_qPCR_16S = anno_barplot(subject_ak_metadata.df[colnames(subject_ak.m),c("Mean_qPCR_16S")]),
#                   which = "column",
#                   col = list(Cohort = c("immunocompetent" = "blue", 
#                                            "organ transplant recipient" = "red")
#                              )
#                   )
# ha <- HeatmapAnnotation(subject_ak_metadata.df[,"Cohort",drop = F],
#                         which = "column",
#                         col = list(Cohort = c("immunocompetent" = "blue", 
#                                               "organ transplant recipient" = "red")
#                         )
# )
# draw(ha %v% temp$heatmap)
# draw(ha %v% hm)
# graphics.off()
# subject_scc_metadata.df
# columnAnnotation()
# 
# 
# 
# temp_meta <- subject_scc_metadata.df[,c("Cohort"),drop = F]
# named_colour_list <- setNames(as.character(subject_scc_pl_metadata.df[, c("Cohort_colour")]), as.character(subject_scc_pl_metadata.df[,"Cohort"]))
# colour_lists <- list()
# colour_lists[["Cohort"]] <- named_colour_list
# 
# draw(anno_barplot(subject_scc_metadata.df$Mean_qPCR_16S))
# temp_meta <- subject_scc_metadata.df[,c("Cohort", "Mean_qPCR_16S")]
# temp_meta <- temp_meta[order(temp_meta$Cohort),]
# ha <- columnAnnotation(df = temp_meta[,c("Cohort", "Mean_qPCR_16S")])
# draw(ha %v% temp$heatmap)
# draw(ha)
# graphics.off()
# metadata.df
# 
# meta_test <- metadata_unfiltered.df %>% 
#   dplyr::group_by(Subject, Sample_type) %>%
#   dplyr::summarise(Mean_qPCR_16S = mean(qPCR_16S)) %>%
#   as.data.frame()
# ha <- HeatmapAnnotation(foo = anno_barplot(meta_test$Mean_qPCR_16S))
# 
# draw(ha)
# draw(temp$heatmap %v% ha)
# graphics.off()


  

# make_heatmap(log(immunosuppressed_genus_top_genus_rel.m*100+0.0000001,2),
#              # row_labels.df =genus_row_labels.df,
#              # col_labels = genus_col_labels.df,
#              metadata.df = metadata.df %>% filter(Cohort == "organ transplant recipient"),
#              annotation_variables = c("Sample_type", "Subject"),
#              annotation_palette = NULL,
#              
#              
#              grid_colour = "grey20",
#              grid_thickness = .1,
#              
#              row_title = "Genus",
#              plot_height = 7,
#              plot_width = 30,
#              my_padding = unit(c(2,2,2,2),"cm"),
#              cluster_columns = F,
#              cluster_rows = T,
#              show_column_dend = F,
#              show_row_dend = F,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.3,.1))*100, ">= 40"),
#              my_breaks = unlist(lapply(c(0.000001, 0.001, 0.005,0.05, seq(.1,.4,.1))*100,log,2)),
#              discrete_legend = T,
#              
#              legend_title = "Relative abundance %",
#              default_palette_choice = 'blue',
#              row_dend_width = unit(4, "cm"),
#              column_dend_height= unit(1, "cm"),
#              simple_anno_size = unit(.25, "cm"),
#              column_names_rot = 45,
#              col_name_size = 5,
#              
#              filename = "Result_figures/heatmaps/IS_genus_relative_abundance.pdf",
# )
# 
# 
# make_heatmap(log(immunocompetent_genus_top_genus_rel.m*100+0.0000001,2),
#              # row_labels.df =genus_row_labels.df,
#              # col_labels = genus_col_labels.df,
#              metadata.df = metadata.df %>% filter(Cohort == "immunocompetent"),
#              annotation_variables = c("Sample_type", "Subject"),
#              annotation_palette = NULL,
#              
#              
#              grid_colour = "grey20",
#              grid_thickness = .1,
#              
#              row_title = "Genus",
#              plot_height = 7,
#              plot_width = 30,
#              my_padding = unit(c(2,2,2,2),"cm"),
#              cluster_columns = F,
#              cluster_rows = T,
#              show_column_dend = F,
#              show_row_dend = F,
#              
#              legend_labels = c(c(0, 0.001, 0.005,0.05, seq(.1,.3,.1))*100, ">= 40"),
#              my_breaks = unlist(lapply(c(0.000001, 0.001, 0.005,0.05, seq(.1,.4,.1))*100,log,2)),
#              discrete_legend = T,
#              
#              legend_title = "Relative abundance %",
#              default_palette_choice = 'blue',
#              row_dend_width = unit(4, "cm"),
#              column_dend_height= unit(1, "cm"),
#              simple_anno_size = unit(.25, "cm"),
#              column_names_rot = 45,
#              col_name_size = 5,
#              
#              filename = "Result_figures/heatmaps/IC_genus_relative_abundance.pdf",
# )



