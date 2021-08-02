library(vegan)
library(reshape2)


source("code/utility.R")

# Load the processed metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", sep =",", header = T)
rownames(metadata.df) <- metadata.df$Index

# Factorise discrete columns
discrete_variables <- c("Sample_type","Gender","Subject", "Cohort", 
                        "Length_of_immunosuppression_group_1", "Length_of_immunosuppression_group_2")

metadata.df$Sample_type <- factor(metadata.df$Sample_type, levels = c("NS","PDS", "AK", "SCC_PL", "SCC"))
metadata.df$Cohort <- factor(metadata.df$Cohort, levels = c("organ transplant recipient", "immunocompetent"))
metadata.df$Gender <- factor(metadata.df$Gender)

# Need to factorise the colour columns as well
# colour_columns <- names(metadata.df)[grepl("colour", names(metadata.df))]
# metadata.df[colour_columns] <- lapply(metadata.df[colour_columns], factor)

# Load the ASV - taxonomy mapping file
asv_taxonomy_map.df <- read.csv("Result_tables/asv_taxonomy_map.csv", header = T)

# Load count data
asv.m <- df2m(read.csv("Result_tables/count_tables/ASV_counts.csv", header = T))
genus.m <- df2m(read.csv("Result_tables/count_tables/Genus_counts.csv", header = T))

# Order the matrices and metadata to be the same order
metadata.df <- metadata.df[order(rownames(metadata.df)),]
asv.m <- asv.m[,order(rownames(metadata.df))]
genus.m <- genus.m[,order(rownames(metadata.df))]

# asv.m <- asv.m[,colSums(asv.m) > 5000]
# genus.m <- genus.m[,colSums(genus.m) > 5000]
# metadata.df <- metadata.df[colnames(genus.m),]

# CLR transform the count matrices
asv_clr.m <- clr(asv.m)
genus_clr.m <- clr(genus.m)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

asv_pca <- rda(t(asv_clr.m))
genus_pca <- rda(t(genus_clr.m))
  
# Immunocompetent, all sample types
immunocompetent_samples <- as.character(metadata.df$Index[metadata.df$Cohort == "immunocompetent"])
immunocompetent_asv_pca <- rda(t(asv_clr.m[,immunocompetent_samples]))
immunocompetent_genus_pca <- rda(t(genus_clr.m[,immunocompetent_samples]))

# Immunosuppressed, all sample types
immunosuppressed_samples <- as.character(metadata.df$Index[metadata.df$Cohort == "organ transplant recipient"])
immunosuppressed_asv_pca <- rda(t(asv_clr.m[,immunosuppressed_samples]))
immunosuppressed_genus_pca <- rda(t(genus_clr.m[,immunosuppressed_samples]))

# ------------------------------------------------
# Testing

temp <- capscale(t(genus.m[,immunosuppressed_samples])~1, data = metadata.df[immunosuppressed_samples,],distance = "jaccard") # ~1 makes it unconstrained
generate_pca(pca_object = temp,
             metadata.df = metadata.df[immunosuppressed_samples,],
             variable_to_plot = "Sample_type",
             variable_colours_available = T,
             variable_shapes_available = F,
             plot_spiders =T,
             spider_alpha = .5,
             plot_hulls =F,
             plot_ellipses = F)

generate_pca(pca_object = immunosuppressed_genus_pca,
             metadata.df = metadata.df[immunosuppressed_samples,],
             variable_to_plot = "Sample_type",
             variable_colours_available = T,
             variable_shapes_available = F,
             plot_spiders = T,
             spider_alpha = .5,
             plot_hulls =F,
             plot_ellipses = F)

temp <- with(metadata.df[immunocompetent_samples,], 
             betadisper(vegdist(t(genus.m[,immunocompetent_samples]), method = "jaccard"), 
                        group = Sample_type))
plot(temp,hull = F, label.cex = .4)
temp <- with(metadata.df[immunocompetent_samples,], 
             betadisper(vegdist(t(genus_clr.m[,immunocompetent_samples]), method = "euclidean"), 
                        group = Sample_type))
plot(temp,hull = F, label.cex = .4)

temp <- with(metadata.df[immunosuppressed_samples,], 
             betadisper(vegdist(t(genus.m[,immunosuppressed_samples]), method = "jaccard"), 
                        group = Sample_type))
plot(temp,hull = F, label.cex = .4)
temp <- with(metadata.df[immunosuppressed_samples,], 
             betadisper(vegdist(t(genus_clr.m[,immunosuppressed_samples]), method = "euclidean"), 
                        group = Sample_type))
plot(temp,hull = F, label.cex = .4)
boxplot(temp, main = "", xlab = "")
vegan::permutest(temp, permutations = 9999, parallel = 2,pairwise = T)
vegan::permutest(temp, permutations = 999, parallel = 2)

permanova_result <- adonis(t(genus.m[,immunosuppressed_samples])~Sample_type,
                 method = "jaccard",
                 data = metadata.df[immunosuppressed_samples,], 
                 permutations = 999)
dist_matrix <- vegdist(t(genus.m[,immunosuppressed_samples]), method = "jaccard")
betadisper_object <- with(metadata.df[immunosuppressed_samples,], betadisper(dist_matrix, group = Sample_type))
permdisp_results <- permutest(betadisper_object, permutations = 999, parallel = 2,pairwise = T)
permanova_result
permdisp_results

permanova_result <- adonis(t(genus_clr.m[,immunosuppressed_samples])~Sample_type,
                           method = "euclidean",
                           data = metadata.df[immunosuppressed_samples,], 
                           permutations = 999)
dist_matrix <- vegdist(t(genus_clr.m[,immunosuppressed_samples]), method = "euclidean")
betadisper_object <- with(metadata.df[immunosuppressed_samples,], betadisper(dist_matrix, group = Sample_type))
permdisp_results <- permutest(betadisper_object, permutations = 999, parallel = 2,pairwise = T)
permanova_result
permdisp_results
# ------------------------------------------------


file_type <- "pdf"
# ------------------------------------------------------------------------------------
# All samples, immunocompetent

generate_pca(pca_object = immunocompetent_genus_pca,
             metadata.df = subset(metadata.df, Cohort == "immunocompetent"),
             variable_to_plot = "Sample_type",
             variable_colours_available = T,
             variable_shapes_available = T,
             legend_title = "Sample type",
             legend_columns = 1,
             legend_cex = .7,
             legend_x = -4, legend_y = 5,
             legend_key_text_distance = 1,
             
             use_fill_shapes = T,
             point_alpha = 1,
             point_line_thickness = .5,
             point_size = 0.8,
             point_line_is_darker_fill = T,
             
             plot_height = 5,
             plot_width = 5,
             plot_hulls = F,
             plot_spiders = F,
             plot_ellipses = F,
             plot_arrows = F,
             arrow_thickness = .5,
             arrow_colour = "grey60",
             arrow_label_size = .3,
             arrow_scalar = 2.5,
             arrow_label_font_type = 4,
             # arrow_label_colour = "red",
             arrow_label_colour = "royalblue4",
             num_top_species = 2,
             axis_limits = c(-4,6,-4,5.1),
             specie_labeller_function = first_resolved_taxonomy,
             file_type = file_type,
             filename = paste0("Result_figures/ordination_plots/genus/immunocompetent_sample_type.", file_type)
)

# Subject
generate_pca(pca_object = immunocompetent_genus_pca,
             metadata.df = subset(metadata.df, Cohort == "immunocompetent"),
             variable_to_plot = "Subject",
             variable_colours_available = T,
             include_legend = F,
             legend_title = "Subject",
             legend_columns = 1,
             legend_cex = .7,
             legend_x = -4, legend_y = 5,
             
             use_fill_shapes = T,
             point_alpha = 1,
             point_line_thickness = .5,
             point_size = 0.8,
             point_line_is_darker_fill = T,
             
             plot_height = 5,
             plot_width = 5,
             plot_hulls = F,
             plot_spiders = F,
             plot_ellipses = T,
             plot_arrows = F,
             arrow_thickness = .5,
             arrow_colour = "grey60",
             arrow_label_size = .3,
             arrow_scalar = 2.5,
             arrow_label_font_type = 4,
             # arrow_label_colour = "red",
             arrow_label_colour = "royalblue4",
             num_top_species = 2,
             axis_limits = c(-4,6,-4,5.1),
             specie_labeller_function = first_resolved_taxonomy,
             file_type = file_type,
             filename = paste0("Result_figures/ordination_plots/genus/immunocompetent_subject.", file_type)
)

# ------------------------------------------------------------------------------------
# All samples, organ transplant recipient
# Sample type
source("code/utility.R")
generate_pca(pca_object = immunosuppressed_genus_pca,
             metadata.df = subset(metadata.df, Cohort == "organ transplant recipient"),
             variable_to_plot = "Sample_type",
             variable_colours_available = T,
             variable_shapes_available = T,
             legend_title = "Sample type",
             legend_columns = 1,
             legend_cex = .7,
             legend_x = -5.5, legend_y = 6,
             legend_key_text_distance = 1,
             
             use_fill_shapes = T,
             point_alpha = 1,
             point_line_thickness = .5,
             point_size = 0.8,
             point_line_is_darker_fill = T,
             
             plot_height = 5,
             plot_width = 5,
             plot_hulls = F,
             plot_spiders = F,
             spider_alpha = .5,
             plot_ellipses = F,
             plot_arrows = F,
             arrow_thickness = .5,
             arrow_colour = "grey60",
             arrow_label_size = .3,
             arrow_scalar = 4,
             arrow_label_font_type = 4,
             # arrow_label_colour = "red",
             arrow_label_colour = "royalblue4",
             num_top_species = 2,
             axis_limits = c(-5.5,4,-5,6),
             specie_labeller_function = first_resolved_taxonomy,
             file_type = file_type,
             filename = paste0("Result_figures/ordination_plots/genus/immunosuppressed_sample_type.", file_type)
)

# Subject
generate_pca(pca_object = immunosuppressed_genus_pca,
             metadata.df = subset(metadata.df, Cohort == "organ transplant recipient"),
             variable_to_plot = "Subject",
             variable_colours_available = T,
             variable_shapes_available = F,
             include_legend = F,
             legend_title = "Subject",
             legend_columns = 1,
             legend_cex = .7,
             legend_x = -4, legend_y = 5,
             
             use_fill_shapes = T,
             point_alpha = 1,
             point_line_thickness = .5,
             point_size = 0.8,
             point_line_is_darker_fill = T,
             
             plot_height = 5,
             plot_width = 5,
             plot_hulls = F,
             plot_spiders = F,
             plot_ellipses = T,
             ellipse_alpha = 1,
             plot_arrows = F,
             arrow_thickness = .5,
             arrow_colour = "grey60",
             arrow_label_size = .3,
             arrow_scalar = 2.5,
             arrow_label_font_type = 4,
             # arrow_label_colour = "red",
             arrow_label_colour = "royalblue4",
             num_top_species = 2,
             axis_limits = c(-5.5,4,-5,6),
             specie_labeller_function = first_resolved_taxonomy,
             file_type = file_type,
             filename = paste0("Result_figures/ordination_plots/genus/immunosuppressed_subject.", file_type)
)


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Publication figure

make_publication_plot <- function(filetype = "pdf"){
  source("code/utility.R")
  graphics.off()
  if (filetype == "svg"){
    svg(filename = "Result_figures/ordination_plots/sample_type_subject_genus_PCA_for_publication.svg",width = 8,height = 7)
  }
  else{
    pdf(file = "Result_figures/ordination_plots/sample_type_subject_genus_PCA_for_publication.pdf",width = 8,height = 7)  
  }
  # par(mfrow = c(2,2))
  # par(mfrow = c(2,2),
  #     mar = c(1,10,5,1))
  par(mfrow = c(2,2),
      oma = c(1,4,0,0))
  # plot.new()
  par(mar = c(1,4,3,1))
  # par(bg = rgb(61, 55, 72, maxColorValue = 255))
  
  generate_pca(pca_object = immunosuppressed_genus_pca,
               metadata.df = subset(metadata.df, Cohort == "organ transplant recipient"),
               variable_to_plot = "Sample_type",
               variable_colours_available = T,
               variable_shapes_available = T,
               legend_title = "Sample type",
               legend_columns = 1,
               legend_cex = .7,
               legend_x = -5.5, legend_y = 6,
               legend_key_text_distance = 1,
               
               use_fill_shapes = T,
               point_alpha = 1,
               point_line_thickness = .5,
               point_size = 0.8,
               point_line_is_darker_fill = T,
               
               plot_x_tick_labels = F, 
               plot_height = 5,
               plot_width = 5,
               plot_hulls = F,
               plot_spiders = F,
               spider_alpha = .5,
               
               plot_ellipses = F,
               plot_arrows = F,
               arrow_thickness = .5,
               arrow_colour = "grey60",
               arrow_label_size = .3,
               arrow_scalar = 4,
               arrow_label_font_type = 4,
               # arrow_label_colour = "red",
               arrow_label_colour = "royalblue4",
               num_top_species = 2,
               axis_limits = c(-5.5,4,-5,6),
               specie_labeller_function = first_resolved_taxonomy,
               file_type = file_type,
  )
  
  title(outer=F,adj=.5,main="Organ transplant recipient",cex.main=4,col.main="black",font.main=2,line=1, cex.main =1)
  title(outer=T,adj=.77,ylab = "By sample type",font.lab=2, line = 2,cex.lab = 1)
  
  par(mar = c(1,4,3,1))
  
  
  generate_pca(pca_object = immunocompetent_genus_pca,
               metadata.df = subset(metadata.df, Cohort == "immunocompetent"),
               variable_to_plot = "Sample_type",
               variable_colours_available = T,
               variable_shapes_available = T,
               legend_title = "Sample type",
               legend_columns = 1,
               legend_cex = .7,
               legend_x = -4, legend_y = 5,
               legend_key_text_distance = 1,
               
               use_fill_shapes = T,
               point_alpha = 1,
               point_line_thickness = .5,
               point_size = 0.8,
               point_line_is_darker_fill = T,
               
               plot_x_tick_labels = F, 
               plot_height = 5,
               plot_width = 5,
               plot_hulls = F,
               plot_spiders = F,
               plot_ellipses = F,
               plot_arrows = F,
      
               arrow_thickness = .5,
               arrow_colour = "grey60",
               arrow_label_size = .3,
               arrow_scalar = 2.5,
               arrow_label_font_type = 4,
               # arrow_label_colour = "red",
               arrow_label_colour = "royalblue4",
               num_top_species = 2,
               axis_limits = c(-4,6,-4,5.1),
               specie_labeller_function = first_resolved_taxonomy,
               file_type = file_type,
  )
 
  title(outer=F,adj=.5,main="Immunocompetent",cex.main=4,col.main="black",font.main=2,line=1, cex.main =1)
  par(mar = c(4,4,0,1))
  generate_pca(pca_object = immunosuppressed_genus_pca,
               metadata.df = subset(metadata.df, Cohort == "organ transplant recipient"),
               variable_to_plot = "Subject",
               variable_colours_available = T,
               variable_shapes_available = F,
               include_legend = F,
               legend_title = "Subject",
               legend_columns = 1,
               legend_cex = .7,
               legend_x = -4, legend_y = 5,
               
               use_fill_shapes = T,
               point_alpha = 1,
               point_line_thickness = .5,
               point_size = 0.8,
               point_line_is_darker_fill = T,
               
               plot_height = 5,
               plot_width = 5,
               plot_hulls = F,
               plot_spiders = F,
               plot_ellipses = F,
               ellipse_alpha = 1,
               plot_arrows = F,
               arrow_thickness = .5,
               arrow_colour = "grey60",
               arrow_label_size = .3,
               arrow_scalar = 2.5,
               arrow_label_font_type = 4,
               # arrow_label_colour = "red",
               arrow_label_colour = "royalblue4",
               num_top_species = 2,
               axis_limits = c(-5.5,4,-5,6),
               specie_labeller_function = first_resolved_taxonomy,
               file_type = file_type,
  )
  
  title(outer=T,adj = 0.26, ylab = "By subject",font.lab=2, line = 2,cex.lab = 1)
  
  par(mar = c(4,4,0,1))
  generate_pca(pca_object = immunocompetent_genus_pca,
               metadata.df = subset(metadata.df, Cohort == "immunocompetent"),
               variable_to_plot = "Subject",
               variable_colours_available = T,
               include_legend = F,
               legend_title = "Subject",
               legend_columns = 1,
               legend_cex = .7,
               legend_x = -4, legend_y = 5,
               
               use_fill_shapes = T,
               point_alpha = 1,
               point_line_thickness = .5,
               point_size = 0.8,
               point_line_is_darker_fill = T,
               
               plot_height = 5,
               plot_width = 5,
               plot_hulls = F,
               plot_spiders = F,
               plot_ellipses = F,
               plot_arrows = F,
               arrow_thickness = .5,
               arrow_colour = "grey60",
               arrow_label_size = .3,
               arrow_scalar = 2.5,
               arrow_label_font_type = 4,
               # arrow_label_colour = "red",
               arrow_label_colour = "royalblue4",
               num_top_species = 2,
               axis_limits = c(-4,6,-4,5.1),
               specie_labeller_function = first_resolved_taxonomy,
               file_type = file_type,
  )
  
  # title(ylab = "test",outer = T)
  
  # title(main = "",outer = T,ylab = "Test")
  dev.off()
}

make_publication_plot("pdf")
make_publication_plot("svg")

