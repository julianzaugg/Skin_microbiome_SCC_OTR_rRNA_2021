# Utility functions for plotting, processing and re-formatting data

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(gplots)
library(grid)
library(phyloseq)
library(seqinr) # For writing fasta files
library(corrplot)
library(svglite)
library(FSA)
library(ggsignif)

# Function to convert a dataframe to a matrix. First column becomes row names
df2m <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}

# Function to convert a matrix to a dataframe. Rows become first column. 
# Specify 'column_name' to set the name of the new column (defaults to "Row_variable")
m2df <- function(mymatrix, column_name = "Row_variable"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, column_name] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(column_name,cur_names)]
  return(mydf)
}

# Make row names of a dataframe the value of a defined column (default 1st) and then remove the column
df_column_to_rownames <- function(mydf, rowname_col = 1){
  my_clean.df <- mydf
  rownames(my_clean.df) <- my_clean.df[,rowname_col]
  my_clean.df[,1] <- NULL
  return(my_clean.df)
}

# Assign colours to a dataframe for specified discrete columns.
assign_colours_to_df <- function(dataframe.df, 
                                 columns, 
                                 my_palette = NULL, 
                                 auto_assign =T, 
                                 my_default_palette=NULL){
  # TODO - make option for default continuous palette
  # Palettes
  colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
  colour_palette_10 <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
  colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba",
                         "#b67249","#9b4a6f","#df8398")
  colour_palette_20 <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782",
                         "#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
  colour_palette_30 <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8",
                         "#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a",
                         "#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
  colour_palette_206 <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff",
                          "#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d",
                          "#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c",
                          "#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b",
                          "#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff",
                          "#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4",
                          "#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00",
                          "#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8",
                          "#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100",
                          "#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8",
                          "#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f",
                          "#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff",
                          "#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3",
                          "#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff",
                          "#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614",
                          "#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec",
                          "#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd",
                          "#93f6e7","#01b4a4")
  # Create palette list for each column
  palette_list <- list()
  for (col in columns){
    column_values <- as.character(unique(dataframe.df[,col]))
    column_values <- column_values[!is.na(column_values)]
    n_values <- length(column_values)
    if (is.null(my_default_palette)){
      if (n_values <= 10){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_10)(n_values), column_values)
        default_colour_palette <- colour_palette_10
      } else if (n_values > 10 & n_values <= 15){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_15)(n_values), column_values)
        default_colour_palette <- colour_palette_15
      } else if (n_values > 15 & n_values <= 20){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_20)(n_values), column_values)
        default_colour_palette <- colour_palette_20
      } else if (n_values > 20 & n_values <= 30){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_30)(n_values), column_values)
        default_colour_palette <- colour_palette_30
      } else{
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_206)(n_values), column_values)
        default_colour_palette <- colour_palette_206
      }
    } else{
      default_colour_palette <- my_default_palette
    }
    
    if (typeof(my_palette) == "character"){
      if (length(my_palette) >= n_values){
        palette_list[[col]] <- my_palette  
      } else{
        stop(paste0("Provided palette has ", length(my_palette), " colours, but column ", col, " has ", n_values, ". Provide a longer palette."))
      }
    } 
    else if (is.null(my_palette) & auto_assign == T){
      palette_list[[col]] <- default_colour_palette
    }
    else if (typeof(my_palette) == "list"){
      if (! col %in% names(my_palette)){
        if (auto_assign == T){
          warning(paste0("Column ", col, " does not have a specified palette. Assigning a default palette"))
          palette_list[[col]] <- default_colour_palette
        } else{
          warning(paste0("Column ", col, " does not have a specified palette. Skipping."))          
        }
      } else{
        palette_list[[col]] <- my_palette[[col]]
      }
    }
    else {
      stop("Provided palette should be a character vector of colours or a list of separate character vectors, or set auto_assign=T")
    }
    # print(col)
    # print(palette_list)
    # print(palette_list[[col]])
    # print(typeof(palette_list[[col]]))
    # print(palette_list[col])
    # print(typeof(palette_list[col]))
    # print("*******")
    if (!is.null(names(palette_list[[col]]))){
      # print(col)
      col_colours <- palette_list[[col]]
    } else{
      # print(col)
      col_colours <- setNames(palette_list[[col]][1:n_values], column_values)
    }
    
    dataframe.df[,paste0(col, "_colour")] <- as.character(lapply(as.character(dataframe.df[,col]), function(x) as.character(col_colours[x])))      
  }
  dataframe.df
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


generate_tax_level_data <- function(feature_count_taxonomy_data.df, sample_ids, tax_string_levels, remove_zero_row_entries = F){
  # This function will take a dataframe with the feature counts (first column ID and remaining columns the sample counts and taxonomy labels),
  # a list of the samples and the taxa labels you wish to sum over
  # The output is a list with the counts and relative abundances at each specified taxonomy level
  
  df2matrix <- function(mydataframe){
    mymatrix <- mydataframe
    rownames(mymatrix) <- mydataframe[,1]
    mymatrix[,1] <- NULL
    mymatrix <- as.matrix(mymatrix)
    mymatrix
  }
  
  m2df <- function(mymatrix, column_name = "Row_variable"){
    mydf <- as.data.frame(mymatrix)
    cur_names <- names(mydf)
    mydf[, column_name] <- rownames(mydf)
    rownames(mydf) <- NULL
    mydf <- mydf[,c(column_name,cur_names)]
    return(mydf)
  }
  output <- list()
  for (tax_string_level in tax_string_levels){
    counts.df <- melt(feature_count_taxonomy_data.df[c(tax_string_level, sample_ids)], measure.vars = c(sample_ids))
    counts.df <- dcast(counts.df, get(tax_string_level) ~ variable, sum)
    names(counts.df)[1] <- tax_string_level
    
    # Now create the relative abundance matrix at the current taxa level
    abundances.m <- counts.df
    
    rownames(abundances.m) <- abundances.m[[tax_string_level]]
    abundances.m[tax_string_level] <- NULL
    abundances.m <- as.matrix(abundances.m)
    abundances.m <- t(t(abundances.m) / colSums(abundances.m))
    abundances.m[is.nan(abundances.m)] <- 0
    
    counts.m <- df2matrix(counts.df)
    if (remove_zero_row_entries == T){
      abundances.m <- abundances.m[apply(abundances.m, 1, max) > 0,]
      counts.m <- counts.m[apply(counts.m, 1, max) > 0,]
    }
    output[[tax_string_level]] <- list("counts" = counts.m, "abundances" = abundances.m)
  }
  output
}

abundance_summary <- function(abundance.df, top_n = NULL){
  # Generate summary from abundance dataframe. 
  # Assumes first column is the unique taxa/feature ID and other columns are abundances for each sample
  # If top_n specified, will return the top n taxa by mean abundance
  abundance_melt.df <- melt(abundance.df, variable.name = "Sample", value.name = "Relative_abundance")
  abundance_summary.df <- 
    abundance_melt.df %>% 
    dplyr::mutate(Samples_total = n_distinct(Sample)) %>% # number of unique samples/index
    dplyr::group_by_at(1,.drop = F) %>% 
    dplyr::mutate(Samples_present = sum(Relative_abundance > 0, na.rm = TRUE)) %>%
    dplyr::summarise(
      Samples_present = max(Samples_present),
      Samples_total = max(Samples_total),
      # Mean_relative_abundance2 = round(mean(sum(Relative_abundance)/max(Samples_total)), 5),
      Min_relative_abundance = round(min(Relative_abundance),5),
      Max_relative_abundance = round(max(Relative_abundance),5),
      Mean_relative_abundance = round(mean(Relative_abundance), 5),
      Median_relative_abundance = round(median(Relative_abundance), 5),
    ) %>%
    arrange(dplyr::desc(Mean_relative_abundance))
  if (!is.null(top_n)){
    abundance_summary.df <- 
      abundance_summary.df %>% 
      arrange(dplyr::desc(Mean_relative_abundance)) %>%
      head(top_n)
  }
  as.data.frame(abundance_summary.df)
}



collapse_abundance_dataframe <- function(abundance.df, filtering_taxa.v = NULL){
  # Collapse abundance dataframe. 
  # Assumes first column is the unique taxa/feature ID and other columns are abundances for each sample
  # If a list of filtering_taxa.v is provided, taxa not in list are assigned a label 'Other' and their abundance
  # values summed together
  abundance.df <- melt(abundance.df, variable.name = "Sample", value.name = "Relative_abundance")
  if (!is.null(filtering_taxa.v)){
    abundance.df[!abundance.df[,1] %in% filtering_taxa.v,][,1] <- "Other"
  }
  
  abundance.df <- 
    abundance.df %>% 
    group_by_at(c(2,1)) %>% 
    summarise_all(list(sum)) %>%
    # summarise_all(funs(sum)) %>%
    as.data.frame()
  abundance.df
}


first_resolved_taxonomy <- function(x) {
  # Get the lowest taxonomy level (or 'first') that has been resolved from a taxonomy string
  ranks.v <- unlist(strsplit(x, split = ';'))
  # print(x)
  # print(ranks.v)
  for (i in length(ranks.v):1) {
    split.v <- unlist(strsplit(ranks.v[i], split = '__'))
    if (!is.na(split.v[2]) & split.v[2] != "") {
      if (!grepl(split.v[2], pattern = "Unassigned|uncultured")) {
        return(ranks.v[i])
      }
    }
  }
  return(ranks.v[1])
}


calculate_PC_taxa_contributions <- function(pca_object){
  # Calculate the percentage contribution from each taxa for PC1-3. Requires unscaled values that are squared
  pc1_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,1]^2, 3),value.name = "PC1_contribution_percentage")
  pc2_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC2_contribution_percentage")
  pc3_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC3_contribution_percentage")
  
  data.frame(pc1_contribution, pc2_contribution, pc3_contribution)
}

generate_pca <- function(pca_object, # PCA object generated by the rda() function from vegan package
                         metadata.df, # Dataframe containing metadata
                         variable_to_plot, # Variable (column) that we are plotting
                         variable_colours_available = F, 
                         my_colour_palette = NULL,
                         my_levels = NULL,
                         
                         axis_limits = NULL, # list of axes limits c(xmin, xmax,ymin,ymax)
                         hide_grid = F, # Hide the grid
                         highlight_origin = T, # Highlight the origin line
                         show_x_label = T, # Show x-axis label
                         show_y_label = T, # Show y-axis label
                         plot_x_ticks = T, # Show x-axis tick marks
                         plot_y_ticks = T, # Show y-axis tick marks
                         plot_x_tick_labels = T, # Show x-axis tick labels
                         plot_y_tick_labels = T, # Show x-axis tick labels
                         plot_title = NULL, # Show plot title
                         title_cex = 1, # Size of title
                         
                         component_choices = c(1,2),
                         
                         plot_arrows = F, # Show contributing species/OTUs vectors
                         num_top_species = 5, # How many species/OTUs to show with the largest contributions to the variance
                         arrow_colour = "black", # Color of arrows
                         arrow_alpha = 1, # Transparency of arrows, 0 = invisible, 1= fully, visible, 0.5 = 50% transparent
                         arrow_thickness = .2, # Thickness of arrow
                         label_arrows=T, # Label the arrows
                         arrow_label_size = .5, # Size of arrow label
                         arrow_scalar = 1, # Scale of the arrows
                         arrow_label_colour = "black", # Colour of arrow label
                         arrow_label_font_type = 1, # Font type of arrow label : 1 Normal, 2 bold, 3 italic, 4 bold italic
                         arrow_label_offset = 0, # How much to offset the arrow label
                         arrow_label_alpha = 1, # The transparency of the arrow labels
                         
                         seed = NULL, # Set random seed
                         
                         plot_hulls = F, # Plot hulls
                         hull_alpha = 1, # The transparency of the hulls
                         
                         plot_spiders = F, # Plot spiders
                         label_spider = F, # Label the spider centroid
                         spider_label_size = 0.5, # Size of spider label
                         spider_alpha = 1, # The transparency of the spiders
                         
                         plot_ellipses = F, # Plot ellipse
                         label_ellipse = F, # Label ellipse
                         ellipse_label_size = 0.5, # Size of ellipse label
                         ellipse_border_width = 1, # Width of ellipse border
                         ellipse_alpha = 1, # The transparency of the ellipses
                         
                         
                         label_sites = F, # Label the individual sites/samples
                         label_species = F, # Label the individual species/OTUs
                         
                         point_alpha = 1, # Transparency of points, 0 = invisible, 1= fully, visible, 0.5 = 50% transparent
                         point_size = 0.8, # Size of points
                         point_line_thickness = 1, # Thickness of point outline
                         point_line_is_darker_fill = F, # Use a darkened version of the point fill colour as outline colour
                         use_shapes = T,
                         variable_shapes_available = F,
                         use_fill_shapes = T, # Whether to use fill shapes or hollow shapes
                         
                         include_legend = T, # Show the legend
                         legend_x = NULL, # x axis position of legend
                         legend_y = NULL, # y axis position of legend
                         legend_x_offset = 0, # How much to offset the legend on the x axis (assumes legend_x = NULL)
                         legend_y_offset = 0, # How much to offset the legend on the y axis (assumes legend_y = NULL)
                         legend_cex = 0.6, # Size of text in legend
                         legend_columns = 2, # Number of columns used for legend key
                         legend_key_text_distance = 0.5, # Distance between legend keys and text
                         legend_column_spacing = 0.5, # Distance between legend columns
                         legend_title = NULL, # Title of the legend. If NULL, just uses the variable
                         legend_fill_colour = NULL,
                         specie_labeller_function = NULL, # Function to re-label the species/OTUs,
                         
                         is_constrained = F, # This is a constrained ordination object
                         
                         file_type = "pdf",
                         filename = NULL,
                         plot_width = 10,
                         plot_height= 10
){
  # Function to generate a PCA plot
  if (!is.null(seed)){
    set.seed(seed)
  }
  pca.scores <- try(scores(pca_object, choices=c(1,2,3)))
  if(inherits(pca.scores, "try-error")) {
    return()
  }
  # Get component x,y coordinates
  pca_site_scores <- scores(pca_object, display = "sites", choices = component_choices)
  pca_specie_scores <- scores(pca_object, display = "species", choices = component_choices)
  
  # Check all entries in the PCA object are also in the metadata
  if (!all(rownames(pca_site_scores) %in% rownames(metadata.df))){
    error_message <- paste0("There are rows in the PCA object that are not defined in the metadata.
                             Rownames in the metadata need to match those used in the PCA object")
    stop(error_message)
  }
  # Remove NA entries from the metadata and from the PCA
  metadata.df <- metadata.df[!is.na(metadata.df[[variable_to_plot]]),]
  pca_site_scores <- pca_site_scores[rownames(pca_site_scores) %in% rownames(metadata.df),]
  
  if (is_constrained){
    pca_percentages <- (pca_object$CCA$eig/sum(pca_object$CCA$eig)) * 100  
  } else{
    pca_percentages <- (pca_object$CA$eig/sum(pca_object$CA$eig)) * 100  
  }
  
  # Assign a percentage of zero to NA values
  pca_percentages[is.na(pca_percentages)] <- 0
  if (length(pca_percentages) == 1){
    pca_percentages[2] <- 0
  }
  
  if (!is.null(axis_limits)){
    x_min <- axis_limits[1]
    x_max <- axis_limits[2]
    y_min <- axis_limits[3]
    y_max <- axis_limits[4]
  }else {
    x_min <- round(lapply(min(pca_site_scores[,1]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    x_max <- round(lapply(max(pca_site_scores[,1]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    y_min <- round(lapply(min(pca_site_scores[,2]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    y_max <- round(lapply(max(pca_site_scores[,2]), function(x) ifelse(x > 0, x + 1, x - 1))[[1]])
    
  }
  my_xlab <- ""
  my_ylab <- ""
  if (show_x_label){
    my_xlab = paste("PC",component_choices[1],"(", round(pca_percentages[component_choices[1]],1), "%)", sep = "")  
  }
  if (show_y_label){
    my_ylab = paste("PC", component_choices[2],"(", round(pca_percentages[component_choices[2]],1), "%)", sep = "")
  }
  
  metadata.df <- metadata.df[order(rownames(metadata.df)),]
  metadata.df <- metadata.df[order(metadata.df[[variable_to_plot]]),]
  
  # ------------------------------------------------------------------------------------
  # Ensure outcome variable is factored
  # Refactor the variable column so that the levels are consistent
  if (!is.null(my_levels)){
    metadata.df[,variable_to_plot] <- factor(metadata.df[,variable_to_plot], levels = my_levels)
  } else{
    # Uncomment to factorise and order levels alphabetically
    metadata.df[,variable_to_plot] <- factor(metadata.df[,variable_to_plot], 
                                             levels = sort(unique(as.character(metadata.df[,variable_to_plot]))))
    # Uncomment to factorise and inherit the current ordering
    # metadata.df[,variable_to_plot] <- factor(metadata.df[,variable_to_plot])
  }
  # ------------------------------------------------------------------------------------
  
  
  if (!is.null(filename)){
    if (file_type == "pdf"){
      pdf(file = filename,height = plot_height, width = plot_width)
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      svg(filename = filename,height = plot_height, width = plot_width)
      # svglite(file = filename,height = plot_height, width = plot_width)
    }
  }
  
  
  # dataframe.df <- dataframe.df[order(rownames(dataframe.df)),,drop = F]
  # dataframe.df <- dataframe.df[order(dataframe.df[[variable_to_plot]]),,drop = F]
  # dataframe.df[,variable_to_plot] <- factor(dataframe.df[,variable_to_plot], 
  #                                           levels = sort(unique(as.character(dataframe.df[,variable_to_plot]))))
  plot(0,
       type='n',
       # x = 0, y=0,
       xlim = c(x_min,x_max),
       ylim = c(y_min,y_max),
       xlab = my_xlab,
       ylab = my_ylab,
       xaxt = "n",
       yaxt = "n",
       # frame.plot = F,
       frame.plot = T
  )
  # Make grid
  if (hide_grid == F){
    grid(NULL,NULL, lty = 2, lwd = 1, col = "grey80")
  }
  if (highlight_origin == T){
    abline(v = 0, h = 0, col="grey40", lty = 2,lwd = 1)
  }
  
  
  # Add axes 
  axis(side = 1, labels = ifelse(plot_x_tick_labels, T, F), tck = -0.01,tick = ifelse(plot_x_ticks,T,F))
  axis(side = 2, labels = ifelse(plot_y_tick_labels, T, F), tck = -0.01, tick = ifelse(plot_x_ticks,T,F))
  
  # Assign (unique) colours and shapes for each grouping variable
  variable_values <- levels(metadata.df[[variable_to_plot]])
  
  # If variable colour column "variable_colour" in metadata, use colours from there
  if (variable_colours_available == T){
    colour_col_name <- paste0(variable_to_plot, "_colour")
    variable_colours <- setNames(as.character(unique(metadata.df[[colour_col_name]])), 
                                 as.character(unique(metadata.df[[variable_to_plot]])))
  } else if (!is.null(my_colour_palette)) {
    if (typeof(my_colour_palette) == "list" & all(variable_values %in% names(my_colour_palette)) & length(names(my_colour_palette)) == length(variable_values)){
      variable_colours <- my_colour_palette
    } else if (typeof(my_colour_palette) == "character" & length(my_colour_palette) == length(variable_values)) {
      colour_palette_distinct <- my_colour_palette
      variable_colours <- setNames(rep(colour_palette_distinct,length(variable_values))[1:length(variable_values)], variable_values)      
    } else{
      warning("Provided palette must be a list(variable_group_value = colour...) or character vector c(colour1, colour2...), with entries for each variable value.
              Using default palette.")
      colour_palette_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
      variable_colours <- setNames(rep(colour_palette_distinct,length(variable_values))[1:length(variable_values)], variable_values)
    }
  }
  else{ # Otherwise use default palette
    colour_palette_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
    variable_colours <- setNames(rep(colour_palette_distinct,length(variable_values))[1:length(variable_values)], variable_values)
  }
  
  # Whether to use shapes
  if (use_shapes == T){
    # If variable shape column "variable_shape" in metadata, use shapes from there
    if (variable_shapes_available == T){
      shape_col_name <- paste0(variable_to_plot, "_shape")
      variable_shapes <- setNames(unique(metadata.df[[shape_col_name]]), 
                                  as.character(unique(metadata.df[[variable_to_plot]])))
      
    } else{ # Use default pre-defined shapes
      # variable_shapes <- setNames(c(25,24,23,22,21,8,6,5,4,3,2,1)[1:length(variable_values)],variable_values)
      if (use_fill_shapes){
        variable_shapes <- setNames(rep(c(25,24,23,22,21),
                                        length(variable_values))[1:length(variable_values)],variable_values)
      }else{
        variable_shapes <- setNames(rep(c(6,5,2,1,0),
                                        length(variable_values))[1:length(variable_values)],variable_values)      
      }
    }
  } else{
    if (use_fill_shapes){
      variable_shapes <- setNames(rep(c(21),length(variable_values))[1:length(variable_values)],variable_values)    
    } else{
      variable_shapes <- setNames(rep(c(1),length(variable_values))[1:length(variable_values)],variable_values)    
    }
  }
  
  # variable_colours <- colours.l[[variable_to_plot]] # Assumes colours.l a global variable
  # variable_shapes <- pchs.l[[variable_to_plot]]  # Assumes pchs.l a global variable, also assumes shapes are only 25,24,23,22,21 (I think)
  # variable_shapes <- setNames(rep(c(25,24,23,22,21),length(variable_values))[1:length(variable_values)],variable_values)      
  
  # print(variable_values)
  # print(variable_colours)
  # print(variable_shapes)
  annotation_dataframe <- data.frame(variable_colours, variable_shapes, stringsAsFactors = F)
  annotation_dataframe$variable_outline_colours <- as.character(annotation_dataframe$variable_colours)
  if (point_line_thickness != 0){
    darken <- function(color, factor=1.4){
      col <- col2rgb(color)
      col <- col/factor
      col <- rgb(t(col), maxColorValue=255)
      col
    }
    
    lighten <- function(color, factor=1.4){
      col <- col2rgb(color)
      col <- col*factor
      col <- rgb(t(col), maxColorValue=255)
      col
    }
    if (point_line_is_darker_fill == T){
      annotation_dataframe[annotation_dataframe$variable_shapes > 15,"variable_outline_colours"] <- 
        as.character(lapply(annotation_dataframe$variable_outline_colours, darken))
    } else{
      annotation_dataframe[annotation_dataframe$variable_shapes > 15,"variable_outline_colours"] <- "black"    
    }
  }
  
  # Filter annotation dataframe to levels of interest
  if (!is.null(my_levels)){
    annotation_dataframe <- annotation_dataframe[my_levels,]
  }
  
  pca_site_scores <- pca_site_scores[rownames(metadata.df),]
  
  # Get the colours for all samples
  all_sample_colours <- as.character(
    lapply(
      as.character(metadata.df[rownames(pca_site_scores),variable_to_plot]), 
      function(x) as.character(annotation_dataframe[x,"variable_colours"])
    )
  )
  
  # Get the shapes for all samples
  all_sample_shapes <- as.numeric(
    lapply(
      as.character(metadata.df[rownames(pca_site_scores),variable_to_plot]), 
      function(x) annotation_dataframe[x,"variable_shapes"][[1]]
    )
  )
  
  # Set the outline colours for all samples based on the sample colours and refering to the annotation dataframe created above
  all_sample_outline_colours <- as.character(unlist(lapply(all_sample_colours, 
                                                           function(x) annotation_dataframe[annotation_dataframe$variable_colours == x, "variable_outline_colours"])))
  
  points(pca_site_scores, 
         cex = point_size,
         lwd = point_line_thickness,
         pch = all_sample_shapes,
         col = alpha(all_sample_outline_colours,point_alpha),
         # col = alpha("black",point_alpha),
         bg = alpha(all_sample_colours, point_alpha)
  )
  
  # Plot arrows for species / variables
  plot_arrows_func <- function(){
    
    left_pc1.v <- rownames(pca_specie_scores[order(pca_specie_scores[,1]),][1:num_top_species,])
    right_pc1.v <- rownames(pca_specie_scores[order(pca_specie_scores[,1]),][(length(pca_specie_scores[,1]) - num_top_species):length(pca_specie_scores[,1]),])
    
    left_pc2.v <- rownames(pca_specie_scores[order(pca_specie_scores[,2]),][1:num_top_species,])
    right_pc2.v <- rownames(pca_specie_scores[order(pca_specie_scores[,2]),][(length(pca_specie_scores[,2]) - num_top_species):length(pca_specie_scores[,2]),])
    
    top_vars.v <- unique(c(left_pc1.v, right_pc1.v, left_pc2.v, right_pc2.v))
    arrows(0,0, 
           arrow_scalar * pca_specie_scores[top_vars.v,1], 
           arrow_scalar * pca_specie_scores[top_vars.v,2], 
           length =0.05, 
           col = alpha(arrow_colour, arrow_alpha),
           lwd = arrow_thickness,
           lty = 1)
    
    if (label_arrows){
      if (!is.null(specie_labeller_function)){
        for (tv in top_vars.v){
          text(x = pca_specie_scores[tv,1]* arrow_scalar,
               y = pca_specie_scores[tv,2]* arrow_scalar,
               labels = specie_labeller_function(tv),
               cex = arrow_label_size,
               # Values of 1, 2, 3 and 4, respectively indicate positions below, 
               # to the left of, above and to the right of the specified (x,y) coordinates.
               pos = sample(c(1,2,3,4),1),
               # pos = 2,
               offset = arrow_label_offset,
               col = alpha(arrow_label_colour,alpha = arrow_label_alpha),
               font = arrow_label_font_type)
        }
      } else{
        for (tv in top_vars.v){
          text(x = pca_specie_scores[tv,1] * arrow_scalar,
               y = pca_specie_scores[tv,2] * arrow_scalar,
               labels = tv,
               cex = arrow_label_size,
               # Values of 1, 2, 3 and 4, respectively indicate positions below, 
               # to the left of, above and to the right of the specified (x,y) coordinates.
               pos = sample(c(1,2,3,4),1),
               # pos = 2,
               offset = arrow_label_offset,
               col = alpha(arrow_label_colour,alpha = arrow_label_alpha),
               font = arrow_label_font_type)
        }
      }
    }
  }
  
  if (plot_arrows == T){
    plot_arrows_func()
  }
  
  # Plot ellipses that are filled
  plot_ellipses_func <- function () {
    for (member in variable_values) {
      if (nrow(metadata.df[metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse
        ordiellipse(pca_site_scores,
                    groups = metadata.df[[variable_to_plot]],
                    kind = "ehull",
                    lwd = ellipse_border_width,
                    # border = variable_colours[member][[1]],
                    border = alpha(variable_colours[member][[1]],ellipse_alpha),
                    # col = variable_colours[member][[1]],
                    col = alpha(variable_colours[member][[1]],ellipse_alpha),
                    show.groups = member,
                    alpha = .05,
                    draw = "polygon",
                    label = F,
                    cex = .5)
      }
    }
  }
  
  # Plot hulls that are filled
  plot_hulls_func <- function () {
    for (member in variable_values){
      if (nrow(metadata.df[metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse}
        ordihull(pca_site_scores,
                 groups = metadata.df[[variable_to_plot]],
                 lwd = ellipse_border_width,
                 # border = variable_colours[member][[1]],
                 border = alpha(variable_colours[member][[1]],hull_alpha),
                 # col = variable_colours[member][[1]],
                 col = alpha(variable_colours[member][[1]],hull_alpha),
                 show.groups = member,
                 alpha = .05,
                 draw = "polygon",
                 label = F,
                 cex = .5)
      }
    }
  }
  if (hasArg(plot_hulls)){
    if (plot_hulls == T){
      plot_hulls_func()    
    }
  }
  
  #Plot spiders
  plot_spiders_func <- function (label_spider = F) {
    for (member in variable_values){
      if (nrow(metadata.df[metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse
        ordispider(pca_site_scores,
                   groups = metadata.df[[variable_to_plot]],
                   # col = variable_colours[member][[1]],
                   col = alpha(variable_colours[member][[1]],spider_alpha),
                   show.groups = member,
                   #alpha = .05,
                   label = label_spider,
                   cex = spider_label_size)
      }
    }
  }
  if (hasArg(plot_spiders)){
    if (plot_spiders == T){
      plot_spiders_func(label_spider = label_spider)    
    }
  }
  
  plot_ellipses_labels_func <- function(label_ellipse = F){
    # Repeat to have labels clearly on top of all ellipses
    for (member in variable_values){
      if (nrow(metadata.df[metadata.df[[variable_to_plot]] == member,,drop = F]) > 2){ # if too few samples, skip plotting ellipse
        ordiellipse(pca_site_scores,
                    groups = metadata.df[[variable_to_plot]],
                    kind = "ehull",
                    # border = variable_colours[member][[1]],
                    border = NA,
                    # col = variable_colours[member][[1]],
                    col = NA,
                    show.groups = member,
                    alpha = 0,
                    draw = "polygon",
                    label = label_ellipse,
                    cex = ellipse_label_size)
      }
    }
  }
  
  if (hasArg(plot_ellipses) | hasArg(label_ellipse)){
    if (plot_ellipses == T){
      plot_ellipses_func()    
      plot_ellipses_labels_func(label_ellipse = label_ellipse)
    } else if (label_ellipse == T){
      plot_ellipses_labels_func(label_ellipse = label_ellipse)
    }
  } 
  
  if (label_sites == T){
    text(x = pca_site_scores[,1],
         y = pca_site_scores[,2],
         labels = rownames(pca_site_scores),
         cex = .5,
         pos = 2)
  }
  if (label_species == T){
    text(x = pca_specie_scores[,1],
         y = pca_specie_scores[,2],
         labels = rownames(pca_specie_scores),
         cex = .5,
         pos = 2)
  }
  if (is.null(legend_x)){
    legend_x <- x_min + legend_x_offset
  }
  if (is.null(legend_y)){
    legend_y <- y_max + legend_y_offset
  }
  if (is.null(legend_title)){
    legend_title <- variable_to_plot
  }
  
  if (!is.null(plot_title)){
    title(main = plot_title, cex.main = title_cex)
  } 
  
  legend_bty = "n"
  if (!is.null(legend_fill_colour)){
    legend_bty <- "o"
    legend_fill <- legend_fill_colour
    legend_colour <- legend_fill_colour
  }
  
  if (include_legend){
    legend(
      # title = bold(variable_to_plot),
      title = as.expression(bquote(bold(.(legend_title)))),
      # title.adj = 0.5,
      title.col="black",
      # x = x_min-4,
      # y = y_max-6,
      x = legend_x,
      y = legend_y,
      # legend= variable_values,
      # pch= unique(all_sample_shapes),
      # col= legend_point_outline_colours,
      # pt.bg = unique(all_sample_colours),
      legend= rownames(annotation_dataframe),
      pch= annotation_dataframe$variable_shapes,
      col= as.character(annotation_dataframe$variable_outline_colours),
      pt.bg = as.character(annotation_dataframe$variable_colours),
      #bg = "white",
      bty = legend_bty,
      bg = legend_fill,
      box.col = legend_colour,
      ncol = legend_columns,
      cex = legend_cex,
      # pt.cex = 0.6,
      pt.lwd = point_line_thickness,
      y.intersp =1,
      x.intersp = legend_key_text_distance, # Distance between legend keys and text
      xjust = 0,
      # yjust = 1,
      title.adj = 0.5,
      text.width = legend_column_spacing # Distance between legend columns
    )
  }
  if (!is.null(filename)){
    dev.off()
  }
}


# Function to create heatmap
make_heatmap <- function(heatmap.m,
                         metadata.df = NULL, # rownames must match the columns of the heatmap.m
                         filename= NULL, # Location to save heatmap figure
                         
                         # Dataframe with two columns. First must match row entry, second the new label
                         row_labels.df = NULL, 
                         col_labels.df = NULL, # Same as row_labels, though with columns
                         height = 10, # Not currently used
                         width = 10, # Not currently used
                         heatmap_height = 10, # Not currently used
                         heatmap_width = 10, # Not currently used
                         plot_height =10,
                         plot_width =10,
                         
                         # e.g. unit(c(2, 2, 2, 20), "mm")) 
                         # bottom, left, top, right paddings
                         my_padding = NULL, # Padding around plot
                         
                         column_title_size = 10,
                         row_title_size = 10,
                         annotation_bar_name_size = 10,
                         annotation_variables = NULL, # Annotations
                         annotation_palette = NULL,
                         simple_anno_size = unit(.5, "cm"), # size of annotations
                         col_annotation_label_size = 6,
                         col_annotation_title_size = 6,
                         col_annotation_legend_grid_height = .2,
                         col_annotation_legend_grid_width = .2,
                         plot_title = NULL,
                         
                         cluster_columns = T,
                         cluster_rows = T,
                         my_breaks = NULL,
                         heatmap_palette = NULL,
                         default_palette_choice = NULL, # red, blue, purple, bluered, dark_bluered
                         legend_title = NULL,
                         legend_labels = NULL,
                         scale_legend_label_size = 6,
                         scale_legend_title_size = 6,
                         discrete_legend = FALSE, # Whether or not to display continuous legend as discrete
                         
                         column_title = "Sample",
                         row_title = "Taxa",
                         
                         show_column_dend = F,
                         show_row_dend = F,
                         do_not_order = F, # Do not order the heatmap rows by the row labels names
                         show_cell_values = F,
                         # If show_cell_values =T, cells less than this will have a black font colour
                         # and above white
                         cell_fun_value_col_threshold = 15,
                         my_cell_fun = NULL, # function to apply to the cell values
                         show_legend = T,
                         show_top_annotation = T,
                         row_name_size = 6,
                         col_name_size = 6,
                         
                         grid_thickness = 1,
                         grid_colour = "white",
                         
                         draw_plot = T, # Draw the plot
                         ...){
  # print(list(...))
  argList<-list(...) # argument list for checking unspecified optional parameters
  # print(argList$cell_fun)
  # return(1)
  
  colour_palette_206 <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff",
                          "#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d",
                          "#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c",
                          "#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b",
                          "#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff",
                          "#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4",
                          "#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00",
                          "#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8",
                          "#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100",
                          "#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8",
                          "#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f",
                          "#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff",
                          "#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3",
                          "#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff",
                          "#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614",
                          "#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec",
                          "#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd",
                          "#93f6e7","#01b4a4")
  
  # Assign internal objects
  internal_heatmap_matrix.m <- heatmap.m
  internal_metadata.df <- metadata.df
  
  ha <- NULL
  
  if (!is.null(internal_metadata.df)){
    # Order/filter the heatmap matrix to order/entries of metadata
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,rownames(internal_metadata.df),drop = F]
    if (!is.null(annotation_variables)){
      # Order the heatmap matrix by the annotation_variables
      internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,do.call(order, internal_metadata.df[,annotation_variables,drop=F]),drop =F]
      # Order the metadata by the annotation_variables
      internal_metadata.df <- internal_metadata.df[do.call(order, internal_metadata.df[,annotation_variables,drop=F]),,drop=F]
      # Create metadata just containing the annotation_variables
      metadata_just_variables <- internal_metadata.df[,annotation_variables, drop = F]    
    }
    # Check that rownames match colnames
    if (!all(rownames(internal_metadata.df) == colnames(internal_heatmap_matrix.m))){
      stop("Row names in metadata do not match column names in matrix")
    }
    
    # Create annotations
    colour_lists <- list()
    for (myvar in annotation_variables){
      var_colour_name <- paste0(myvar, "_colour")
      # Assumes there is a colour column for each variable in the metadata
      # If there is no colour column, create one and assign from palette
      # internal_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
      # internal_colour_palette_10_distinct <- my_colour_palette_20_distinct
      if (is.null(annotation_palette)){
        annotation_palette <- colour_palette_206
      }
      if (!var_colour_name %in% names(internal_metadata.df)){
        myvar_values <- factor(as.character(sort(unique(internal_metadata.df[,myvar]))))
        myvar_colours <- setNames(annotation_palette[1:length(myvar_values)], myvar_values)
        all_variable_colours <- as.character(lapply(as.character(internal_metadata.df[,myvar]), function(x) myvar_colours[x]))
        internal_metadata.df[,paste0(myvar,"_colour")] <- all_variable_colours
      }
      
      metadata_subset <- unique(internal_metadata.df[,c(myvar, var_colour_name)])
      # Order by the variable column
      metadata_subset <- metadata_subset[order(metadata_subset[,myvar]),]
      # Factorise the variable column
      metadata_subset[,myvar] <- factor(metadata_subset[,myvar])
      metadata_subset <- metadata_subset[!is.na(metadata_subset[,myvar]),]
      named_colour_list <- setNames(as.character(metadata_subset[, var_colour_name]), as.character(metadata_subset[,myvar]))
      colour_lists[[myvar]] <- named_colour_list
    }
    
    # Appearance of the column annotations
    #HeatmapAnnotation
    if (show_top_annotation == T){
      if(is.null(annotation_variables)){
        print("No annotation variables specified, cannot add annotations")
        show_top_annotation <- F
      } else{
        # TODO make flexible for different types of annotations, e.g. bar
        # Or allow user to provide annotations that are added in, 'provided_ha'
        ha <- columnAnnotation(df = metadata_just_variables,
                               
                               # which = "column",
                               col = colour_lists,
                               gp = gpar(col = "black",lwd =.2),
                               gap = unit(.1,"cm"),
                               show_annotation_name = T,
                               # annotation_legend_param, # ?color_mapping_legend for options
                               show_legend = show_legend,
                               simple_anno_size = simple_anno_size,
                               annotation_legend_param = list(labels_gp = gpar(fontsize = col_annotation_label_size),
                                                              title_gp = gpar(fontsize = col_annotation_title_size,fontface = "bold"),
                                                              grid_height = unit(col_annotation_legend_grid_height,"cm"),
                                                              grid_width = unit(col_annotation_legend_grid_width,"cm")),
                               annotation_name_gp = gpar(fontsize = annotation_bar_name_size))      
      }
    }
  }

  # TODO - add option for row annotation
  if (is.null(heatmap_palette)){
    if (is.null(default_palette_choice)) {default_palette_choice <- "blue"}
    if (!default_palette_choice %in% c("blue", "purple","red","dark_bluered","bluered")) { default_palette_choice <- "blue"}
    if (default_palette_choice == "blue"){
      heatmap_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))
    } 
    else if (default_palette_choice == "purple"){
      heatmap_palette <- colorRampPalette(c("white", "#f9cdac","#f3aca2", "#ee8b97", "#e96a8d","#db5087","#b8428c", "#973490", "#742796","#5e1f88", "#4d1a70", "#3d1459","#2d0f41"))
    } else if (default_palette_choice == "red"){
      heatmap_palette <- colorRampPalette(c("white", "#fded86","#fde86e", "#f9d063", "#f5b857","#f0a04b","#eb8a40", "#e77235","#e35b2c", "#c74e29","#9d4429","#753c2c","#4c3430"))
    } else if (default_palette_choice == "dark_bluered"){
      heatmap_palette <- colorRampPalette(c("#08306B","#FFD92F","#67001F"))
    } else if (default_palette_choice == "bluered"){
      heatmap_palette <- colorRampPalette(c("#17468a","#ffdd47","#99113a"))
    }
  } else{
    heatmap_palette <- colorRampPalette(heatmap_palette)
  }
  
  if (!is.null(my_breaks)){
    internal_breaks <- my_breaks
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = heatmap_palette(length(internal_breaks)))
  } else{
    internal_breaks <- seq(min(internal_heatmap_matrix.m), max(internal_heatmap_matrix.m), length.out = 6)
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = heatmap_palette(length(internal_breaks)))
  }
  
  row_labels.v = rownames(internal_heatmap_matrix.m)
  if (!is.null(row_labels.df)){
    row_labels.v <- as.character(lapply(row_labels.v, function(x) as.character(row_labels.df[row_labels.df[,1] == x,][,2])))
  }
  col_labels.v = colnames(internal_heatmap_matrix.m)
  if (!is.null(col_labels.df)){
    col_labels.v <- as.character(lapply(col_labels.v, function(x) as.character(col_labels.df[col_labels.df[,1] == x,][,2])))
  }
  
  if (do_not_order != T){
    # Order the heatmap rows by the row labels names
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[order(row_labels.v),]
    row_labels.v <- row_labels.v[order(row_labels.v)]    
  }
  
  # if show values and no function provided
  if (show_cell_values == T & is.null(my_cell_fun)){ 
    my_cell_fun <- function(j, i, x, y, width, height, fill) {
      # if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold & internal_heatmap_matrix.m[i, j] != 0){
      if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold){
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "black"))}
      else if(internal_heatmap_matrix.m[i, j] >= cell_fun_value_col_threshold ) {
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "white"))
      }
    }
  }
  
  # Legend appearance
  if (is.null(legend_labels)){
    my_labels <- internal_breaks
  } else{
    my_labels <- legend_labels
  }
  if (discrete_legend == TRUE){
    hm_legend <- Legend(
      labels = rev(my_labels),
      at = internal_breaks,
      labels_gp = gpar(fontsize = scale_legend_label_size),
      legend_gp = gpar(fill = rev(col_fun(internal_breaks))), # For discrete
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = scale_legend_title_size,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  } else{
    hm_legend <- Legend(
      col_fun = col_fun, # For continuous
      labels = my_labels,
      at = internal_breaks,
      labels_gp = gpar(fontsize = scale_legend_label_size),
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = scale_legend_title_size,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  }
  
  hm <- Heatmap(matrix = internal_heatmap_matrix.m,

                top_annotation = ha,
                
                # Colours
                col = col_fun,
                na_col = "grey",
                
                # Sizing
                show_heatmap_legend = F,
                heatmap_legend_param = list(hm_legend),
                row_names_max_width = unit(35,"cm"),
                row_labels = row_labels.v,
                column_labels = col_labels.v,
                # row_names_side = "left",
                # height = unit(height,"cm"),
                # width = unit(width,"cm"),
                # heatmap_height = unit(heatmap_height,"cm"),
                # heatmap_width = unit(heatmap_width,"cm"),
                # heatmap_width = unit(15,"cm"),
                
                # Titles
                column_title = column_title,
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize = column_title_size),
                row_title = row_title,
                row_title_side = "left",
                row_title_gp = gpar(fontsize = row_title_size),
                
                # Clustering
                cluster_columns = cluster_columns,
                cluster_rows = cluster_rows,
                clustering_method_columns = "average",
                clustering_method_rows = "average",
                show_column_dend = show_column_dend, 
                show_row_dend = show_row_dend,
                # column_dend_height = unit(2, "cm"),
                # row_dend_width = unit(3, "cm"),
                
                # Borders
                border = F,
                rect_gp = gpar(col = grid_colour, lwd = grid_thickness),
                
                # Text appearance
                row_names_gp = gpar(fontsize = row_name_size),
                column_names_gp = gpar(fontsize = col_name_size),
                cell_fun = my_cell_fun,
                ...
  )
  
  # if (show_top_annotation == T & exists("ha")){
  #   hm <- ha %v% hm
  # }
  padding <- unit(c(2, 2, 2, 2), "mm")
  if (!is.null(my_padding)){
    padding <- my_padding 
  }
  
  if (!is.null(plot_title)){
    column_title <- plot_title
  } else{
    column_title <- character(0)
  }
  if (!is.null(filename)){
    pdf(filename,height=plot_height,width=plot_width)
    draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T,padding = padding, column_title = column_title,
         column_title_gp = gpar(hjust = 0))
    dev.off()
  }
  if (draw_plot == T){
    draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T,padding = padding,column_title = column_title)
  }
  return(list("heatmap" = hm, "legend" = hm_legend))
  
}


gm_mean <- function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

# Center log ratio transform
clr <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}


read_counts_and_unique_features <- function(count_matrix, sample_list){
  # Number of features and read counts
  temp <- count_matrix
  temp[temp > 0] <- 1
  read_sums <- colSums(count_matrix)
  feature_sums <- colSums(temp)
  sample_read_counts <- data.frame(sample_list, unlist(lapply(sample_list, 
                                                              function(x) ifelse(x %in% names(read_sums), read_sums[x],0))))
  sample_feature_counts <- data.frame(sample_list, unlist(lapply(sample_list, 
                                                                 function(x) ifelse(x %in% names(feature_sums), feature_sums[x],0))))
  out <- list("sample_read_counts" = sample_read_counts,
              "sample_feature_counts" = sample_feature_counts)
  out
}



run_permanova_custom <- function(my_metadata, my_formula, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  result <- adonis(my_formula,data = my_metadata, permu=permutations,method= my_method)
  # result <- adonis(my_formula,data = my_metadata, permu=999,method="bray")
  for (r in rownames(result$aov.tab)){
    variable <- r
    Degress_of_freedom <- result$aov.tab[r,]$Df[1]
    SumOfSqs <- round(result$aov.tab[r,]$SumsOfSqs[1], 3)
    meanSqs <- round(result$aov.tab[r,]$MeanSqs[1], 3)
    F.model <- round(result$aov.tab[r,]$F.Model[1], 3)
    R2 <- round(result$aov.tab[r,]$R2[1], 3)
    p_value <- round(result$aov.tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degress_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       R2,
                                                       p_value))
  }
  # my_formula_string <- paste0(as.character(my_formula)[2], as.character(my_formula)[1], as.character(my_formula)[3])
  my_formula_string <- paste0(as.character(my_formula)[1], as.character(my_formula)[3])
  print(paste0("FORMULA: ", my_formula_string))
  print(result)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Formula <- my_formula_string
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}



run_permdisp_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  dist_matrix <- vegdist(t(my_data), method = my_method)
  betadisper_object <- with(my_metadata, betadisper(dist_matrix, group = get(my_group)))
  permutest_results <- permutest(betadisper_object, permutations = permutations, parallel = 2)
  
  for (r in rownames(permutest_results$tab)){
    variable <- r
    Degrees_of_freedom <- permutest_results$tab[r,]$Df[1]
    SumOfSqs <- round(permutest_results$tab[r,]$`Sum Sq`[1],3)
    meanSqs <- round(permutest_results$tab[r,]$`Mean Sq`[1], 3)
    F.model <- round(permutest_results$tab[r,]$F[1], 3)
    N_permutations <- permutest_results$tab[r,]$N.Perm[1]
    p_value <- round(permutest_results$tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degrees_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       N_permutations,
                                                       p_value))
  }
  print(permutest_results)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","Permutations","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Group <- my_group
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}

run_anosim_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  anosim_object <- with(my_metadata, anosim(x = t(my_data), 
                                            grouping = get(my_group),
                                            permutations = permutations,
                                            distance = my_method,
                                            parallel = 2))
  
  stat_sig_table <- rbind(stat_sig_table, data.frame(my_group,
                                                     anosim_object$statistic,
                                                     anosim_object$signif,
                                                     anosim_object$permutations))
  names(stat_sig_table) <- c("Variable","R_statistic", "Significance","Permutations")
  stat_sig_table$Method <- my_method
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}






calculate_correlation_matrix <- function(mydata.m, method = "pearson", adjust = "BH"){
  
  # Remove row entries that don't vary across all samples
  internal_data.m <- mydata.m
  zv <- apply(internal_data.m, 1, function(x) length(unique(x)) == 1)
  internal_data.m <- internal_data.m[!zv, ]
  
  # Take a two vectors and perform a significance/correlation test
  calculate_stats <- function(x, y, dist.name) {
    k <- cor.test(x, y, method=dist.name)
    c(k$estimate, k$stat, k$p.value)
  }
  
  # cor_result <- corr.test(t(internal_data.m), adjust = adjust, method = method)
  r <- apply(t(internal_data.m), 2,function(col) 
    t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,1,drop =F])
  p <- apply(t(internal_data.m), 2,function(col) 
    t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,3,drop =F])
  padj <- apply(t(internal_data.m), 2,function(col) 
    p.adjust(t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,3,drop =F],
             method = adjust))
  
  rownames(r) <- colnames(r)
  rownames(p) <- colnames(p)
  rownames(padj) <- colnames(padj)
  r <- signif(r, digits = 5)
  p <- signif(p, digits = 5)
  padj <- signif(padj, digits = 5)
  
  # r	The matrix of correlations
  # n	Number of cases per correlation
  # t	value of t-test for each correlation
  # p	two tailed probability of t for each correlation. For symmetric matrices, p values adjusted for multiple tests are reported above the diagonal.
  # se	standard error of the correlation
  # ci	the alpha/2 lower and upper values, as well as the (Holm or Bonferroni) adjusted confidence intervals.
  # list(cor_matrix = cor_result$r, cor_pval_matrix = cor_result$p)
  list(cor_matrix = r, cor_pval_matrix = p, cor_padj_matrix = padj)
}

# Generate a dot correlation plot generated by corrplot
plot_corrplot <- function(correlation_matrix, 
                          p_value_matrix = NULL, 
                          plot_title = "", 
                          plot_title_size = .6,
                          plot_height = 10, plot_width = 10,
                          p_value_threshold = 0.05, 
                          label_size = 1,
                          relabeller_function = NULL, 
                          insig = "blank", insig_pch = 4, insig_pch_cex = 1, insig_pch_col = "black",
                          make_insig_na = F,
                          to_exclude = NULL, 
                          method = "circle", 
                          outline = T,
                          label_colour = "black", 
                          colour_label_size = 1,
                          grid_colour = "black",
                          pairs_to_na = NULL, # Should be two column dataframe (row column / column row)
                          order = "hclust", 
                          col = NULL, 
                          filename = NULL,
                          file_type = "pdf",
                          ...
                          ){
  
  cor.m <- correlation_matrix
  cor_pval.m <- p_value_matrix
  
  # Entries to remove completely
  if (!is.null(to_exclude)){ 
    cor.m <- cor.m[!rownames(cor.m) %in% to_exclude, !colnames(cor.m) %in% to_exclude]
    if (!is.null(p_value_matrix)){
      cor_pval.m <- cor_pval.m[!rownames(cor_pval.m) %in% to_exclude, !colnames(cor_pval.m) %in% to_exclude]  
    }
    
  }
  # Entries to make NA. Should be a two column dataframe
  if (!is.null(pairs_to_na)){
    if (order == "hclust"){
      # stop("Cannot use hclust ordering with NA values in matrix")
      # So before inserting NA values, first order the matrix
      # print(summary(hclust(dist(cor.m),method = "average")))
      ord <- hclust(dist(cor.m),method = "average")$order
      cor.m <- cor.m[ord,ord]
      if (!is.null(p_value_matrix)){
        cor_pval.m <- cor_pval.m[ord,ord]
      }
      order <- "original"
    }
    for (row in 1:nrow(pairs_to_na)){
      a <- as.character(pairs_to_na[row,1])
      b <- as.character(pairs_to_na[row,2])
      if (a %in% rownames(cor.m) & b %in% rownames(cor.m)){
        cor.m[a, b] <- NA
        cor.m[b, a] <- NA
        if (!is.null(p_value_matrix)){
          cor_pval.m[a, b] <- NA
          cor_pval.m[b, a] <- NA
        }
      }
    }
  }
  
  if (make_insig_na == T){
    cor.m[which(cor_pval.m >= p_value_threshold)] <- NA
    cor_pval.m[cor_pval.m >= p_value_threshold] <- NA
  }
  
  if (!is.null(relabeller_function)){
    colnames(cor.m) <- unlist(lapply(colnames(cor.m), relabeller_function))
    rownames(cor.m) <- unlist(lapply(colnames(cor.m), relabeller_function))
    if (!is.null(p_value_matrix)){
      colnames(cor_pval.m) <- unlist(lapply(colnames(cor_pval.m), relabeller_function))
      rownames(cor_pval.m) <- unlist(lapply(rownames(cor_pval.m), relabeller_function))
      # colnames(cor_pval.m) <- relabeller_function(colnames(cor_pval.m))
      # rownames(cor_pval.m) <- relabeller_function(rownames(cor_pval.m))
    }
  }
  if (is.null(col)){
    col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                  "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061")))(200)
  }
  
  if (!is.null(filename)){
    # pdf(file = filename,height = plot_height, width = plot_width)
    if (file_type == "pdf"){
      pdf(file = filename, height = plot_height, width = plot_width)
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      # svg(filename = filename,height = plot_height, width = plot_width)
      svglite(file = filename,height = plot_height, width = plot_width)
    }
  }
  
  corrplot(corr = cor.m,
           method = method,
           outline = outline,
           tl.col = label_colour,
           tl.cex = label_size,
           addgrid.col = grid_colour,
           # tl.srt = 45,
           # title = plot_title,
           col = col,
           # col = brewer.pal(n = 8, name = "RdYlBu"),
           type = "lower",
           diag = F,
           na.label = "square",
           na.label.col = "grey",
           order = order,
           hclust.method = "average",
           p.mat = cor_pval.m,
           sig.level = p_value_threshold,
           # insig = "blank",
           insig = insig,
           pch = insig_pch,
           pch.cex = insig_pch_cex,
           pch.col = insig_pch_col,
           cl.pos = 'r',
           cl.cex = colour_label_size,
           mar=c(1,0,3,1),
           ...)
  title(main = plot_title,cex.main = plot_title_size)
  if (!is.null(filename)){
    dev.off()
  }
}
