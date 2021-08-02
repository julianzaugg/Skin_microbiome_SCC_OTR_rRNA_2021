# Diversity calculations Alpha (Shannon, Chao1, Simpson) for each sample.

library(vegan)
library(phyloseq)
library(ggnewscale)
library(tidyverse)
library(FSA)

source("code/utility.R")


summarise_alpha_diversities <- function(mydata, group_by_columns){
  summary.df <- mydata %>% 
    dplyr::group_by_(.dots = c(group_by_columns)) %>%
    dplyr::summarise(Shannon_Mean =mean(Shannon),
                     Shannon_Stdev=sd(Shannon),
                     Shannon_Max=max(Shannon), 
                     Shannon_Min=min(Shannon), 
                     Shannon_Median=median(Shannon), 
                     
                     Simpson_Mean=mean(Simpson), 
                     Simpson_Stdev=sd(Simpson),
                     Simpson_Max=max(Simpson), 
                     Simpson_Min=min(Simpson), 
                     Simpson_Median=median(Simpson), 
                     
                     Chao1_Mean=mean(Chao1), 
                     Chao1_Stdev=sd(Chao1),
                     Chao1_Max=max(Chao1), 
                     Chao1_Min=min(Chao1), 
                     Chao1_Median=median(Chao1),
                     
                     N_patients=n_distinct(Subject),
                     N_samples = n_distinct(Index)
    ) %>% as.data.frame()
  summary.df
}


summarise_diversities_each_variable <- function(mydata, variables){
  mydata_summary.df <- NULL
  for (myvar in variables){
    if (is.null(mydata_summary.df)){
      mydata_summary.df <- summarise_alpha_diversities(mydata, myvar)
      # Melt because we are combining variables into a single table
      mydata_summary.df <- melt(mydata_summary.df, measure.vars = myvar, value.name = "Group", variable.name = "Variable")
    } else{
      mydata_summary.df <- rbind(mydata_summary.df, melt(summarise_alpha_diversities(mydata, myvar),
                                                         measure.vars = myvar, value.name = "Group", variable.name = "Variable"))      
    }
    
  }
  mydata_summary.df
}

calculate_alpha_diversity_significance_multiple <- function(mydata.df, variable){
  # Assumes there are Shannon, Chao1 and Simpson columns
  n_groups = length(as.character(unique(mydata.df[,variable])))
  if (any(is.na(mydata.df[,variable]))){
    return()
  }
  if (n_groups > 2){
    kw_shannon <- kruskal.test(Shannon~get(variable), data = mydata.df)
    kw_simpson <- kruskal.test(Simpson~get(variable), data = mydata.df)
    kw_chao1 <- kruskal.test(Shannon~get(variable), data = mydata.df)
    dunn_shannon <- dunnTest(x = Shannon~get(variable), data = mydata.df, method = "bh", alpha = 0.05)
    dunn_simpson <- dunnTest(x = Simpson~get(variable), data = mydata.df, method = "bh", alpha = 0.05)
    dunn_chao1 <- dunnTest(x = Chao1~get(variable), data = mydata.df, method = "bh", alpha = 0.05)
    
    dunn_shannon <- separate(dunn_shannon$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
    names(dunn_shannon) <- c("Group_1","Group_2","Shannon_Dunn_pvalue","Shannon_Dunn_padj")
    
    dunn_simpson <- separate(dunn_simpson$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
    names(dunn_simpson) <- c("Group_1","Group_2","Simpson_Dunn_pvalue","Simpson_Dunn_padj")
    
    dunn_chao1 <- separate(dunn_chao1$res, Comparison, into = c("Group_1", "Group_2"), sep = " - ")[,c("Group_1","Group_2","P.unadj","P.adj")]
    names(dunn_chao1) <- c("Group_1","Group_2","Chao1_Dunn_pvalue","Chao1_Dunn_padj")
    
    multiple_group_comparison.df <- merge(merge(x = dunn_shannon, y = dunn_simpson, by = c("Group_1", "Group_2")), y = dunn_chao1,  by = c("Group_1", "Group_2"))
    multiple_group_comparison.df$Shannon_KrusW_pvalue <- kw_shannon$p.value
    multiple_group_comparison.df$Simpson_KrusW_pvalue <- kw_simpson$p.value
    multiple_group_comparison.df$Chao1_KrusW_pvalue <- kw_chao1$p.value
    multiple_group_comparison.df$Variable = variable
    return(multiple_group_comparison.df)
  }
  multiple_group_comparison.df <- data.frame("Variable" = character(),
                                             "Group_1" = character(),
                                             "Group_2" = character(),
                                             "Shannon_Dunn_pvalue" = character(),
                                             "Shannon_Dunn_padj" = character(),
                                             "Simpson_Dunn_pvalue" = character(),
                                             "Simpson_Dunn_padj" = character(),
                                             "Chao1_Dunn_pvalue" = character(),
                                             "Chao1_Dunn_padj" = character(),
                                             "Shannon_KrusW_pvalue" = character(),
                                             "Simpson_KrusW_pvalue" = character(),
                                             "Chao1_KrusW_pvalue" = character())
  return(multiple_group_comparison.df)
  
}

generate_p_labels <- function(sig_table){
  for (sig_column in c("Chao1_Dunn_padj", "Shannon_Dunn_padj", "Simpson_Dunn_padj")){
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
  axis.text = element_text(size = 9, colour = "black"),
  axis.title = element_text(size = 10,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))



# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Load the processed metadata
metadata.df <- read.csv("Result_tables/other/processed_metadata.csv", sep =",", header = T)

# Set the Index to be the rowname
rownames(metadata.df) <- metadata.df$Index

# Define the discrete variables
discrete_variables <- c("Sample_type", "Cohort", "Subject")

# Factorise
metadata.df$Sample_type <- factor(metadata.df$Sample_type, levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))

# Load the counts
genus.m <-  as.matrix(read.csv("Result_tables/count_tables/Genus_counts.csv", header =T, row.names = 1))

# Order the matrices and metadata to be the same order
genus.m <- genus.m[,rownames(metadata.df)]

# Create the rarefied matrices.
# First check the counts for each sample and use that to decide the threshold
sort(colSums(genus.m))


# Rarefaction curve
# rarecurve(t(genus.m[,colSums(genus.m) > 1]),step = 500, label = F,xlim = c(0,30000))

# Below we set the threshold and any samples with less than this amount are discarded
rarefy_threshold <- 5000

# How many rows/columns (taxa/samples) before we filter?
nrow(genus.m)
ncol(genus.m)

genus_rare.m <- t(rrarefy(t(genus.m[,colSums(genus.m) >= rarefy_threshold]), rarefy_threshold))

# How many rows/columns (taxa/samples) remaining?
nrow(genus_rare.m)
ncol(genus_rare.m)

# Filter metadata to match
metadata.df <- metadata.df[colnames(genus_rare.m),]

# Now that we have rarefied, create phyloseq objects with the rarefied matrices
genus_rare_phyloseq <- otu_table(genus_rare.m, taxa_are_rows=TRUE)

# Estimate alpha diversities
both_cohorts_genus_rare_alpha.df <- estimate_richness(genus_rare_phyloseq, measures = c("Chao1", "Simpson","Shannon"))
# both_cohorts_genus_rare_alpha.df <- estimate_richness(genus_rare_phyloseq, measures = c("Observed", "Chao1", "Simpson","Shannon"))
# both_cohorts_genus_rare_alpha.df <- both_cohorts_genus_rare_alpha.df[rownames(metadata.df),]

# Combine with metadata
both_cohorts_genus_rare_alpha.df$Sample <- rownames(both_cohorts_genus_rare_alpha.df)
both_cohorts_genus_rare_alpha.df <- left_join(both_cohorts_genus_rare_alpha.df, metadata.df, by = c("Sample" = "Index"))


# Create Cohort specific datasets
immunosuppressed_genus_rare_alpha.df <- subset(both_cohorts_genus_rare_alpha.df, Cohort == "organ transplant recipient")
immunocompetent_genus_rare_alpha.df <- subset(both_cohorts_genus_rare_alpha.df, Cohort == "immunocompetent")


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ANOVA tests on diversity

anova_genus_chao1 <- with(both_cohorts_genus_rare_alpha.df, aov(Chao1~Sample_type))
anova_genus_shannon <- with(both_cohorts_genus_rare_alpha.df, aov(Shannon~Sample_type))
anova_genus_simpson <- with(both_cohorts_genus_rare_alpha.df, aov(Simpson~Sample_type))
anova_immunosuppressed_genus_chao1 <- with(immunosuppressed_genus_rare_alpha.df, aov(Chao1~Sample_type))
anova_immunosuppressed_genus_shannon <- with(immunosuppressed_genus_rare_alpha.df, aov(Shannon~Sample_type))
anova_immunosuppressed_genus_simpson <- with(immunosuppressed_genus_rare_alpha.df, aov(Simpson~Sample_type))
anova_immunocompetent_genus_chao1 <- with(immunocompetent_genus_rare_alpha.df, aov(Chao1~Sample_type))
anova_immunocompetent_genus_shannon <- with(immunocompetent_genus_rare_alpha.df, aov(Shannon~Sample_type))
anova_immunocompetent_genus_simpson <- with(immunocompetent_genus_rare_alpha.df, aov(Simpson~Sample_type))

hist(residuals(anova_genus_chao1))
hist(residuals(anova_genus_shannon))
hist(residuals(anova_genus_simpson))

anova_genus_chao1 %>% summary()
anova_genus_shannon %>% summary()
anova_genus_simpson %>% summary()
anova_immunosuppressed_genus_chao1 %>% summary()
anova_immunosuppressed_genus_shannon %>% summary()
anova_immunosuppressed_genus_simpson %>% summary()
anova_immunocompetent_genus_chao1 %>% summary()
anova_immunocompetent_genus_shannon %>% summary()
anova_immunocompetent_genus_simpson %>% summary()


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Kruskal-wallis and Dunn-tests


immunosuppressed_genus_alpha_dunn_significances.df <- calculate_alpha_diversity_significance_multiple(immunosuppressed_genus_rare_alpha.df,variable = "Sample_type")
immunosuppressed_genus_alpha_dunn_significances.df <- immunosuppressed_genus_alpha_dunn_significances.df[with(immunosuppressed_genus_alpha_dunn_significances.df, which(Simpson_KrusW_pvalue <= 0.05 | Shannon_KrusW_pvalue <= 0.05 | Chao1_KrusW_pvalue <= 0.05)),]
immunosuppressed_genus_alpha_dunn_significances.df <- immunosuppressed_genus_alpha_dunn_significances.df[with(immunosuppressed_genus_alpha_dunn_significances.df, which(Shannon_Dunn_padj <= 0.05 | Simpson_Dunn_padj <= 0.05 | Chao1_Dunn_padj <= 0.05)),]
write.csv(temp,"Result_tables/dunn_tests/immunosuppressed_dunn_diversity_genus_sample_type.csv",quote = F, row.names = F)

immunocompetent_genus_alpha_dunn_significances.df <- calculate_alpha_diversity_significance_multiple(immunocompetent_genus_rare_alpha.df,variable = "Sample_type")
immunocompetent_genus_alpha_dunn_significances.df <- immunocompetent_genus_alpha_dunn_significances.df[with(immunocompetent_genus_alpha_dunn_significances.df, which(Simpson_KrusW_pvalue <= 0.05 | Shannon_KrusW_pvalue <= 0.05 | Chao1_KrusW_pvalue <= 0.05)),]
immunocompetent_genus_alpha_dunn_significances.df <- immunocompetent_genus_alpha_dunn_significances.df[with(immunocompetent_genus_alpha_dunn_significances.df, which(Shannon_Dunn_padj <= 0.05 | Simpson_Dunn_padj <= 0.05 | Chao1_Dunn_padj <= 0.05)),]
write.csv(temp,"Result_tables/dunn_tests/immunocompetent_dunn_diversity_genus_sample_type.csv",quote = F, row.names = F)

immunosuppressed_genus_alpha_diversity_summary.df <- summarise_diversities_each_variable(immunosuppressed_genus_rare_alpha.df, variables = discrete_variables)
immunocompetent_genus_alpha_diversity_summary.df <- summarise_diversities_each_variable(immunocompetent_genus_rare_alpha.df, variables = discrete_variables)


immunosuppressed_group_count.df <- immunosuppressed_genus_rare_alpha.df %>% 
  filter(!is.na(Shannon)) %>%
  group_by(Sample_type) %>% 
  dplyr::summarise(N_Group_1 = n())
immunosuppressed_group_count.l <- setNames(immunosuppressed_group_count.df$N_Group_1,immunosuppressed_group_count.df$Sample_type)

immunocompetent_group_count.df <- immunocompetent_genus_rare_alpha.df %>% 
  filter(!is.na(Shannon)) %>%
  group_by(Sample_type) %>% 
  dplyr::summarise(N_Group_1 = n())
immunocompetent_group_count.l <- setNames(immunocompetent_group_count.df$N_Group_1,immunocompetent_group_count.df$Sample_type)

# Add group sizes
immunosuppressed_genus_alpha_dunn_significances.df$N_Group_1 <- 
  unlist(lapply(immunosuppressed_genus_alpha_dunn_significances.df$Group_1, function(x)immunosuppressed_group_count.l[[x]]))
immunosuppressed_genus_alpha_dunn_significances.df$N_Group_2 <- 
  unlist(lapply(immunosuppressed_genus_alpha_dunn_significances.df$Group_2, function(x)immunosuppressed_group_count.l[[x]]))

immunocompetent_genus_alpha_dunn_significances.df$N_Group_1 <- 
  unlist(lapply(immunocompetent_genus_alpha_dunn_significances.df$Group_1, function(x)immunocompetent_group_count.l[[x]]))
immunocompetent_genus_alpha_dunn_significances.df$N_Group_2 <- 
  unlist(lapply(immunocompetent_genus_alpha_dunn_significances.df$Group_2, function(x)immunocompetent_group_count.l[[x]]))

# Add rarefied depth
immunosuppressed_genus_alpha_dunn_significances.df$Rarefied_depth <- rarefy_threshold
immunocompetent_genus_alpha_dunn_significances.df$Rarefied_depth <- rarefy_threshold

# Add p-value labels
immunosuppressed_genus_alpha_dunn_significances.df <- generate_p_labels(immunosuppressed_genus_alpha_dunn_significances.df)
immunocompetent_genus_alpha_dunn_significances.df <- generate_p_labels(immunocompetent_genus_alpha_dunn_significances.df)

# Add Cohort information
immunosuppressed_genus_alpha_dunn_significances.df$Cohort  <- "organ transplant recipient"
immunocompetent_genus_alpha_dunn_significances.df$Cohort  <- "immunocompetent"

# Write cohort significance values
write.csv(x = rbind(immunosuppressed_genus_alpha_dunn_significances.df,immunocompetent_genus_alpha_dunn_significances.df),
          file = paste0("Result_tables/diversity/IS_IC_diversity_signficances_",rarefy_threshold,".csv"),
          row.names =F,
          quote = F)

# ----------------------------------------
# Plotting
sample_type_colours <- unique(metadata.df[,c("Sample_type", "Sample_type_colour")])
sample_type_colours <- setNames(as.character(sample_type_colours[,"Sample_type_colour"]), sample_type_colours[,"Sample_type"])
sample_type_shapes <- unique(metadata.df[,c("Sample_type", "Sample_type_shape")])
sample_type_shapes <- setNames(as.numeric(sample_type_shapes[,"Sample_type_shape"]), sample_type_shapes[,"Sample_type"])

sample_type_outline_colours <- unique(metadata.df[,c("Sample_type", "Sample_type_colour")])
sample_type_outline_colours <- setNames(darken(as.character(sample_type_outline_colours[,"Sample_type_colour"])), sample_type_outline_colours[,"Sample_type"])

immunosuppressed_genus_shannon_boxplot <-
  ggplot(immunosuppressed_genus_rare_alpha.df,aes(x = Sample_type, 
                                                  y = Shannon, 
                                                  fill = Sample_type, 
                                                  # colour = Sample_type,
                                                  shape = Sample_type)) +
  # geom_violin() +
  geom_boxplot(position = position_dodge(width =.75), outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .75,
                                              dodge.width = .75))  +
  scale_shape_manual(values = sample_type_shapes,
                     name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, 
                    name = "Sample type") +
  scale_colour_manual(values = sample_type_outline_colours,
                    name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  xlab("Sample type") +
  ylab("Shannon") +
  common_theme +
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5, .5))
immunosuppressed_genus_shannon_boxplot


immunosuppressed_genus_chao1_boxplot <- 
  ggplot(immunosuppressed_genus_rare_alpha.df,aes(x = Sample_type, 
                                                  y = Chao1, 
                                                  fill = Sample_type, 
                                                  shape = Sample_type)) +
  # geom_violin() +
  geom_boxplot(position = position_dodge(width =.75), outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .75,
                                              dodge.width = .75))  +
  scale_shape_manual(values = sample_type_shapes,
                     name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, 
                    name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  xlab("Sample type") +
  ylab("Chao1") +
  common_theme +
  scale_y_continuous(limits = c(0,300), breaks = seq(0,350, 50))
immunosuppressed_genus_chao1_boxplot

immunosuppressed_genus_simpson_boxplot <- 
  ggplot(immunosuppressed_genus_rare_alpha.df,aes(x = Sample_type, 
                                                  y = Simpson, 
                                                  fill = Sample_type, 
                                                  shape = Sample_type)) +
  # geom_violin() +
  geom_boxplot(position = position_dodge(width =.75), outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .75,
                                              dodge.width = .75))  +
  scale_shape_manual(values = sample_type_shapes,
                     name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, 
                    name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  xlab("Sample type") +
  ylab("Simpson") +
  common_theme +
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1, .1))
immunosuppressed_genus_simpson_boxplot




immunocompetent_genus_shannon_boxplot <- 
  ggplot(immunocompetent_genus_rare_alpha.df,aes(x = Sample_type, 
                                                 y = Shannon, 
                                                 fill = Sample_type, 
                                                 shape = Sample_type)) +
  # geom_violin() +
  geom_boxplot(position = position_dodge(width =.75), outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .75,
                                              dodge.width = .75))  +
  scale_shape_manual(values = sample_type_shapes,
                     name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, 
                    name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  xlab("Sample type") +
  ylab("Shannon") +
  common_theme +
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5, 0.5))

immunocompetent_genus_chao1_boxplot <- 
  ggplot(immunocompetent_genus_rare_alpha.df,aes(x = Sample_type, 
                                                 y = Chao1, 
                                                 fill = Sample_type, 
                                                 shape = Sample_type)) +
  # geom_violin() +
  geom_boxplot(position = position_dodge(width =.75), outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .75,
                                              dodge.width = .75))  +
  scale_shape_manual(values = sample_type_shapes,
                     name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, 
                    name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  xlab("Sample type") +
  ylab("Chao1") +
  common_theme +
  scale_y_continuous(limits = c(0,300), breaks = seq(0,350, 50))
immunocompetent_genus_chao1_boxplot

immunocompetent_genus_simpson_boxplot <- 
  ggplot(immunocompetent_genus_rare_alpha.df,aes(x = Sample_type, 
                                                 y = Simpson, 
                                                 fill = Sample_type, 
                                                 shape = Sample_type)) +
  # geom_violin() +
  geom_boxplot(position = position_dodge(width =.75), outlier.shape = NA, width=.5,lwd =.3) +
  geom_jitter(size = .6,stroke =.1,
              position = position_jitterdodge(jitter.width = .75,
                                              dodge.width = .75))  +
  scale_shape_manual(values = sample_type_shapes,
                     name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, 
                    name = "Sample type") +
  stat_summary(fun = "mean", colour = "grey2", geom = "point",
               shape = 16,size = 1, position = position_dodge(width = .75),show.legend = F) +
  xlab("Sample type") +
  ylab("Simpson") +
  common_theme +
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1, .1))
immunocompetent_genus_simpson_boxplot

# Extract the legend
my_legend <- cowplot::get_legend(immunosuppressed_genus_chao1_boxplot + 
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


library(cowplot)
chao1 <- plot_grid(plotlist = list(immunosuppressed_genus_chao1_boxplot + 
                                     # labs(title = "Organ transplant recipient") + 
                                     theme(axis.title.x = element_blank(),
                                           axis.text.x = element_blank(),
                                           legend.position = "none",
                                           plot.title = element_text(face = "bold", hjust = 0.5,size = 10))
                                   , 
                                   immunocompetent_genus_chao1_boxplot + 
                                     # labs(title = "Immunocompetent") + 
                                     theme(axis.title.x = element_blank(),
                                           axis.title.y = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.text.y = element_blank(),
                                           legend.position = "none",
                                           plot.title = element_text(face = "bold", hjust = 0.5,size = 10))
),
rel_widths = c(1,.8))

shannon <- plot_grid(plotlist = list(immunosuppressed_genus_shannon_boxplot + 
                                       theme(
                                         # axis.title.x = element_blank(),
                                         # axis.title.y = element_blank(),
                                         # axis.text.x = element_blank(),
                                         # axis.text.y = element_blank(),
                                         legend.position = "none",
                                         plot.margin = unit(c(5.5,5.5,5.5,6),"pt")
                                       ), 
                                     immunocompetent_genus_shannon_boxplot + 
                                       theme(
                                         # axis.title.x = element_blank(),
                                         axis.title.y = element_blank(),
                                         # axis.text.x = element_blank(),
                                         axis.text.y = element_blank(),
                                         legend.position = "none")
),
rel_widths = c(1,.8))
simpson <- plot_grid(plotlist = list(immunosuppressed_genus_simpson_boxplot + 
                                       theme(legend.position = "none",
                                             plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
                                             axis.title.x = element_blank(),
                                             # axis.title.y = element_blank(),
                                             axis.text.x = element_blank(),
                                             # axis.text.y = element_blank(),
                                       )+
                                       labs(title = "Organ transplant recipient"),
                                     immunocompetent_genus_simpson_boxplot + 
                                       theme(
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_blank(),
                                         axis.text.x = element_blank(),
                                         axis.text.y = element_blank(),
                                         legend.position = "none",
                                         plot.title = element_text(face = "bold", hjust = 0.5,size = 10)
                                       ) +
                                       labs(title = "Immunocompetent")
                                     
                                     
),
rel_widths = c(1,.8))

grid_plot <- plot_grid(plotlist = list(simpson,chao1,shannon),
                       ncol = 1,nrow = 3,
                       rel_heights = c(1,1,1), 
                       align = "hv")

# grid_plot <- plot_grid(plotlist = list(grid_plot, my_legend), 
# nrow = 2, rel_heights = c(1,.1))


# grid_plot
ggsave(filename = "Result_figures/diversity/IS_IC_diversity_boxplots.pdf", 
       plot = grid_plot, width = 13, 
       height = 20, units = "cm", device = "pdf")

ggsave(filename = "Result_figures/diversity/IS_IC_diversity_boxplots.svg", 
       plot = grid_plot, width = 13, 
       height = 20, units = "cm", device = "svg")


