library(vegan)
library(tidyverse)
library(reshape2)


source("code/utility.R")


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

# Load metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", header = T)

# Set rownames to match the sample ID
rownames(metadata.df) <- metadata.df$Index

# Order the metadata.df by the sample value
# metadata.df <- metadata.df[order(metadata.df$Sample),]

# Load count data
asv.m <- df2m(read.csv("Result_tables/count_tables/ASV_counts.csv", header = T))
genus.m <- df2m(read.csv("Result_tables/count_tables/Genus_counts.csv", header = T))

# CLR transform the count matrices
asv_clr.m <- clr(asv.m)
genus_clr.m <- clr(genus.m)


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# PERMANOVA tests whether distance differ between groups.

# Permutational Multivariate Analysis of Variance (PERMANOVA) can be used to 
# determine if the structure of the microbial communities is significantly different between
# environmental variables. This is done using the adonis function from Vegan with a distance metric, e.g. Bray-Curtis, and a
# specified number of permutations, e.g. 9,999.
# The analysis measures the degree each environmental variable affects the community composition and indicates 
# the significance of that effect on beta diversity (described by p-values and R2 values). 
# The R2 value corresponds to the proportion of variability observed in the dissimilarity.

# Genus, clr euclidean
print("Centred-log ratio transformed counts - Euclidean distance")

asv_permanova_results <- data.frame()
genus_permanova_results <- data.frame()

asv_interaction_subject_permanova_results <- data.frame()
genus_interaction_subject_permanova_results <- data.frame()

asv_within_cohort_permanova_results <- data.frame()
genus_within_cohort_permanova_results <- data.frame()

asv_interaction_subject_within_cohort_permanova_results <- data.frame()
genus_interaction_subject_within_cohort_permanova_results <- data.frame()

# run_permanova_custom(my_metadata = metadata.df,
#                      my_formula = as.formula(paste0("t(genus_clr_subset.m)~Subject+Cohort+Sample_type+Subject:Sample_type + Cohort:Sample_type")),
#                      my_method = "euclidean",label = "CLR",permutations = 999)
discrete_variables <- c("Sample_type", "Cohort", "Subject")
# Compare groups for each variable
for (myvar in discrete_variables){
  print(myvar)
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  
  asv_clr_subset.m <- asv_clr.m[,rownames(metadata_subset.df)]
  genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
  number_of_samples <- nrow(metadata_subset.df)
  
  # asv
  temp <- run_permanova_custom(my_metadata = metadata_subset.df,
                               my_formula = as.formula(paste0("t(asv_clr_subset.m)~", myvar)),
                               my_method = "euclidean",label = "CLR",permutations = 999)
  temp$Number_of_samples <- number_of_samples
  asv_permanova_results <- rbind(asv_permanova_results,temp)
  
  # Genus
  temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                               my_formula = as.formula(paste0("t(genus_clr_subset.m)~", myvar)),
                               my_method = "euclidean",label = "CLR",permutations = 999)
  temp$Number_of_samples <- number_of_samples
  genus_permanova_results <- rbind(genus_permanova_results, temp)
  
  # Include an interaction with Subject in the model
  if (myvar != "Subject"){
    temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                                 my_formula = as.formula(paste0("t(asv_clr_subset.m)~Subject+", myvar, "+Subject:",myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 999)
    temp$Number_of_samples <- number_of_samples
    asv_interaction_subject_permanova_results <- rbind(asv_interaction_subject_permanova_results, temp)
    
    temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                                 my_formula = as.formula(paste0("t(genus_clr_subset.m)~Subject+", myvar, "+Subject:",myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 999)
    temp$Number_of_samples <- number_of_samples
    genus_interaction_subject_permanova_results <- rbind(genus_interaction_subject_permanova_results, temp)
  }
  # Also compare the groups after accounting for the cohort
  if (myvar == "Cohort") {next}
  for (cohort in unique(metadata.df$Cohort)){
    print(paste0("Processing cohort: ", cohort))
    metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
    metadata_subset.df <- subset(metadata_subset.df, Cohort == cohort)
    if (length(unique(metadata_subset.df[,myvar])) < 2){
      next
    }
    
    asv_clr_subset.m <- asv_clr.m[,rownames(metadata_subset.df)]
    genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
    
    number_of_samples <- nrow(metadata_subset.df)
    
    # Include an interaction with Subject in the model
    if (myvar != "Subject"){
      temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                                   my_formula = as.formula(paste0("t(asv_clr_subset.m)~Subject+", myvar, "+Subject:",myvar)),
                                   my_method = "euclidean",label = "CLR",permutations = 999)
      temp$Cohort <- cohort
      temp$Number_of_samples <- number_of_samples
      asv_interaction_subject_within_cohort_permanova_results <- rbind(asv_interaction_subject_within_cohort_permanova_results, temp)
      
      temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                                   my_formula = as.formula(paste0("t(genus_clr_subset.m)~Subject+", myvar, "+Subject:",myvar)),
                                   my_method = "euclidean",label = "CLR",permutations = 999)
      temp$Cohort <- cohort
      temp$Number_of_samples <- number_of_samples
      genus_interaction_subject_within_cohort_permanova_results <- rbind(genus_interaction_subject_within_cohort_permanova_results, temp)
    }
    
    temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                                 my_formula = as.formula(paste0("t(asv_clr_subset.m)~", myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 999)
    temp$Cohort <- cohort
    temp$Number_of_samples <- number_of_samples
    asv_within_cohort_permanova_results <- rbind(asv_within_cohort_permanova_results, temp)
    
    temp <- run_permanova_custom(my_metadata = metadata_subset.df, 
                                 my_formula = as.formula(paste0("t(genus_clr_subset.m)~", myvar)),
                                 my_method = "euclidean",label = "CLR",permutations = 999)
    temp$Cohort <- cohort
    temp$Number_of_samples <- number_of_samples
    genus_within_cohort_permanova_results <- rbind(genus_within_cohort_permanova_results,temp)
    
  }
}
# genus_within_cohort_permanova_results

write.csv(asv_permanova_results, file = "Result_tables/permanova/asv_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_permanova_results, file = "Result_tables/permanova/genus_PERMANOVA.csv", row.names = F, quote = F)

write.csv(asv_within_cohort_permanova_results, file = "Result_tables/permanova/asv_within_cohort_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_within_cohort_permanova_results, file = "Result_tables/permanova/genus_within_cohort_PERMANOVA.csv", row.names = F, quote = F)

write.csv(asv_interaction_subject_permanova_results, file = "Result_tables/permanova/asv_interaction_subject_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_interaction_subject_permanova_results, file = "Result_tables/permanova/genus_interaction_subject_PERMANOVA.csv", row.names = F, quote = F)

write.csv(asv_interaction_subject_within_cohort_permanova_results, file = "Result_tables/permanova/asv_interaction_subject_within_cohort_PERMANOVA.csv", row.names = F, quote = F)
write.csv(genus_interaction_subject_within_cohort_permanova_results, file = "Result_tables/permanova/genus_interaction_subject_within_cohort_PERMANOVA.csv", row.names = F, quote = F)

# ---------------------------------------------
# PERMDISP (betadisper)

asv_permdisp_results <- data.frame()
genus_permdisp_results <- data.frame()

asv_within_cohort_permdisp_results <- data.frame()
genus_within_cohort_permdisp_results <- data.frame()


for (myvar in discrete_variables){
  print(myvar)
  metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
  
  asv_clr_subset.m <- asv_clr.m[,rownames(metadata_subset.df)]
  genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
  
  number_of_samples <- nrow(metadata_subset.df)
  
  temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                              my_data = asv_clr_subset.m,
                              my_group = myvar,
                              my_method = "euclidean",
                              permutations = 999,
                              label = "CLR")
  temp$Number_of_samples <- number_of_samples
  asv_permdisp_results <- rbind(asv_permdisp_results, temp)
  
  temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                              my_data = genus_clr_subset.m,
                              my_group = myvar,
                              my_method = "euclidean",
                              permutations = 999,
                              label = "CLR")
  temp$Number_of_samples <- number_of_samples
  genus_permdisp_results <- rbind(genus_permdisp_results, temp)
  
  if (myvar == "Cohort") {next}
  for (cohort in unique(metadata.df$Cohort)){
    print(paste0("Processing cohort: ", cohort))
    metadata_subset.df <- metadata.df[!is.na(metadata.df[,myvar]),]
    metadata_subset.df <- subset(metadata_subset.df, Cohort == cohort)
    if (length(unique(metadata_subset.df[,myvar])) < 2){
      next
    }
    number_of_samples <- nrow(metadata_subset.df)
    
    asv_clr_subset.m <- asv_clr.m[,rownames(metadata_subset.df)]
    genus_clr_subset.m <- genus_clr.m[,rownames(metadata_subset.df)]
    
    temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                                my_data = asv_clr_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 999,
                                label = "CLR")
    temp$Cohort <- cohort
    temp$Number_of_samples <- number_of_samples
    asv_within_cohort_permdisp_results <- rbind(asv_within_cohort_permdisp_results, temp)
    
    temp <- run_permdisp_custom(my_metadata = metadata_subset.df, 
                                my_data = genus_clr_subset.m,
                                my_group = myvar,
                                my_method = "euclidean",
                                permutations = 999,
                                label = "CLR")
    temp$Cohort <- cohort
    temp$Number_of_samples <- number_of_samples
    genus_within_cohort_permdisp_results <- rbind(genus_within_cohort_permdisp_results, temp)
  }
}

write.csv(asv_permdisp_results, file = "Result_tables/permdisp/asv_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_permdisp_results, file = "Result_tables/permdisp/genus_PERMDISP.csv", row.names = F, quote = F)

write.csv(asv_within_cohort_permdisp_results, file = "Result_tables/permdisp/asv_within_cohort_PERMDISP.csv", row.names = F, quote = F)
write.csv(genus_within_cohort_permdisp_results, file = "Result_tables/permdisp/genus_within_cohort_PERMDISP.csv", row.names = F, quote = F)
