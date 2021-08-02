library(DESeq2)
# library(BiocParallel)




filter_matrix_rows <- function(my_matrix, row_max){
  rows_before <- dim(my_matrix)[1]
  filtered_matrix <- my_matrix[apply(my_matrix,1,max) >= row_max,]
  rows_after <- dim(filtered_matrix)[1]
  print(paste0("Rows before = ", rows_before))
  print(paste0("Rows after = ", rows_after))
  print(paste0("Lost % = ", round((rows_before-rows_after)/rows_before*100, 2), "%"))
  return(filtered_matrix)
}



# Load utility functions
source("code/utility.R")


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Load the ASV - taxonomy mapping file
asv_taxonomy_map.df <- read.csv("Result_tables/asv_taxonomy_map.csv", header = T)

# Load the processed metadata
metadata.df <- read.csv("Result_tables/processed_metadata.csv", sep =",", header = T)
rownames(metadata.df) <- metadata.df$Index

# Order the metadata.df by the index value
metadata.df <- metadata.df[order(metadata.df$Index),]

# Define the variables of interest
discrete_variables <- c("Sample_type","Subject", "Cohort", "Length_of_immunosuppression_group_1", "Length_of_immunosuppression_group_2")

# Load count data
asv.m <- df2m(read.csv("Result_tables/count_tables/ASV_counts.csv", header = T))
genus.m <- df2m(read.csv("Result_tables/count_tables/Genus_counts.csv", header = T))

# Order the matrices and metadata to be the same order
asv.m <- asv.m[,rownames(metadata.df)]
genus.m <- genus.m[,rownames(metadata.df)]

# Filter out features/taxa that do not have at # reads in at least one sample
head(melt(sort(colSums(asv.m))))
asv.m <- filter_matrix_rows(asv.m,30)
genus.m <- filter_matrix_rows(genus.m,30)
# dim(otu.m[apply(otu.m, 1, max) == 0,])
head(melt(sort(colSums(asv.m))))

# Ensure names of the otu / genus count matrices match the order of the metadata.df!
# Assumes number of samples in metadata.df and count data are the same
all(colnames(asv.m) == metadata.df$Index) # Should be 'True'
all(colnames(asv.m) == rownames(metadata.df)) # Should be 'True'

# Convert variables to factors
metadata.df[discrete_variables] <- lapply(metadata.df[discrete_variables], factor)

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Perform differential abundance calculations at the ASV level and genus level, 
# comparing between the groups within variables of interest

# DESeq requires a count matrix ('countData'), a corresponding metadata.df ('colData') and a 'design' formula. The formula expresses
# how the counts for each ASV/genus depend on the variables defined in the 'colData'. See help(DESeqDataSetFromMatrix) for more information.
# The first column of the metadata.df ('colData') must match the ordering of the columns of the countData

# metadata.df %>% group_by(Subject, Sample_type) %>% dplyr::summarise(Count = n()) %>% as.data.frame()


compare_groups_deseq <- function(mydata.m, mymetadata.df, myvariables.v, asv_taxonomy_map.df = NULL, assign_taxonomy = T){
  
  # For each rowname (ASV_ID), get the corresponding taxonomy_species
  # Assumes "ASV_ID" and "taxonomy_species" columns in the provided mapping dataframe
  assign_taxonomy_to_asv <- function(ASV_table.df, taxon_map.df){
    taxonomies <- c()
    for (asv_id in rownames(ASV_table.df)){
      taxonomies <- c(taxonomies, as.character(taxon_map.df[taxon_map.df$ASV_ID == asv_id,]$taxonomy_species))
    }
    return(taxonomies)
  }
  
  # Filter and sort DESeq result tables
  filter_and_sort_dds_results <- function(x, p_value_threshold = 0.05){
    filtered_table <- x
    filtered_table <- filtered_table[!is.na(filtered_table$padj),]
    filtered_table <- filtered_table[filtered_table$padj <= p_value_threshold,]
    filtered_table <- filtered_table[order(filtered_table$padj),]
    return(filtered_table)
  }
  
  # Compare groups for all variables
  combined_results_ordered.df <- data.frame()
  
  for (myvar in myvariables.v){
    print(paste0("Processing ", myvar))
    
    # Get all non-NA entries in the metadata
    mymetadata_filtered.df <- mymetadata.df[!is.na(mymetadata.df[,myvar]),]
    
    # Ensure factored variable
    mymetadata_filtered.df[,myvar] <- factor(mymetadata_filtered.df[,myvar])
    
    # Extract corresponding entries from data
    mydata_filtered.m <- mydata.m[,rownames(mymetadata_filtered.df)]
    
    # If the number of samples is 1 or there is only one unique variable
    if (dim(mymetadata_filtered.df)[2] == 1 | length(unique(mymetadata_filtered.df[,myvar])) == 1){
      print("Only one sample or only one unique group")
      break
    }
    if (dim(mymetadata_filtered.df)[1] == 0 | dim(mydata_filtered.m)[2] == 0){
      print("No samples after filtering")
      break
    }
    
    # If the column and rownames do not match, entries are missing
    if (!all(rownames(mymetadata_filtered.df) == colnames(mydata_filtered.m))){
      print("Colnames and metadata names don't match!!!")
      break
    }
    
    # Run DESeq
    dds <- DESeqDataSetFromMatrix(countData = mydata_filtered.m, colData = mymetadata_filtered.df, design = as.formula(paste0("~", myvar)))
    
    geoMeans <- apply(counts(dds), 1, gm_mean)
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    dds <- try(DESeq(dds, test = "Wald", fitType = "parametric", parallel = F))
    if(inherits(dds, "try-error")) {
      next
    }
    group_combinations <- combn(sort(unique(mymetadata_filtered.df[,myvar])),2)
    
    for (i in 1:ncol(group_combinations)){
      group_1 <- as.character(group_combinations[1,i])
      group_2 <- as.character(group_combinations[2,i])
      
      # group_1_meta <- subset(full, get(myvar) == group_1)
      # group_2_meta <- subset(full, get(myvar) == group_2)
      # n_group_1 <- dim(group_1_meta)[1]
      # n_group_2 <- dim(group_2_meta)[1]
      
      n_group_1 <- dim(subset(mymetadata_filtered.df, get(myvar) == group_1))[1]
      n_group_2 <- dim(subset(mymetadata_filtered.df, get(myvar) == group_2))[1]
      
      # Extract results for contrasted groups
      print(paste0(myvar, ": ", group_1, " vs ", group_2))
      resMFSource <- results(dds, contrast = c(myvar,group_1,group_2), alpha=0.05, independentFiltering = F, cooksCutoff = F,parallel = T)
      # print(resMFSource)
      resMFSource$Group_1 <- group_1
      resMFSource$Group_2 <- group_2
      resMFSource$Variable <- myvar
      resMFSource$N_Group_1 <- n_group_1
      resMFSource$N_Group_2 <- n_group_2
      
      # Assign the taxonomy to the results. Assumes feature.
      if (assign_taxonomy == T){
        if (! is.null(asv_taxonomy_map.df)){
          resMFSource$Taxonomy <- assign_taxonomy_to_asv(resMFSource, asv_taxonomy_map.df)     
        } else{
          resMFSource$Taxonomy <- ""
        }
        
        # Convert to dataframe
        resMFSource <- m2df(resMFSource, "ASV_ID")
      } else{
        # Convert to dataframe
        resMFSource <- m2df(resMFSource, "Taxonomy")
      }
      # print(resMFSource)
      resMFSource <- filter_and_sort_dds_results(resMFSource, 0.05)
      combined_results_ordered.df <- rbind(combined_results_ordered.df, resMFSource)
    }
  }
  combined_results_ordered.df
}

compare_groups_deseq_within_group <- function(mydata.m, mymetadata.df, myvariables.v, within_group_variable, asv_taxonomy_map.df = NULL, assign_taxonomy = F){
  combined_results.df <- data.frame()
  reduced_variables <- myvariables.v[which(!myvariables.v == within_group_variable)]
  for (myvar_value in unique(metadata.df[,within_group_variable])){
    print(paste0("Processing ", myvar_value))
    temp <- compare_groups_deseq(mydata.m = mydata.m, 
                                 mymetadata.df = subset(mymetadata.df, get(within_group_variable) == myvar_value), 
                                 myvariables.v = reduced_variables, 
                                 asv_taxonomy_map.df,
                                 assign_taxonomy = assign_taxonomy)
    if (dim(temp)[1] == 0){
      next
    }
    temp[,within_group_variable] <- myvar_value
    combined_results.df <- rbind(combined_results.df, temp)
  }
  combined_results.df
}


# May be commented out to avoid re-running (very slow)

# Compare all lesion types within each cohort
asv_group_comparison_within_cohort.df <- compare_groups_deseq_within_group(mydata.m = asv.m,
                                                                           mymetadata.df = metadata.df,
                                                                           myvariables = c("Sample_type"),
                                                                           within_group_variable = "Cohort",
                                                                           asv_taxonomy_map.df = asv_taxonomy_map.df,
                                                                           assign_taxonomy = T)
write.csv(x =asv_group_comparison_within_cohort.df,file ="Result_tables/DESeq_results/ASV_within_cohort_deseq.csv",quote = F, row.names =F)

genus_group_comparison_within_cohort.df <- compare_groups_deseq_within_group(mydata.m = genus.m,
                                                                             mymetadata.df = metadata.df,
                                                                             myvariables = c("Sample_type"),
                                                                             within_group_variable = "Cohort",
                                                                             asv_taxonomy_map.df = asv_taxonomy_map.df,
                                                                             assign_taxonomy = F)
write.csv(x =genus_group_comparison_within_cohort.df,file ="Result_tables/DESeq_results/Genus_within_cohort_deseq.csv",quote = F, row.names =F)

# Compare each length of immunosuppresion group within each sample type
genus_lois_comparison_within_sample_type.df <- compare_groups_deseq_within_group(mydata.m = asv.m,
                                                                                 mymetadata.df = metadata.df,
                                                                                 myvariables = c("Length_of_immunosuppression_group_1","Length_of_immunosuppression_group_2"),
                                                                                 within_group_variable = "Sample_type",
                                                                                 asv_taxonomy_map.df = asv_taxonomy_map.df,
                                                                                 assign_taxonomy = T)
write.csv(x =genus_lois_comparison_within_sample_type.df,file ="Result_tables/DESeq_results/ASV_length_of_immunosuppresion_within_sample_type_deseq.csv",quote = F, row.names =F)


genus_lois_comparison_within_sample_type.df <- compare_groups_deseq_within_group(mydata.m = genus.m,
                                                                             mymetadata.df = metadata.df,
                                                                             myvariables = c("Length_of_immunosuppression_group_1","Length_of_immunosuppression_group_2"),
                                                                             within_group_variable = "Sample_type",
                                                                             asv_taxonomy_map.df = asv_taxonomy_map.df,
                                                                             assign_taxonomy = F)
write.csv(x =genus_lois_comparison_within_sample_type.df,file ="Result_tables/DESeq_results/Genus_length_of_immunosuppresion_within_sample_type_deseq.csv",quote = F, row.names =F)


# Comparing the same lesion types between cohorts. Always compare suppressed vs competent, e.g. suppressed AK vs competent AK
# The trick is to group by the lesion type and then only compare groups within the Cohort variable
# asv_cohort_comparison_within_sample_type.df <- compare_groups_deseq_within_group(mydata.m = asv.m,
#                                                                             mymetadata.df = metadata.df,
#                                                                             myvariables = c("Cohort"),
#                                                                             within_group_variable = "Sample_type",
#                                                                             assign_taxonomy = T)
# write.csv(x =asv_cohort_comparison_within_sample_type.df,file ="Result_tables/DESeq_results/ASV_cohort_within_sample_type_deseq.csv",quote = F, row.names =F)

# genus_cohort_comparison_within_sample_type.df <- compare_groups_deseq_within_group(mydata.m = genus.m,
#                                                                               mymetadata.df = metadata.df,
#                                                                               myvariables = c("Cohort"),
#                                                                               within_group_variable = "Sample_type",
#                                                                               assign_taxonomy = F)
# write.csv(x =genus_cohort_comparison_within_sample_type.df,file ="Result_tables/DESeq_results/Genus_cohort_within_sample_type_deseq.csv",quote = F, row.names =F)


# Compare groups
# asv_group_comparison.df <- compare_groups_deseq(mydata.m = asv.m, mymetadata.df = metadata.df, myvariables = c("Sample_type"), assign_taxonomy = T)
# write.csv(x =asv_group_comparison.df,file ="Result_tables/DESeq_results/ASV_deseq.csv",quote = F, row.names =F)

# genus_group_comparison.df <- compare_groups_deseq(mydata.m = genus.m, mymetadata.df = metadata.df, myvariables = c("Sample_type"), assign_taxonomy = F)
# write.csv(x =genus_group_comparison.df,file ="Result_tables/DESeq_results/Genus_deseq.csv",quote = F, row.names =F)

# Compare all lesion types within each patient
# asv_group_comparison_within_patient.df <- compare_groups_deseq_within_group(mydata.m = asv.m,
#                                                                             mymetadata.df = metadata.df,
#                                                                             myvariables = c("Sample_type"),
#                                                                             within_group_variable = "Subject",
#                                                                             assign_taxonomy = T)
# write.csv(x =asv_group_comparison_within_patient.df,file ="Result_tables/DESeq_results/ASV_within_patient_deseq.csv",quote = F, row.names =F)

# genus_group_comparison_within_patient.df <- compare_groups_deseq_within_group(mydata.m = genus.m,
                                                                              # mymetadata.df = metadata.df,
                                                                              # myvariables = c("Sample_type"),
                                                                              # within_group_variable = "Subject",
                                                                              # assign_taxonomy = F)
# write.csv(x =genus_group_comparison_within_patient.df,file ="Result_tables/DESeq_results/Genus_within_patient_deseq.csv",quote = F, row.names =F)



