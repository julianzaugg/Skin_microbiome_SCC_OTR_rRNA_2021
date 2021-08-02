#************************************
# Main data preparation script for publication:
# "Changes in the skin microbiome associated with squamous cell carcinoma
# in transplant recipients"
# 
# Data is formatted for use in companion scripts.
#************************************


# ------------------------------------------------------------------------------------------------------------
#                                                LOAD AND INSTALL PACKAGES

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# install.packages("phyloseq")
# install.packages("reshape2")
# "The tidyverse is an opinionated collection of R packages designed for data science.", https://www.tidyverse.org
# install.packages("tidyverse")
# install.packages("vegan")
# install.packages("ggplot2")
# install.packages("openxlsx")


library(phyloseq)
library(reshape2)
library(tidyverse)
library(vegan)
library(ggplot2)
library(seqinr)
library(openxlsx)

# BiocManager::install("decontam")
library(decontam)

# install.packages("devtools") #Installs devtools (if not already installed)
# devtools::install_github("donaldtmcknight/microDecon") #Installs microDecon
library(microDecon)

# Load utility functions
source("code/utility.R")

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#                                           COLOUR PALETTES

colour_palette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072", "#80b1d3", "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5", "#cc0000")
# From http://tools.medialab.sciences-po.fr/iwanthue/
colour_palette_10 <- c("#4bafd0","#d14c38","#52a874","#b05cc6","#77b241","#7178ca","#d19b42","#c75a8e","#837f38","#c06a5a")
colour_palette_20 <- c("#66bd79","#a35bcf","#5bb643","#d14ea6","#a2b239","#5c6bcc","#dc892e","#5e93cd","#d64737","#49b6a8","#dc3c6e","#4f7e3c","#bd8cd5","#caab55","#914c88","#867230","#df82a2","#a65429","#ab4a5a","#e0896a")
colour_palette_20_distinct <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782","#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
colour_palette_25_distinct <- c("#005e27","#34adff","#352138","#f10039","#92bb00","#99e0ff","#373300","#680015","#fca000","#d0dbb1","#507400","#9900b7","#3a32cd","#ff6f21","#ff9794","#e700a6","#54032d","#ffc74c","#dfda84","#8d2700","#271472","#eecbf4","#013e77","#d59bff","#ff81ce")
colour_palette_25_distinct_b <- c("#82ff7a","#f1beff","#26324c","#c82be3","#006f91","#7a003e","#efe5ff","#ff5bbc","#9bff42","#004e1c","#ce8c00","#3a246d","#c6b700","#7c1100","#2affbc","#92008e","#9a9900","#c60024","#ff753e","#009b8c","#ffd68a","#009008","#489aff","#6e19cc","#ecffc8")
colour_palette_30_distinct <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8","#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a","#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
colour_palette_206_distinct <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff","#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d","#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c","#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b","#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff","#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4","#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00","#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8","#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100","#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8","#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f","#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff","#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3","#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff","#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614","#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec","#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd","#93f6e7","#01b4a4")
colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba","#b67249","#9b4a6f","#df8398")
colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
# colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#a85bd2","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")

# Subject
subject_palette_45 <- c("#d64530","#585fb1","#795d97","#9e4773","#3f6921","#71692c","#a2b93c","#d571cc","#9b3e97","#33947a","#98ad66","#448a4e","#869ae0","#5ce7af","#e085a3","#dfdc87","#d19be2","#5cb735","#e38269","#3db6c0","#50b565","#50902c","#a98a2c","#dde84a","#db3d76","#5fe485","#7c8329","#b3e791","#6fe965","#5ebce9","#3c86c1","#2a6a45","#65b688","#6651d1","#af4ed3","#df872f","#56e4db","#737cea","#ac464b","#dd37b5","#995b2b","#daac6f","#92e2be","#a2e24b","#e0be3a")

# Lesion / Sampletype
lesion_palette_10 <- c("#d4a33e","#5ca876","#687fc9","#ce5944","#51b2d0","#9b62c8","#d14a8e","#79b041","#bc759a","#9c7f45")

# Ulceration
ulceration_palette <- c(no = "#55a3d8", yes = "#d64a4a")

# Forearm
forearm_palette <- c(no = "#55a3d8", yes = "#d64a4a")

# Cohort
#335fa5 - blue
#c12a2a - red
#61b4c1 blue lighter
#c93434 red lighter
# cohort_palette_2 <- c("#61b4c1", "#c93434")
cohort_palette <- c("organ transplant recipient" = "#c93434", "immunocompetent" = "#61b4c1")

# Gender
gender_palette <- c("#61b4c1", "#c93434")

# Sample_type
sampletype_negative_col <- rgb(red=154,green=154,blue=154,maxColorValue = 255)
sampletype_hs_col <- rgb(red=121,green=168,blue=122,maxColorValue = 255)
sampletype_pds_col <- rgb(red=121,green=175,blue=201,maxColorValue = 255)
sampletype_ak_col <- rgb(red=255,green=196,blue=0,maxColorValue = 255)
sampletype_sccpl_col <- rgb(red=240,green=110,blue=123,maxColorValue = 255)
sampletype_scc_col <- rgb(red=179,green=19,blue=19,maxColorValue = 255)
sample_type_palette_final <- setNames(c(sampletype_negative_col, sampletype_hs_col, sampletype_pds_col, sampletype_ak_col, sampletype_sccpl_col,sampletype_scc_col),
                                      c("negative", "NS", "PDS", "AK", "SCC_PL", "SCC"))


# Length of suppression palette
length_of_suppression_palette <- c("#78a34a","#a464c4","#ce624c")

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#                                             CREATE DIRECTORIES 

dir.create(file.path(".", "Result_figures"), showWarnings = FALSE)
dir.create(file.path("./Result_figures", "various"), showWarnings = FALSE,recursive = T)

dir.create(file.path("./Result_figures", "ordination_plots"), showWarnings = FALSE)
dir.create(file.path("./Result_figures/ordination_plots", "asv"), showWarnings = FALSE)
dir.create(file.path("./Result_figures/ordination_plots", "genus"), showWarnings = FALSE)
dir.create(file.path("./Result_figures", "diversity"), showWarnings = FALSE)
dir.create(file.path("./Result_figures", "abundance_analysis_plots"), showWarnings = FALSE)
dir.create(file.path("./Result_figures", "heatmaps"), showWarnings = FALSE)

dir.create(file.path(".", "Result_tables"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "count_tables"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "relative_abundance_tables"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "combined_counts_abundances_and_metadata_tables"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "other"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "DESeq_results"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "permanova"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "permdisp"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "dunn_tests"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "diversity"), showWarnings = FALSE)
dir.create(file.path("./Result_tables", "abundance_analysis"), showWarnings = FALSE)



dir.create(file.path(".", "Result_objects"), showWarnings = FALSE)

dir.create(file.path(".", "Result_other"), showWarnings = FALSE)
dir.create(file.path("./Result_other", "sequences"), showWarnings = FALSE,recursive = T)


# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# Load the and process metadata
metadata.df <- read.table("data/skin_SCC_metadata.tsv", header = T, sep = "\t")
names(metadata.df)

# Re-assign the index column to the internal name + Job ID.
metadata.df$Index <- with(metadata.df,paste0(Sample_name, "_", Job_ID))

# Make the index the rowname
rownames(metadata.df) <- metadata.df$Index

# Make empty cells NA
metadata.df[metadata.df == ''] <- NA

# For some columns, we want a value of 'no' instead of NA
metadata.df$Forearm_sample[is.na(metadata.df$Forearm_sample)] <- "no"
metadata.df$Ulceration[is.na(metadata.df$Ulceration)] <- "no"

# Change lesion types that are labelled SwabCo to negative (as they are negative)
metadata.df$Sample_type_original <- as.character(metadata.df$Sample_type_original)
metadata.df[metadata.df$Sample_type_original == "SwabCo",]$Sample_type_original <- "negative"
metadata.df$Sample_type_original <- factor(metadata.df$Sample_type_original)

# Get the sample ids
sample_ids <- metadata.df$Index
sample_ids_original <- metadata.df$Index # List of samples prior to filtering (if any) occurs

# Get the negative sample IDs
negative_sample_ids <- as.character(metadata.df[metadata.df$Sample_type == "negative",]$Index)
not_negative_sample_ids <- sample_ids[!sample_ids %in% negative_sample_ids]

# Count of sample types per cohort
metadata.df %>% group_by(Cohort) %>% tally()
metadata.df %>% group_by(Cohort, Sample_type) %>% tally()
metadata.df %>% group_by(Cohort) %>% dplyr::summarise(Mean_age = mean(Age))

# Set the discrete variables
discrete_variables <- c("Sample_type","Cohort", "Subject", 
                        "Length_of_immunosuppression_group_1","Length_of_immunosuppression_group_2",
                        "Gender", "Ulceration","Forearm_sample")

# Assign unique colours for each discrete state for variables
metadata.df <- assign_colours_to_df(metadata.df,
                                    discrete_variables,
                                    auto_assign = T,
                                    my_palette = 
                                      list(Sample_type = sample_type_palette_final,
                                           Cohort = cohort_palette,
                                           Subject = subject_palette_45,
                                           Length_of_immunosuppression_group_1 = length_of_suppression_palette,
                                           Length_of_immunosuppression_group_2 = length_of_suppression_palette,
                                           Gender = gender_palette,
                                           Ulceration = ulceration_palette,
                                           Forearm_sample = forearm_palette
                                           ),
                                    # my_default_palette = colour_palette_10_distinct
                                    my_default_palette = colour_palette_soft_8
)
# unique(metadata.df[c("Ulceration","Ulceration_colour")])
# unique(metadata.df[c("Forearm_sample","Forearm_sample_colour")])
# unique(metadata.df[c("Gender","Gender_colour")])
# unique(metadata.df[c("Cohort","Cohort_colour")])
# unique(metadata.df[c("Sample_type","Sample_type_colour")])
# unique(metadata.df[c("Length_of_immunosuppression_group_1","Length_of_immunosuppression_group_1_colour")])

# Assign unique shapes for each discrete state for variables
shapes_for_sampletype <- setNames(c(25,24,23,22,21), c("NS","PDS", "AK", "SCC_PL", "SCC"))
metadata.df$Sample_type_shape <- unlist(lapply(metadata.df$Sample_type, function(x) shapes_for_sampletype[x][[1]]))

# Retain copy of metadata prior to filtering samples
metadata_unfiltered.df <- metadata.df
write.csv(metadata.df,file = "Result_tables/processed_unfiltered_metadata.csv", row.names = F, quote = F)
# ------------------------------------------------------------------------------------------------------------

# Load the ASV/feature count table into a dataframe object
mydata.df <- read.csv(file = "data/skin_SCC_feature_counts.csv", header = T)

mydata.df <- mydata.df %>% select(-Confidence, -Frequency)

# Get the sample ids from the feature table
sample_ids <- names(mydata.df)[!names(mydata.df) %in% c("OTU.ID","ASV_MD5","ASV","ASV_ID","Frequency", "Taxon", "Confidence", "RepSeq", "Sequence") ]
summary(sample_ids %in% rownames(metadata.df))
summary( rownames(metadata.df) %in% sample_ids )

# (optional) Sort the dataframe by the ASV hash 
mydata.df <- mydata.df[order(mydata.df$ASV_MD5),]

# Assign short form ID for each ASV
mydata.df$ASV_ID <- paste0("ASV_", 1:length(mydata.df$ASV))

# Remove spaces after semicolon in taxonomy string
mydata.df$Taxon <- gsub("; ", ";", mydata.df$Taxon)


# Split the Taxon string into Domain, Phylum...Species columns
mydata.df <- separate(mydata.df, "Taxon", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"), remove =F, sep = ";")

# Splitting taxonomy strings  that are not specified at certain taxonomy levels will produce 'NA' entries at those levels.
# NA entries should be changed to "Unassigned"
mydata.df[is.na(mydata.df)] <- "Unassigned"

# Create full taxonomy strings for phylum, class, order, family, genus and specie levels
mydata.df$taxonomy_phylum <- with(mydata.df, paste(Domain, Phylum, sep =";"))
mydata.df$taxonomy_class <- with(mydata.df, paste(Domain, Phylum, Class, sep =";"))
mydata.df$taxonomy_order <- with(mydata.df, paste(Domain, Phylum, Class, Order, sep =";"))
mydata.df$taxonomy_family <- with(mydata.df, paste(Domain, Phylum, Class, Order, Family, sep =";"))
mydata.df$taxonomy_genus <- with(mydata.df, paste(Domain, Phylum, Class, Order, Family, Genus, sep =";"))
mydata.df$taxonomy_species <- with(mydata.df, paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep =";"))


# Store a version of the unfiltered data
mydata_unfiltered.df <- mydata.df

# (optional) Remove entries in the project (count) table that have no counts
# zero_count_samples <- sample_ids[colSums(mydata.df[,sample_ids]) == 0]
# mydata.df <- mydata.df[,!names(mydata.df) %in% zero_count_samples]
# sample_ids <- sample_ids[!sample_ids %in% zero_count_samples]
# ---------------------------------------------------------------------------------------------------------------- 
#  ----------------------------------------------------------------------------------------------------------------
#  ------------------------------------------- Remove unwanted lineages ------------------------------------------- 

fungal_phyla <- c("Basidiomycota","Ascomycota","Mucoromycota","Cryptomycota",
                  "Peronosporomycetes","Blastocladiomycota","Chytridiomycota",
                  "Zoopagomycota","Neocallimastigomycota")
#"Aphelidea","LKM15" - unclear if these are fungal
fungal_pattern <- paste0(fungal_phyla,collapse = "|p__")
fungal_pattern <- gsub("^", "p__", fungal_pattern)

# Discard anything not Bacterial or Archaeal or fungal
mydata.df <- mydata.df[grepl(paste0("d__Bacteria|d__Archaea|",fungal_pattern), mydata.df$Taxon),]

# Discard anything that is Unassigned at the Phylum level
mydata.df <- mydata.df[!mydata.df$Phylum == "Unassigned",]

# Discard chloroplast features
mydata.df <- mydata.df[!grepl("o__Chloroplast", mydata.df$Taxon,ignore.case = T),]

# Discard mitochondria features
mydata.df <- mydata.df[!grepl("f__Mitochondria", mydata.df$Taxon,ignore.case = T),]

# ---------------------------------------------------------------------------------------------------------------- 
#  ----------------------------------------------------------------------------------------------------------------
# Remove original Taxon column
mydata.df$Taxon <- NULL
mydata_unfiltered.df$Taxon <- NULL

# Store the ASVs, corresponding taxonomy information and sequence in a separate dataframe
asv_taxonomy_map.df <- mydata.df[c("ASV_MD5",
                                   "ASV_ID",
                                   "Domain", 
                                   "Phylum", 
                                   "Class", 
                                   "Order", 
                                   "Family",
                                   "Genus",
                                   "Species",
                                   "taxonomy_phylum",
                                   "taxonomy_class",
                                   "taxonomy_order",
                                   "taxonomy_family",
                                   "taxonomy_genus",
                                   "taxonomy_species", 
                                   "RepSeq")]

# Save this ASV taxonomy map for later reference
write.table(asv_taxonomy_map.df, file = "Result_tables/asv_taxonomy_map.csv", sep = ",", quote = F, row.names = F)

# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
#                                         GENERATE TAXONOMY LEVEL TABLES

# Now we can generate the tables that we will need for different analyses at both the ASV and higher taxa levels

# First at the ASV level
asv.df <- mydata.df[c("ASV_ID", sample_ids)]
asv_unfiltered.df <- mydata_unfiltered.df[c("ASV_ID", sample_ids)]

asv.m <- df2m(asv.df) # count matrix
asv_unfiltered.m <- df2m(asv_unfiltered.df) # count matrix
temp <- read_counts_and_unique_features(asv.m, sample_ids)
df2m(temp$sample_read_counts)
df2m(temp$sample_feature_counts)

# Create relative abundance matrix from counts matrix
asv_rel.m <- t(t(asv.m)/ colSums(asv.m))
asv_unfiltered_rel.m <- t(t(asv_unfiltered.m)/ colSums(asv_unfiltered.m))

# Change nans to 0. Occurs when a sample has no hits at this point.
asv_rel.m[is.nan(asv_rel.m)] <- 0
asv_unfiltered_rel.m[is.nan(asv_unfiltered_rel.m)] <- 0

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
#                                     CONTAMINANT IDENTIFICATION AND REMOVAL

# For each cohort, identify contaminants and remove them.

# One simple approach to removing contaminants, or 'noise', is to simply remove low abundance features or
# remove features that are present / more prevalent in the negative controls.
# Alternatively (or in addition), packages like microdecon and decontam can be used to identity putative contaminants.

# Since we are processing each cohort separately, create feature tables for each cohort
immunosuppressed_asv.m <- asv.m[,rownames(metadata.df[metadata.df$Cohort == "organ transplant recipient",])]
immunosuppressed_asv_rel.m <- asv_rel.m[,rownames(metadata.df[metadata.df$Cohort == "organ transplant recipient",])]

immunocompetent_asv.m <- asv.m[,rownames(metadata.df[metadata.df$Cohort == "immunocompetent",])]
immunocompetent_asv_rel.m <- asv_rel.m[,rownames(metadata.df[metadata.df$Cohort == "immunocompetent",])]

# Get the sample ids for each cohort
immunosuppressed_sample_ids <- as.character(subset(metadata.df, Cohort == "organ transplant recipient")$Index)
immunocompetent_sample_ids <- as.character(subset(metadata.df, Cohort == "organ transplant recipient")$Index)

# Get the negative sample IDs for each cohort
immunosuppressed_negative_sample_ids <- as.character(subset(metadata.df, Sample_type_original == "negative" & Cohort == "organ transplant recipient")$Index)
immunocompetent_negative_sample_ids <- as.character(subset(metadata.df, Sample_type_original == "negative" & Cohort == "immunocompetent")$Index)

# Create feature table for negative samples
negative_asv.m <- asv.m[,negative_sample_ids]
negative_asv_rel.m <- asv_rel.m[,negative_sample_ids]
immunosuppressed_negative_asv.m <- asv.m[,immunosuppressed_negative_sample_ids]
immunosuppressed_negative_asv_rel.m <- asv_rel.m[,immunosuppressed_negative_sample_ids]
immunocompetent_negative_asv.m <- asv.m[,immunocompetent_negative_sample_ids]
immunocompetent_negative_asv_rel.m <- asv_rel.m[,immunocompetent_negative_sample_ids]

# ------------------------------------------------------------------------
# 1. Use Decontam to identify contaminants
decontam_contaminants_immunosuppressed.df <- isContaminant(t(immunosuppressed_asv.m), 
                                                           method = "prevalence", 
                                                           neg = rownames(t(immunosuppressed_asv.m)) %in% immunosuppressed_negative_sample_ids, 
                                                           threshold = 0.5)
decontam_contaminants_immunocompetent.df <- isContaminant(t(immunocompetent_asv.m), 
                                                          method = "prevalence", 
                                                          neg = rownames(t(immunocompetent_asv.m)) %in% immunocompetent_negative_sample_ids, 
                                                          threshold = 0.5)

hist(decontam_contaminants_immunosuppressed.df$p,100, ylim = c(0,7000), xlim = c(0,1))
abline(v = 0.5, col = 'red')
hist(decontam_contaminants_immunocompetent.df$p,100, ylim = c(0,7000), xlim = c(0,1))
abline(v = 0.5, col = 'red')

decontam_contaminants_features_immunosuppressed <- rownames(subset(decontam_contaminants_immunosuppressed.df, contaminant == T))
decontam_contaminants_features_immunocompetent <- rownames(subset(decontam_contaminants_immunocompetent.df, contaminant == T))
# ------------------------------------------------------------------------
# 2. Use microdecon
# SKIP

# asv Neg1 ....Sample 1....taxa (optional)
# microdecon_data_immunosuppressed.df <- m2df(immunosuppressed_asv.m[,colnames(immunosuppressed_asv.m) %in% immunosuppressed_negative_sample_ids], column_name = "ASV_ID")
# microdecon_data_immunosuppressed.df <- cbind(microdecon_data_immunosuppressed.df, immunosuppressed_asv.m[,!colnames(immunosuppressed_asv.m) %in% immunosuppressed_negative_sample_ids])
# microdecon_data_immunosuppressed.df$taxonomy_species <- subset(asv_taxonomy_map.df, ASV_ID %in% microdecon_data_immunosuppressed.df$ASV_ID)$taxonomy_species
# microdecon_contaminants_immunosuppressed.df <- decon(microdecon_data_immunosuppressed.df,
#                                     numb.blanks = length(immunosuppressed_negative_sample_ids),
#                                     numb.ind = length(immunosuppressed_sample_ids) - length(immunosuppressed_negative_sample_ids),
#                                     taxa = T,
#                                     runs = 2,regression = 1)
# microdecon_contaminants_features_immunosuppressed <- microdecon_contaminants_immunosuppressed.df$OTUs.removed$ASV_ID
# 
# microdecon_data_immunocompetent.df <- m2df(immunocompetent_asv.m[,colnames(immunocompetent_asv.m) %in% immunocompetent_negative_sample_ids], column_name = "ASV_ID")
# microdecon_data_immunocompetent.df <- cbind(microdecon_data_immunocompetent.df, immunocompetent_asv.m[,!colnames(immunocompetent_asv.m) %in% immunocompetent_negative_sample_ids])
# microdecon_data_immunocompetent.df$taxonomy_species <- subset(asv_taxonomy_map.df, ASV_ID %in% microdecon_data_immunocompetent.df$ASV_ID)$taxonomy_species
# microdecon_contaminants_immunocompetent.df <- decon(microdecon_data_immunocompetent.df,
#                                                     numb.blanks = length(immunocompetent_negative_sample_ids),
#                                                     numb.ind = length(immunocompetent_sample_ids) - length(immunocompetent_negative_sample_ids),
#                                                     taxa = T,
#                                                     runs = 2,regression = 1)
# microdecon_contaminants_features_immunocompetent <- microdecon_contaminants_immunocompetent.df$OTUs.removed$ASV_ID
# ------------------------------------------------------------------------
# 3. Identify features that are present in negative controls (extreme approach) 
negative_present_features <- rownames(negative_asv.m[rowSums(negative_asv.m) > 0,])
negative_present_features_immunosuppressed <- rownames(immunosuppressed_negative_asv.m[rowSums(immunosuppressed_negative_asv.m) > 0,])
negative_present_features_immunocompetent <- rownames(immunocompetent_negative_asv.m[rowSums(immunocompetent_negative_asv.m) > 0,])

# ------------------------------------------------------------------------
# 4. Identify features that are more prevalent in negative controls.
# Use 'rarefied' (capped) data to reduce false positives
# e.g. asv_rare_count.m <- t(rrarefy(x = t(asv_count_raw.m), sample = 2000))
# Large differences in read depth between samples will make it harder to determine the true-positive contaminants from false-positive.
asv_rare_count.m <- t(rrarefy(x = t(asv.m), sample=2000))
negative_asv_rare_count.m <- t(rrarefy(x = t(negative_asv.m), sample=2000))
immunosuppressed_asv_rare.m <-  t(rrarefy(x = t(immunosuppressed_asv.m), sample=2000))
immunocompetent_asv_rare.m <-  t(rrarefy(x = t(immunocompetent_asv.m), sample=2000))
immunosuppressed_negative_asv_rare.m <- t(rrarefy(x = t(immunosuppressed_negative_asv.m), sample=2000))
immunocompetent_negative_asv_rare.m <- t(rrarefy(x = t(immunocompetent_negative_asv.m), sample=2000))

asv_negative_sample_prevalences <- apply(negative_asv_rare_count.m, 1,
                                         function(x) {length(which(x > 0))}) / length(negative_sample_ids)
immunosuppressed_negative_asv_sample_prevalences <- apply(immunosuppressed_negative_asv_rare.m, 1, 
                                                          function(x) {length(which(x > 0))}) / length(immunosuppressed_negative_sample_ids)

immunocompetent_negative_asv_sample_prevalences <- apply(immunocompetent_negative_asv_rare.m, 1, 
                                                         function(x) {length(which(x > 0))}) / length(immunocompetent_negative_sample_ids)



asv_not_negative_sample_prevalences <- apply(asv_rare_count.m[,not_negative_sample_ids], 1, function(x) {length(which(x > 0))}) /length(not_negative_sample_ids)
contaminating_asvs_from_prevalences <- names(asv_negative_sample_prevalences[asv_negative_sample_prevalences > asv_not_negative_sample_prevalences])

immunosuppressed_not_negative_sample_ids <- colnames(immunosuppressed_asv_rare.m)[colnames(immunosuppressed_asv_rare.m) %in% not_negative_sample_ids]
immunosuppressed_asv_not_negative_sample_prevalences <-
  apply(immunosuppressed_asv_rare.m[,immunosuppressed_not_negative_sample_ids], 1,
        function(x) {length(which(x > 0))}) / length(immunosuppressed_not_negative_sample_ids)
immunosuppressed_contaminating_asvs_from_prevalences <- names(immunosuppressed_negative_asv_sample_prevalences[immunosuppressed_negative_asv_sample_prevalences > immunosuppressed_asv_not_negative_sample_prevalences])

immunocompetent_not_negative_sample_ids <- colnames(immunocompetent_asv_rare.m)[colnames(immunocompetent_asv_rare.m) %in% not_negative_sample_ids]
immunocompetent_asv_not_negative_sample_prevalences <-
  apply(immunocompetent_asv_rare.m[,immunocompetent_not_negative_sample_ids], 1,
        function(x) {length(which(x > 0))}) / length(immunocompetent_not_negative_sample_ids)
immunocompetent_contaminating_asvs_from_prevalences <- names(immunocompetent_negative_asv_sample_prevalences[immunocompetent_negative_asv_sample_prevalences > immunocompetent_asv_not_negative_sample_prevalences])


# ------------------------------------------------------------------------
# Final contaminating features to remove
contaminating_features <- unique(c(decontam_contaminants_features_immunosuppressed, decontam_contaminants_features_immunocompetent))
contaminating_features_immunosuppressed <- unique(c(decontam_contaminants_features_immunosuppressed))
contaminating_features_immunocompetent <- unique(c(decontam_contaminants_features_immunocompetent))

# Remove contaminants from data (all samples).
# asv_decontaminated.m <- asv.m[!rownames(asv.m) %in% unique(c(contaminating_features)),]
# asv_rel_decontaminated.m <- asv_rel.m[!rownames(asv_rel.m) %in% unique(c(contaminating_features)),]

# Remove contaminants from data for each cohort
immunosuppressed_asv.m <- immunosuppressed_asv.m[!rownames(immunosuppressed_asv.m) %in% unique(c(contaminating_features_immunosuppressed)),]
immunosuppressed_asv_rel.m <- immunosuppressed_asv_rel.m[!rownames(immunosuppressed_asv_rel.m) %in% unique(c(contaminating_features_immunosuppressed)),]
immunocompetent_asv.m <- immunocompetent_asv.m[!rownames(immunocompetent_asv.m) %in% unique(c(contaminating_features_immunocompetent)),]
immunocompetent_asv_rel.m <- immunocompetent_asv_rel.m[!rownames(immunocompetent_asv_rel.m) %in% unique(c(contaminating_features_immunocompetent)),]

# Combine cohort tables together to create new primary otu.m that contains decontaminated counts
# DO NOT re-normalise the abundances as we need them to summarise the contaminates
temp1 <- melt(immunosuppressed_asv.m)
temp2 <- melt(immunocompetent_asv.m)
temp1 <- temp1[temp1$value != 0,]
temp2 <- temp2[temp2$value != 0,]
temp3 <- rbind(temp1,temp2)

temp4 <- spread(temp3,key = "Var2", value = "value",fill = 0)

# If samples had a zero count after removing contaminants, add them back
missing_samples <- c(colnames(immunosuppressed_asv.m), colnames(immunocompetent_asv.m))[!c(colnames(immunosuppressed_asv.m), colnames(immunocompetent_asv.m))  %in% colnames(temp4)]
temp4[,missing_samples] <- 0
asv_decontaminated.m <- df2m(temp4)

# temp <- read_counts_and_unique_features(asv_decontaminated.m, sample_ids_original)
# df2m(temp$sample_read_counts)
# df2m(temp$sample_feature_counts)

# Re-generate relative abundances
asv_rel_decontaminated.m <- t(t(asv_decontaminated.m) / colSums(asv_decontaminated.m))
asv_rel_decontaminated.m[is.nan(asv_rel_decontaminated.m)] <- 0

# (Optional) Use de-contaminated data going forward
asv_rel.m <- asv_rel_decontaminated.m
asv.m <- asv_decontaminated.m

# colSums(asv.m)["SB4911_J1426"]

# (Optional) Remove negative samples as they are no longer needed
# asv_rel.m <- asv_rel.m[,!colnames(asv_rel.m) %in% negative_sample_ids]
# asv.m <- asv.m[,!colnames(asv.m) %in% negative_sample_ids]
# negative_metadata.df <- metadata.df[rownames(metadata.df) %in% negative_sample_ids,]
# metadata.df <- metadata.df[!rownames(metadata.df) %in% negative_sample_ids,]

# Reassign the sample ids
sample_ids <- colnames(asv.m)

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------

# Filter those ASVs that are low abundance in all samples
# Check how many ASVs would be removed if we filtered any whose abundance is less than 0.05% (0.0005)
filter_percentage = 0.0005
asv_rel_low_abundance_asvs.m <- asv_rel.m[apply(asv_rel.m[,sample_ids],1,function(z) all(z < filter_percentage)),]
print(paste("There are a total of", dim(asv_rel.m)[1], "OTUs before filtering"))
print(paste("A total of", dim(asv_rel_low_abundance_asvs.m)[1], 
            "ASVs will be filtered at a threshold of", 
            filter_percentage * 100, "percent"))

# If you are happy with the filtering threshold, apply it.
asv_rel.m <- asv_rel.m[apply(asv_rel.m[,sample_ids],1,function(z) any(z>=filter_percentage)),]

# Re-normalise the abundance matrix after filtering
asv_rel.m <- t(t(asv_rel.m) / colSums(asv_rel.m))

# Change nans to 0. Occurs when a sample has no hits at this point (if filtered out earlier)
asv_rel.m[is.nan(asv_rel.m)] <- 0
asv_rel.df <- m2df(asv_rel.m, "ASV_ID")

# Also remove low abundance ASVs from the ASV count matrix
asv.m  <- asv.m[rownames(asv_rel.m),]

# (optional) Discard samples with less than # reads.
count_threshold <- 2000
asv_pre_low_count_filter.m <- asv.m
asv_pre_low_count_filter_rel.m <- asv_rel.m

summary(colSums(asv.m) < count_threshold)
dim(asv.m)
asv.m <- asv.m[,colSums(asv.m) >= count_threshold]
dim(asv.m)

# There might be rows whose maximum is 0 at this point. Remove them.
dim(asv.m)
asv.m <- asv.m[apply(asv.m, 1, max) != 0,]
dim(asv.m)

# Reassign the sample ids
sample_ids <- colnames(asv.m)



# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# GET ABUNDANT UNASSIGNED FEATURES

# Original unfiltered data with just unassigned features
unassigned_mydata_unfiltered.df <- mydata_unfiltered.df[mydata_unfiltered.df$Domain == "Unassigned",]
rownames(unassigned_mydata_unfiltered.df) <- unassigned_mydata_unfiltered.df$ASV_ID

# Get the abundances for the unassigned features and convert to dataframe
unassigned_asv_unfiltered_rel.df <- m2df(asv_unfiltered_rel.m[unassigned_mydata_unfiltered.df$ASV_ID,], "ASV_ID")

# Melt dataframe
unassigned_asv_unfiltered_rel.df <- melt(unassigned_asv_unfiltered_rel.df, 
                                         variable.name = "Sample", value.name = "Relative_abundance")

if (max(unassigned_asv_unfiltered_rel.df$Relative_abundance) != 0){
  
  # Get the # top most abundant unassigned features per sample
  most_abundant_unassigned.df <- unassigned_asv_unfiltered_rel.df %>%
    group_by(Sample) %>%
    filter(Relative_abundance > 0) %>%
    top_n(n = 10, wt = Relative_abundance) %>%
    mutate(Relative_abundance = round(Relative_abundance*100,3)) %>%
    as.data.frame()

  temp <- melt(round(colSums(unassigned_mydata_unfiltered.df[unique(most_abundant_unassigned.df$ASV_ID),sample_ids])/colSums(mydata_unfiltered.df[,sample_ids]) * 100,4), value.name = "Relative_abundance")
  length(unique(most_abundant_unassigned.df$ASV_ID)) 
  length(unique(most_abundant_unassigned.df$ASV_ID)) / length(mydata_unfiltered.df$ASV_ID) *100
  temp <- m2df(temp,column_name = "Index")
  temp <- merge(temp,metadata.df,by = "Index") %>% select(Index, Sample_type, Cohort, Relative_abundance)
  temp <- temp[order(temp$Relative_abundance,decreasing = T),]
  # Note — these stats are based on only the top 10, not all
  print(max(temp$Relative_abundance))
  print(min(temp$Relative_abundance))
  print(mean(temp$Relative_abundance))
  print(median(temp$Relative_abundance))
  
  # Get sequences corresponding to top features
  most_abundant_unassigned.df$RepSeq <- unlist(lapply(most_abundant_unassigned.df$ASV_ID,
                                                      function(x) as.character(unassigned_mydata_unfiltered.df[unassigned_mydata_unfiltered.df$ASV_ID == x,]$RepSeq)))
  write.csv(x = most_abundant_unassigned.df, file = paste0("Result_tables/other/most_abundant_unassigned.csv"), row.names = F)
  
  # Get unique set of features
  unique_most_abundant_unassigned.df <- unique(most_abundant_unassigned.df[c("ASV_ID", "RepSeq")])
  
  # Write fasta file
  write.fasta(sequences = as.list(unique_most_abundant_unassigned.df$RepSeq),open = "w",
              names = as.character(unique_most_abundant_unassigned.df$ASV_ID),
              file.out = paste0("Result_other/sequences/most_abundant_unassigned_features.fasta"))
}

round(max(colSums(asv_unfiltered_rel.m[unique(unique_most_abundant_unassigned.df$ASV_ID),]))*100,3)
round(min(colSums(asv_unfiltered_rel.m[unique(unique_most_abundant_unassigned.df$ASV_ID),]))*100,3)
round(mean(colSums(asv_unfiltered_rel.m[unique(unique_most_abundant_unassigned.df$ASV_ID),]))*100,3)
round(median(colSums(asv_unfiltered_rel.m[unique(unique_most_abundant_unassigned.df$ASV_ID),]))*100,3)

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
#         (optional) RAREFYING (capping maximum reads)

# Note - For samples with a read count lower then the sample=# parameter, 
# the rrarefy function in vegan will return the this sample with its existing read count.
# Sample with counts higher than the sample parameter will be rarefied as normal and 
# their counts will be capped at the sample parameter value


# Counts for each Sample
column_sums <- colSums(asv.m)
column_sums.df <- m2df(melt(column_sums[order(column_sums)]),"sample")
column_sums.df <- column_sums.df[c("sample", "value")]
column_sums.df$sample <- factor(column_sums.df$sample, levels = column_sums.df$sample)


myplot <- ggplot(column_sums.df, aes(x = sample, y = value)) + 
  geom_histogram(stat = "identity", fill = "grey", colour = "black", lwd = .3, width = 0.75) +
  # geom_hline(yintercept = 30000, color = 'red')+
  # geom_hline(yintercept = 20000, color = 'red')+
  geom_hline(yintercept = mean(column_sums.df$value), color = 'blue')+
  geom_hline(yintercept = median(column_sums.df$value), color = 'purple')+
  scale_y_continuous(breaks = seq(0,max(column_sums.df$value), 5000),expand = c(0,0)) +
  xlab("Sample") +
  ylab("Read count") +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust= 1,size= 6),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white", size = 1),
        axis.line = element_line(colour = "black", size = 0.5),
  )
myplot
ggsave(plot = myplot, filename = "./Result_figures/various/reads_per_sample_hist.pdf", width=20, height=6)

myplot <- ggplot(column_sums.df, aes(x = value)) + 
  xlab("Read count") +
  ylab("Number of samples") +
  scale_x_continuous(breaks = seq(500,max(column_sums.df$value)+2000, 2000)) +
  scale_y_continuous(breaks = seq(0,10, 1)) +
  geom_histogram(stat = "bin", bins = 50, colour = "black",fill = "grey") +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust= 1,size= 6),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white", size = 1),
        axis.line = element_line(colour = "black", size = 0.5),
  )
myplot
ggsave(plot = myplot, filename = "./Result_figures/various/samples_per_read_hist.pdf", width=10, height=6)


# Rarefaction curve
rarecurve(t(asv.m[,colSums(asv.m) > 1]),step = 500, label = F,xlim = c(0,30000))

# Generate a rarefied ASV count matrix. The rrarefy function just caps samples to a maximum of threshold. 
# Samples with less will still have less!!
set.seed(1234)
asv_rare_count.m <- t(rrarefy(x = t(asv.m), sample=30000))

# How many samples have reads removed
summary(colSums(asv_rare_count.m - asv.m) < 0)

# (optional) only use rrarefied capped counts
asv.m <- asv_rare_count.m
asv.df <- m2df(asv.m, "ASV_ID")
# temp <- read_counts_and_unique_features(asv.m, sample_ids)
# df2matrix(temp$sample_read_counts)
# df2matrix(temp$sample_feature_counts)

# -----------------------------------------------------
# -----------------------------------------------------

# And re-calculate the abundances after filtering/rarefying (if done)
asv_rel.m <- t(t(asv.m)/ colSums(asv.m))
asv_rel.m[is.nan(asv_rel.m)] <- 0

# Remove negative samples from metadata and feature tables
not_negative_sample_ids <- colnames(asv.m)[!colnames(asv.m) %in% negative_sample_ids]
asv_rel.m <- asv_rel.m[,not_negative_sample_ids]
asv.m <- asv.m[,not_negative_sample_ids]
asv.df <- m2df(asv.m, "ASV_ID")
asv_rel.df <- m2df(asv_rel.m, "ASV_ID")
metadata.df <- metadata.df[metadata.df$Index %in% colnames(asv.m),]

# Reassign sample IDs
sample_ids <- colnames(asv_rel.m)

# Write the final ASV counts and abundances to file
write.table(asv.df, file = "Result_tables/count_tables/ASV_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(asv_rel.df, file = "Result_tables/relative_abundance_tables/ASV_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)

# Write the metadata for the remaining samples to file
write.csv(metadata.df, "Result_tables/processed_metadata.csv", row.names = F, quote = F)

# Above we processed the frequencies for each ASV to calculate the relative abundances.
# However, we also want the abundances at not just the ASV level, but also different taxonomy levels.
# Now we will generate the abundance tables at each taxonomy level from Phylum, Class, Order, Family, Genus and Species

# First merge the ASV counts back with the taxonomy data
asv_taxonomy_merged.df <- left_join(asv.df, asv_taxonomy_map.df, by.x = "ASV_ID", by.y = "ASV_ID")

# Run the taxa summarising function at each taxonomy level 
tax_data.l <- generate_tax_level_data(asv_taxonomy_merged.df, 
                                      tax_string_levels = c("taxonomy_species", "taxonomy_genus", "taxonomy_family", "taxonomy_class", "taxonomy_order", "taxonomy_phylum"),
                                      sample_ids = sample_ids,
                                      remove_zero_row_entries = T)

# Extract the counts and abundances into separate dataframe objects
species.df <- m2df(tax_data.l$taxonomy_species$counts,"taxonomy_species")
species_rel.df <-  m2df(tax_data.l$taxonomy_species$abundances,"taxonomy_species")
genus.df <- m2df(tax_data.l$taxonomy_genus$counts,"taxonomy_genus")
genus_rel.df <-  m2df(tax_data.l$taxonomy_genus$abundances,"taxonomy_genus")
family.df <- m2df(tax_data.l$taxonomy_family$counts,"taxonomy_family")
family_rel.df <-  m2df(tax_data.l$taxonomy_family$abundances,"taxonomy_family")
order.df <- m2df(tax_data.l$taxonomy_order$counts,"taxonomy_order")
order_rel.df <-  m2df(tax_data.l$taxonomy_order$abundances,"taxonomy_order")
class.df <- m2df(tax_data.l$taxonomy_class$counts,"taxonomy_class")
class_rel.df <-  m2df(tax_data.l$taxonomy_class$abundances,"taxonomy_class")
phylum.df <- m2df(tax_data.l$taxonomy_phylum$counts,"taxonomy_phylum")
phylum_rel.df <-  m2df(tax_data.l$taxonomy_phylum$abundances,"taxonomy_phylum")

# Write the final counts and abundances for each taxonomy level to file
write.table(species.df, file = "Result_tables/count_tables/Specie_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(species_rel.df, file = "Result_tables/relative_abundance_tables/Specie_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus.df, file = "Result_tables/count_tables/Genus_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus_rel.df, file = "Result_tables/relative_abundance_tables/Genus_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family.df, file = "Result_tables/count_tables/Family_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family_rel.df, file = "Result_tables/relative_abundance_tables/Family_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order.df, file = "Result_tables/count_tables/Order_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order_rel.df, file = "Result_tables/relative_abundance_tables/Order_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class.df, file = "Result_tables/count_tables/Class_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class_rel.df, file = "Result_tables/relative_abundance_tables/Class_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum.df, file = "Result_tables/count_tables/Phylum_counts.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum_rel.df, file = "Result_tables/relative_abundance_tables/Phylum_relative_abundances.csv", sep = ",", quote = F, col.names = T, row.names = F)


# Finally create and save a dataframe, separately for each Phylum, Class, Order, Family, Genus, Species and ASV,
# containing the abundances/counts in each sample, metadata and taxonomy information.
create_combined_dataframe <- function(counts.df, abundances.df, metadata.df, mylevel = "ASV_ID", sample_id_column = "Index", asv_taxonomy_map.df = NULL){
  counts <- df_column_to_rownames(counts.df)
  rel_abundances <- df_column_to_rownames(abundances.df)
  
  # Ensure ordering is the same
  rel_abundances <- rel_abundances[rownames(counts),,drop = F]
  
  # Combine the datasets. Passing as.matrix(counts) captures the rownames as a column. This can be renamed after
  combined_data <- cbind(melt(as.matrix(counts), variable.name = "sample", value.name = "Read_count"),
                         melt(rel_abundances, value.name = "Relative_abundance")[,2, drop = F])
  
  # Remove samples with a read count of zero
  combined_data <- combined_data[combined_data$Read_count > 0,]
  
  # Calculate logged read counts
  combined_data$Read_count_log10 <- log(combined_data$Read_count, 10)
  
  # Fix the Var2 column
  names(combined_data)[2] <- "Sample"
  # Merge with metadata. Assumes an Index column matching Sample
  combined_data <- merge(combined_data, metadata.df, by.x = "Sample", by.y = sample_id_column)
  if (mylevel == "ASV_ID"){
    names(combined_data)[names(combined_data) == "Var1"] <- "ASV_ID"
    asv_taxonomy_map_reduced.df <- unique(asv_taxonomy_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species", "taxonomy_species", "ASV_MD5","ASV_ID")])
    combined_data <- merge(combined_data, asv_taxonomy_map_reduced.df, by.x = "ASV_ID", by.y = "ASV_ID")
  }
  else if (mylevel == "Species"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_species"
    asv_taxonomy_map_reduced.df <- unique(asv_taxonomy_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species", "taxonomy_species")])
    combined_data <- merge(combined_data, asv_taxonomy_map_reduced.df, by.x = "taxonomy_species", by.y = "taxonomy_species")
  }
  else if (mylevel == "Genus"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_genus"
    asv_taxonomy_map_reduced.df <- unique(asv_taxonomy_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "Genus", "taxonomy_genus")])
    combined_data <- merge(combined_data, asv_taxonomy_map_reduced.df, by.x = "taxonomy_genus", by.y = "taxonomy_genus")
  }
  else if (mylevel == "Family"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_family"
    asv_taxonomy_map_reduced.df <- unique(asv_taxonomy_map.df[,c("Domain","Phylum", "Class", "Order", "Family", "taxonomy_family")])
    combined_data <- merge(combined_data, asv_taxonomy_map_reduced.df, by.x = "taxonomy_family", by.y = "taxonomy_family")
  }
  else if (mylevel == "Order"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_order"
    asv_taxonomy_map_reduced.df <- unique(asv_taxonomy_map.df[,c("Domain","Phylum", "Class", "Order", "taxonomy_order")])
    combined_data <- merge(combined_data, asv_taxonomy_map_reduced.df, by.x = "taxonomy_order", by.y = "taxonomy_order")
  }
  else if (mylevel == "Class"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_class"
    asv_taxonomy_map_reduced.df <- unique(asv_taxonomy_map.df[,c("Domain","Phylum", "Class", "taxonomy_class")])
    combined_data <- merge(combined_data, asv_taxonomy_map_reduced.df, by.x = "taxonomy_class", by.y = "taxonomy_class")
  }
  else if (mylevel == "Phylum"){
    names(combined_data)[names(combined_data) == "Var1"] <- "taxonomy_phylum"
    asv_taxonomy_map_reduced.df <- unique(asv_taxonomy_map.df[,c("Domain","Phylum", "taxonomy_phylum")])
    combined_data <- merge(combined_data, asv_taxonomy_map_reduced.df, by.x = "taxonomy_phylum", by.y = "taxonomy_phylum")
  }
  return(combined_data)
}


# Generate the combined datasets
phylum_combined <- create_combined_dataframe(counts.df = phylum.df,
                                             abundances.df = phylum_rel.df,
                                             metadata.df = metadata.df,
                                             mylevel = "Phylum",
                                             sample_id_column = "Index",
                                             asv_taxonomy_map.df = asv_taxonomy_map.df)

class_combined <- create_combined_dataframe(counts.df = class.df,
                                            abundances.df = class_rel.df,
                                            metadata.df = metadata.df,
                                            mylevel = "Class",
                                            sample_id_column = "Index",
                                            asv_taxonomy_map.df = asv_taxonomy_map.df)

order_combined <- create_combined_dataframe(counts.df = order.df,
                                            abundances.df = order_rel.df,
                                            metadata.df = metadata.df,
                                            mylevel = "Order",
                                            sample_id_column = "Index",
                                            asv_taxonomy_map.df = asv_taxonomy_map.df)

family_combined <- create_combined_dataframe(counts.df = family.df,
                                             abundances.df = family_rel.df,
                                             metadata.df = metadata.df,
                                             mylevel = "Family",
                                             sample_id_column = "Index",
                                             asv_taxonomy_map.df = asv_taxonomy_map.df)

genus_combined <- create_combined_dataframe(counts.df = genus.df,
                                            abundances.df = genus_rel.df,
                                            metadata.df = metadata.df,
                                            mylevel = "Genus",
                                            sample_id_column = "Index",
                                            asv_taxonomy_map.df = asv_taxonomy_map.df)

species_combined <- create_combined_dataframe(counts.df = species.df,
                                              abundances.df = species_rel.df,
                                              metadata.df = metadata.df,
                                              mylevel = "Species",
                                              sample_id_column = "Index",
                                              asv_taxonomy_map.df = asv_taxonomy_map.df)

asv_combined <- create_combined_dataframe(counts.df = asv.df,
                                          abundances.df = asv_rel.df,
                                          metadata.df = metadata.df,
                                          mylevel = "ASV_ID",
                                          sample_id_column = "Index",
                                          asv_taxonomy_map.df = asv_taxonomy_map.df)

write.table(asv_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/ASV_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(species_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Specie_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(genus_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Genus_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(family_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Family_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(order_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Order_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(class_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Class_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)
write.table(phylum_combined, file = "Result_tables/combined_counts_abundances_and_metadata_tables/Phylum_counts_abundances_and_metadata.csv", sep = ",", quote = F, col.names = T, row.names = F)


