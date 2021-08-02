library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(lemon)

common_theme <- theme(
  panel.border = element_blank(), 
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
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
  legend.key.height= unit(.3,"cm"),
  legend.key.width = unit(.3,"cm"),
  legend.text = element_text(size = 8),
  axis.line = element_line(colour = "black", size = 0.5),
  axis.text = element_text(size = 6, colour = "black"),
  axis.title = element_text(size = 7,face = "bold"),
  complete = F,
  plot.title = element_text(size = 8))


source("code/utility.R")


# Load the processed metadata that is unfiltered. This is for the qPCR results.
metadata_unfiltered.df <- read.csv("Result_tables/processed_unfiltered_metadata.csv", sep =",", header = T)
rownames(metadata_unfiltered.df) <- metadata_unfiltered.df$Index
metadata_unfiltered.df <- metadata_unfiltered.df[is.na(metadata_unfiltered.df$qPCR_exclude),]
metadata_unfiltered.df <- metadata_unfiltered.df[,c("Index","Subject" ,"Sample_type","Swab_ID", "Cohort","Sample_type_colour","Sample_type_shape","S_aureus_qPCR","Staph_spp_qPCR", "qPCR_16S")]

# Remove negative samples
metadata_unfiltered.df <- metadata_unfiltered.df[! metadata_unfiltered.df$Sample_type == "negative",]

# Set sample type colours and shapes
sample_type_colours <- unique(metadata_unfiltered.df[,c("Sample_type", "Sample_type_colour")])
sample_type_colours <- setNames(as.character(sample_type_colours[,"Sample_type_colour"]), sample_type_colours[,"Sample_type"])
sample_type_shapes <- unique(metadata_unfiltered.df[,c("Sample_type", "Sample_type_shape")])
sample_type_shapes <- setNames(as.numeric(sample_type_shapes[,"Sample_type_shape"]), sample_type_shapes[,"Sample_type"])

sample_type_outline_colours <- unique(metadata_unfiltered.df[,c("Sample_type", "Sample_type_colour")])
sample_type_outline_colours <- setNames(darken(as.character(sample_type_outline_colours[,"Sample_type_colour"]),factor = 2), 
                                        sample_type_outline_colours[,"Sample_type"])

# Sort by subject and qPCR_16S value
metadata_unfiltered.df <- metadata_unfiltered.df %>% 
  dplyr::arrange_at(c("Subject","qPCR_16S"))

# Factorise sample type
metadata_unfiltered.df$Sample_type <- factor(metadata_unfiltered.df$Sample_type, 
                                             levels = c("NS", "PDS", "AK", "SCC_PL","SCC"))

metadata_unfiltered.df %>% 
  group_by(Cohort, Subject, Sample_type) %>%
  dplyr::summarise(Mean_qPCR_16S = mean(qPCR_16S))

immunosuppresed_barplot <-
  ggplot(metadata_unfiltered.df %>% filter(Cohort == "organ transplant recipient"), 
       aes(x = Sample_type, y = qPCR_16S, fill = Sample_type, colour = Sample_type)) +
  geom_bar(position = position_dodge2(preserve = "single",padding = .2),  
           stat = "identity",
           # colour = "black",
           lwd = .1) +
  labs(title = "Organ transplant recipients") +
  xlab("Sample type") +
  ylab("SSU rRNA equivalents per sampling area") +
  scale_shape_manual(values = sample_type_shapes,name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, name = "Sample type") +
  scale_colour_manual(values = sample_type_outline_colours, name = "Sample type")+
  scale_y_log10(limits=c(10^-1, 10^5),
                breaks=10^(-1:5),
                labels = trans_format('log10',math_format(10^.x)),
                expand = c(0,0)) +
  coord_cartesian(ylim = c(10^0, 10^5))+
  lemon::facet_rep_wrap(~Subject, scales = "free_x", nrow = 7,repeat.tick.labels = F) +
  common_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5,size = 5),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 6),
        panel.spacing = unit(0, "cm")
        )


immunocompetent_barplot <-
  ggplot(metadata_unfiltered.df %>% filter(Cohort == "immunocompetent"), 
         aes(x = Sample_type, y = qPCR_16S, fill = Sample_type, colour = Sample_type)) +
  geom_bar(position = position_dodge2(preserve = "single",padding = .3),  
           stat = "identity",
           # colour = "black",
           lwd = .1) +
  labs(title = "Immunocompetent") +
  xlab("Sample type") +
  ylab("SSU rRNA equivalents per sampling area") +
  scale_shape_manual(values = sample_type_shapes,name = "Sample type") +
  scale_fill_manual(values = sample_type_colours, name = "Sample type") +
  scale_colour_manual(values = sample_type_outline_colours, name = "Sample type")+
  scale_y_log10(limits=c(10^-1, 10^5),
                breaks=10^(-1:5),
                labels = trans_format('log10',math_format(10^.x)),
                expand = c(0,0)) +
  coord_cartesian(ylim = c(10^0, 10^5))+
  lemon::facet_rep_wrap(~Subject, scales = "free_x", nrow = 5,repeat.tick.labels = F) +
  common_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5,size = 5),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 6),
        panel.spacing = unit(0, "cm")
        ) 

ggsave(filename = "Result_figures/various/immunosuppressed_microbial_loads_per_subject.pdf",
       plot = immunosuppresed_barplot,
       device = "pdf",
       width = 20,height = 20,units = "cm")
ggsave(filename = "Result_figures/various/immunocompetent_microbial_loads_per_subject.pdf",
       plot = immunocompetent_barplot,
       device = "pdf",
       width = 15,height = 18,units = "cm")



# ggplot(metadata_unfiltered.df %>% filter(Cohort == "organ transplant recipient"), 
#        aes(x = Sample_type, y = qPCR_16S, fill = Sample_type, colour = Sample_type, shape = Sample_type)) +
#   stat_summary(fun = "mean", geom = "bar",
#                shape = 16,size = 1,
#                position = position_dodge2(preserve = "single",padding = .2)) +
#                # position = position_dodge(width = .75),show.legend = F)+
#   
#   geom_jitter(size = .6,stroke =.1,
#               position = position_jitterdodge(jitter.width = .2,
#                                               dodge.width = .75)) +
#   labs(title = "Organ transplant recipients") +
#   xlab("Sample type") +
#   ylab("SSU rRNA equivalents per sampling area") +
#   scale_shape_manual(values = sample_type_shapes,name = "Sample type") +
#   scale_fill_manual(values = sample_type_colours, name = "Sample type") +
#   scale_colour_manual(values = sample_type_outline_colours, name = "Sample type")+
#   scale_y_log10(limits=c(10^-1, 10^5),
#                 breaks=10^(-1:5),
#                 labels = trans_format('log10',math_format(10^.x)),
#                 expand = c(0,0)) +
#   coord_cartesian(ylim = c(10^0, 10^5))+
#   lemon::facet_rep_wrap(~Subject, scales = "free_x", nrow = 7,repeat.tick.labels = F) +
#   common_theme +
#   theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5,size = 5),
#         strip.text = element_text(face = "bold", size = 8),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#         plot.subtitle = element_text(hjust = 0.5, size = 6),
#         panel.spacing = unit(0, "cm")
#   )
# graphics.off()

