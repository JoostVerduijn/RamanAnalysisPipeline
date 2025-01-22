############################################################################################
#### Composed by Joost Verduijn 01/22/2025                                              ####
#### Galluzzi Lab                                                                       ####
#### Cancer Signaling and Microenvironment Program, FCCC, Philadelphia, PA, USA         ####
#### Johannes.Verduijn@fccc.edu                                                         ####
############################################################################################
library(minpack.lm)             # Load the minpack.lm package
library(nls2)
library(hyperSpec)
library(baseline)
rm(list=ls())  
#Load dataset
setwd("%Filelocation%")
df<-readRDS("PeakFittedData_R.Rds") #Created by "3) Raman_Analysis_Joost Verduijn_Peakfit_Gaussians_AtPresetCenters.R"

library(dplyr)
presetCenters<-c(717,777,820,851,885,938,1000,1040,1085,1120,1210,1250,1300,1333,1445,1569,1610,1657) #plusminus 10
tolerance <- 15
grouped_data <- df %>%
  group_by(grp = cut(Center_restricted, breaks = sort(c(presetCenters - tolerance, (presetCenters) + tolerance))))

grouped_data$RCD_type<-sub("_.*","",dirname(grouped_data$Filepath)) #RCD not set proper
grouped_data<-grouped_data[complete.cases(grouped_data), ]#clean out NA data
pred_order<-c("Control","Apoptosis", "Ferroptosis","Necroptosis")
grouped_data$RCD_type<-factor(grouped_data$RCD_type, levels = pred_order)
###### BOXPLOTS SAVE IMAGES
library(ggplot2)
library(patchwork)
library(ggsignif)
# Create a list to store the plots
plots <- list()

group.colors=c("Apoptosis"="#729161", "Control"="#f89948","Ferroptosis"="#e18181","Necroptosis"="#82a2bb")
# Iterate over each group in grouped_data
for (grp_val in unique(grouped_data$grp)) {
  # Subset the data for the current group
  subset_data <- filter(grouped_data, grp == grp_val)
  
  # Create the boxplot for the current group
  plot <- ggplot(subset_data, aes(x = RCD_type, y = AUC_restricted,fill=RCD_type)) +
    geom_violin(outlier.colour="red") +
    scale_fill_manual(values=group.colors)+
    labs(title = paste("Group:", round(mean(subset_data$Center_restricted),2)),
         x = "RCD_type",
         y = "AUC_restricted")+
    theme(legend.position = "none")+ 
    ylim(min(subset_data$AUC_restricted),max(subset_data$AUC_restricted)*1.4)+
    geom_signif(y_position=max(subset_data$AUC_restricted)*1.0, comparisons = list(c(levels(subset_data$RCD_type)[1],levels(subset_data$RCD_type)[2])),map_signif_level=TRUE)+
    geom_signif(y_position=max(subset_data$AUC_restricted)*1.3, comparisons = list(c(levels(subset_data$RCD_type)[1],levels(subset_data$RCD_type)[3])),map_signif_level=TRUE)+
    geom_signif(y_position=max(subset_data$AUC_restricted)*1.2, comparisons = list(c(levels(subset_data$RCD_type)[1],levels(subset_data$RCD_type)[4])),map_signif_level=TRUE)+
    geom_signif(y_position=max(subset_data$AUC_restricted)*1.0, comparisons = list(c(levels(subset_data$RCD_type)[2],levels(subset_data$RCD_type)[3])),map_signif_level=TRUE)+
    geom_signif(y_position=max(subset_data$AUC_restricted)*1.1, comparisons = list(c(levels(subset_data$RCD_type)[2],levels(subset_data$RCD_type)[4])),map_signif_level=TRUE)+
    geom_signif(y_position=max(subset_data$AUC_restricted)*1.0, comparisons = list(c(levels(subset_data$RCD_type)[3],levels(subset_data$RCD_type)[4])),map_signif_level=TRUE)
  # Add the plot to the list
  plots[[as.character(grp_val)]] <- plot
}

# Save the plots as separate files
for (grp_val in names(plots)) {
  plot <- plots[[grp_val]]
  ggsave(filename = paste0("boxplot_", grp_val, ".png"), plot, width = 6, height = 4)
}

# Combine the plots into a single image using patchwork (optional)
combined_plot <- wrap_plots(plots) #, ncol=3

# Save the combined plot as a file (optional)
ggsave(filename = "combined_boxplots.png", combined_plot, width = 12, height = 8)







