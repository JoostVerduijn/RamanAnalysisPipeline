############################################################################################
#### Composed by Joost Verduijn 15/07/2021                                              ####
#### Laboratory of Nano-Biotechnology, Department of Biotechnology, Ghent University    ####
#### Johannes.Verduijn@ugent.be                                                         ####
############################################################################################

# this script creates a dataset starting from .spc files, performs baseline correction,
# normalisation and double k-means clustering. It also plots the heatmaps, the 
# results of clustering and the spectra of selected groups of data.


# remove old dataset
rm(list=ls())

# set folder where to save data
week.label = "Joost"



##### LOAD PACKAGES #####
library(vegan)
library("stringr")
library("hyperSpec") # hyperSpec package
library("baseline") # Contains all kinds of functions for baseline removal
#library(MALDIquant)
library(gridExtra )
library(hyperSpec.utils)
library(readxl)
library(writexl)



##### LOAD FILES #####
setwd("C:/Users/Joost Verduijn/OneDrive - UGent/_ Research/RAMAN/Raman data Alive Cells LAURA")
files=list.files(pattern=".spc",recursive = TRUE) #Loads all .spc files from subfolders


##### REMOVE UNWANTED EXPERIMENTS #####
# 
 files <- files[-grep("STS", files)]               # STS Apoptosis experiments
 files <- files[-grep("Control_Laura_5", files)]   # suspicious control experiment
 files <- files[-grep("Control_Laura_3", files)] 
 files <- files[-grep("Control_Eugenia_14/cell03-control_44x44.spc", files)] #spectra at x=32,33 and y=11,11 show weird behavior
 files <- files[-grep("Necroptosis_Eugenia", files)]
 files <- files[-grep("Necroptosis_Laura_mTNF", files)]
# 



##### SET LOADING FUNCTION AND PARAMETERS #####
xList=sub("x.*", "",sub(".spc","",sub(".*_","",files))) #get x size from file name "File-name_20x25.spc" gives 20
yList=sub(".*x", "",sub(".spc","",sub(".*_","",files))) #get y size from file name "File-name_20x25.spc" gives 25

#### Choose normalisation type
# "SNV" -> snv-normalisation
# "AREA" -> area-normalisation
# "BAND" -> band-normalisation (over phenylalanine peak)
norm.type = "SNV"

#### Function to analyse each cell, with f=file, x=xSize, y=ySize, ntype=type of normalisation
AnalyzeSingleCell <- function(f, x, y, ntype){
  hs=read.spc(f)
  #select the area of interest. It can be either 500-1800 or 2500-3200
  hs=hs[,,500~1800]
  wavelengths <- hs@wavelength #get wavelengths
  if (nrow(hs)!=as.numeric(x)*as.numeric(y)){print("ERROR in size")}
  else{
    ppl <- x
    lines <- y
    pxl.Âµm <- 1
    XY <- data.frame(cbind(rep((1:ppl)*1/pxl.Âµm,lines),rep((1:lines)*1/pxl.Âµm,each=ppl))) #Makes dataframe giving X and Y coordinates to each spectra
    colnames(XY) <- c("X","Y")
    hs$x <- XY$X
    hs$y <- sort(XY$Y, decreasing=TRUE)
    labels(hs,"x") <- as.character("x / µm")
    labels(hs,"y") <- as.character("y / µm")
    
    # Set group (experiment type) and specific cell (in file_spc)
    hs$group<-as.factor(dirname(hs$filename)) #Experiment name
    hs$file_spc<-as.factor(basename(hs$filename)) #File name
    hs$RCD<-as.factor(sub("_.*","",hs$group)) #RCD type name (since folder name is "Ferroptosis_Laura 1", gives Ferroptosis)
    
    
    hs.raw = hs #hs.raw is only used in export (so is obsolete)
    # ### Create baseline per LAS, polynomial fit with 4th degree gives best background reduction
    # b<-baseline(hs$spc,method="modpolyfit",degree=4)#Old =method="als",lambda = 4, p = 0.001, maxit = 50) 
    # #plot(b)
    # hs$spc=getCorrected(b)
    # # hs$spc[hs$spc<0]=0 #After baseline corrections some spc values might be below zero --> set to zero
    # hs.base<-hs
    # # NORMALISATION
    # if (ntype == "SNV") {
    #   hs.norm<-snv_normalization(hs)
    # } else if (ntype == "AREA") {
    #   hs.norm<-area_normalization(hs)
    # } else if (ntype == "BAND") {
    #   bnd = 999~1009 # phenylalanine peak @ 1004 +- 5
    #   hs.norm<-band_normalization(hs, bnd)
    # } 
    # hs.norm*1000000
    # #outputs the Hyperspec(HS), the baseline corrected HS and the normalised HS
    # #Add the new factors to the Global environment
    # newList <- list("hs.base" = hs.base, "hs.norm" = hs.norm,"hs.raw"=hs.raw)
    newList <- list("hs.raw"=hs.raw)
    list2env(newList ,.GlobalEnv)
  }
}





##### USE FUNCTION FOR ALL FILES #####
#Create an empty HysperSpec object to append to
EmptyHs=read.spc(files[1])
EmptyHs=EmptyHs[,,500~1800]                   # choose wavenumber range
Cell.normalised<-hyperSpec::empty(EmptyHs)
Cell.base<-hyperSpec::empty(EmptyHs)
Cell.raw<-hyperSpec::empty(EmptyHs)

#For all files, run function, and add them together (rbind)
for (i in 1:length(files)){ #1:length(files)
  AnalyzeSingleCell(files[i],xList[i], yList[i], norm.type) # load, baseline correct, normalise
  #Cell.normalised<-rbind2(Cell.normalised, hs.norm, wl.tolerance=0.3 ) #0.3
  #Cell.base<-rbind2(Cell.base, hs.base, wl.tolerance=0.3 ) #0.3
  Cell.raw<-rbind2(Cell.raw, hs.raw, wl.tolerance=0.3 ) #0.3
  print(files[i])
  #print(nrow(Cell.normalised))
}
##BASELINE
BaselineCorrected<-Cell.raw
BaselineCorrected$spc<-getCorrected(baseline(BaselineCorrected$spc,method="modpolyfit",degree=4))
## CRR
CRRCorrected<-crr(BaselineCorrected,threshold = 10) #Correct using Hyperspec.Utils, threshold is optimized
#Cell.normalised<-Cell.normalised[-which(as.numeric(summary(Cell.normalised$crr)[,1])>0)] #Removes all spectra that have been edited by CRR
## NORMALIZATION
HS_Normalized<-snv_normalization(CRRCorrected)

#HS_band<-band_normalization(CRRCorrected, 999~1009)# phenylalanine peak @ 1004 +- 5


# #Test for cell specific cluster?
# plotmap(HS_Normalized[HS_Normalized$file_spc==unique(HS_Normalized$file_spc)[120]][,,1000~1010],spc~x*y, col.regions = topo.colors(100))
# plotmap(HS_Normalized[HS_Normalized$file_spc==unique(HS_Normalized$file_spc)[120]][,,1000~1010][rowMeans(HS_Normalized[HS_Normalized$file_spc==unique(HS_Normalized$file_spc)[120]][,,1000~1010]$spc)>-0.1],spc~x*y, col.regions = topo.colors(100))
# 
# plotmap(HS_Normalized[HS_Normalized$file_spc==unique(HS_Normalized$file_spc)[120]],cluster_kmeans~x*y, col.regions = topo.colors(100))
# # Set group (experiment type) and specific cell (in file_spc)
# Cell.normalised$group<-as.factor(dirname(Cell.normalised$filename)) #Experiment name
# Cell.normalised$file_spc<-as.factor(basename(Cell.normalised$filename)) #File name
# Cell.normalised$RCD<-as.factor(sub("_.*","",Cell.normalised$group)) #RCD type name (since folder name is "Ferroptosis_Laura 1", gives Ferroptosis)



Cell.normalised<-HS_Normalized
# Amount of cells per each RCD type
table(as.factor(sub("_.*","",unique(Cell.normalised$filename))))






# ##### COSMIC RAY REMOVAL #####
# tmp3 <- Cell.normalised
# # Cell.normalised <- tmp
# 
# ### SET THRESHOLD CRR
# Cell.normalised<-crr(Cell.normalised,threshold = 10) #Correct using Hyperspec.Utils, threshold is optimized
# Cell.normalised<-Cell.normalised[-which(as.numeric(summary(Cell.normalised$crr)[,1])>0)] #Removes all spectra that have been edited by CRR
# #Cell.normalised[-row(Cell.normalised$spc)[which(is.na(Cell.normalised$spc))],]
# #Have to look into this, instead of removing CR it levels them to their respectable max.
# #which(as.numeric(summary(Cell.normalised$crr)[,1])>0) #the spectra with an correction by CRR
# #which(is.na(Cell.normalised$spc)) #should be shorter but less spectra??
# #Cell.normalised<-Cell.normalised[-row(Cell.normalised$spc)[which(is.na(Cell.normalised$spc)))],]
# 
# #Cell.normalised<-Cell.normalised[-row(Cell.normalised$spc)[which(Cell.normalised$spc>7)],]
# 
# 
# 


##### SAVE NORMALISED DATASET #####
date.label = format(Sys.time(), "%y-%m-%d-%H%M%S")   # unique label with date and time
# folder where to save data

dir.create(paste0("C:/Users/Joost Verduijn/OneDrive - UGent/_ Research/RAMAN/",week.label,"/databases"),showWarnings= FALSE)
setwd(paste0("C:/Users/Joost Verduijn/OneDrive - UGent/_ Research/RAMAN/",week.label,"/databases"))
save(Cell.normalised, file = paste0("full_dataset_", date.label, ".R"))





##### CLUSTERING FOR WHOLE DATASET #####
nclu = 3   # number of clusters
# k-means clustering
Cell.normalised$cluster_kmeans<-as.factor(kmeans(Cell.normalised$spc,
                                                 nclu,
                                                 iter.max = 100,
                                                 nstart = 10)$cluster)





##### SAVE DATASET WITH CLUSTERS #####
save(Cell.normalised, file = paste0("full_dataset_", date.label, "_full_clustered.R"))





##### CHOOSE CELL CLUSTER #####
#Select each cell only the inner (30-70% both x and y) part to find the correct cluster
clusterlist_total<-factor()
for (i in 1:length(unique(Cell.normalised$filename))){
  
  xmax=max(Cell.normalised$x[Cell.normalised$filename==unique(Cell.normalised$filename)[i]])
  ymax=max(Cell.normalised$y[Cell.normalised$filename==unique(Cell.normalised$filename)[i]])
  scalefactor=0.7 #Select inner 70%
  #Extract K-means cluster for in inner level (0.3*xmax till 0.7xmax) same for y
  #This line looks insane but it is not, so select all cluster_kmeans values for the pixels with x AND Y between 0.3-0.7 of its max
  clusterlist<-Cell.normalised$cluster_kmeans[Cell.normalised$filename==unique(Cell.normalised$filename)[i]&Cell.normalised$x<ceiling(xmax*scalefactor)&Cell.normalised$x>floor(xmax*(1-scalefactor))&Cell.normalised$y<ceiling(ymax*scalefactor)&Cell.normalised$y>floor(ymax*(1-scalefactor))]
  clusterlist_total<-c(clusterlist_total,clusterlist) #c of combine
  #print(unique(Cell.normalised$filename)[i])
}

# choose biggest cluster overall
max_cluster_kmeans <- order(-tabulate(clusterlist_total))  # Change to which.max(tabulate( ))?





##### REMOVE CELLS WHOSE BIGGEST CLUSTER IS < 20% or > 85% #####
tmp2 <- Cell.normalised
# set low and high thresholds
lowTh = 0.20
highTh = 0.85
perc = NULL
delIdx = NULL
for (j in 1:length(unique(Cell.normalised$filename))) {
  hs <- Cell.normalised[Cell.normalised$filename==unique(Cell.normalised$filename)[j]]
  #nmax = table(hs$cluster_kmeans == max_cluster_kmeans[1])
  #perc1 = nmax["TRUE"]/sum(nmax)
  perc1<-sum(hs$cluster_kmeans == max_cluster_kmeans[1])/length(hs$cluster_kmeans)
  perc = c(perc, perc1)
# if ((is.na(perc1)) | (perc1<lowTh) | (perc1>highTh)) {
#    delIdx = c(delIdx, j)
#  }
}
perc
delIdx1<-is.na(perc) | (perc<lowTh) | (perc>highTh)

length(unique(Cell.normalised$filename))
sum(delIdx1)
Cell.normalised2 <- subset(Cell.normalised, subset = filename %in% unique(Cell.normalised$filename)[-which(delIdx1)]) #should be indexes
Cell.normalised2$file_spc<-droplevels(Cell.normalised2$file_spc) # drop unused levels for clarity
length(unique(Cell.normalised2$filename))

##### REPEAT CLUSTERING WITH NEW DATASET #####
Cell.normalised2$cluster_kmeans2<-as.factor(kmeans(Cell.normalised2$spc,
                                                 nclu,
                                                 iter.max = 100,
                                                 nstart = 10)$cluster)

### SAVE DATASET WITH CLUSTERS
save(Cell.normalised2, file = paste0("full_dataset_", date.label, "_full_clustered2.R"))

### CHOOSE CELL CLUSTER
#Select each cell only the inner (30-70% both x and y) part to find the correct cluster
clusterlist_total<-factor()
for (i in 1:length(unique(Cell.normalised2$filename))){
  
  xmax=max(Cell.normalised2$x[Cell.normalised2$filename==unique(Cell.normalised2$filename)[i]])
  ymax=max(Cell.normalised2$y[Cell.normalised2$filename==unique(Cell.normalised2$filename)[i]])
  scalefactor=0.7 #Select inner 70%
  #Extract K-means cluster for in inner level (0.3*xmax till 0.7xmax) same for y
  #This line looks insane but it is not, so select all cluster_kmeans2 values for the pixels with x AND Y between 0.3-0.7 of its max
  clusterlist<-Cell.normalised2$cluster_kmeans2[Cell.normalised2$filename==unique(Cell.normalised2$filename)[i]&Cell.normalised2$x<ceiling(xmax*scalefactor)&Cell.normalised2$x>floor(xmax*(1-scalefactor))&Cell.normalised2$y<ceiling(ymax*scalefactor)&Cell.normalised2$y>floor(ymax*(1-scalefactor))]
  clusterlist_total<-c(clusterlist_total,clusterlist) #c of combine
  #print(unique(Cell.normalised$filename)[i])
}

# choose biggest cluster overall
max_cluster_kmeans2 <- order(-tabulate(clusterlist_total)) 


##### SELECT CLUSTERED DATASET #####
#Select biggest cluster(s), in number of pixels. Selects 1 cluster if k=3, 2 clusters if k>3
Cell.normalised_filter<-Cell.normalised2[Cell.normalised2$cluster_kmeans2==max_cluster_kmeans2[1]] #Select all clusters equal to max central cluster
Cell.normalised_background<-Cell.normalised2[!Cell.normalised2$cluster_kmeans2==max_cluster_kmeans2[1]] #Select all clusters equal THAT IS NOT max central cluster
unique(Cell.normalised_filter$file_spc)

##### SAVE FILTERED DATASET #####
save(Cell.normalised_filter, file = paste0("full_dataset_", date.label, "_filtered.R"))

# for (z in 1:10){     
#   #win.graph()
#   print(unique(Cell.normalised$file_spc)[z])
#   
#   
#   p1<-plotmap(Cell.raw[Cell.raw$file_spc==unique(Cell.raw$file_spc)[z]],spc~x*y, col.regions = topo.colors(100))
#   p2<-plotmap(Cell.base[Cell.base$file_spc==unique(Cell.base$file_spc)[z]],spc~x*y, col.regions = topo.colors(100))
#   p3<-plotmap(tmp3[tmp3$file_spc==unique(tmp3$file_spc)[z]],spc~x*y, col.regions = topo.colors(100))
#   p4<-plotmap(Cell.normalised[Cell.normalised$file_spc==unique(Cell.normalised$file_spc)[z]],spc~x*y, col.regions = topo.colors(100))
#   p5<-plotmap(Cell.normalised_filter[Cell.normalised_filter$file_spc==unique(Cell.normalised_filter$file_spc)[z]],spc~x*y, col.regions = topo.colors(100))
#   
#   
#   q1<-qplotspc(Cell.raw[Cell.raw$file_spc==unique(Cell.raw$file_spc)[z]],spc.nmax = 10000)
#   q2<-qplotspc(Cell.base[Cell.base$file_spc==unique(Cell.base$file_spc)[z]],spc.nmax = 10000)
#   q3<-qplotspc(tmp3[tmp3$file_spc==unique(tmp3$file_spc)[z]],spc.nmax = 10000)
#   q4<-qplotspc(Cell.normalised[Cell.normalised$file_spc==unique(Cell.normalised$file_spc)[z]],spc.nmax = 10000)
#   q5<-qplotspc(Cell.normalised_filter[Cell.normalised_filter$file_spc==unique(Cell.normalised_filter$file_spc)[z]],spc.nmax = 10000)
#   
#   qplotspc(Cell.raw[Cell.raw$file_spc==unique(Cell.raw$file_spc)[z]],spc.nmax = 10000)
#   #p1<-plotmap(Cell.HCA_total[Cell.HCA_total$filename==unique(Cell.HCA_total$filename)[i]], 
#   #            clusters~x*y,col.regions=topo.colors(50))
#   grid.arrange(p1,q1,p2,q2,p3,q3,p4,q4,p5,q5, ncol = 2, top =unique(Cell.normalised$filename)[i] )
# }



#### Outlier Removal ####
RCDs <- unique(Cell.normalised_filter$RCD)
library('mt')

Delnames<-c()
temp<-hyperSpec::empty(Cell.normalised_filter)
for (x in 1:length(RCDs)){
  temp<-Cell.normalised_filter[(Cell.normalised_filter$RCD %in% RCDs[x])] #select only one RCD at the time
  temp$group<-droplevels(temp$group)     # drop levels for clarity
  unique(temp$group)
  # takes median per cell
  temp<-aggregate(temp,temp$file_spc,median) #Median per cell
  
  
  DelIdx2<-pca.outlier(temp$spc)$outlier #Give index number of outliers
  
  Delnames<-c(Delnames,temp[DelIdx2]$filename) #Gives filename (used for delete later)
}
Raman_Removed_outliers<-Cell.normalised_filter[-which(Cell.normalised_filter$filename %in% Delnames)]

save(Raman_Removed_outliers, file = paste0("full_dataset_", date.label, "_filtered_Removed_Outliers.R"))



 # ##### CLUSTERING PLOT ----- TO CHECK #####
# # load function
# 
# 
# source("C:/Users/Joost Verduijn/OneDrive - UGent/_ Research/RAMAN/final/general_plot_EDR.R")
# 
# # plot and save clustering for all cells
# FigureFolder<-paste0("C:/Users/Joost Verduijn/OneDrive - UGent/_ Research/RAMAN/",week.label,"/Figures/")
# png(paste0(FigureFolder,"full_dataset_", date.label, ".png"), width = 2500, height = 2850)
# generalPlot(Cell.normalised, Cell.normalised$cluster_kmeans)
# dev.off()
# 
# # plot and save second clustering for all cells
# png(paste0(FigureFolder,"full_dataset_2ND_CLU", date.label, ".png"), width = 2500, height = 2850)
# generalPlot(Cell.normalised2, Cell.normalised2$cluster_kmeans2)
# dev.off()
# 
# ##### COMPARE HEAT MAPS ----- to check #####
# # load function
# source("C:/Users/edremigi/OneDrive - UGent/RAMAN/HeatMaps_EDR.R")
# 
# # select parameters
# # cnr = 59 # cell number
# wn = 1150 # central wavenumber
# wnr = 650  # range (doubled)
# 
# # whole data
# png(paste0("heatmaps_", norm.type, "_k", nclu, "_wn", wn, "_", date.label, ".png"), width = 1800, height = 2850)
# heatMap(Cell.normalised, wn, wnr, norm.type)
# dev.off()
# 
# # only cell cluster
# png(paste0("heatmaps_", norm.type, "_k", nclu, "_wn", wn, "_", date.label, "_full.png"), width = 1800, height = 2850)
# heatMap(Cell.normalised_filter, wn, wnr, norm.type)
# dev.off()
# 
# # #### compare to raw dataset
# # # load raw dataset
# # # load("C:/Users/edremigi/OneDrive - UGent/RAMAN/raw_dataset_bc.R")
# # 
# # png(paste0("heatmaps_raw_k_", norm.type, "_k", nclu, "_wn", wn, "_", date.label, ".png"), width = 1800, height = 2850)
# # heatMap(Cell.HCA_total, wn, wnr, "raw")
# # dev.off()
# 
# ##### SPC PLOTS #####
# # for more spc plots see "useful functions" folder
# 
# # calculate median per experiment
# hs <- aggregate(Cell.normalised_background, Cell.normalised_background$group, median)
# 
# # select a certain group (e.g. RCD)
# RCDi = "Apoptosis"
# hs <-subset(hs, subset = RCD %in% unique(hs$RCD)[which(startsWith(as.character(unique(hs$RCD)),RCDi))])
# 
# # plot selected group
# qplotspc(hs, spc.nmax = 100000)+
#   geom_line(mapping = aes(x = .wavelength, y = spc, colour = group))
#   # theme(legend.position="none")
# ggsave(filename=paste0("SPC_median_", RCD, ".png"), plot=last_plot(), width = 3600, height = 2000, units = "px")

