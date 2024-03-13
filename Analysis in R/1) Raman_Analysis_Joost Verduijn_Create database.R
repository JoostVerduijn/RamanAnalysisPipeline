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
week.label = "Joost Verduijn"



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
setwd("%YOURFOLDER%")
files=list.files(pattern=".spc",recursive = TRUE) #Loads all .spc files from subfolders


##### SET LOADING FUNCTION AND PARAMETERS #####
xList=sub("x.*", "",sub(".spc","",sub(".*_","",files))) #get x size from file name "File-name_20x25.spc" gives 20
yList=sub(".*x", "",sub(".spc","",sub(".*_","",files))) #get y size from file name "File-name_20x25.spc" gives 25

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
       
    hs.raw = hs 
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
  AnalyzeSingleCell(files[i],xList[i], yList[i], norm.type) # load, baseline correct, normalize
  Cell.raw<-rbind2(Cell.raw, hs.raw, wl.tolerance=0.3 ) #0.3
  print(files[i])
 }
##BASELINE
BaselineCorrected<-Cell.raw
BaselineCorrected$spc<-getCorrected(baseline(BaselineCorrected$spc,method="modpolyfit",degree=4))
## Cosmic ray removal
CRRCorrected<-crr(BaselineCorrected,threshold = 10) #Correct using Hyperspec.Utils, threshold is optimized

## NORMALIZATION
HS_Normalized<-snv_normalization(CRRCorrected)

Cell.normalised<-HS_Normalized
# Amount of cells per each RCD type
table(as.factor(sub("_.*","",unique(Cell.normalised$filename))))

##### SAVE NORMALISED DATASET #####
date.label = format(Sys.time(), "%y-%m-%d-%H%M%S")   # unique label with date and time
# folder where to save data

dir.create(paste0("%YOURFOLDER%/",week.label,"/databases"),showWarnings= FALSE)
setwd(paste0("%YOURFOLDER%/",week.label,"/databases"))
save(Cell.normalised, file = paste0("full_dataset_", date.label, ".R"))





##### CLUSTERING FOR WHOLE DATASET #####
nclu = 3   # number of clusters
# k-means clustering
Cell.normalised$cluster_kmeans<-as.factor(kmeans(Cell.normalised$spc,nclu,iter.max = 100,nstart = 10)$cluster)

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
