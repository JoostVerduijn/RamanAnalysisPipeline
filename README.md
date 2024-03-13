# Raman Analysis Pipeline (RAP)

## Overview
This repository sets out the skeleton of an organizational structure used for scientific research. It loosely follows what I have used for several of my research projects and I hope it inspires you to conduct your research in an open, reproducible, and honest manner.

## Layout

The repository is split into seven main directories, many of which have subdirectories. This structure has been designed to be easily navigable by humans and computers alike, allowing for rapid location of specific files and instructions. Within each directory is a `README.md` file which summarizes the purpose of that directory as well as some examples where necessary. This structure may not be perfect for your intended us and may need to be modified. Each section is briefly described below. 

### **`Data Raman scans`** 
Folder containing all raw data from Raman scans as a .spc file. Each folder is a biological replicate, which each .spc being a different cell. Size of a Large Area Scan (LAS) are given in the file names name.

## **`Analysis in R`**

In this folder some of the scripts essential for data analysis are deposited.

* "Raman_Analysis_Joost Verduijn_Create database.R" Pulls all .spc files from a folder. Adds them together to one HyperSpec variable "Cell.raw". This is normalized and baseline corrected in "Cell.normalised". This is clustered using k-means clustering and filtered to select only cell specific clusters (normally in the middle of the scan). If the cell specific cluster is <20% or >85% cells are removed. End result is "Cell.normalised_filter" and is saved accordingly
* 2) Raman_Analysis_Joost Verduijn_Create SVM_Balanced.R creates SVM's on different Holdout and Training sets. Making the selection of experiments less important since this is averaged. 
* 3) Raman_Analysis_Joost Verduijn_Peakfit_Gaussians_AtPresetCenters.R. Fits multiple gaussian peaks are given wavenumbers.
* 4) Raman_Analysis_Joost Verduijn_Make_Boxplots_on_Peakfitted data.R. Makes boxplots from the area under the curve (AUC) of fitted peaks. Including statistics

