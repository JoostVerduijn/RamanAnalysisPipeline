## `Analysis in R

In this folder some of the scripts essential for data analysis are deposited.

* "Raman_Analysis_Joost Verduijn_Create database.R" Pulls all .spc files from a folder. Adds them together to one HyperSpec variable "Cell.raw". This is normalized and baseline corrected in "Cell.normalised". This is clustered using k-means clustering and filtered to select only cell specific clusters (normally in the middle of the scan). If the cell specific cluster is <20% or >85% cells are removed. End result is "Cell.normalised_filter" and is saved accordingly
* 2) Raman_Analysis_Joost Verduijn_Create SVM_Balanced.R creates SVM's on different Holdout and Training sets. Making the selection of experiments less important since this is averaged. 
* 3) Raman_Analysis_Joost Verduijn_Peakfit_Gaussians_AtPresetCenters.R. Fits multiple gaussian peaks are given wavenumbers.
* 4) Raman_Analysis_Joost Verduijn_Make_Boxplots_on_Peakfitted data.R. Makes boxplots from the area under the curve (AUC) of fitted peaks. Including statistics

