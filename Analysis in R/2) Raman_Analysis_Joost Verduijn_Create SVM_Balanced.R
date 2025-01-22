############################################################################################
#### Composed by Joost Verduijn 01/22/2025                                              ####
#### Galluzzi Lab                                                                       ####
#### Cancer Signaling and Microenvironment Program, FCCC, Philadelphia, PA, USA         ####
#### Johannes.Verduijn@fccc.edu                                                         ####
############################################################################################

##### LOADING #####
rm(list=ls())
library(caret)
library(writexl)
library(e1071)
library("hyperSpec")
library(rminer)

assign('All_Raman_data', get(load("%/_GitHub/Joost Verduijn/databases/full_dataset_XX-XX-XX-XXXXXX_full_clustered2.R"))) #From 1) Raman_Analysis_Joost Verduijn_Create database.R

RCDs <- unique(All_Raman_data$RCD)
##### BALANCED TS #####
# chooses n experiments for each RCD
n = 4
week.label = "Joost Verduijn"
date.folder<-format(Sys.time(), "%y-%m-%d")
dir.create(paste0("%/_GitHub/",week.label,"/SVM/"),showWarnings= FALSE)
dir.create(paste0("%/_GitHub/",week.label,"/SVM/",date.folder),showWarnings= FALSE)
setwd(paste0("%/_GitHub/",week.label,"/SVM/",date.folder))
all_stats = NULL
for (k in 1:100) { #100 SVM's creation
  idx = NULL
  for (i in 1:length(RCDs)) { #Select n=4 experiments from each RCD type
    idx <- c(idx, sample(which(startsWith(as.character(unique(All_Raman_data$group)),
                                          as.character(RCDs[i]))), n))
  }

  TrainingSet<-subset(All_Raman_data, subset = group %in% unique(All_Raman_data$group)[-idx]) 
  HoldoutSet<-subset(All_Raman_data, subset = group %in% unique(All_Raman_data$group)[idx])
  idx3 = NULL
  for (i in 1:length(RCDs)) {
    idx3 <- c(idx3, sample(which(startsWith(as.character(unique(TrainingSet$group)),
                                          as.character(RCDs[i]))), n))
  }
  TrainingSet<-subset(TrainingSet, subset = group %in% unique(TrainingSet$group)[idx3])
  
  ##### balance # of spectra in training set ######
  minSpc = min(table(TrainingSet$RCD))
  minSpc = min(minSpc, floor(length(TrainingSet$RCD)*0.05)) #Take 5% of dataset for each group more is not ness.
  idx2 = NULL
  for (i in 1:length(RCDs)) {
    idx2 <- c(idx2, sample(which(startsWith(as.character(TrainingSet$group),
                                          as.character(RCDs[i]))),minSpc))
  }
  TrainingSet<-TrainingSet[idx2]
  
  
  ##### SVM #####
  svm.model<-best.svm( x=TrainingSet$spc, y=TrainingSet$RCD, kernel='radial', cost=1,probability = TRUE)
  
  svm.model
  svm.pred <- predict(svm.model, HoldoutSet,type = "class")
  cm <- confusionMatrix(svm.pred,HoldoutSet$RCD)
  cm
    
  ##### RECORD STATS IN VECTOR #####
  date.label = format(Sys.time(), "%y-%m-%d-%H%M%S")
  st <- c(cm$overall[1], cm$byClass[,1],as.double(format(Sys.time(), "%H%M%S")),minSpc) #Converts date label into double so the rest is not changed to character, a tibble would be better.
  all_stats = rbind(all_stats, st)
  
  ##### SAVE DATA #####
  #date.label = format(Sys.time(), "%y-%m-%d-%H%M%S")
  write_xlsx(list("Stats" = data.frame(cm$byClass)[c(1,2,6)],
                  "Confusion" = as.data.frame(matrix(as.numeric(cm$table), nrow=4)),
                  "Training" = data.frame(unique(TrainingSet$group)),
                  "N spectra" = data.frame(minSpc),
                  "Holdout" = data.frame(unique(HoldoutSet$group))),
             paste0("confusionMatrix_", date.label,".xlsx"))
  
  save(HoldoutSet, file=paste0("HoldoutSet_", date.label, ".R"))
  save(TrainingSet, file=paste0("TrainingSet_", date.label, ".R"))
  save(svm.model, file=paste0("SVM_", date.label, ".R"))
}

write_xlsx(data.frame(all_stats),
           paste0("overall_stats_", date.label,".xlsx"))

 ### Plot results

HoldoutSet$prediction<-predict(svm.model, HoldoutSet,type = "class")
