
library(minpack.lm)
library(nls2)
library(alsace)

library(hyperSpec)
library(baseline)
rm(list=ls())
assign('All_Raman_data', get(load("C:/Users/Joost Verduijn/OneDrive - UGent/_ Research/RAMAN/Joost/databases/full_dataset_23-05-10-135126_filtered_Removed_Outliers.R")))
All_Raman_data$cellnumber<-as.numeric(gsub("\\D","",sub("_.*","",sub("-.*","",All_Raman_data$file_spc))))

gaussian <- function(x, amplitude, center, width) {
  #amplitude * exp(-(x - center)^2 / (2 * width^2)) #gauss?
  
  amplitude * exp(((-4) * log(2) * (x - center)^2) / (width ^ 2)) #gauss? with adjusted width (half max full width)
}

# Define the sum of Lorentzian functions
sum_gaussian <- function(x, params) {
  n_peaks <- length(params) / 3
  result <- rep(0, length(x))
  for (i in 1:n_peaks) {
    result <- result + gaussian(x, params[i], params[(n_peaks+i)], params[(n_peaks*2+i)])
  }
  return(result)
}
test<-All_Raman_data
general_ALS_baseline<-baseline::baseline(spectra = test$spc, method = "als", lambda = 4, p = 0.004)
base_test<-test #Copy metadata
base_test$spc<-getCorrected(general_ALS_baseline)

presetCenters<-c(717,777,820,851,885,938,1000,1040,1085,1120,1210,1250,1300,1333,1445,1569,1610,1657) #plusminus 10

df<-data.frame()
total_iterations <- 10
# Preallocate the params_df data frame
params_df <- data.frame(
  RCD_type = character(total_iterations * n_peaks), 
  Peaknumber = numeric(total_iterations * n_peaks), 
  Filepath = character(total_iterations * n_peaks),
  x_coordinate = numeric(total_iterations * n_peaks),
  y_coordinate = numeric(total_iterations * n_peaks),
  Cell_number = integer(total_iterations * n_peaks),
  Amplitude_restricted = numeric(total_iterations * n_peaks),
  Center_restricted = numeric(total_iterations * n_peaks),
  Width_restricted = numeric(total_iterations * n_peaks),
  AUC_restricted = numeric(total_iterations * n_peaks),
  Amplitude_init = numeric(total_iterations * n_peaks),
  Center_init = numeric(total_iterations * n_peaks),
  Width_init = numeric(total_iterations * n_peaks),
  row.names = character(total_iterations * n_peaks)
)
for (q in 1:length(base_test)){
  x<-base_test@wavelength[]
  y<-base_test$spc[q,]
  RCD_type<-base_test$RCD[q] 
  maxima<-data.frame(xmax=(base_test[,,presetCenters][q,]@wavelength),ymax=t(base_test[,,presetCenters][q,]$spc))
  fwhm_test<-numeric(nrow(maxima))
  for (qq in 1:nrow(maxima)) {
    half_max <- maxima$ymax[qq] / 2
    x_right<-min(x[y>half_max&x>maxima$xmax[qq]]) #right
    x_left<-max(x[y>half_max&x<maxima$xmax[qq]]) #left
    fwhm_test[qq] <- (2*min(maxima$xmax[qq]-x_left, x_right-maxima$xmax[qq]))
}
  width_peaks<-fwhm_test#[maxima$ymax>0.05]
  center_peaks<-maxima$xmax#[maxima$ymax>0.05]
  height_peaks<-maxima$ymax#[maxima$ymax>0.05]
  n_peaks <- length(center_peaks)  
  # Generate initial parameters
  init_params <- data.frame(
    amplitude = height_peaks,#runif(n_peaks, 0.5, 2),
    center = center_peaks,#runif(n_peaks, min(x), max(x)),
    width = width_peaks#runif(n_peaks, 15, 60)
  )
  init_params<-unlist(init_params)
  fit_data <- data.frame(x = x, y = y)
  fit_formula_gaussian <- as.formula(paste("y ~", paste(paste0("gaussian(x, amplitude", 1:n_peaks, ", center", 1:n_peaks, ", width", 1:n_peaks, ")"), collapse = " + ")))
  lower_limit <- unlist(data.frame(
    amplitude = 0.3*height_peaks,
    center = center_peaks-10,
    width = 0.8*width_peaks
  ))
  upper_limit <-  unlist(data.frame(
    amplitude = 1.1*height_peaks,
    center = center_peaks+10,
    width = width_peaks^2
  ))
  nlcLM<-nls.lm.control(maxiter = 1000,nprint = 0)
  fit_result_gaussian_restricted <-nlsLM(formula = fit_formula_gaussian, start = init_params,control=nlcLM, data = fit_data,upper = upper_limit,lower = lower_limit)
  fitted_params_gaussian_restricted <- coef(fit_result_gaussian_restricted)
  params_df <- data.frame(RCD_type=rep(RCD_type,n_peaks), #RCD identifier
                          Peaknumber= 1:n_peaks, 
                          Filepath=base_test[q,]$filename,
                          x_coordinate=base_test[q,]$x,
                          y_coordinate=base_test[q,]$y,
                          
                           #
                          Cell_number=rep(base_test$cellnumber[q] ,n_peaks),#indication of time
                          #Amplitude = fitted_params_gaussian[1:n_peaks],
                          #Center = fitted_params_gaussian[(n_peaks+1):(2*n_peaks)],
                          #Width = fitted_params_gaussian[(2*n_peaks+1):(3*n_peaks)],
                          
                          Amplitude_restricted = fitted_params_gaussian_restricted[1:n_peaks],
                          Center_restricted = fitted_params_gaussian_restricted[(n_peaks+1):(2*n_peaks)],
                          Width_restricted = fitted_params_gaussian_restricted[(2*n_peaks+1):(3*n_peaks)],
                          AUC_restricted = unlist(lapply(1:n_peaks, function(i) sum(gaussian(x, fitted_params_gaussian_restricted[i], fitted_params_gaussian_restricted[n_peaks + i], fitted_params_gaussian_restricted[2 * n_peaks + i])))),       
                          
                          Amplitude_init = init_params[1:n_peaks],
                          Center_init = init_params[(n_peaks+1):(2*n_peaks)],
                          Width_init = init_params[(2*n_peaks+1):(3*n_peaks)],
                          row.names=paste(RCD_type,'counter',q,"peak", 1:n_peaks,sep = '_')) #add an counter for multiple spectra
  df<-rbind(df,params_df)
  print(q)
}

df$Number<-1:dim(df)[1]
df$Amplitude_Ratio<-df$Amplitude_restricted/df$Amplitude_init
df$Center_Ratio<-df$Center_restricted/df$Center_init
df$Center_Offset<-abs(df$Center_restricted-df$Center_init)
df$Width_Ratio<-df$Width_restricted/df$Width_init

pred_order<-c("Control","Apoptosis", "Ferroptosis","Necroptosis")
df$RCD_type<-factor(df$RCD_type, levels = pred_order)
#df_dump<-df
library(dplyr)
tolerance <- 15
grouped_data <- df %>%
  group_by(grp = cut(Center_restricted, breaks = sort(c(presetCenters - tolerance, (presetCenters) + tolerance))))

#Write to file
library(xlsx)
library(openxlsx)

#xlsx::write.xlsx(df,"test_restricted_presetcentra4.xlsx")
openxlsx::write.xlsx(x = df, file = "dataframe6.xlsx") #Seems to work better with larger files
saveRDS(df, file = "FullDataset_V2.Rds")

