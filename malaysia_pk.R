#####################################
# Set up and running Spatial NetRate
#####################################

library(VGAM)
library(tensorflow)
library(reticulate)
library(scales)
library(tidyverse)
library(lubridate)

# set working directory to directory where the function and data are stored
setwd('C:/Users/htopazia/OneDrive - Imperial College London/Github/spatial_rcs')

# source function to generate RC
source("function_genRC.R")

# set up anaconda environment and download tensorflow and tensorflow_probability 
# (old versions for compatibility with code)
tensorflow::install_tensorflow(version = "1.14",
                               envname = 'r-reticulate',
                               extra_packages = c("tensorflow-probability==0.7", "scipy==1.5.0"),
                               conda_python_version = "3.6")

# tell R which python to use
reticulate::use_python('r-reticulate')
Sys.setenv(RETICULATE_PYTHON="C:/Users/htopazia/Anaconda3/envs/r-reticulate")
options(reticulate.conda_binary = "C:/Users/htopazia/Anaconda3/envs/r-reticulate")

# check that python and conda versions are correct
reticulate::py_config()
reticulate::conda_version()

# anaconda command prompt arguments to use as alternatives:
# https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf
# conda create -n r-reticulate tensorflow
# conda env remove --name r-reticulate
# conda activate r-reticulate
# import tensorflow-probability

# tell R which environment and package to use
use_condaenv("r-reticulate", required = TRUE)
tfp <- reticulate::import("tensorflow_probability",convert=FALSE)
tf_config() 

# simulated data test ----------------------------------------------------------
res <- read.csv('res.csv')

n_imp <- 100 # imported cases are the first n_seed rows

# make time and distance matrices and format for function  
dd <- list(n=nrow(res), # n observations
           I= c(rep(1,length.out = n_imp), rep(0, length.out = nrow(res)-n_imp)), # imported cases
           t=res$inf_times, # time of infection
           d=res$inf_dist) # spatial point

# empty matrices for time and distance for each pair of observations
tmat <- dmat <- matrix(0, nrow=dd$n, ncol=dd$n) 

# shift time difference minimum serial interval 
SI <- 15 # serial interval = 15 days

for (i in 1:dd$n) { # for infectees
  for(j in 1:(dd$n)){ # for all cases/infectors
    tmat[i,j]=dd$t[i]-dd$t[j]-SI # find difference between case times - serial interval
    
  }
}

tmat[tmat<0]=0 # anything less than zero is in the past so set time difference to zero
tmat <- tmat[-which(dd$I==1),] # imported cases cannot be infectees, remove these

# calculate distance in meters
for (i in 1:dd$n) { # for infectees
  for(j in 1:(dd$n)){ # for all cases/infectors
    dmat[i,j]= abs(dd$d[i]-dd$d[j]) # find difference between spatial points		
    
  }
}

dmat <- dmat/1000 # convert meters to kilometers
dmat <- dmat[-which(dd$I==1),] # imported cases cant be infectees, remove these

# RUN FUNCTION
sd_mid <- spatialnetrate(tp=tmat, 
                         dp= dmat, 
                         DataType = euclidian,
                         fixed = "epsilon", 
                         alpha = c(0.002, 0.001), # transmission rates
                         delta=c(0.01,0.001), # spatial parameters
                         SpatialKernel = "exponential", 
                         epsilon = 1e-20) 

# inspect results 

# distribution of Rc
hist(sd_mid[[1]])

# delta parameter
sd_mid[[2]]

# alpha parameter 
hist(sd_mid[[3]])

# epsilon
sd_mid[[4]]

# AIC
sd_mid[[5]]



# Pk example -------------------------------------------------------------------
MYSdata <- read_csv("MYSdata.csv")
MYSdata <- MYSdata %>%          # fill in with variable names
   mutate(case_number = case_number,
          date_onset =  dmy(date_onset), # DD-MM-YYY
          import =      import,          # if cases are imported (suspected from forest? no parent case in dataset)
          longitude =   longitude,       # decimal degrees WGS84
          latitude =    latitude) %>%    # decimal degrees WGS84

   mutate(date = as.numeric(difftime(date_onset, dmy('01/01/2000'), unit='days')), # days since 2000
          date = ifelse(duplicated(date), jitter(date, 1), date)) # if date is a duplicate, jitter point

MYSdata <- st_as_sf(MYSdata, coords = c('longitude', 'latitude'), crs = 4326) # WGS84


# make time and distance matrices
dd <- list(n = nrow(MYSdata), # n observations
           I = MYSdata$import, # imported cases
           t = MYSdata$date, # time of infection
           d = MYSdata$geometry) # spatial point

# empty matrices for time and distance for each pair of observations
tmat <- dmat <- matrix(0, nrow=dd$n, ncol=dd$n) 
      
      # shift time difference minimum serial interval 
      SI <- 15 # set serial interval = 15 days
      
      for (i in 1:dd$n) { # for infectees
        for(j in 1:(dd$n)){ # for all cases/infectors
          tmat[i,j]=dd$t[i]-dd$t[j]-SI # find difference between case times - serial interval
          
        }
      }
      
      tmat[tmat<0]=0 # anything less than zero is in the past so set time difference to zero
      tmat <- tmat[-which(dd$I==1),] # imported cases cannot be infectees, remove these
      
      # calculate distance in decimal degrees
      dmat = st_distance(MYSdata, MYSdata, which="Euclidean", by_element=F) 
      dmat <- dmat[-which(dd$I==1),] # imported cases cannot be infectees, remove these
      
      # RUN FUNCTION
      sd_mid <- spatialnetrate(tp = tmat, 
                               dp = dmat, 
                               DataType = euclidian,
                               fixed = "epsilon", # fixed variables
                               alpha = c(0.002, 0.001), # transmission rates
                               delta = c(0.01,0.001), # spatial parameters
                               SpatialKernel = "exponential", 
                               epsilon = 1e-20) # high = lots of unobserved infection, 
                                                # low = little missing data
      
      # inspect results 
      
      # distribution of Rc
      hist(sd_mid[[1]])
      
      # delta parameter
      sd_mid[[2]]
      
      # alpha parameter 
      hist(sd_mid[[3]])
      
      # epsilon
      sd_mid[[4]]
      
      # AIC
      sd_mid[[5]]

