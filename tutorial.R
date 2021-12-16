#####################################
## Set up and running Spatial NetRate
#####################################

library(VGAM)
library(tensorflow)
library(reticulate)
library(scales)
# anaconda prompt cheatsheet
# https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf


# set working directory to directory where the function and simulated data are stored
setwd('C:/Users/htopazia/OneDrive - Imperial College London/Github/spatial_rcs')

# source function to generate RC
source("function_genRC.R")

# Set up miniconda environment and download tensorflow and tensorflow probability (old versions for compatibility with code)
# Note: need to change paths 

# tensorflow::install_tensorflow(envname = 'r-reticulate',
#                                extra_packages = c("tensorflow-probability", "scipy"))

tensorflow::install_tensorflow(version = "1.14",
                               envname = 'r-reticulate',
                               extra_packages = c("tensorflow-probability==0.7", "scipy==1.5.0"),
                               conda_python_version = "3.6")

reticulate::use_python('r-reticulate')
Sys.setenv(RETICULATE_PYTHON="C:/Users/htopazia/Anaconda3/envs/r-reticulate")
options(reticulate.conda_binary = "C:/Users/htopazia/Anaconda3/envs/r-reticulate")

reticulate::py_config()
reticulate::conda_version()

# already installed these with install_tensorflow above
# reticulate::conda_install('r-reticulate', 'tensorflow')
# reticulate::conda_install('r-reticulate', 'tensorflow-probability')
# reticulate::conda_install('r-reticulate', 'scipy')

# anaconda command prompt arguments:
# conda create -n r-reticulate tensorflow
# conda env remove --name r-reticulate
# conda activate r-reticulate
# python
# import tensorflow-probability

# reticulate::conda_install('r-reticulate', 'tensorflow==1.14')
# reticulate::conda_install('r-reticulate', 'tensorflow-probability==0.7')
# reticulate::conda_install('r-reticulate', 'scipy==1.5.0')

# reticulate = a way to use python in R (import = library)
use_condaenv("r-reticulate", required = TRUE)
tfp <- reticulate::import("tensorflow_probability",convert=FALSE)
tf_config() 

# tfp <- reticulate::import_from_path("tensorflow-probability",path = "C:/Users/htopazia/Anaconda3/pkgs/tensorflow-probability-0.7-py_3", convert=FALSE)

# Load simulated data
res<-read.csv('res.csv')

n_seed <- 100 # imported cases are the first n_seed rows

# make matrices and reformat for function  
dd <- list(n=nrow(res), # n observations
           I= c(rep(1,length.out = n_seed), rep(0, length.out = nrow(res)-n_seed)), # imported cases
           t=res$inf_times, # time of infection
           d=res$inf_dist) # spatial coordinates

  # empty matrices for time and distance for each pair of observations
tmat <- dmat <- matrix(0, nrow=dd$n, ncol=dd$n) 

  # shift for 15 day minimum serial interval 
for (i in 1:dd$n) { # for infectees
  for(j in 1:(dd$n)){ # for all cases/infectors
    tmat[i,j]=dd$t[i]-dd$t[j]-15 #otherwise calculate time difference		
    
  }
}


tmat[tmat<0]=0 # anything less than zero is in the past so set time difference to zero
tmat <- tmat[-which(dd$I==1),] # imported cases cant be infectees and so remove these

# calculate distance in meters
for (i in 1:dd$n) { # for infectees
  for(j in 1:(dd$n)){ # for all cases/infectors
    dmat[i,j]= abs(dd$d[i]-dd$d[j]) # otherwise calculate spatial absolute difference		
    
  }
}

dmat <- dmat/1000 # convert to kilometers
dmat <- dmat[-which(dd$I==1),] # imported cases cant be infectees and so remove these

# RUN FUNCTION
sd_mid <- spatialnetrate(tp=tmat, 
                         dp= dmat, 
                         DataType = euclidian,
                         fixed = "epsilon", 
                         alpha = c(0.002, 0.001), # transmission rates
                         delta=c(0.01,0.001), # spatial parameters
                         SpatialKernel = "exponential", 
                         epsilon = 1e-20) # high = lots of unobserved infection, low = little missing data

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

