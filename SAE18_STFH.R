



## This R code is a set of routines that correspond to chapter 18 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, install and/or load the package sae and load the data set grapes

if(!require(sae)){
  install.packages("sae")
  library(sae)
}
data(grapes); dim(grapes)            # Load data set
data(grapesprox); dim(grapesprox)    # Load proximity matrix
apply(grapesprox, 1, sum)            # Check that rows sum up to one



## Fit the spatial Fay-Herriot model by using the ML method

resultML <- eblupSFH(grapehect ~ area + workdays - 1, var, grapesprox,
                     method="ML", data=grapes)
head(resultML$eblup, 10)             # EBPLUPs of first 10 domains
resultML$fit$estcoef                 # Regression coefficients
resultML$fit$refvar                  # Variance of random effect
resultML$fit$spatialcorr             # Spatial correlation parameters
resultML$fit$goodness                # AIC



## Fit the spatial Fay-Herriot model using REML method

resultREML <- eblupSFH(grapehect ~ area + workdays - 1, var, grapesprox,
                       data=grapes)
head(resultREML$eblup, 10)           # EBPLUPs of first 10 domains
resultREML$fit$estcoef               # Regression coefficients
resultREML$fit$refvar                # Variance of random effect
resultREML$fit$spatialcorr           # Spatial correlation parameters
resultREML$fit$goodness              # AIC



# Load the data set spacetime

data(spacetime); dim(spacetime)            # Load data set
data(spacetimeprox); dim(spacetimeprox)    # Load proximity matrix
D <- nrow(spacetimeprox)                   # number of domains
T <- length(unique(spacetime$Time))        # number of time instants



# Fit the spatio-temporal Fay-Herriot model with uncorrelated time effects for each domain

resultS <- eblupSTFH(Y ~ X1 + X2, D, T, Var, spacetimeprox,
                     model="S", data=spacetime)
resultS$fit$estcoef                        # Regression coefficients
resultS$fit$estvarcomp                     # Variance components
# EPLUPs for the last time domain
subset(resultS$eblup, spacetime$Time==T)



# Fit the spatio-temporal Fay-Herriot model with AR(1)-correlated time effects for each domain

resultST <- eblupSTFH(Y ~ X1 + X2, D, T, Var, spacetimeprox,
                      model="ST", data=spacetime)
resultST$fit$estcoef                       # Regression coefficients
resultST$fit$estvarcomp                    # Variance components
# EPLUPs for the last time domain
subset(resultST$eblup, spacetime$Time==T)



## The R code to save the results is

output <- data.frame(Area=subset(spacetime$Area, spacetime$Time==T),
                     Indep=subset(resultS$eblup, spacetime$Time==T),
                     AR1=subset(resultST$eblup, spacetime$Time==T))
head(output, 10)





###
