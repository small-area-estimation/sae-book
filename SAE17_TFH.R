



## This R code is a set of routines that correspond to chapter 17 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, install and/or load the package saery and load the data set

if(!require(saery)){
  install.packages("saery")
  library(saery)
}
data(datos)



## Redefine variables and calculate auxiliary information

sigma2edi <- datos[,6]            # Sampling variance
X <- as.matrix(datos[,5])         # Auxiliary variable
ydi <- datos[,3]                  # Target variable
D <- length(unique(datos[,1]))    # Number of domains
md <- table(datos[,1])            # Number of time instants per domain



## Fit the non-correlated area-level temporal linear mixed model

output.fit.indep <- fit.saery(X, ydi, D, md, sigma2edi, model="INDEP",
                              conf.level=0.9)
output.fit.indep                  # Summary of results
output.fit.indep$Regression       # Regression parameters
output.fit.indep$SIGMA            # Variance components



## Calculate the EBLUP and the explicit MSE estimator.

eblup.output.indep <- eblup.saery(X, ydi, D, md, sigma2edi,
                                  model="INDEP", plot=TRUE)                                  
head(round(eblup.output.indep, 4), 10)      # Summary of results



## Fit the AR(1)-correlated area-level temporal linear mixed model.

output.fit.ar1 <- fit.saery(X, ydi, D, md, sigma2edi, model="AR1",
                            conf.level=0.9)
output.fit.ar1                    # Summary of results
output.fit.ar1$Regression         # Regression parameters
output.fit.ar1$SIGMA              # Variance components



## Calculate the EBLUP and and the explicit MSE estimator under the
## AR(1)-correlated area-level temporal linear mixed model.

eblup.output.ar1 <- eblup.saery(X, ydi, D, md, sigma2edi, "AR1", type="III")
head(round(eblup.output.ar1, 4), 10)        # Summary of results





###
