



## This R code is a set of routines that correspond to chapter 16 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we load the dir2 function from the second chapter

source('SAE2_DesignDirect_Functions.R')



## We install and/or load the R package sae

if(!require(sae)){
  install.packages("sae")
  library(sae)
}



## The following code reads the data files and rename some variables.

# Read dat
dat <- read.table("datLCS.txt",header=TRUE,sep="\t",dec=",")    
# Read aux
aux <- read.table("auxLCS.txt",header=TRUE,sep="\t",dec=",")    
# Sort aux by dom
aux <- aux[order(aux$dom),]
# Rename variables
income <- dat$income; w <- dat$w; dom <- dat$dom; D <- nrow(aux)    
Mwork <- aux$Mwork; Mnowork <- aux$Mnowork
Minact <- aux$Minact; Nd <- aux$TOT



## We calculate direct estimators of domain average incomes, their variance estimates and the population sizes by domain.

income.dir <- dir2(data=income, w, domain=list(dom=dom))
diry <- income.dir$mean; vardiry <- income.dir$var.mean
CVdir <- round(100*sqrt(vardiry)/diry, 2)
hatNd <- income.dir$Nd.hat; nd <- income.dir$nd
log.var <- log(vardiry)



## We apply the Generalized Variance Function (GVF) approach. We fit the model

model.A <- lm(log.var ~ diry + nd + diry:nd)    # GVF Model A
summary(model.A)                                # Model coefficients
res <- resid(model.A)                           # Residuals
p <- predict(model.A)                           # Fitted values
# Error variance of GVF model
v <- deviance(model.A)/df.residual(model.A)
Vgvf <- exp(p)*exp(v/2)                         # GVF variance estimates	
CVgvf <- round(100*sqrt(Vgvf)/diry, 2)          # GVF CV estimates



## We select a Fay-Herriot model. First we fit the model with the three auxiliary variables,
## but employed is not significative. We select the model with unemployed and inactive as regressors.

# Model 123
M123 <- eblupFH(diry ~ Mwork + Mnowork + Minact, Vgvf)
M123$fit$estcoef; M123$fit$goodness[2]
# Take Mwork out. It has the highest p-value
M23 <- eblupFH(diry ~ Mnowork + Minact, Vgvf)
M23$fit$estcoef; M23$fit$goodness[2]
# Select model M23



## We calculate the EBLUPs of income means by domain and the corresponding MSE estimators.

# EBLUPs
eblup <- M23$eblup
mseM23 <- mseFH(diry ~ Mnowork + Minact, Vgvf)
MSEeblup <- mseM23$mse
# CV with respect to direct estimator
CVeblup <- round(100*sqrt(MSEeblup)/diry, 2)



## Saving the results

output <- data.frame(nd, DIR=round(diry), Vdir=round(vardiry),
                     Vgvf=round(Vgvf), CVgvf, 
                     EB=round(eblup), MSEeb=round(MSEeblup), CVeb=CVeblup)
head(output, 10)





###
