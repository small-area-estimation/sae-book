



## This R code is a set of routines that correspond to chapter 19 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## In this chapter the code has been separated in two parts.
## This part contains the main program. Then, it is neccessary
## to load the other part previously

source('SAE19_BFH_Functions.R')



## We install the R packages sae and Matrix.

if(!require(sae)){
  install.packages("sae")
  library(sae)
}
if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}



## It reads the data files and renames some variables.

dat <- read.table("datLCS.txt", header=TRUE, sep="\t", dec=",")    # Read dat
# poverty threshold fixed at 7280
z0 <- 7280                                                         
# variable poor
poor <- as.numeric(dat$income<z0)
# variable gap
gap <- (z0-dat$income)*poor/z0                                     
# Read auxiliary data
aux <- read.table("auxLCS2.txt", header=TRUE, sep="\t", dec=",")
aux <- aux[order(aux$dom),]     # Sort aux by dom
D <- nrow(aux)                  # Number of domains
ss <- aux$ss                    # Select ss as aux variable



## We calculate the direct estimators of domain poverty proportions and gaps,
## their variance and covariance estimates and the population sizes by domain.
## We look for those cases where the estimates of the variances are equal to zero.

poor.dir <- direct.bfh(data1=poor, data2=gap, w=dat$w,
                       domain=list(domain=dat$dom))
y <- lapply(1:D,
            function(d) matrix(c(poor.dir$y1.mean[d], poor.dir$y2.mean[d])))
sel.pprop <- poor.dir$var.y1.mean>0
sel.pgap <- poor.dir$var.y2.mean>0
poor.dir$var.y1.mean[!sel.pprop] <- min(poor.dir$var.y1.mean[sel.pprop])
poor.dir$var.y2.mean[!sel.pgap] <- min(poor.dir$var.y2.mean[sel.pgap])



## We define the covariance matrices Ved of the vectors of sampling errors.

Ved <- list()
for(d in 1:D){
  Ved[[d]] <- matrix(NA, nrow=2, ncol=2)
  sup <- as.numeric(poor.dir$cov.y12.mean[d])
  Ved[[d]][upper.tri(Ved[[d]])] <- sup
  Ved[[d]][lower.tri(Ved[[d]])] <- sup
  diag(Ved[[d]]) <- as.numeric(c(poor.dir$var.y1.mean[d],
                                 poor.dir$var.y2.mean[d]))
}



## We use the function mseFH of the R package sae
## for fitting the FH models, that will be used as seeds.

# Target variables
directy1 <- poor.dir$y1.mean; directy2 <- poor.dir$y2.mean
# Sampling error variances
vardiry1 <- poor.dir$var.y1.mean; vardiry2 <- poor.dir$var.y2.mean
fit.FH.sae1 <- mseFH(directy1 ~ ss, vardiry1)    # FH model for y1
fit.FH.sae2 <- mseFH(directy2 ~ ss, vardiry2)    # FH model for y2



## The following code calculates the function inputs.

# Construct the matrix of auxiliary variables.
X <- lapply(1:D, 
            function(d) as.matrix(t(bdiag(as.numeric(c(1,ss[d])),
                                          as.numeric(c(1,ss[d]))))))
# Number of regression parameters
p <- ncol(X[[1]])
# Seeds for variance component parameters
thetas.0 <- c(fit.FH.sae1$est$fit$refvar, fit.FH.sae2$est$fit$refvar, 0)
Vud <- UveU(thetas.0)



## Apply the R function REML.BFH to fit a bivariate Fay-Herriot model.

fit <- try(REML.BFH(X, y, D, Ved, Vud), TRUE)
cat("BFH model parameters", fit[[1]], "\n")
cat("Number of iterations", fit[[3]], "\n")
cat("Errors: ", fit[[4]], "\n")



## R functions to calculate estimated parameters and CIs of BFH model

beta.u.hat <- BETA.U.BFH(X, y, D, Ved , fit[[1]])
beta.hat <- beta.u.hat[[1]]        # Estimates of regression parameters
theta.hat <- fit[[1]]              # Estimates of variance components
u <- beta.u.hat[[2]]               # Predictions of random effects
pv <- pvalue(beta.hat, fit)        # Apply function pvalue
# Anova table with standard errors, z-values and p-values
betas <- data.frame(beta.hat, pv, Sig=pv[[3]]<0.05)
row.names(betas) <- c("Intercept.p", "ss.p", "Intercept.g",  "ss.g")
# Apply function CI for confidence intervals
C.I <- CI(beta.hat, fit)



## In the variable eblup.bfh we calculate the EBLUPs of poverty proportions and gaps.

eblup.bfh <- mapply(FUN="+", lapply(X, FUN="%*%", beta.hat), u)
eblup.bfh.1 <- eblup.bfh[1,]    # EBLUPs of poverty proportions
eblup.bfh.2 <- eblup.bfh[2,]    # EBLUPs of poverty gaps



## The MSE.BFH function calculates analytic estimators of the MSEs.

mse.bfh <- MSE.BFH(X, D, Ved, fit)    # Apply function mse.bfh
# MSE estimates
mse.bfh.1 <- mse.bfh[[1]]; mse.bfh.2 <- mse.bfh[[2]]
# coefficients of variation
CVdir1 <- round(100*sqrt(poor.dir$var.y1.mean/abs(poor.dir$y1.mean)), 2)
CVdir2 <- round(100*sqrt(poor.dir$var.y2.mean/abs(poor.dir$y2.mean)), 2)
CV.bfh.1 <- round(100*sqrt(mse.bfh.1/abs(eblup.bfh.1)), 2)
CV.bfh.2 <- round(100*sqrt(mse.bfh.2/abs(eblup.bfh.2)), 2)



## Saving the results

output1 <- data.frame(nd=poor.dir$nd, DIR=round(poor.dir$y1.mean,5),
                      Vdir=round(poor.dir$var.y1.mean,5), CVdir=CVdir1, 
                      EB=round(eblup.bfh.1,5), MSEeb=round(mse.bfh.1,5),
                      CVeb=CV.bfh.1)
head(output1, 10)
output2 <- data.frame(nd=poor.dir$nd, DIR=round(poor.dir$y2.mean,5),
                      Vdir=round(poor.dir$var.y2.mean,5), CVdir=CVdir2, 
                      EB=round(eblup.bfh.2,5), MSEeb=round(mse.bfh.2,5),
                      CVeb=CV.bfh.2)
head(output2, 10)





###
