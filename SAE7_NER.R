



## This R code is a set of routines that correspond to chapter 7 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## In this chapter the code has been separated in two parts.
## This part contains the main program. Then, it is neccessary
## to load the other part previously

source('SAE7_NER_Functions.R')



## First we read the sample data file and and make some calculations

dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
dat$EDUCATION <- as.factor(dat$EDUCATION)  # EDUCATION as factor
y <- 10^(-4)*dat$INCOME
age <- dat$ageG
edu2 <- as.numeric(dat$EDUCATION==2)
edu3 <- as.numeric(dat$EDUCATION==3)
n <- length(dat$AREA)                      # Global sample size
domain <- sort(unique(dat$AREA))           # Domains
D <- length(domain)                        # Number of domains
# Domain sizes
nd <- tapply(rep(1, n), list(dat$AREA), sum, simplify = FALSE)
yd <- Xd <- Vd <- list()
for (d in 1:D) {
  condition <- dat$AREA==domain[d]
  # Auxiliary variables
  Xd[[d]] <- cbind(1, edu2[condition],edu3[condition])
  # Target variable
  yd[[d]] <- y[condition]
}



## We load/install the R package lme4

if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}
lmm.1 <- lmer(formula=y~EDUCATION+(1|AREA), data=dat, REML=FALSE)
summary(lmm.1)



## We calculate the seeds for the model parameters

start <- seed(Xd, yd, n, nd)



## The function FisherScoring applies the Fisher-scoring algorithm
## described in the book. This function already has R codes for calculating the seeds

ML.fit <- FisherScoring(Xd, yd, n, nd, n.iter=200)



## The R function nlm also calculates the MLE estimators of the parameters
## of a NER model by minimizing the log-likelihood function log_lik.

ML.nlm <- nlm(log_lik, p=start, Xd=Xd, yd=yd, n=n, nd=nd, hessian=TRUE)



## We use the function optimx of the R package optimx
## for fitting the FH models, that will be used as seeds.

if(!require(optimx)){
  install.packages("optimx")
  library(optimx)
}
ML.optimx <- optimx(par=start, log_lik, Xd=Xd, yd=yd, n=n, nd=nd,
                    method="BFGS")



## The function nleqslv of the R package nleqslv can be applied
## to solve a system of nonlinear equations.

if(!require(nleqslv)){
  install.packages("nleqslv")
  library(nleqslv)
}
ML.nleqslv <- nleqslv(start, scores, Xd=Xd, yd=yd, n=n, nd=nd, jacobian=TRUE)



## The functions dfsane and BBsolve of the R package BB also solve
## the system of log-likelihood equations.

if(!require(BB)){
  install.packages("BB")
  library(BB)
}
ML.dfsane <- dfsane(start, scores, Xd=Xd, yd=yd, n=n, nd=nd, method=1)
# Alternative for dfsane
ML.BBsolve <- BBsolve(start, scores, Xd=Xd, yd=yd, n=n, nd=nd, method=1)



## Saving the results

var <- as.data.frame(VarCorr(lmm.1))
output <- data.frame(seeds=start,
                     lmer=c(lmm.1@beta, var$vcov[1], var$vcov[2]),
                     FisherScoring=ML.fit[[5]],
                     nlm=ML.nlm$estimate,
                     optimx=as.numeric(ML.optimx[1,1:5]),
                     nleqslv=ML.nleqslv$x, dfsane=ML.dfsane$par)
round(output, 5)





###
