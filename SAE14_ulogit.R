



## This R code is a set of routines that correspond to chapter 14 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we load the dir2 function from the second chapter

source('SAE2_DesignDirect_Functions.R')



## Reads the survey data file and calculates new variables

dat <- read.table("datLCS.txt", header=TRUE, sep="\t", dec=",")
z0 <- 7280                          # Poverty threshold.
# Poverty variable: 1 if income<z0, 0 if income>z0
poor <- as.numeric(dat$income<z0)
# Labor situation: 0 if < 16 years, 1 if employed,
# 2 if unemployed, 3 if inactive.
work <- as.numeric(dat$lab=="1")    # Employed
nowork <- as.numeric(dat$lab=="2")  # Unemployed
inact <- as.numeric(dat$lab=="3")   # Inactive




## We rename some variables

income <- dat$income; w <- dat$w; dom <- dat$dom
one <- rep(1,nrow(dat))             # Auxiliary 1-variable
domains <- sort(unique(dat$dom))    # Domains sorted in ascending order
ndom <- length(domains)             # Number of domains




## The following code reads the file with the aggregated auxiliary variables and sort the file by domain.

aux <- read.table("auxLCS.txt", header=TRUE, sep="\t", dec=",")
# Totals of employed people
aux$Twork <- round(aux$TOT*aux$Mwork, 0)
# Totals of unemployed people
aux$Tnowork <- round(aux$TOT*aux$Mnowork, 0)    
# Totals in the categories "inactive" or "<16 years"
aux$Tres <- aux$TOT - aux$Twork - aux$Tnowork   
aux <- aux[order(aux$dom), ]   # Sort aux by dom (in ascending order)



## We calculate sample sizes

# sample sizes by domains and othres
nd <- tapply(X=one, INDEX=dom, FUN=sum)
ndworking <- tapply(one, list(dom,work), sum, default=0)
# sample sizes by domains and employment category
ndwork <- ndworking[,2]            
ndnoworking <- tapply(one, list(dom,nowork), sum, default=0)
# sample sizes by domains and unemployment category
ndnowork <- ndnoworking[,2]        
# sample sizes by domains and inactive or <16 categories
ndres <- nd - ndwork - ndnowork    



## We calculate direct estimators and estimated sizes by domains, by using {\textsf dir2} function described in section \ref{Rfunction.dir}.

# Direct estimates of poverty proportions
poor.dir <- dir2(data=poor, w=w, domain=list(dom))
dirp <- poor.dir$mean
# Estimated sizes
hatNd <- poor.dir$Nd.hat   
# Direct estimates of totals of employed people
hatNwork <- tapply(work*w, dom, sum)
# Direct estimates of totals of unemployed people
hatNnowork <- tapply(nowork*w, dom, sum)    
# Direct estimates of totals of remaining people
hatNdres <- hatNd - hatNwork - hatNnowork   



## We install and/or load the R package \textsf{lme4}.

if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}



## We fit a unit-level logit mixed model to the variable poor with the R library lme4

glmm <- glmer(formula=poor ~ work + nowork + (1|dom), data=dat,
              family=binomial)
summary(glmm)             # Summary of model results
pihat <- fitted(glmm)     # Predicted probabilities
beta <- fixef(glmm)       # Beta parameters
# Fixed effect for x0=1, x1=1, x2=0 (employed)
bwork <- beta[1] + beta[2]    
# Fixed effect for x0=1, x1=0, x2=1 (unemployed)
bnowork <- beta[1] + beta[3]
# Fixed effect for x0=1, x1=0, x2=0 (remaining categories)
binact <- beta[1]
# Variance components
var <- as.data.frame(VarCorr(glmm)) 
# Standard deviation of the random effect
phi <- var$sdcor[1]
ud <- ranef(glmm)$dom     # Random effects



## We calculate the plug-in predictors of poverty proportions by domains

# Poverty proportions by domains and employment category
pdworking <- tapply(pihat, list(dom,work), mean, default=0)
pdwork <- pdworking[,2]
# Poverty proportions by domains and unemployment category
pdnoworking <- tapply(pihat, list(dom,nowork), mean, default=0)
pdnowork <- pdnoworking[,2]
# Poverty proportions by domains and remaining categories
pdinactive <- tapply(pihat, list(dom,inact), mean, default=0)
pdinact <- pdinactive[,2]
# Poverty sample totals
totp <- tapply(poor, dom, sum)
# Plug-in estimates of poverty proportions by domains
plug.p <- (totp + (aux$Twork-ndwork)*pdwork + 
             (aux$Tnowork-ndnowork)*pdnowork +
             (aux$Tres-ndres)*pdinact)/aux$TOT



## We calculate the EBPs of poverty proportions by domains

S <- 1000             # Auxiliary terms
set.seed(123)         # Set seed of random number generator
# Generate the random values
v <- rnorm(n=S*ndom, mean=0, sd=1)
vv <- Map(split(v, rep(1:ndom,each=S)),
          split(-v, rep(1:ndom, each=S)), f=c)
vvphi <- lapply(vv, phi, FUN="*")
# Define some auxiliary variables
xbeta <- split(beta[1]+beta[2]*work+beta[3]*nowork, dat$dom)
# Calculate some terms of the EBP expression
a1 <- Map(phi*totp, vv, f="*")
a2 <- Dd <- vector()
term2 <- list()
for(d in 1:ndom){
  for(s in 1:(2*S)){
    a2[s] <- sum(log(1+exp(xbeta[[d]]+vvphi[[d]][s])))
  }
  term2[[d]] <- exp(a1[[d]] - a2)
  Dd[d] <- mean(term2[[d]])
}
expwork <- lapply(lapply(vvphi, bwork, FUN="+"), FUN=exp)
expnowork <- lapply(lapply(vvphi, bnowork, FUN="+"), FUN=exp)
expinact <- lapply(lapply(vvphi, binact, FUN="+"), FUN=exp)
expworkplus1 <- lapply(expwork, 1, FUN="+")
expnoworkplus1 <- lapply(expnowork, 1, FUN="+")
expinactplus1 <- lapply(expinact, 1, FUN="+")
Ads.work <- Map(Map(expwork, expworkplus1, f="/"), term2, f="*")
Ads.nowork <- Map(Map(expnowork, expnoworkplus1, f="/"), term2, f="*")
Ads.inact <- Map(Map(expinact, expinactplus1, f="/"), term2, f="*")
Ad.work <- sapply(Ads.work, FUN=mean)
Ad.nowork <- sapply(Ads.nowork, FUN=mean)
Ad.inact <- sapply(Ads.inact, FUN=mean)
# Calculating the EBPs by domains
ebp.work <- Ad.work/Dd
ebp.nowork <- Ad.nowork/Dd
ebp.inact <- Ad.inact/Dd
# Calculation of poverty proportion EBPs
ebp.p <- (totp + (aux$Twork-ndwork)*ebp.work +
            (aux$Tnowork-ndnowork)*ebp.nowork +
            (aux$Tres-ndres)*ebp.inact)/aux$TOT



## Saving the results

output <- data.frame(nd, dirp=round(dirp,3),
                     ebpp=round(ebp.p,3), plugp=round(plug.p,3))
head(output, 10)





###
