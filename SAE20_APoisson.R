



## This R code is a set of routines that correspond to chapter 20 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we load the dir2 function from the second chapter

source('SAE2_DesignDirect.R')



## Install and/or load the package lme4.

if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}



## Read the unit-level data file and defines the variable poor.

dat <- read.table("datLCS.txt", header=TRUE, sep="\t", dec=",")
# Poverty variable: 1 if yes (income < z0), 0 if no (income > z0).
z0 <- 7280                             # poverty threshold.
poor <- as.numeric(dat$income < z0)    # variable poor
one <- rep(1, nrow(dat))               # variable one
domains <- sort(unique(dat$dom))       # domains
ndom <- length(domains)                # number of domains



## Read the auxiliary data file.

aux <- read.table("auxLCS.txt", header=TRUE, sep="\t", dec=",")
aux <- aux[order(aux$dom),]            # sort aux by dom



## Calculate sample sizes and sample counts.

nd <- tapply(X=one, INDEX=dat$dom, FUN=sum)     # sample sizes by domains
# sample counts of poor people by domains
ndpoor <- tapply(poor, dat$dom, sum)



## Calculate direct estimators of sizes and poverty proportions by domains, by using dir2

# Direct estimates of poverty proportions
poor.dir <- dir2(data=poor, w=dat$w, domain=list(dat$dom))
dirp <- poor.dir$mean                         # poverty proportions
hatNd <- poor.dir$Nd.hat                      # estimated sizes
# root-MSE estimates of direct estimators
rmsedirp <- sqrt(poor.dir$var.mean)



## Fit an area-level Poisson mixed model to the sample counts of poor people

glmm <- glmer(formula=ndpoor ~ Minact + (1|dom), family=poisson,
              offset=log(nd), data=aux)
summary(glmm)
beta <- fixef(glmm)                           # regression parameters
var <- as.data.frame(VarCorr(glmm))           # variance parameters
phi <- var$sdcor                              # standard deviation parameter
# ML-Laplace predictions of random effects
ud <- ranef(glmm)$dom



## Calculate the plug-in predictors of sample counts of poor people and of poverty proportions.

# plug-in estimators of sample counts of poor people by domains
plugmu <- fitted(glmm)
# plug-in estimators of poverty proportions by domains
plugp <- plugmu/nd



## Calculate the EBPs of sample counts of poor people and of poverty proportions.

S <- 1000
set.seed(123)         # Set seed of random number generator
# Generate the random values
v <- rnorm(n=S*ndom, mean=0, sd=1)
vv <- Map(split(v, rep(1:ndom,each=S)), split(-v, rep(1:ndom, each=S)), f=c)
vvphi <- lapply(vv, phi, FUN="*")
Minact <- aux$Minact
# Define some auxiliary variables
xbeta.vvphi <- Map(beta[1]+beta[2]*aux$Minact, vvphi, f="+")
sum2 <- mapply(lapply(xbeta.vvphi, FUN=exp), nd, FUN="*")
Ad <- apply(exp(mapply(xbeta.vvphi, ndpoor+1, FUN="*") - sum2), 2, FUN=mean)
Dd <- apply(exp(mapply(xbeta.vvphi, ndpoor, FUN="*") - sum2), 2, FUN=mean)
ebp.p <- Ad/Dd        # EBPs of poverty proportions
ebp.mu <- ebp.p*nd    # EBPs of sample counts of poor people



## Calculate the parametric bootstrap estimates of the root-MSEs
## of the plug-in predictors and EBPs of poverty proportions.

set.seed(123); B <- 200
mseplugpcum.star <- mseebppcum.star <- 0
for (b in 1:B){
  cat("iteration = ", b, "\n")
  v.star <- rnorm(n=ndom, mean=0, sd=1)
  # True proportions of bootsrap model
  p.star <- exp(beta[1] + beta[2]*Minact + phi*v.star)
  # Bootstrap counts of poor people
  y.star <- rpois(n=ndom, lambda=nd*p.star)            
  mod.star <- glmer(formula=y.star ~ Minact + (1|dom), family=poisson,
                    offset=log(nd), data=aux)
  # regression parameters of bootstrap model
  b.star <- fixef(mod.star)
  # standard deviation of bootstrap model
  phi.star <- as.data.frame(VarCorr(mod.star))$sdcor
  # Estimated proportions of bootstrap model
  plugp.star <- fitted(mod.star)/nd                    
  # Calculation of EBP with bootstrap data
  v.star <- rnorm(n=S*ndom, mean=0, sd=1)   # Standard normal random numbers
  vv.star <- Map(split(v.star, rep(1:ndom,each=S)),
                 split(-v.star, rep(1:ndom, each=S)), f=c)
  vvphi.star <- lapply(vv.star, phi, FUN="*")
  # Define some auxiliary variables
  xbeta.vvphi.star <- Map(beta[1]+beta[2]*Minact, vvphi.star, f="+")
  sum2.star <- mapply(lapply(xbeta.vvphi.star, FUN=exp), nd, FUN="*")
  Ad.star <- apply(exp(mapply(xbeta.vvphi.star, y.star+1, FUN="*") -
                         sum2.star), 2, FUN=mean)
  Dd.star <- apply(exp(mapply(xbeta.vvphi.star, y.star, FUN="*") -
                         sum2.star), 2, FUN=mean)
  ebp.p.star <- Ad.star/Dd.star       # EBPs of poverty proportions
  mseplugpcum.star <- mseplugpcum.star + (plugp.star - p.star)^2
  mseebppcum.star <- mseebppcum.star + (ebp.p.star - p.star)^2
}
# root-MSE bootstrap estimates of plug-in
rmseplugp.star <- sqrt(mseplugpcum.star/B)
# root-MSE bootstrap estimates of EBPs
rmseebpp.star <- sqrt(mseebppcum.star/B)



## Saving the results

output <- data.frame(nd, ndmu=ndpoor, ebpmu=round(ebp.mu,2),
                     plugmu=round(plugmu,2), dirp, ebp.p, plugp,
                     rdirp=rmsedirp, rplugp=rmseplugp.star,
                     rebpp=rmseebpp.star)
head(round(output,4), 10)





###
