



## This R code is a set of routines that correspond to chapter 15 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we load the dir2 function from the second chapter

source('SAE2_DesignDirect_Functions.R')



## The following R code reads the aggregated auxiliary data in Ndsa20.txt (true population sizes) and calculates the corresponding sizes of covariate classes by area and age group.

aux <- read.table("Ndsa20.txt", header=TRUE, sep="\t", dec=".")
# Totals by subdomains
Ndt <- tapply(X=aux$N, INDEX=list(aux$area,aux$age), FUN=sum)    
# Totals by subdomains-edu1
Ndt.edu1 <- tapply(aux$edu1, list(aux$area,aux$age), sum)        
# Totals by subdomains-edu2
Ndt.edu2 <- tapply(aux$edu2, list(aux$area,aux$age), sum)
# Totals by subdomains-edu3
Ndt.edu3 <- tapply(aux$edu3, list(aux$area,aux$age), sum)



## The following R code reads the unit-level data file and defines the age groups and variables indicating the \textsc{education} categories.

dat <- read.table("LFS20.txt", header=TRUE, sep="\t", dec=".")
ns <- nrow(dat)                            # Global sample size
n.areas <- length(unique(dat$AREA))        # Number of areas
z0 <- 36500                                # Poverty threshold
# Poverty variable: 1 if INCOME<z0, 0 if INCOME>z0
poor <- as.numeric(dat$INCOME<z0)
Ga <- cut(dat$AGE, breaks=c(0,25,54,max(dat$AGE)), labels=1:3, right=TRUE)
ageG <- as.numeric(Ga)
edu2 <- as.numeric(dat$EDUCATION==2)
edu3 <- as.numeric(dat$EDUCATION==3)
one <- rep(1, ns)



## We calculate sample sizes, counts and means in subdomains

# Sizes by subdomains
ndt <- tapply(X=one, INDEX=list(dat$AREA,ageG), FUN=sum)  
# Sample counts of edu3
ndt.edu3 <- tapply(edu3, list(dat$AREA,ageG), sum)        
# Sample counts of edu2
ndt.edu2 <- tapply(edu2, list(dat$AREA,ageG), sum)        
# Sample counts of edu1
ndt.edu1 <- ndt - ndt.edu3 - ndt.edu2                     
# Sample counts of poor
ndtpoor <- tapply(poor, list(dat$AREA,ageG), sum)         
# Sample means of poor
mdtpoor <- ndtpoor/ndt                                    



## We calculate direct estimators and estimated sizes by subdomains, by using  dir2  function described in the second chapter of the book

dir.poor <- dir2(data=poor, w=dat$WEIGHT,
                 domain=list(area=dat$AREA, ageG=ageG))
hatNdt <- dir.poor$Nd.hat    # Estimated sizes
dir.p <- dir.poor$mean       # Direct estimates of poverty proportions



## We install and/or load the R package \textsf{lme4}.

if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}



## We fit the unit-level two-fold logit mixed model to the variable y=poor by applying  the glmer function of the R library lme4.

glmm <- glmer(formula=poor ~ edu3 + edu2 + (1|AREA/ageG), data=dat,
              family=binomial)
summary(glmm)                           # Summary of model results
pihat <- fitted(glmm)                   # Predicted probabilities
beta <- fixef(glmm)                     # Beta parameters
bedu3 <- beta[1]+beta[2]                # Fixed effect for x1=1, x2=1, x3=0
bedu2 <- beta[1]+beta[3]                # Fixed effect for x1=1, x2=0, x3=1
bedu1 <- beta[1]                        # Fixed effect for x1=1, x2=0, x3=0
var <- as.data.frame(VarCorr(glmm))     # Variance parameters
# Standard deviation of area:age random effects
phi2.e <- var$sdcor[1]
# Standard deviation of area random effects
phi1.e <- var$sdcor[2]
r.effects <- ranef(glmm)                # Random effetcs
# Modes of the random effects for area
ud <- as.matrix(r.effects[[2]])         
# Modes of the random effects for area:age
re.dt <- as.matrix(r.effects[[1]])
udt <- matrix(re.dt, nrow=n.areas)



## Using the fitted model, we first calculate the plug-in estimators of poverty proportions by subdomains.

uud <- sweep(udt, 1, ud, "+")     # Random effects for subdomains
etadt.edu1 <- bedu1 + uud
# pihat by subdomains and education level 1
pd.edu1 <- exp(etadt.edu1)/(1+exp(etadt.edu1))
etadt.edu2 <- bedu2 + uud
# pihat by subdomains and education level 2
pd.edu2 <- exp(etadt.edu2)/(1+exp(etadt.edu2))   
etadt.edu3 <- bedu3 + uud
# pihat by subdomains and education level 3
pd.edu3 <- exp(etadt.edu3)/(1+exp(etadt.edu3))   
# Poverty proportion plug-in estimators
plug.p <- (ndtpoor + (Ndt.edu1-ndt.edu1)*pd.edu1 +
             (Ndt.edu2-ndt.edu2)*pd.edu2 +
             (Ndt.edu3-ndt.edu3)*pd.edu3)/Ndt



## We prepare the data in the format used by the R function calc.EBP.y

D <- n.areas
TT <- 3          # Number of levels of aggregation (or time instants)
K <- 3           # Number of covariate classes
# z_k for edu1, edu2, edu3
class.vec <- rbind(c(1,0,0), c(1,0,1), c(1,1,0)) 
# Calculating the unobserved covariate classes sizes
N.dtk.r <- array(0, dim=c(D,TT,K))
N.dtk.r[,,1] <- Ndt.edu1 - ndt.edu1
N.dtk.r[,,2] <- Ndt.edu2 - ndt.edu2
N.dtk.r[,,3] <- Ndt.edu3 - ndt.edu3
S1 <- S2 <- 30              # Number of points for Monte-Carlo approximation
beta.e <- as.matrix(beta)   # Estimated parameters beta in a column format
X <- cbind(one, edu3, edu2) # Matrix X of our model
y <- poor                   # Vector of the target variable y



## The R code of the function calc.EBP.y that calculates 
## the EBPs of subdomain poverty proportions is

calc.EBP.y <- function(y, X, beta.e, phi1.e, phi2.e,
                       S1, S2, D, TT, K, class.vec, w.dtk.r){
  # Generate the random effects
  v1 <- rnorm(D * S1,mean = 0,sd = 1)
  v1 <- matrix(v1, nrow = D)
  v1 <- cbind(v1, -v1)
  v2 <- rnorm(D*TT*S2, mean=0, sd=1)
  v2 <- array(v2, dim=c(D,TT,S2))
  v2 <- abind(v2, -v2)
  # Initialize variables
  y.dtau <- matrix(0, D, TT)
  q.dtk <- rep(0,K)  #dim=c(D, TT, K))
  mu.dt <- matrix(0, D, TT)
  # Begin the loop for calculating the EBPs by domains
  for(d in 1:D) {
    # Begin the calculation of inner sum-products over s1, tau, s2
    sum.s1.D <- 0
    sum.s1.N <- matrix(0, TT, K)
    # Begin the loop s1
    for(s1 in 1:(2*S1)){
      prod.tau.D <- 1
      prod.tau.N <- matrix(1, TT, K)
      # Begin the loop tau
      for(tau in 1:TT){
        # Specify the matrix xx_dt and sample vector y.dt
        x.dt <- X[dat$AREA==d & dat$ageG==tau,]
        y.dt <- y[dat$AREA==d & dat$ageG==tau]
        # Calculate the sums of y.dt
        y.dtau[d,tau] <- sum(y.dt)
        sum.s2.D <- 0
        sum.s2.N <- matrix(0,TT,K)
        for(s2 in 1:(2*S2)){
          ran.term <- phi1.e*v1[d,s1] + phi2.e*v2[d,tau,s2]
          # Calculate the inner sum over i
          eta <- x.dt%*%beta.e + ran.term
          vec <- log(1+exp(eta))
          sum.i <- sum(vec)
          # End of the calculation of the inner sum over i
          exponent.D <- y.dtau[d,tau]*ran.term - sum.i
          sum.s2.D <- sum.s2.D + exp(exponent.D)
          # Calculate the sum over s2 for N.dtk
          for(t in 1:TT){
            for(k in 1:K){
              if(t==tau){
                exponent <- class.vec[k,]%*%beta.e +
                  (y.dtau[d,tau]+1)*ran.term -
                  log(1+exp(class.vec[k,]%*%beta.e+ran.term)) -
                  sum.i
              }
              else{
                exponent <- exponent.D
              }
              sum.s2.N[t,k] <- sum.s2.N[t,k] + exp(exponent)
            }
          }
        }
        # End of the loop in s2
        prod.tau.D <- prod.tau.D * sum.s2.D
        prod.tau.N <- prod.tau.N * sum.s2.N
      }
      # End of the loop in tau
      sum.s1.D <- sum.s1.D + prod.tau.D
      sum.s1.N <- sum.s1.N + prod.tau.N
      cat("d=", d, s1, "\n")
    }
    # End of the loop in s1
    ## End of the calculation of the inner sum-products over s1, tau, s2
    D.term <- sum.s1.D
    N.dtk <- sum.s1.N
    q.dtk <- N.dtk/D.term
    for(t in 1:TT){
      mu.dt[d,t] <- as.numeric(w.dtk.r[d,t,] %*% q.dtk[t,] + y.dtau[d,t] )
    }
  }
  ## End of the loop for calculating the EBPs by domains
  ebpp <- mu.dt/Ndt
  return(ebpp)
}



## First of all, we must load the abind package and the function calc.EBP.y
## code that appears below, then set the seed for the random number generator
## and calculate the EBPs of poverty proportions.
if(!require(abind)){
  install.packages("abind")
  library(abind)
}
set.seed(123)
ebp.p <- calc.EBP.y(y, X, beta.e, phi1.e, phi2.e,
                    S1, S2, D, TT, K, class.vec, N.dtk.r)



## Saving the results

dir.p <- round(data.frame(dir.p[1:20], dir.p[21:40], dir.p[41:60]),4)
output <- data.frame(ndt, round(ebp.p,4), round(plug.p,4), dir.p)
names(output) <- c("n1", "n2", "n3", "ebpp1", "ebpp2", "ebpp3",
                   "plugp1", "plugp2", "plugp3",
                   "dirp1", "dirp2", "dirp3")
head(output, 10)





###
