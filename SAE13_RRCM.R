



## This R code is a set of routines that correspond to chapter 13 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we install/load load some packages and load the dir2 function
## from the second chapter

source('SAE2_DesignDirect_Functions.R')

if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}
if(!require(sae)){
  install.packages("sae")
  library(sae)
}



## We read the data files and calculate some variables.

# Read unit-level data
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
# Education level 2
edu2 <- as.numeric(dat$EDUCATION==2)
# Education level 3
edu3 <- as.numeric(dat$EDUCATION==3)
# Read domain-level data
aux <- read.table("Nds20.txt", header=TRUE, sep = "\t", dec = ".")
# Prop. of registered people
aux$mreg <- aux$reg/aux$N
# Proportion of edu2 people
aux$medu2 <- aux$edu2/aux$N
# Proportion of edu3 people
aux$medu3 <- aux$edu3/aux$N



## We also define some new variables.

income.dir <- dir2(data=dat$INCOME, w=dat$WEIGHT, domain=list(sex=dat$SEX,
                                                              area=dat$AREA))
diry <- income.dir$mean         # Direct estimates of domain means
hatNd <- income.dir$Nd          # Direct estimates of population sizes
nd <- income.dir$nd             # Sample sizes
fd <- nd/aux$N                  # Sample fractions



## We calculate sample means by domains.

dat2 <- data.frame(income=dat$INCOME, edu2, edu3, reg=dat$REGISTERED)
smeans <- aggregate(dat2, by=list(sex=dat$SEX,area=dat$AREA), mean)
meany <- smeans$income            # Sample means of income
meanedu2 <- smeans$edu2           # Sample means of edu2
meanedu3 <- smeans$edu3           # Sample means of edu3
meanreg <- smeans$reg             # Sample means of registered



## We fit a random regression coefficient model with income as dependent variable

dat$EDUCATION <- as.factor(dat$EDUCATION)
rrc <- lmer(formula=INCOME ~ REGISTERED + EDUCATION + (EDUCATION|AREA:SEX),
            data=dat, REML=FALSE)
summary(rrc)                              # Summary of the fitting procedure
anova(rrc)                                # Analysis of Variance Table
beta <- fixef(rrc); beta                  # Regression parameters
var <- as.data.frame(VarCorr(rrc))        # Variance parameters
ref <- ranef(rrc)[[1]]                    # Modes of the random effects
head(fitted(rrc))                         # Predicted values
residuals <- resid(rrc)                   # Residuals
p.values <- 2*pnorm(abs(coef(summary(rrc))[,3]), low=F)
p.values                                  # p values



## We calculate the EBLUPs of income means by domain and the corresponding MSE estimators.

Xbeta <- beta[1] + beta[2]*aux$mreg + beta[3]*aux$medu2 + beta[4]*aux$medu3
Xubeta <- ref[,1] + aux$medu2*ref[,2] + aux$medu3*ref[,3]
mu <- Xbeta + Xubeta              # Projective estimates of income means
xbeta <- beta[1] + beta[2]*meanreg + beta[3]*meanedu2 + beta[4]*meanedu3
xubeta <- ref[,1] + meanedu2*ref[,2] + meanedu3*ref[,3]
mu.s <- meany - xbeta - xubeta
eb <- mu + fd*mu.s                # EBLUPs of income means



## Saving the results

output <- data.frame(Nd=aux$N[c(T,F)], hatNd=hatNd[c(T,F)] , nd=nd[c(T,F)],
                     meany=round(meany[c(T,F)],0), dir=round(diry[c(T,F)],0),
                     hatmu=round(mu[c(T,F)],0), eblup=round(eb[c(T,F)],0))
head(output, 10)





###
