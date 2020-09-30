



## This R code is a set of routines that correspond to chapter 11 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we load the dir2 function from the second chapter

source('SAE2_DesignDirect_Functions.R')



## We read the auxiliary variables

aux <- read.table("Ndsa20.txt", header=TRUE, sep = "\t", dec = ".")
# sort auxiliary data by area (ascending), sex (ascending), age (ascending)
aux <- aux[order(aux$sex, aux$age, aux$area), ]



## We read the sample data files and rename some variables

dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
ns <- nrow(dat)    # global sample size
one <- rep(1,ns)   # auxiliary 1-variable
# Age groups
ageG <- as.numeric(cut(dat$AGE, breaks=c(0,25,54,max(dat$AGE)),
                       labels=c(1,2,3), right=TRUE))
# EDUCATION categories
edu2 <- as.numeric(dat$EDUCATION==2)
edu3 <- as.numeric(dat$EDUCATION==3)



## We calculate direct estimates of means

dir.income <- dir2(data=dat$INCOME, w=dat$WEIGHT, domain=list(area=dat$AREA,
                                                              age=ageG, sex=dat$SEX))
dir <- dir.income$mean
n <- dir.income$nd
py <- dir2(data=dat$INCOME, w=one, domain=list(area=dat$AREA, age=ageG,
                                               sex=dat$SEX))$mean
pedu2 <- dir2(data=edu2, w=one, domain=list(area=dat$AREA, age=ageG,
                                            sex=dat$SEX))$mean
pedu3 <- dir2(data=edu3, w=one, domain=list(area=dat$AREA, age=ageG,
                                            sex=dat$SEX))$mean



## We install and/or load some R packages.

if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}



## We fit a NER2 model with the R library lme4.

# Define new variables
edu <- as.factor(dat$EDUCATION)
sexage <- paste(dat$SEX, ageG, sep="")
# Apply lmer function
lmm <- lmer(formula=INCOME ~ edu + (1|AREA/sexage), data=dat, REML=TRUE)
# Summary of the fitting procedure
summary(lmm)
# Analysis of Variance Table
anova <- anova(lmm)
# Regression parameters
beta <- fixef(lmm); beta
# Variance parameters
var <- as.data.frame(VarCorr(lmm))
# Standard deviation of u2dt
sigmau2 <- var$sdcor[1]
# Standard deviation of u1d
sigmau1 <- var$sdcor[2]
# Standard deviation of edtj
sigmae <- var$sdcor[3]
# Modes of the random effects
ranef(lmm)
# Modes of subdomain random effects
udt <- ranef(lmm)[[1]]
# Modes of domain random effects
ud <- ranef(lmm)[[2]]
# Predicted values
ypred <- fitted(lmm)
# p values
p.values <- 2*pnorm(abs(coef(summary(lmm))[,3]), low=F); p.values



## We calculate the EBLUPs of the average incomes by areas and sex-age groups.

xbeta <- beta[1] + beta[2]*(aux$edu2/aux$N) + beta[3]*(aux$edu3/aux$N)
zu <- ud[,1] + udt[,1]
e1 <- xbeta + zu
fdt <- n/aux$N
xsbeta <- beta[1] + beta[2]*pedu2 + beta[3]*pedu3
e2 <- fdt*(py-xsbeta)
eblup <- e1+e2



## Saving the results

output <- data.frame(aux[1:4], n, dir, eblup)
head(output, 10)





###
