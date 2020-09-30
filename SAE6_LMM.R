



## This R code is a set of functions that correspond to chapter 6 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First we read the sample data file and load libraries

# Auxiliary data
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
Ga <- cut(dat$AGE, breaks=c(0,25,54,1000), labels=c(1,2,3), right=TRUE)
dat$ageG <- as.numeric(Ga)                 # Age groups
dat$EDUCATION <- as.factor(dat$EDUCATION)  # EDUCATION as factor
# Install and load libraries
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}
if(!require(nlme)){
  install.packages("nlme")
  library(nlme)
}



## We fit the model LM1 to the data

lm.1 <- lm(formula=INCOME~AGE+EDUCATION, data=dat)
summary(lm.1)



## We fit the model LMM1 to the data

lmm.1 <- lmer(formula=INCOME~AGE+EDUCATION+(1|AREA), data=dat, REML=TRUE)
summary(lmm.1)                         # Summary of the fitting procedure
anova(lmm.1)                           # Analysis of Variance Table
fixef(lmm.1)                           # Regression parameters
var <- as.data.frame(VarCorr(lmm.1))   # Variance parameters
sigmau <- var$sdcor[1]; sigmau         # Random effect standard deviation
sigmae <- var$sdcor[2]; sigmae         # Residual standard deviation
ranef(lmm.1)                           # Modes of the random effects
lmm.1.fit <- fitted(lmm.1)             # Predicted values



## We fit the model LMM3 to the data

lmm.3 <- lmer(formula=INCOME~AGE+EDUCATION+(1|AREA:ageG), data=dat,
              REML=TRUE)
summary(lmm.3)                # Summary of the fitting procedure



## We fit the model LMM4 to the data

lmm.4 <- lmer(formula=INCOME~AGE+EDUCATION+(1|AREA/ageG), data=dat,
              REML=TRUE)
summary(lmm.4)      # Summary of the fitting procedure
anova(lmm.4)        # Analysis of Variance Table
fixef(lmm.4)        # Regression parameters
var <- as.data.frame(VarCorr(lmm.4)); var   # Variance parameters
sigmau2 <- var$sdcor[1]; sigmau2            # u_{2,dt} standard deviation
sigmau1 <- var$sdcor[2]; sigmau1            # u_{1,d} standard deviation
sigmae <- var$sdcor[3]; sigmae              # Residual standard deviation



## Now we are going to fit the models LMM1_LMM4 with the lme function of nlme package

# We fit the model LMM1 to the data

LMM.1 <- lme(INCOME~AGE+EDUCATION, random=~1|AREA, data=dat, method="REML")
summary(LMM.1)



## We fit the model LMM2 to the data

LMM.2 <- lme(INCOME~AGE+EDUCATION, random=~AGE|AREA, data=dat, method="ML")
summary(LMM.2)



## We fit the model LMM3 to the data

inter <- interaction(dat$AREA, dat$ageG)
LMM.3 <- lme(INCOME~AGE+EDUCATION, random=~1|inter, data=dat, method="REML")
summary(LMM.3)



## We fit the model LMM4 to the data

LMM.4 <- lme(INCOME~AGE+EDUCATION, random=~1|AREA/ageG, data=dat,
             method="REML")
summary(LMM.4)





###
