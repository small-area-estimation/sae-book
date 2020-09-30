



## This R code is a set of routines that correspond to chapter 5 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First read the auxiliary and sample data files and rename some variables

dataux <- read.table("Nds20.txt", header=TRUE, sep = "\t", dec = ".")
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
# Sample size
n <- nrow(dat); n
# Rename variables
y <- dat$INCOME; x <- dat$REGISTERED; w <- dat$WEIGHT; area <- dat$AREA



## Calculate the estimates of average incomes by area-sex by a regression
## synthetic estimator

# Assumed model
mod1 <- lm(y~x-1, weights=w)
# Model error variance
sigma12 <- anova(mod1)[2,3]
# Regression parameter
beta1 <- as.numeric(mod1$coefficients)
# Population size and sampling fraction
Npop <- sum(dataux$N); f <- n/Npop
# Totals of x by area-sex
Xtot.ds <- tapply(dataux$reg,list(dataux$area,dataux$sex),mean)
# Sizes by area sex
Ntot.ds <- tapply(dataux$N,list(dataux$area,dataux$sex),mean)
# Means of x by area-sex
Xmean.ds <- Xtot.ds/Ntot.ds
# Sample means of y by area-sex
ymean.ds <- tapply(y,list(dat$AREA,dat$SEX),mean)
# Sample means of x by area-sex
xmean.ds <- tapply(x,list(dat$AREA,dat$SEX),mean)
# Regression synthetic estimators of y-means by area-sex
yreg1.ds <- beta1*Xmean.ds
# BLUP estimators of y-means by area-sex
yblup1.ds <- (1-f)*yreg1.ds+f*(ymean.ds+(Xmean.ds-xmean.ds)*beta1)
yblup1.ds



## Calculate the estimates of average incomes by area-sex by using estimators
## without domain auxiliary information

# Assumed model
mod2 <- lm(y~x,weights=w)
# Model error variance
sigma22 <- anova(mod2)[2,3]
# Regression parameters
beta2 <- as.numeric(mod2$coefficients)
# Regression parameter for x
beta2x <- beta2[2]; beta2x
# Population size and sampling fraction
Npop <- sum(dataux$N); f <- n/Npop
# Totals of x by area-sex
Xtot.ds <- tapply(dataux$reg,list(dataux$area,dataux$sex),mean)
# Sizes by area sex
Ntot.ds <- tapply(dataux$N,list(dataux$area,dataux$sex),mean)
# Means of x by area-sex
Xmean.ds <- Xtot.ds/Ntot.ds
# Sample means of y by area-sex
ymean.ds <- tapply(y,list(dat$AREA,dat$SEX),mean)
# Sample means of x by area-sex
xmean.ds <- tapply(x,list(dat$AREA,dat$SEX),mean)
# Auxiliary terms for avoiding overflows
Ndir <- sum(w*10^(-3)); ydir <- sum(y*w*10^(-3)); xdir <- sum(x*w*10^(-3))
# Direct estimators of global means
ymeandir <- ydir/Ndir; xmeandir <- xdir/Ndir
# Projective estimator without domain auxiliary information
yreg2.ds <- ymeandir+(Xmean.ds-xmeandir)*beta2x; yreg2.ds
# BLUP estimators of domain means
yblup2.ds <- (1-f)*yreg2.ds+f*(ymean.ds+(Xmean.ds-xmean.ds)*beta2x)
yblup2.ds



## Calculate the estimates of average incomes by area-sex by using estimators
## with domain auxiliary information

# Assumed model
mod3 <- lm(y~as.factor(area)+x,weights=w)
# Model error variance
sigma32 <- anova(mod3)[3,3]
# Regression parameters
beta3 <- as.numeric(mod3$coefficients)
# Regression parameter for x
beta3x <- beta3[length(beta3)]
# Population size and sampling fraction
Npop <- sum(dataux$N); f <- n/Npop
# Totals of x by area-sex
Xtot.ds <- tapply(dataux$reg,list(dataux$area,dataux$sex),mean)
# Sizes by area sex
Ntot.ds <- tapply(dataux$N,list(dataux$area,dataux$sex),mean)
# Means of x by area-sex
Xmean.ds <- Xtot.ds/Ntot.ds
# Direct estimators of sizes
Ndir.ds <- tapply(w,list(dat$AREA,dat$SEX),sum)
# Direct estimators of y-totals
ydir.ds <- tapply(y*w,list(dat$AREA,dat$SEX),sum)
# Direct estimators of x-totals
xdir.ds <- tapply(x*w,list(dat$AREA,dat$SEX),sum)
# Direct estimators of y-means
ymeandir.ds <- ydir.ds/Ndir.ds
# Direct estimators of x-means
xmeandir.ds <- xdir.ds/Ndir.ds
# Projective estimator with domain auxiliary information
yreg3.ds <- ymeandir.ds+(Xmean.ds-xmeandir.ds)*beta3x; yreg3.ds
# BLUP estimators of domain means
yblup3.ds <- (1-f)*yreg3.ds+f*(ymean.ds+(Xmean.ds-xmean.ds)*beta3x)
yblup3.ds



## Save the results

output <- data.frame(reg1=yreg1.ds[,1], blup1=yblup1.ds[,1],
                     reg2=yreg2.ds[,1], blup2=yblup2.ds[,1],
                     reg3=yreg3.ds[,1], blup3=yblup3.ds[,1],
                     reg1=yreg1.ds[,2], blup1=yblup1.ds[,2],
                     reg2=yreg2.ds[,2], blup2=yblup2.ds[,2],
                     reg3=yreg3.ds[,2], blup3=yblup3.ds[,2]
                     )
head(output, 10)





###
