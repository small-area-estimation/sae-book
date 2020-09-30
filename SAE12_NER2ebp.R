



## This R code is a set of routines that correspond to chapter 12 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we load the dir2 function from the second chapter

source('SAE2_DesignDirect_Functions.R')



## We read the unit-level data file and defines some variables

if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
z0 <- 36500                                 # poverty threshold
ns <- nrow(dat)                             # global sample size
poor <- as.numeric(dat$INCOME<z0)           # variable poor
gap <- (z0-dat$INCOME)*poor/z0              # variable gap
one <- rep(1,nrow(dat))                     # variable one
Ga <- cut(dat$AGE, breaks=c(0,25,54,max(dat$AGE)), labels=c(1,2,3),
          right=TRUE)
ageG <- as.numeric(Ga)                      # age group
edu2 <- as.numeric(dat$EDUCATION==2)        # secondary education
edu3 <- as.numeric(dat$EDUCATION==3)        # superior education
y <- log(dat$INCOME)                        # variable y=log(income)



## We read the auxiliary variables

aux <- read.table("Ndsa20.txt", header=TRUE, sep = "\t", dec = ".")
# Sort aux by sex, age and area
aux <- aux[order(aux$sex, aux$age, aux$area), ]



## We calculate sample sizes, counts and means by subdomains

# Sizes
ndt <- tapply(X=one, INDEX=list(dat$AREA,ageG), FUN=sum)
# Sample counts of edu3
ndtedu3 <- tapply(edu3, list(dat$AREA,ageG), sum)
# Sample counts of edu2
ndtedu2 <- tapply(edu2, list(dat$AREA,ageG), sum)
# Sample counts of edu1
ndtedu1 <- ndt - ndtedu3 - ndtedu2
# Sample counts of poor people
ndtpoor <- tapply(poor, list(dat$AREA,ageG), sum)
# Sample poverty proportion
mdtpoor <- ndtpoor/ndt
# Sample sum of gap variable
ndtgap <- tapply(gap, list(dat$AREA,ageG), sum)
# Sample poverty gap
mdtgap <- ndtgap/ndt
# Sample sum of log-income
ndty <- tapply(y, list(dat$AREA,ageG), sum)
# Sample log-income mean
mdty <- ndty/ndt



## We calculate direct estimators

dir.poor <- dir2(data=poor, w=dat$WEIGHT, domain=list(area=dat$AREA,
                                                      ageG=ageG))
dir.gap <- dir2(data=gap, w=dat$WEIGHT, domain=list(area=dat$AREA,
                                                    ageG=ageG))
dir.income <- dir2(data=dat$INCOME, w=dat$WEIGHT, domain=list(area=dat$AREA,
                                                              ageG=ageG))
hatNdt <- dir.poor$Nd.hat   # sizes
dirp <- dir.poor$mean       # poverty proportions
dirg <- dir.gap$mean        # poverty gaps
diri <- dir.income$mean     # income means



## We fit a two-fold nested error regression model by applying REML method

lmm <- lmer(formula=y ~ edu3 + edu2 + (1|AREA/ageG), data=dat, REML=TRUE)
summary(lmm)                             # summary of the fitting procedure
anova(lmm)                               # analysis of variance table
beta <- fixef(lmm)                       # regression parameters
bedu3 <- beta[1] + beta[2]               # beta for x1=1, x2=1, x3=0 (edu3)
bedu2 <- beta[1] + beta[3]               # beta for x1=1, x2=0, x3=1 (edu2)
bedu1 <- beta[1]                         # beta for x1=1, x2=0, x3=0 (edu0)
var <- as.data.frame(VarCorr(lmm))       # variance parameters
sigmau2 <- var$sdcor[1]                  # standard deviation of u_{2,dt}
sigmau1 <- var$sdcor[2]                  # standard deviation of u_{1,d}
sigmae <- var$sdcor[3]; sigmae           # residual standard deviation
ranef(lmm)                               # modes of the random effects
ypred <- fitted(lmm)                     # predictions
residuals <- resid(lmm)                  # residuals
p.values <- 2*pnorm(abs(coef(summary(lmm))[,3]), low=F); p.values # p values



## We calculate the means and the variances of unobserved values 

# Calculation of gammadt, gammadt, by subdomains
gammadt <- sigmau2^2*ndt/(sigmau2^2*ndt+sigmae^2)
# Calculation of deltad
gammad <- apply((1-gammadt)*ndt, 1, sum)
# phid by domains
phid <- sigmau1^2/(sigmae^2+sigmau1^2*gammad)
# Calculation of the conditioned variances, vdt, by subdomains
vdt <- sigmae^2*(1+phid*(1+gammadt*(gammadt-2))) + sigmau2^2*(1-gammadt)



## We do some calculations previous to EBPs

# Preliminary calculations
mdtxbeta <- (bedu3*ndtedu3 + bedu2*ndtedu2 + bedu1*ndtedu1)/ndt
gammayxd <- apply(gammadt*(mdty-mdtxbeta), 1, sum)
uudt <- gammadt*(mdty-mdtxbeta+(sigmae^2/sigmau2^2)^2*(phid/ndt)*gammayxd)
# Calculation of the conditioned means
muedu3 <- bedu3 + uudt; muedu2 <- bedu2 + uudt; muedu1 <- bedu1 + uudt
# alphadt
y0 <- log(z0)
alphadedu3 <- vdt^(-1/2)*(y0-muedu3)
alphadedu2 <- vdt^(-1/2)*(y0-muedu2)
alphadedu1 <- vdt^(-1/2)*(y0-muedu1)



## We calculate the EBPs of the poverty proportions by subdomains

# Normal CDF
noredu3 <- pnorm(alphadedu3, mean=0, sd=1)
noredu2 <- pnorm(alphadedu2, mean=0, sd=1)
noredu1 <- pnorm(alphadedu1, mean=0, sd=1)
# Populations sizes by subdomains
Ndt <- tapply(X=aux$N, INDEX=list(aux$area,aux$age), FUN=sum)   # global
Nedu3 <- tapply(aux$edu3, list(aux$area,aux$age), sum)          # edu3
Nedu2 <- tapply(aux$edu2, list(aux$area,aux$age), sum)          # edu2
Nedu1 <- Ndt - Nedu3 - Nedu2                                    # edu1
# Poverty proportion EBPs by subdomains
ebpptot <- ndtpoor + (Nedu3-ndtedu3)*noredu3 + (Nedu2-ndtedu2)*noredu2 +
  (Nedu1-ndtedu1)*noredu1
ebp.poor <- ebpptot/Ndt; ebp.poor



## We calculate the poverty gaps of subdomains

gap3 <- noredu3 - exp(vdt/2+muedu3)*pnorm(alphadedu3-vdt^{1/2})/z0
gap2 <- noredu2 - exp(vdt/2+muedu2)*pnorm(alphadedu2-vdt^{1/2})/z0
gap1 <- noredu1 - exp(vdt/2+muedu1)*pnorm(alphadedu1-vdt^{1/2})/z0
# Poverty gap EBPs by subdomains
ebpgtot <- ndtgap + (Nedu3-ndtedu3)*gap3 + (Nedu2-ndtedu2)*gap2 +
  (Nedu1-ndtedu1)*gap1
ebp.gap <- ebpgtot/Ndt; ebp.gap



## The following R code calculates the income means by subdomains.

inc3 <- exp(vdt/2+muedu3)
inc2 <- exp(vdt/2+muedu2)
inc1 <- exp(vdt/2+muedu1)
# EBPs of income means by subdomains
ebpitot <- ndty + (Nedu3-ndtedu3)*inc3 + (Nedu2-ndtedu2)*inc2 +
  (Nedu1-ndtedu1)*inc1
ebpi <- ebpitot/Ndt; ebpi



## Saving the results

# Poverty proportions
output1 <- data.frame(n1=ndt[,1], n2=ndt[,2], n3=ndt[,3], 
                      ebpp1=round(ebp.poor[,1],4), ebpp2=round(ebp.poor[,2],4),
                      ebpp3=round(ebp.poor[,3],4),
                      dirp1=round(subset(dir.poor, dom.ageG=="1")$mean,4),
                      dirp2=round(subset(dir.poor, dom.ageG=="2")$mean,4),
                      dirp3=round(subset(dir.poor, dom.ageG=="3")$mean,4))
head(output1, 10)

# Poverty gaps
output2 <- data.frame(n1=ndt[,1], n2=ndt[,2], n3=ndt[,3],
                      ebpg1=round(ebp.gap[,1],4), ebpg2=round(ebp.gap[,2],4),
                      ebpg3=round(ebp.gap[,3],4),
                      dirg1=round(subset(dir.gap, dom.ageG=="1")$mean,4),
                      dirg2=round(subset(dir.gap, dom.ageG=="2")$mean,4),
                      dirg3=round(subset(dir.gap, dom.ageG=="3")$mean,4))
head(output1, 10)

# Results for income means

output3 <- data.frame(n1=ndt[,1], n2=ndt[,2], n3=ndt[,3],
                      ebpi1=round(ebpi[,1],0), ebpi2=round(ebpi[,2],0),
                      ebpi3=round(ebpi[,3],0),
                      diri1=round(subset(dir.income, dom.ageG=="1")$mean,0),
                      diri2=round(subset(dir.income, dom.ageG=="2")$mean,0),
                      diri3=round(subset(dir.income, dom.ageG=="3")$mean,0))
head(output3, 10)





###
