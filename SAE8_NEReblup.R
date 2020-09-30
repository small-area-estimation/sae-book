



## This R code is a set of functions that correspond to chapter 8 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First we read the sample data file and load libraries

if(!require(sae)){
  install.packages("sae")
  library(sae)
}
dataux <- read.table("Nds20.txt", header=TRUE, sep="\t", dec=".")
dat <- read.table("LFS20.txt", header=TRUE, sep="\t", dec=".")
dat$AREASEX <- paste(dat$AREA, dat$SEX, sep="s") # Domains
n <- dim(dat)[1]                        # Global sample size
dat <- dat[order(dat$AREA, dat$SEX),]   # Sort dat by area and sex
narea <- length(unique(dataux$area))    # Number of areas
nsex <- length(unique(dataux$sex))      # Number of sex categories



## We rename some variables.

y <- dat$INCOME; x1 <- dat$REGISTERED; x2 <- dat$EDUCATION
area <- dat$AREA; sex<-dat$SEX; AREASEX <- dat$AREASEX
edu2 <- as.numeric(x2==2)    # Category 2 of EDUCATION
edu3 <- as.numeric(x2==3)    # Category 3 of EDUCATION



## We construct a data frame with the auxiliary variables and domain sizes.

areasex <- paste(dataux$area, dataux$sex, sep="s")
Xmean <- data.frame(areasex, reg=dataux$reg/dataux$N,
                    edu2=dataux$edu2/dataux$N,
                    edu3=dataux$edu3/dataux$N)
Popn <- data.frame(areasex, N=dataux$N)



## We fit a NER model and we calculate the EBLUPs of domain-sex average incomes.

result <- eblupBHF(y~x1+edu2+edu3, dom=dat$AREASEX, meanxpop=Xmean,
                   popnsize=Popn, method="ML")
result$eblup[1:10,]             # Domains, EBLUPs and sample sizes
eblup <- result$eblup[,2]       # EBLUPS
# Additional results
result$fit$summary
# Regression parameters of the fitted NER models
result$fit$fixed
result$fit$random               # Domain random effects
result$fit$errorvar             # Error variance
result$fit$refvar               # Random effect variance
result$fit$loglike              # Log-likelihood estimated value



## We calculate the parametric bootstrap MSEs of the EBLUPs of domain-sex average incomes.

set.seed(123)
msey <- pbmseBHF(y~x1+edu2+edu3, dom=dat$AREASEX, meanxpop=Xmean,
                 popnsize=Popn, B=500, method="ML")
msey$mse                                       # Domains and MSEs
mseEBLUP <- msey$mse[,2]                       # MSEs
cvEBLUP <- round(100*sqrt(mseEBLUP)/eblup,2)   # Coefficients of variation



## Saving the results

output1 <- data.frame(area=result$eblup[1:10,1], eblup=eblup, mse=msey$mse[,2],
                      cv=cvEBLUP)
output1



## First we read the sample data file and load libraries

if(!require(sae)){
  install.packages("sae")
  library(sae)
}
source('SAE2_DesignDirect_Functions.R')
aux <- read.table("auxLCS.txt", header=TRUE, sep = "\t", dec = ",")
dat <- read.table("datLCS.txt", header=TRUE, sep = "\t", dec = ",")
aux <- aux[order(aux$dom),]                  # sort the file by domains
D <- nrow(aux)                               # number of domains



## We define some variables

income <- dat$income; w <- dat$w; dom <- dat$dom
work <- as.numeric(dat$lab==1)
nowork <- as.numeric(dat$lab==2)
inact <- as.numeric(dat$lab==3)



## We calculate the direct Hajek estimators of domain average incomes

income.dir <- dir2(data=income, w, domain=list(dom=dom))
diry <- income.dir$mean; hatNd <- income.dir$Nd.hat; nd <- income.dir$nd



## We calculate the direct estimators of proportions

work.dir <- dir2(data=work, w, domain=list(dom=dom))
dirwork <- work.dir$mean
nowork.dir <- dir2(data=nowork, w, domain=list(dom=dom))
dirnowork <- nowork.dir$mean



## We calculate domain sizes and $x$-means.

Nd <- aux$TOT                                # domain sizes
round(hatNd[1:10],0)                         # estimated domain sizes
# Dataframes needed for applying the function eblupBHF of the library sae
Xmean <- data.frame(dom=aux$dom, work=aux$Mwork, nowork=aux$Mnowork)
PopN <- data.frame(dom=aux$dom, Nd)



## We fit a NER model and we calculate the EBLUPs of income means by domains.

EB <- eblupBHF(income~work+nowork, dom=dom, meanxpop=Xmean, popnsize=PopN,
               method="REML", data=dat)
anova <- EB$fit$summary                      # ANOVA Table
2*pnorm(abs(anova$coefficients[,3]), low=F)  # p-values
beta <- EB$fit$fixed                         # betas
sigmae2 <- EB$fit$errorvar                   # Error variance
sigmau2 <- EB$fit$refvar                     # Random effect variance
head(EB$fit$random, 10)                      # Random effects
Residuals <- EB$fit$residuals                # Residuals
EBLUP <- EB$eblup[order(EB$eblup$domain),]   # Sort EB$EBLUP by dom
# Inserts EB$fit$random in EBLUP
EBLUP$random <- EB$fit$random[,1]
eblup <- EBLUP$eblup                         # EBLUPs sorted by domains
# Random effects sorted by domains
ud <- EBLUP$random



## We calculate the parametric bootstrap MSE of the  EBLUP.

set.seed(123)
MSEeblup <- pbmseBHF(income~work+nowork, dom=dom, meanxpop=Xmean,
                     popnsize=PopN, B=500, method="REML", data=dat)
# Sort MSEeblup$mse by dom
EBMSE <- MSEeblup$mse[order(MSEeblup$mse$domain),]
ebmse <- EBMSE$mse                             # Model-based MSEs of EBLUPS
mCVeb <- round(100*sqrt(ebmse)/eblup,2)        # Model-based CVs of EBLUPS



## We calculate the model-assisted counterpart of the EBLUP.

s1 <- beta[1] + beta[2]*Xmean$work + beta[3]*Xmean$nowork + ud; s1
s2 <- hatNd*(diry-beta[1]-beta[2]*dirwork-beta[3]*dirnowork-ud)/Nd; s2
MA <- s1+s2



## Saving the results

output2 <- data.frame(dom=aux$dom, Nd, hatNd=round(hatNd,0), nd,
                      DIR=round(diry,0), MA=round(MA,0), EB=round(eblup,0),
                      mMSEeb=round(ebmse,0), mCVeb)
head(output2, 10)





###
