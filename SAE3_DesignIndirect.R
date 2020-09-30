



## This R code is a set of functions that correspond to chapter 3 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First read the auxiliary and sample data files

# Auxiliary data
dataux <- read.table("Ndsa20.txt", header=TRUE, sep = "\t", dec = ".")
# Sort dataux by sex, area and age
dataux <- dataux[order(dataux$area, dataux$sex, dataux$age),]
head(dataux)
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")



## We build the age and sex-age groups

Ga <- cut(dat$AGE, breaks=c(0,25,54,max(dat$AGE)), labels=FALSE, right=TRUE)
Gsa <- paste0(dat$SEX, Ga)
gsa <- as.numeric(factor(Gsa, levels=c("11","12","13","21","22","23"),
                         labels=c(1,2,3,4,5,6)))



## We fix the number of domains and groups.

D <- length(unique(dat$AREA))     # Number of areas
Ns <- length(unique(dat$SEX))     # Number of sex groups
Nga <- length(unique(Ga))         # Number of age groups
Ngsa <- length(unique(gsa))       # Number of sex-age groups



## We rename some variables.

y1 <- dat$UNEMPLOYED; y2 <- dat$EMPLOYED
w <- dat$WEIGHT; area <- dat$AREA; sex <- dat$SEX



## We put the estimated population sizes by area and sex-age group in a matrix of dimension D$\times$Ngsa.

hatNdg <- matrix(0, D, Ngsa)
for(i in 1:D) {
  for(j in 1:Ngsa) {
    hatNdg[i,j] <- sum(w[area==i & gsa==j])
  }
}
head(hatNdg)         # Some estimated sizes



## We do the same arrangement for the true population sizes.

Ndg <- matrix(0, D, Ngsa)
for(i in 1:D) {
  for(j in 1:Ns) {
    for(k in 1:Nga) {
      Ndg[i,(j-1)*Nga+k] <- dataux$N[dataux$area==i & dataux$sex==j &
                                       dataux$age==k]
    }
  }
}
head(Ndg)            # Some true sizes



## Calculates the direct estimates of means in sex-age groups

# Mean of unemployed people
Y1bar.g <- aggregate(w*y1, by=list(gsa=gsa), sum)[,2]/colSums(hatNdg)
# Mean of employed people
Y2bar.g <- aggregate(w*y2, by=list(gsa=gsa), sum)[,2]/colSums(hatNdg)



## Estimation of the variance of direct estimators of means by sex-age groups.

varY1bar.g <- varY2bar.g <- vector() # Define variance vectors
den <- colSums(hatNdg)^2
for(i in 1:Ngsa) {
  varY1bar.g[i] <- sum(w[gsa==i]*(w[gsa==i]-1)*(y1[gsa==i]-Y1bar.g[i])^2)/
    den[i]
  varY2bar.g[i] <- sum(w[gsa==i]*(w[gsa==i]-1)*(y2[gsa==i]-Y2bar.g[i])^2)/
    den[i]
}
# Variance estimates in sex-age groups
varY1bar.g; varY2bar.g



## We build a data.frame with the obtained results for sex-age groups.

gresults <- data.frame(g=1:Ngsa, Y1bar.g, Y2bar.g, varY1bar.g, varY2bar.g)
gresults



## The synthetic estimators of totals by area is

Synth.Y1 <- Ndg%*%Y1bar.g
Synth.Y2 <- Ndg%*%Y2bar.g



## Estimating the variance of synthetic estimators of totals by area

V.Synth.Y1 <- Ndg^2%*%varY1bar.g
V.Synth.Y2 <- Ndg^2%*%varY2bar.g



## The unemployment rates (in \%) by areas can be calculated as follows.

Synth.R.g <- Synth.Y1*100/(Synth.Y1 + Synth.Y2)



## Saving the results

dresults <- data.frame(d=1:D, Synth.Y1, Synth.Y2, V.Synth.Y1, V.Synth.Y2,
                       Synth.R.g)



## We calculate the synthetic estimators of totals and their estimated variances by area for men.

Synth.Y1.gM <- Ndg[,1:3]%*%Y1bar.g[1:3]
Synth.Y2.gM <- Ndg[,1:3]%*%Y2bar.g[1:3]
V.Synth.Y1.gM <- Ndg[,1:3]^2%*%varY1bar.g[1:3]
V.Synth.Y2.gM <- Ndg[,1:3]^2%*%varY2bar.g[1:3]



## Men unemployment rate (in \%) by area 

Synth.R.gM <- Synth.Y1.gM*100/(Synth.Y1.gM + Synth.Y2.gM)



## We build a data.frame with the obtained results.

SYNTH.M <- data.frame(area=1:D, y1tot=Synth.Y1.gM, y2tot=Synth.Y2.gM,
                        y1var=V.Synth.Y1.gM, y2var=V.Synth.Y2.gM,
                        rate=Synth.R.gM)



## We calculate the synthetic estimators of totals by area and their estimated variances for women.

Synth.Y1.gW <- Ndg[,4:6]%*%Y1bar.g[4:6]
Synth.Y2.gW <- Ndg[,4:6]%*%Y2bar.g[4:6]
V.Synth.Y1.gW <- Ndg[,4:6]^2%*%varY1bar.g[4:6]
V.Synth.Y2.gW <- Ndg[,4:6]^2%*%varY2bar.g[4:6]


## Women unemployment rate (in \%) by area

Synth.R.gW <- Synth.Y1.gW*100/(Synth.Y1.gW + Synth.Y2.gW)



## We build a data.frame with the obtained results.

SYNTH.W <- data.frame(area=1:D, y1tot=Synth.Y1.gW, y2tot=Synth.Y2.gW,
                      y1var=V.Synth.Y1.gW, y2var=V.Synth.Y2.gW,
                      rate=Synth.R.gW)



## Saving the results

output1 <- rbind(SYNTH.M, SYNTH.W)
head(output1, 10)



## We calculate the direct estimates of domain-group means

# Mean of unemployed people
Y1bar.dg <- aggregate(w*y1, by=list(area=area, gsa=gsa), sum)[,3]/hatNdg
# Mean of employed people
Y2bar.dg <- aggregate(w*y2, by=list(area=area, gsa=gsa), sum)[,3]/hatNdg



## We estimate the variances of direct estimators of means in domains crossed by sex-age groups

varY1bar.dg <- varY2bar.dg <- matrix(0, D, Ngsa)
den <- hatNdg^2
for(i in 1:D) {
  for(j in 1:Ngsa) {
    varY1bar.dg[i,j] <- sum(w[area==i&gsa==j]*(w[area==i&gsa==j]-1)
                            *(y1[area==i&gsa==j]-Y1bar.dg[i,j])^2)/den[i,j]
    varY2bar.dg[i,j] <- sum(w[area==i&gsa==j]*(w[area==i&gsa==j]-1)
                            *(y2[area==i&gsa==j]-Y2bar.dg[i,j])^2)/den[i,j]
  }
}



## Then we calculate post-stratified estimators of totals by area

# Post-stratified estimates of totals
Post.Y1 <- rowSums(Ndg*Y1bar.dg)
Post.Y2 <- rowSums(Ndg*Y2bar.dg)



## To estimate its variances, we can do

# Variance estimates in domains crossed by groups
V.Post.Y1 <- rowSums(Ndg^2*varY1bar.dg)
V.Post.Y2 <- rowSums(Ndg^2*varY2bar.dg)



## We calculate unemployment rates (in \%) by area 

Post.R.g <- Post.Y1*100/(Post.Y1 + Post.Y2)



## We build a data.frame with these results

dresults <- data.frame(area=1:D, Post.Y1, Post.Y2, V.Post.Y1, V.Post.Y2,
                       Post.R.g)



## We calculate the post-stratified estimators of totals by area, only for men.

# Post-stratified estimates of totals men
Post.Y1.gM <- rowSums(Ndg[,1:3]*Y1bar.dg[,1:3])
Post.Y2.gM <- rowSums(Ndg[,1:3]*Y2bar.dg[,1:3])



## We calculate the variance estimates only for men

# Variance estimates
V.Post.Y1.gM <- rowSums(Ndg[,1:3]^2*varY1bar.dg[,1:3])
V.Post.Y2.gM <- rowSums(Ndg[,1:3]^2*varY2bar.dg[,1:3])



## We calculate unemployment rates (in \%) by area, only for men

Post.R.gM <- Post.Y1.gM*100/(Post.Y1.gM + Post.Y2.gM)



## We build a data.frame only for men with these results

Post.area.M <- data.frame(area=1:D, Post.Y1.gM, Post.Y2.gM, V.Post.Y1.gM,
                          V.Post.Y2.gM, Post.R.gM)



## We calculate the post-stratified estimators of totals by area, only for women.

# Post-stratified estimates of totals women
Post.Y1.gW <- rowSums(Ndg[,4:6]*Y1bar.dg[,4:6])
Post.Y2.gW <- rowSums(Ndg[,4:6]*Y2bar.dg[,4:6])



## We calculate the variance estimates only for women

# Variance estimates
V.Post.Y1.gW <- rowSums(Ndg[,4:6]^2*varY1bar.dg[,4:6])
V.Post.Y2.gW <- rowSums(Ndg[,4:6]^2*varY2bar.dg[,4:6])



## We calculate unemployment rates (in \%) by area, only for women

Post.R.gW <- Post.Y1.gW*100/(Post.Y1.gW + Post.Y2.gW)



## We build a data.frame only for women with these results

Post.area.W <- data.frame(area=1:D, Post.Y1.gW, Post.Y2.gW, V.Post.Y1.gW,
                          V.Post.Y2.gW, Post.R.gW)



## Saving the results

output2 <- cbind(Post.area.M, Post.area.W)
head(output2, 10)



## For the GREG estimator, first we read the data files

dataux <- read.table("Nds20.txt", header=TRUE, sep = "\t", dec = ".")
# Auxiliary data: population sizes by area and sex
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
n <- nrow(dat)                       # Global sample size
narea <- length(unique(dataux$area)) # Number of areas
nsex <- length(unique(dataux$sex))   # Number of sex categories



## We rename some variables

y <- dat$INCOME; x1 <- dat$REGISTERED; x2 <- dat$EDUCATION;
w <- dat$WEIGH; area <- dat$AREA; sex <- dat$SEX
edu2 <- as.numeric(x2==2); edu3 <- as.numeric(x2==3)



## We calculate the estimator of the regression parameter

X <- cbind(1, x1, edu2, edu3)              # Matrix with auxiliary variables
p <- ncol(X)
Y <- matrix(y, nrow=n)                     # Vector with target variable
W <- diag(w)                               # Matrix with sampling weights
Q <- solve( crossprod(X,W)%*%X )
beta <- tcrossprod(Q,X) %*% crossprod(W,Y) # Estimator of beta parameter



## Alternatively

mod <- lm(y ~ x1 + edu2 + edu3, weights=w)
coef(mod)



## We calculate direct estimates of sizes and means by area-sex.

# Estimated size
hatNds <- aggregate(w, by=list(sex=sex, area=area), sum)
# x1 means
dir.reg <- aggregate(w*x1, by=list(sex=sex, area=area), sum)[,3]/hatNds[,3]
# edu2 means
dir.edu2 <- aggregate(w*edu2, by=list(sex=sex, area=area), sum)[,3]/
  hatNds[,3]
# edu3 means
dir.edu3 <- aggregate(w*edu3, by=list(sex=sex, area=area), sum)[,3]/
  hatNds[,3]
# y means
dir.income <- aggregate(w*y, by=list(sex=sex, area=area), sum)[,3]/hatNds[,3]
hatdir <- data.frame(area=hatNds[,2], sex=hatNds[,1], N=hatNds[,3],
                     reg=dir.reg, edu2=dir.edu2, edu3=dir.edu3,
                     income=dir.income)



## For estimating the GREG estimator of averages by area-sex, we do

# Direct estimates of y-mean
Ymean.dir <- as.matrix(hatdir$income)
# Direct estimates of X-mean
Xmean.dir <- as.matrix(cbind(1, hatdir[,4:6]))
# True X-means
Xmean <- cbind(1, dataux$reg/dataux$N, dataux$edu2/dataux$N,
               dataux$edu3/dataux$N)
# GREG estimates of y-means
Ymean.greg <- Ymean.dir + (Xmean-Xmean.dir)%*%beta



## For estimating the variance of direct estimators of means by area-sex,
## we apply the function dir2 of the second chapter

source('SAE2_DesignDirect_Functions.R')
dir2.est <- dir2(y, w, domain=list(sex=sex, area=area), Nd=dataux$N)
Yvar.dir <- dir2.est$var.mean



## We calgulate de g-weights

g <- vector()
for (k in 1:n) {
  condition <- dataux$area==dat$AREA[k]&dataux$sex==dat$SEX[k]
  Nk <- dataux$N[condition]
  regk <- dataux$reg[condition]/Nk
  edu2k <- dataux$edu2[condition]/Nk
  edu3k <- dataux$edu3[condition]/Nk
  Xk <- c(1, regk, edu2k, edu3k)
  condition2 <- hatdir$area==dat$AREA[k]&hatdir$sex==dat$SEX[k]
  Nhatk <- hatdir$N[condition2]
  reghatk <- hatdir$reg[condition2]
  edu2hatk <- hatdir$edu2[condition2]
  edu3hatk <- hatdir$edu3[condition2]
  Xdirk <- c(1, reghatk, edu2hatk ,edu3hatk)
  g[k] <- Nk/Nhatk + Nk*(Xk-Xdirk)%*%Q%*%X[k,]
}



## For estimating the variance, we do

yyds <- (y-fitted(mod))^2
# (alternatively fitted(mod) can be computed as X%*%beta)
vargreg.ds <- aggregate(w*(w-1)*g^2*yyds, by=list(sex, area), sum)
Yvar.greg <- vargreg.ds[,3]/dataux$N^2



## We calculate the estimated coefficients of variation

cvdir <- round(100*sqrt(Yvar.dir)/as.vector(Ymean.dir),2)
cvgreg <- round(100*sqrt(Yvar.greg)/as.vector(Ymean.greg),2)



## Saving the results

output3 <- data.frame(area=dataux$area, sex=dataux$sex, dir=Ymean.dir,
                      GREG=Ymean.greg, dirvar=Yvar.dir, GREGvar=Yvar.greg)
head(round(output3), 10)





###
