



## This R code is a set of routines that correspond to chapter 10 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## The following code loads the package lme4 and reads the auxiliary data file

if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}
aux <- read.table("auxLCS.txt", header=TRUE, sep = "\t", dec = ",")
# Number of employees by domain
aux$Twork <- round(aux$TOT*aux$Mwork,0)
# Number of unemployees by domain
aux$Tnowork <- round(aux$TOT*aux$Mnowork,0)
# Rest of labour cases by domain
aux$Tres <- aux$TOT - aux$Twork - aux$Tnowork
# Sort aux by dom
aux <- aux[order(aux$dom),]



## The following code reads the sample file, define some variables and calculate sample sizes.

dat <- read.table("datLCS.txt", header=TRUE, sep = "\t", dec = ",")
z0 <- 7280                          # poverty threshold.
income <- dat$income                # income variable
dom <- dat$dom; w <- dat$w          # sampling weight and domain
poor <- as.numeric(income<z0)       # poverty variable: 1 if yes (income<z0),
# 0 otherwise
gap <- (z0-income)*poor/z0          # gap variable
work <- as.numeric(dat$lab==1)      # employed
nowork <- as.numeric(dat$lab==2)    # unemployed
inact <- as.numeric(dat$lab==3)     # inactive
one <- rep(1, nrow(dat))            # variable one



## We calculate sample sizes by domain and labour status.

# Sample sizes by domains
nd <- tapply(X=one, INDEX=dom, FUN=sum)
# Sample sizes by domains and employment category
ndworking <- tapply(X=one, INDEX=list(dom,work), FUN=sum, default=0)
ndwork <- ndworking[,2]
# Sample sizes by domains and unemployment category
ndnoworking <- tapply(X=one, INDEX=list(dom,nowork), FUN=sum, default=0)
ndnowork <- ndnoworking[,2]
# By domains and innactive or <16 categories
ndres <- nd - ndwork - ndnowork



## We calculate some direct estimators and population sizes by domains.

# Sizes
hatNd <- tapply(X=w, INDEX=dom, FUN=sum)
# Totals of employed people
hatNwork <- tapply(X=work*w, INDEX=dom, FUN=sum)
# Totals of unemployed people
hatNnowork <- tapply(X=nowork*w, INDEX=dom, FUN=sum)
# Totals of innactive or <16 people
hatNdres <- hatNd - hatNwork - hatNnowork
dirptot <- tapply(X=poor*w, INDEX=dom, FUN=sum)
# Poverty proportions
dirp <- dirptot/hatNd
dirgtot <- tapply(X=gap*w, INDEX=dom, FUN=sum)
# Poverty gaps
dirg <- dirgtot/hatNd
diritot <- tapply(X=income*w, INDEX=dom, FUN=sum)
# Average incomes
diri <- diritot/hatNd



## We fit a nested error regression model and do some calculations

# Log transformation
c <- 10
y <- log(income+c)
y0 <- log(z0+c)
# Fitting the model
lmm <- lmer(formula=y ~ work+nowork+(1|dom), data = dat, REML = TRUE)
# Beta
beta <- fixef(lmm); beta
bwork <- beta[1] + beta[2]           # x1=1, x2=0
bnowork <- beta[1] + beta[3]         # x1=0, x2=1
bres <- beta[1]                      # x1=0, x2=0
# u (sorted by dom)
ud <- ranef(lmm)$dom; ud; dim(ud)
# mu
muwork <- bwork + ud
munowork <- bnowork + ud
mures <- bres + ud
# sigma
var <- as.data.frame(VarCorr(lmm))
sigmau2 <- var$vcov[1]
sigmae2 <- var$vcov[2]
# gammad
gammad <- sigmau2*nd/(sigmau2*nd+sigmae2)
# vd|s
vd <- sigmau2*(1-gammad) + sigmae2; vd
# alphad
alphadwork <- vd^(-1/2)*(y0-muwork)
alphadnowork <- vd^(-1/2)*(y0-munowork)
alphadres <- vd^(-1/2)*(y0-mures)



## We calculate the EBPs of poverty proportions by domains

# Normal CDF
nor1work <- pnorm(alphadwork[,1], mean=0, sd=1)
nor1nowork <- pnorm(alphadnowork[,1], mean=0, sd=1)
nor1res <- pnorm(alphadres[,1], mean=0, sd=1)
# Poverty sample totals
totp <- tapply(X=poor, INDEX=dom, FUN=sum)
# Poverty proportion EBPs
ebpp <- (totp + (aux$Twork-ndwork)*nor1work + (aux$Tnowork-ndnowork)*
           nor1nowork + (aux$Tres-ndres)*nor1res)/aux$TOT



## We calculate the EBPs of poverty gaps by domains

# Normal CDF
nor2work <- pnorm(alphadwork[,1]-vd^(1/2), mean=0, sd=1)
nor2nowork <- pnorm(alphadnowork[,1]-vd^(1/2), mean=0, sd=1)
nor2res <- pnorm(alphadres[,1]-vd^(1/2), mean=0, sd=1)
# Exponential terms
expwork <- exp(vd/2 + muwork)
expnowork <- exp(vd/2 + munowork)
expres <- exp(vd/2 + mures)
# Poverty gap summands
gapwork <- ((z0+c)/z0)*nor1work - expwork*nor2work/z0
gapnowork <- ((z0+c)/z0)*nor1nowork - expnowork*nor2nowork/z0
gapres <- ((z0+c)/z0)*nor1res - expres*nor2res/z0; gapres
# Poverty gap sample totals
totg <- tapply(X=gap, INDEX=dom, FUN=sum); totg
# Poverty gap EBPs
ebpg <- (totg+(aux$Twork-ndwork)*gapwork+(aux$Tnowork-ndnowork)*gapnowork +
           (aux$Tres-ndres)*gapres)/aux$TOT



## We calculate the EBPs of average incomes by domains.

# Income summands
incomework <- expwork-c
incomenowork <- expnowork-c
incomeres <- expres-c
# Income sample totals
toti <- tapply(X=income, INDEX=dom, FUN=sum); toti
# Average income EBPs
ebpi <- (toti + (aux$Twork-ndwork)*incomework + (aux$Tnowork-ndnowork)*
           incomenowork + (aux$Tres-ndres)*incomeres)/aux$TOT



## Saving the results

output <- data.frame(dom=aux$dom, nd, dirp=round(dirp,5), ebpp=round(ebpp,5),
                     dirg=round(dirg,5), ebpg=round(ebpg,5), diri=round(diri,0),
                     ebpi=round(ebpi,0))
output





###
