



## This R code is a set of routines that correspond to chapter 2 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First we read the auxiliary and sample data files and rename some variables

dataux <- read.table("Nds20.txt", header=TRUE, sep = "\t", dec = ".")
# Sort dataux by sex and area:
dataux <- dataux[order(dataux$sex, dataux$area),]
# Sample data
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
# number of rows (cases) in dat:
n <- nrow(dat)
# Rename some variables
y1 <- dat$UNEMPLOYED; y2 <- dat$EMPLOYED
w <- dat$WEIGHT
area <- dat$AREA; sex <- dat$SEX



## To estimate the totals of unemployed and employed people we run

dir1.ds <- aggregate(w*data.frame(y1,y2), by=list(Area=area,Sex=sex), sum)
# Assign column names
names(dir1.ds) <- c("area", "sex", "y1tot", "y2tot")



## To estimate variances and the coefficients of variation we run

vardir1.ds <- aggregate(w*(w-1)*data.frame(y1^2,y2^2),
                        by=list(Area=area,Sex=sex), sum)
# Assign column names
names(vardir1.ds) <- c("area", "sex", "y1var", "y2var")
# Add columns y1var and y2var
dir1.ds <- cbind(dir1.ds, vardir1.ds$y1var, vardir1.ds$y2var)
# CV for y1
y1cv <- 100*sqrt(vardir1.ds$y1var)/abs(dir1.ds$y1tot)
# CV for y2
y2cv <- 100*sqrt(vardir1.ds$y2var)/abs(dir1.ds$y2tot)
# Add columns y1cv and y2cv
dir1.ds <- cbind(dir1.ds, y1cv, y2cv)
# Change column names for dir1.ds
namesds <- c("area", "sex", "y1tot", "y2tot", "y1var", "y2var", "y1cv",
             "y2cv")
names(dir1.ds) <- namesds



## We calculate the estimators of  means and their variances

# Add column with population sizes
dir1.ds <- cbind(dir1.ds, dataux$N)
# Add columns with HT estimates of means
dir1.ds <- cbind(dir1.ds, dir1.ds$y1tot/dataux$N,
                 dir1.ds$y2tot/dataux$N)
# Variance estimates of HT estimator
dir1.ds <- cbind(dir1.ds, dir1.ds$y1var/dataux$N^2, dir1.ds$y2var/
                   dataux$N^2)
# Change column names for dir1.ds
names(dir1.ds) <- c(namesds, "Nds", "y1mean", "y2mean", "y1meanvar",
                    "y2meanvar")



## To estimate unemployment rates we run

dirrate.ds <- 100*dir1.ds$y1tot/(dir1.ds$y1tot + dir1.ds$y2tot)
dir1.ds <- cbind(dir1.ds, rate=dirrate.ds)



## To estimate the covariances we run

covardir1.ds <- aggregate(w*(w-1)*data.frame(y1*y2),
                          by=list(Area=area,Sex=sex), sum)
# Column names
names(covardir1.ds) <- c("area", "sex", "covar")



## Calculation of variance of unemployment rate estimator

s1.ds <- dir1.ds$y2tot^2*dir1.ds$y1var/(dir1.ds$y1tot+dir1.ds$y2tot)^4
s2.ds <- dir1.ds$y1tot^2*dir1.ds$y2var/(dir1.ds$y1tot+dir1.ds$y2tot)^4
s12.ds <- 2*dir1.ds$y1tot*dir1.ds$y2tot*covardir1.ds$covar/
  (dir1.ds$y1tot + dir1.ds$y2tot)^4
# Estimates of variances and coefficients of variation
dir1.ds$vrate <- 10^4*(s1.ds+s2.ds-s12.ds)
dir1.ds$cvrate <- 100*sqrt(dir1.ds$vrate)/abs(dir1.ds$rate)



## Saving the results

output1 <- data.frame(dir1.ds[,1:6], rate=round(dirrate.ds,2))
head(output1, 10)




## Calculation of estimator of the totals of unemployed and employed people by AREA and SEX

dir <- aggregate(w*data.frame(1/w,1,y1,y2), by=list(Area=area,Sex=sex), sum)
# Column names
names(dir) <- c("area", "sex", "nds", "hatNds", "y1tot", "y2tot")



## We calculate the direct estimates of means

dir2.ds <- data.frame(area=dir$area, sex=dir$sex, nds=dir$nds,
                      hatNds=dir$hatNds)
# Estimates of means of unemployed people
dir2.ds$y1mean <- dir$y1tot/dir$hatNds
# Estimates of means of employed people
dir2.ds$y2mean <- dir$y2tot/dir$hatNds



## To estimate variances and the coefficients of variation we run

difference1 <- difference2 <- numerator1 <- numerator2 <- ww1 <- list()
for(d in 1:nrow(dir2.ds)){
  # Create a logic vector with the indexes of the corresponding domains
  condition <- paste(dat$AREA,dat$SEX,sep="")==paste(dir2.ds$area,
                                                     dir2.ds$sex,sep="")[d]
  # Calculate the difference between data and mean of each domain
  difference1[[d]] <- y1[condition]-dir2.ds$y1mean[d]
  difference2[[d]] <- y2[condition]-dir2.ds$y2mean[d]
  ww1[[d]] <- w[condition]*(w[condition]-1)
  numerator1[[d]] <- ww1[[d]]*difference1[[d]]^2
  numerator2[[d]] <- ww1[[d]]*difference2[[d]]^2
}



## Calculation of variances

dir2.ds$y1meanvar <- sapply(numerator1, sum)/dir2.ds$hatNds^2
dir2.ds$y2meanvar <- sapply(numerator2, sum)/dir2.ds$hatNds^2



## Calculation of CV's (in %)

# cv of y1-mean (in %)
dir2.ds$y1cv <- 100*sqrt(dir2.ds$y1meanvar)/abs(dir2.ds$y1mean)
# cv of y2-mean (in %)
dir2.ds$y2cv <- 100*sqrt(dir2.ds$y2meanvar)/abs(dir2.ds$y2mean)



## Calculation of estimators of means and their variances

dir2.ds$y1tot <- dir2.ds$y1mean*dataux$N
dir2.ds$y2tot <- dir2.ds$y2mean*dataux$N
dir2.ds$y1totvar <- dir2.ds$y1meanvar*dataux$N^2
dir2.ds$y2totvar <- dir2.ds$y2meanvar*dataux$N^2



## To estimate the unemployment rates we run

dir2.ds$rate <- 100*dir2.ds$y1tot/(dir2.ds$y1tot + dir2.ds$y2tot)



## To estimate the covariances we run

ww1s1s2 <- mapply(ww1, mapply(difference1, difference2, FUN="*"),
                  FUN="*")
sumcovardir2 <- sapply(ww1s1s2, sum)
covardir2.ds <- sumcovardir2*dataux$N^2/dir2.ds$hatNds^2



## Calculation of variance of the unemployment rate estimator

# Summands in formula of covariance estimator
s1.ds <- dir2.ds$y2tot^2*dir2.ds$y1totvar/(dir2.ds$y1tot+
                                             dir2.ds$y2tot)^4
s2.ds <- dir2.ds$y1tot^2*dir2.ds$y2totvar/(dir2.ds$y1tot+
                                             dir2.ds$y2tot)^4
s12.ds <- 2*dir2.ds$y1tot*dir2.ds$y2tot*covardir2.ds/
  (dir2.ds$y1tot+dir2.ds$y2tot)^4
# Estimates of variances and coefficients of variation
dir2.ds$vrate <- 10^4*(s1.ds+s2.ds-s12.ds)
dir2.ds$cvrate <- 100*sqrt(dir2.ds$vrate)/abs(dir2.ds$rate)



## Saving the results

output2 <- data.frame(dir2.ds[,1:2], round(dir2.ds[,11:14]),
                      rate=round(dir2.ds[,15],2))
head(output2, 10)




## Calculation of variances by Jackknife estimator
## We calculate some auxiliary parameters of the sample data file LFS20.txt

# Number of domains
D <- length(unique(dat$AREA))
# Domain sample sizes
nd <- tapply(rep(1,n),INDEX=list(dat$AREA),FUN=sum)
# Clusters
nCLUSTER <- unique(dat$CLUSTER)
# Number of clusters
J <- length(unique(dat$CLUSTER))
md <- vector()
# Number of clusters by domains
for (d in 1:D)
  md[d] <- length(unique(dat$CLUSTER[dat$AREA==d]))



## We calculate the direct estimates dir1, of the totals of unemployed and employed people

dir.d <- aggregate(w*data.frame(y1,y2), by=list(dat$AREA), sum)
# Assign column names
names(dir.d) <- c("area", "y1tot", "y2tot")



## Calculation of direct estimators of the variances

vardir.d <- aggregate(w*(w-1)*data.frame(y1^2,y2^2), by=list(dat$AREA), sum)
# Assign column names
names(vardir.d) <- c("area", "y1var", "y2var")



## The direct estimators of the coefficients of variations are

cvdir1 <- round(100*sqrt(vardir.d$y1var)/abs(dir.d$y1tot),2) # CV for y1
cvdir2 <- round(100*sqrt(vardir.d$y2var)/abs(dir.d$y2tot),2) # CV for y2



## For calculating the jackknife estimators we do the loop

jackdir1 <- jackdir2 <- matrix(0, nrow=D, ncol=J)
for (j in 1:J) {
  set <- subset(dat, dat$CLUSTER!=nCLUSTER[j], na.rm=TRUE)
  # Jackknife weights
  if (length(dat$AREA[dat$CLUSTER==j])>0) {
    domjack <- unique(dat$AREA[dat$CLUSTER==j])
    jfactor <- sum(dat$WEIGHT[dat$AREA==domjack])/
      sum(set$WEIGHT[set$AREA==domjack])
    set$WEIGHT[set$AREA==domjack] <- set$WEIGHT[set$AREA==domjack]*
      jfactor
  }
  # Direct estimators
  jdir.d <- aggregate(set$WEIGHT*data.frame(set$UNEMPLOYED,
                                            set$EMPLOYED), by=list(set$AREA), sum)
  # Assign column names
  names(jdir.d) <- c("area","y1tot","y2tot")
  jackdir1[,j] <- jdir.d$y1tot
  jackdir2[,j] <- jdir.d$y2tot
}



## We calculate the jackknife means.

jmeandir1 <- rowMeans(jackdir1)
jmeandir2 <- rowMeans(jackdir2)



## We calculate the jackknife variances and coefficients of variation.

# Number of clusters by jackknife domain
md.J <- list()
for (d in 1:D){
  md.J[[d]] <- md
  md.J[[d]][d] <- md.J[[d]][d]-1
}
factor <- Map(f="/", lapply(md.J,1,FUN="-"), md.J)
# Jackknife variances
diff.cuad.1 <- (jackdir1-jmeandir1)^2
diff.cuad.2 <- (jackdir2-jmeandir2)^2
group <- rep(1:D, md)
jvardir1 <- jvardir2 <- vector()   # declare objects for indexing
for (d in 1:D) {
  jvardir1[d] <- sum(sapply(split(diff.cuad.1[d,],group), sum)*factor[[d]])
  jvardir2[d] <- sum(sapply(split(diff.cuad.2[d,],group), sum)*factor[[d]])
}
# Jackknife coefficients of variation
jcvdir1 <- round(100*sqrt(jvardir1)/jmeandir1,2)
jcvdir2 <- round(100*sqrt(jvardir2)/jmeandir2,2)



## The R code to save the results is

output3 <- data.frame(nd, y1=dir.d$y1tot, v.y1=vardir.d$y1var,
                      vJ.y1=round(jvardir1), cv.y1=cvdir1, cvJ.y1=jcvdir1,
                      y2=dir.d$y2tot, v.y2=vardir.d$y2var,
                      vJ.y2=round(jvardir2), cv.y2=cvdir2, cvJ.y2=jcvdir2)
head(output3, 10)






## Here the code has been separated in two parts.
## This part contains the main program.
## Then, it is neccessary load the other part previously

source('SAE2_DesignDirect_Functions.R')



## We load the files and rename some variables

# Auxiliary data
dataux <- read.table("Nds20.txt", header=TRUE, sep = "\t", dec = ".")
# Sort dataux by sex and area:
dataux <- dataux[order(dataux$sex, dataux$area),]
# Sample data
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
# number of rows (cases) in dat:
n <- nrow(dat)
# Rename some  variables
y1 <- dat$UNEMPLOYED
w <- dat$WEIGHT



## Here we calculate HT and Hajek direct estimators by functions implemented ad-hoc

# Horvitz-Thompson direct estimator for unemployed people
direct1 <- dir1(data=y1, w=dat$WEIGHT, domain=list(area=dat$AREA,
                                                   sex=dat$SEX), Nd=dataux$N)
head(direct1, 10)
# Hajek direct estimator for unemployed people
direct2 <- dir2(data=y1, w=dat$WEIGHT, domain=list(area=dat$AREA,
                                                   sex=dat$SEX), Nd=dataux$N)
head(direct2, 10)





###
