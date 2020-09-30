



# This R code is a set of routines that correspond to chapter 1 of the book 
# A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
# www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we read the data file LFS20.txt and calculate some values

dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")

datcol <- ncol(dat); datcol           # number of columns (variables) in dat
n <- nrow(dat)                        # number of rows (cases) in dat
narea <- length(unique(dat$AREA))     # number of areas
nsex <- length(unique(dat$SEX))       # number of sex categories
output1 <- data.frame(area=dat$AREA, cluster=dat$CLUSTER, age=dat$AGE,
                      sex=dat$SEX, weight=dat$WEIGHT, employed=dat$EMPLOYED,
                      unemployed=dat$UNEMPLOYED, inactive=dat$INACTIVE,
                      registered=dat$REGISTERED, education=dat$EDUCATION,
                      income=dat$INCOME)
head(output1, 10)                     # some survey data



## We calculate some variables and tables

summary(dat$AGE)
AGEG <- cut(dat$AGE, breaks=c(0,25,54,100), labels=c(1,2,3))
AGEG <- as.numeric(AGEG)          # variable age group
table(AGEG, dat$AREA, dat$SEX)    # sample sizes per age group, area and sex
table(dat$SEX, dat$AREA)          # sample sizes per sex and area



## We read the file Ndsa20.txt that contains auxiliary aggregated data by area, sex and age

dataux <- read.table("Ndsa20.txt", header=TRUE, sep = "\t", dec = ".")
# Population sizes by area, sex and age group
output2 <- data.frame(area=dataux$area, sex=dataux$sex, age=dataux$age,
                      N=dataux$N, reg=dataux$reg, edu1=dataux$edu1, edu2=dataux$edu2,
                      edu3=dataux$edu3)
head(output2, 10)
dim(dataux)           # file dimensions



## We read the file Nds20.txt that contains auxiliary aggregated data by area and sex

dataux <- read.table("Nds20.txt", header=TRUE, sep = "\t", dec = ".")
# Population sizes by area and sex
output3 <- data.frame(area=dataux$area, sex=dataux$sex, N=dataux$N,
                      reg=dataux$reg, edu1=dataux$edu1, edu2=dataux$edu2,
                      edu3=dataux$edu3)
head(output3, 10)
dim(dataux)           # file dimensions



## We read the file datLCS.txt that contains unit-level data from a living conditions survey

dat <- read.table("datLCS.txt", header=TRUE, sep = "\t", dec = ",")
dat <- dat[order(dat$dom),]          # sort dat by dom
datcol <- ncol(dat)                  # number of columns (variables) in dat
n <- nrow(dat)                       # number of rows (cases) in dat
ndom <- length(unique(dat$dom))      # number of domains in dat
output4 <- data.frame(dom=dat$dom, sex=dat$sex, house=dat$house, w=dat$w,
                      income=dat$income, lab=dat$lab)
head(output4, 10)                    # some survey data
table(dat$sex, dat$dom)              # sample sizes per dom and sex



## We read the file Nds20.txt that contains auxiliary aggregated data by domain

aux <- read.table("auxLCS.txt", header=TRUE, sep = "\t", dec = ",")
aux <- aux[order(aux$dom),]      # sort aux by dom
output5 <- data.frame(dom=aux$dom, TOT=aux$TOT, Mwork=aux$Mwork,
                      Mnowork=aux$Mnowork, Minact=aux$Minact, ss=aux$ss)
head(output5,10)
dim(aux)              # file dimensions





###
