



## This R code is a set of routines that correspond to chapter 4 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First read the auxiliary and sample data files and rename some variables

# Survey data
dat <- read.table("LFS20.txt", header=TRUE, sep = "\t", dec = ".")
# Sample size
n <- nrow(dat); n
# Rename variables
y <- dat$INCOME; x <- dat$REGISTERED
# Auxiliary data
dataux <- read.table("Nds20.txt", header=TRUE, sep = "\t", dec = ".")



## Expansive estimation of the average income

mod1 <- lm(y~1)                              # Assumed model
sigma12 <- as.numeric(anova(mod1)[3])        # Model error variance
beta1 <- as.numeric(mod1$coefficients)       # Regression parameter
Npop <- sum(dataux$N)                        # Population size
f <- n/Npop; f                               # Sampling fraction
Mincome1 <- beta1                            # Expansive estimator
Mincome1; mean(y)                            # Checking
varMincome1 <- (1-f)*sigma12/n               # Estimator error variance



## Linear regression estimation of the average income
mod2 <- lm(y~x)                               # Assumed model
sigma22 <- anova(mod2)[2,3]                   # Model error variance
beta2 <- mod2$coefficients                    # Regression parameters
Npop <- sum(dataux$N)                         # Population size
f <- n/Npop; f                                # Sampling fraction
ymean <- mean(y); xmean <- mean(x)            # Sample means of y and x
Xmean <- sum(dataux$reg)/sum(dataux$N)        # Population mean of x
Mincome2 <- as.numeric(ymean+(Xmean-xmean)*   # Linear regression estimator
                         beta2[2])
xvar <- (n-1)*var(x)/n; xvar                  # Sample variance of x
varMincome2 <- ((1-f)*sigma22/n)*
  (1+((xmean-Xmean)^2/((1-f)*xvar))) # Estimator error variance



## The R code to save the results is

model1 <- c(beta1, NA, sigma12, Mincome1, varMincome1)
model2 <- c(beta2, sigma22, Mincome2, varMincome2)
labels <- c("intercept", "beta1", "sigma2", "Mincome", "Mincome variance")
output <- data.frame(labels, model1, model2)





###
