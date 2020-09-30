



## This R code is a set of routines that correspond to chapter 21 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, install and/or load the package lme4.

if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}



## The following R codes simulates the area-level temporal count data.

# Number of domains and time instants
D <- 50; T <- 4
d1 <- seq(1:D); dd <- rep(d1,T)              # Domain index
# Regression parameters
beta0 <- -3; beta1 <- 0.8
beta <- c(beta0, beta1)
phi1 <- phi2 <- 0.5                          # Standard deviation parameters
set.seed(123)                                # Set the simulation seed
vv1 <- rnorm(D, 0, 1); v1 <- rep(vv1, T)     # Domain random effects
v2 <- rnorm(D*T, mean=0, sd=1)               # Domain-time random effects
m <- rep(100, D*T)                           # Offset parameter
# Time index and x-variables
x <- tt <- vector()
for (t in 1:T){
  for (d in 1:D){
    x[(t-1)*D+d] <- (d + t/T)/D
    tt[(t-1)*D+d] <- t
  }
}
# Generation of y
eta <- beta0 + beta1*x + phi1*v1 + phi2*v2   # Natural parameter
p <- exp(eta)                                # Proportion parameter
mu <- m*p                                    # Mean parameter
y <- rpois(D*T, mu)                          # Poisson counts



## Fit an area-level Poisson mixed model with independent time effects to the simulated sample counts.

glmm <- glmer(formula=y ~ x + (1|dd/tt), family=poisson, offset=log(m))
summary(glmm)
hbeta <- fixef(glmm)                   # Estimated regression parameters
# Estimated variances and standard deviations
var <- as.data.frame(VarCorr(glmm))
hphi2 <- var$sdcor[1]                  # Estimated standard deviation of udt
hphi1 <- var$sdcor[2]                  # Estimated standard deviation of ud



## Calculate the modes of the random effects, the plug-in predictors
## of the sample counts and the model residuals.

r.effects <- ranef(glmm)     # Modes of random effects
udt <- r.effects[[1]]        # Modes of udt
ud <- r.effects[[2]]         # Modes of ud
mutilde <- fitted(glmm)      # Plug-in predictions of counts
residuals <- resid(glmm)     # Model residuals
# Summary of results
result <- data.frame(d=dd, t=tt, m, x, mu, y, mutilde)



## Calculate calculate the EBPs of sample counts and add the output to result

set.seed(123)
S1 <- S2 <- 100
for (t in 1:T){                                 # Loop for time instants
  for (d in 1:D){                               # Loop for domains
    numdt <- dendt <- 0
    for (s1 in 1:S1) {                          # Loop for outer sum
      v1 <- rnorm(n=1, mean=0, sd=1)
      v1 <- c(v1, -v1)
      numdtau <- dendtau <- 1
      for (tau in 1:T){                         # Loop for product
        condition <- tt==tau & dd==d
        vv2 <- rnorm(n=S2, mean=0, sd=1)
        vv2 <- c(vv2, -vv2)
        ydtau <- y[condition]
        deltau <- as.numeric(t==tau)
        xbetau <- hbeta[1] + hbeta[2]*x[condition] + hphi1*v1 + hphi2*vv2
        sum1num <- (ydtau + deltau)*xbetau
        sum2 <- m[condition] * exp(xbetau)
        # Numerator term
        numdtau <- numdtau * sum(exp(sum1num-sum2))
        sum1den <- ydtau * xbetau
        # Denominator term
        dendtau <- dendtau * sum(exp(sum1den-sum2))
      }
      # EBP numerator
      numdt <- numdt + numdtau
      # EBP denominator
      dendt <- dendt + dendtau
    }
    condition2 <- tt==t&dd==d
    # EBP
    result$ebp[condition2] <- m[condition2]*numdt/dendt
  }
}



## Calculate the parametric bootstrap estimates of the root-MSEs of the EBPs of sample counts.

set.seed(123); B <- 200
mseebppcum.star <- 0
dif.ebp.star <- dif.plugin.star <- list()
ebp.star <- plugin.star <- list()
for (b in 1:B){
  cat("iteration = ", b, "\n")
  # Generation of y
  v1.star <- rnorm(n=D, mean=0, sd=1)
  v2.star <- rnorm(n=D*T, mean=0, sd=1)
  p.star <- exp(hbeta[1] + hbeta[2]*x + hphi1*v1.star + hphi2*v2.star)
  y.star <- rpois(D*T, m*p.star)
  # Fitting boostrap model
  glmm.star <- glmer(formula=y.star ~ x + (1|dd/tt), family=poisson,
                     offset=log(m))
  hbeta.star <- fixef(glmm.star)
  var.star <- as.data.frame(VarCorr(glmm.star))
  hphi2.star <- var.star$sdcor[1]
  hphi1.star <- var.star$sdcor[2]
  # Plug-in bootstrap predictions of counts
  plugin.star[[b]] <- as.numeric(fitted(glmm.star))
  dif.plugin.star[[b]] <- (plugin.star[[b]] - m*p.star)^2
  # Loops to calculate EBPs
  ebp.star[[b]] <- vector()
  for (t in 1:T){
    for (d in 1:D){
      numdt.star <- dendt.star <- 0
      for (s1 in 1:S1) {
        vv1.star <- rnorm(n=1, mean=0, sd=1)
        vv1.star <- c(vv1.star, -vv1.star)
        numdtau.star <- dendtau.star <- 1
        for (tau in 1:T){
          condition <- tt==tau & dd==d
          vv2.star <- rnorm(n=S2, mean=0, sd=1)
          vv2.star <- c(vv2.star, -vv2.star)
          ydtau.star <- y.star[condition]
          deltau.star <- as.numeric(t==tau)
          xbetau.star <- hbeta.star[1] + hbeta.star[2]*x[condition] +
            hphi1.star*vv1.star + hphi2.star*vv2.star
          sum1num.star <- (ydtau.star + deltau.star)*xbetau.star
          sum2.star <- m[condition] * exp(xbetau.star)
          numdtau.star <- numdtau.star * sum(exp(sum1num.star-sum2.star))
          sum1den.star <- ydtau.star * xbetau.star
          dendtau.star <- dendtau.star * sum(exp(sum1den.star-sum2.star))
        }
        numdt.star <- numdt.star + numdtau.star
        dendt.star <- dendt.star + dendtau.star
      }
      condition2 <- tt==t&dd==d
      # EBP
      ebp.star[[b]][condition2] <- m[condition2]*numdt.star/dendt.star
    }
  }
  dif.ebp.star[[b]] <- (ebp.star[[b]] - m*p.star)^2
}
mse.plugin.star <- Reduce(f="+", dif.plugin.star)/B
mse.ebp.star <- Reduce(f="+", dif.ebp.star)/B
# root-MSE bootstrap estimates of plug-in and EBPs
result$rmseplugp.star <- sqrt(mse.plugin.star)
result$rmseebpp.star <- sqrt(mse.ebp.star)
result$cv.plugp.star <- 100*result$rmseplugp.star/mu
result$cv.eb.star <- 100*result$rmseebpp.star/mu



## Saving the results

output <- result[c(1,2,5,7,8,11,12)]
head(round(output,3), 10)





###
