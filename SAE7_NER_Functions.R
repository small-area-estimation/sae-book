



## This R code is a set of functions that correspond to chapter 7 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## In this chapter the code has been separated in two parts.
## This part contains all the functions that it will be used at the main program.

########################################################################################


## The function seed calculates initial estimates of the parameters of a NER model.

seed <- function(Xd, yd, n, nd){
  XtX <- Reduce(Map(Xd, Xd, f=crossprod), f="+")
  Xty <- Reduce(Map(Xd, yd, f=crossprod), f="+")
  beta0 <- solve(XtX)%*%Xty
  sigmau2.0 <- mean((sapply(yd, mean) - sapply(lapply(Xd, colMeans), beta0,
                                               FUN=crossprod))^2)
  yxbeta <- Map(yd, lapply(Xd, beta0, FUN="%*%"), f="-")
  sigmae2.0 <- sum(mapply(yxbeta, yxbeta, FUN=crossprod))/(n-ncol(Xd[[1]]))
  return(c(beta0, sigmau2.0, sigmae2.0))
}



## The function log-lik calculates the log-likelihood function of the NER model,
## with negative sign, for applying minimization algorithms

log_lik <- function (parameters, Xd, yd, n, nd) {
  beta <- parameters[1:3]
  sigmau <- parameters[4]
  sigmae <- parameters[5]
  sum1 <- n/2*(log(2*pi))
  Vd <- lapply(mapply(mapply(sigmau/sigmae,nd,nd, FUN="matrix"),
                      lapply(nd, FUN=diag),FUN="+"), sigmae, FUN="*")
  Ad <- mapply(mapply(sigmau/sigmae,nd,nd, FUN="matrix"),
               lapply(nd, FUN=diag),FUN="+")
  Ad.det <- lapply(Ad, FUN=det)
  sum2 <- 0.5*sum(mapply(lapply(nd, log(sigmae), FUN="*"),
                         lapply(Ad.det, FUN=log), FUN="+"))
  Vd.inv <- lapply(Vd,solve)
  yxbeta <- mapply(yd, lapply(Xd, beta, FUN="%*%"), FUN="-")
  yxbetat <- mapply(mapply(yd, lapply(Xd, beta, FUN="%*%"), FUN="-"),
                    FUN=t)
  yxbetat.v.yxbeta <- mapply(mapply(yxbetat, Vd.inv, FUN="%*%"), yxbeta,
                             FUN="%*%")
  sum3 <- 0.5*sum(mapply(mapply(yxbetat, Vd.inv, FUN="%*%"), yxbeta,
                         FUN="%*%"))
  l <- sum1 + sum2 + sum3
  return(l)
}



## The function scores calculates the first partial derivatives
## of the log-likelihood function of the NER model.

scores <- function(parameters, Xd, yd, n, nd, small=TRUE){
  beta <- parameters[1:3]
  sigmau <- parameters[4]
  sigmae <- parameters[5]
  Vd <- lapply(mapply(mapply(sigmau/sigmae,nd,nd, FUN="matrix"),
                      lapply(nd, FUN=diag),FUN="+"), sigmae, FUN="*")
  Vd.inv <- lapply(Vd,solve)
  S.beta <- apply(mapply(lapply(Xd, FUN=t), mapply(Vd.inv, mapply(yd,
                                                                  lapply(Xd, beta, FUN="%*%"),FUN="-"), FUN="%*%"),
                         FUN="%*%"), 1, sum)
  yxbeta <- mapply(yd, lapply(Xd, beta, FUN="%*%"), FUN="-")
  yxbetat <- mapply(mapply(yd, lapply(Xd, beta, FUN="%*%"), FUN="-"), FUN=t)
  S.sigmau1 <- Reduce(lapply(lapply(mapply(Vd.inv, mapply(1, nd, nd,
                                                          FUN="matrix"), FUN="%*%"), FUN=diag), FUN=sum), f=sum)
  S.sigmau2 <- sum(mapply(mapply(yxbetat, mapply(mapply(Vd.inv,
                                                        mapply(1, nd, nd, FUN="matrix"), FUN="%*%"), Vd.inv,
                                                 FUN="%*%"), FUN="%*%"), yxbeta, FUN="%*%"))
  S.sigmau <- -0.5*S.sigmau1 + 0.5*S.sigmau2
  S.sigmae1 <- Reduce(lapply(lapply(Vd.inv, diag), sum), f=sum)
  S.sigmae2 <- sum(mapply(yxbetat, mapply(Vd.inv, mapply(Vd.inv, yxbeta,
                                                         FUN="%*%"), FUN="%*%"), FUN="%*%"))
  S.sigmae <- -0.5*S.sigmae1 + 0.5*S.sigmae2
  S <- c(S.beta, S.sigmau, S.sigmae)
  if(small==TRUE)
    S <- c(S.beta/D, S.sigmau/n, S.sigmae/n)
  return(S)
}



## The function Inf.Fisher calculates the Fisher information matrix.

Inf.Fisher <- function(parameters, Xd, yd, n, nd){
  beta <- parameters[1:3]
  sigmau <- parameters[4]
  sigmae <- parameters[5]
  Vd <- lapply(mapply(mapply(sigmau/sigmae, nd, nd, FUN="matrix"),
                      lapply(nd,FUN=diag),FUN="+"), sigmae, FUN="*")
  Vd.inv <- lapply(Vd, solve)
  yxbeta <- mapply(yd, lapply(Xd, beta, FUN="%*%"), FUN="-")
  yxbetat <- mapply(mapply(yd, lapply(Xd, beta, FUN="%*%"), FUN="-"), FUN=t)
  F.beta.beta <- Reduce(Map(lapply(Xd, FUN=t), mapply(Vd.inv,  Xd,
                                                      FUN="%*%"), f="%*%"), f="+")
  Vd.inv.Jnd <- mapply(Vd.inv, mapply(1,nd,nd, FUN="matrix"), FUN="%*%")
  F.sigmau.sigmau <- 0.5*sum(sapply(lapply(Map(Vd.inv.Jnd, Vd.inv.Jnd,
                                               f="%*%"), FUN=diag),FUN=sum))
  F.sigmau.sigmae <- 0.5*sum(sapply(lapply(Map(Vd.inv, Vd.inv.Jnd, f="%*%"),
                                           FUN=diag), FUN=sum))
  F.sigmae.sigmae <- 0.5*sum(sapply(lapply(Map(Vd.inv, Vd.inv, f="%*%"),
                                           FUN=diag), FUN=sum))
  F.sigma <- matrix(c(F.sigmau.sigmau, F.sigmau.sigmae, F.sigmau.sigmae,
                      F.sigmae.sigmae), ncol=2)
  return(list(F.beta.beta, F.sigma))
}



## The R function FisherScoring runs the Fisher-Scoring algorithm described in the book.
## This algorithm calculates the ML estimates of the parameters of a NER model

FisherScoring <- function(Xd, yd, n, nd, n.iter=100){
  parameters <- list()
  parameters[[1]] <- seed(Xd, yd, n, nd)
  p <- ncol(Xd[[1]])
  iter <- 1
  while(iter < n.iter){
    S <- scores(parameters[[iter]], Xd, yd, n, nd, small=FALSE)
    F <- Inf.Fisher(parameters[[iter]], Xd, yd, n, nd)
    diff.beta <- solve(F[[1]])%*%S[1:p]
    diff.sigma <- solve(F[[2]])%*%S[(p+1):(p+2)]
    diff <- rbind(diff.beta, diff.sigma)
    parameters[[iter+1]] <- parameters[[iter]] + diff
    RT <- abs(diff)/abs(parameters[[iter+1]])
    out <- RT>0.01
    if(sum(out)==0)
      break
    iter <- iter + 1
  }
  return(parameters)
}





###
