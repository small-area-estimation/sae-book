



## This R code is a set of functions that correspond to chapter 19 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## In this chapter the code has been separated in two parts.
## This part contains all the functions that it will be used at the main program.

########################################################################################


## The R function direct.bfh calculates the direct estimators of two domain totals
## and means and the corresponding direct estimators of their variances and covariances.

direct.bfh <- function(data1, data2, w, domain) {
  if(!inherits(data1,"data.frame")){
    last <- length(domain) + 1
    # Population and sample size
    Nd <- aggregate(w, by=domain, sum)[,last]
    nd <- aggregate(rep(1, length(data1)), by=domain, sum)[,last]
    # y1 direct
    y1 <- aggregate(w*data1, by=domain, sum)
    names(y1) <- c(names(domain), "y1")
    y1.mean <- y1[,last]/Nd
    # y2 direct
    y2 <- aggregate(w*data2, by=domain, sum)
    names(y2) <- c(names(domain), "y2")
    y2.mean <- y2[,last]/Nd
    # Variances and covariances
    difference1 <- difference2 <- difference12 <- list()
    for(d in 1:length(y1.mean)){
      condition <- domain[[1]]==y1[d,1]
      difference1[[d]] <- w[condition]*(w[condition]-1) *
        (data1[condition]-y1.mean[d])^2
      difference2[[d]] <- w[condition]*(w[condition]-1) *
        (data2[condition]-y2.mean[d])^2
      difference12[[d]] <- w[condition]*(w[condition]-1) *
        (data1[condition]-y1.mean[d]) *
        (data2[condition]-y2.mean[d])
    }
    var.y1 <- unlist(lapply(difference1, sum))
    var.y2 <- unlist(lapply(difference2, sum))
    cov.y12 <- unlist(lapply(difference12, sum))
    var.y1.mean <- var.y1/Nd^2
    var.y2.mean <- var.y2/Nd^2
    cov.y12.mean <- cov.y12/Nd^2
    return(cbind(y1, var.y1, y1.mean, var.y1.mean,
                 y2, var.y2, y2.mean, var.y2.mean,
                 cov.y12, cov.y12.mean, nd, Nd))
  }
  else{
    warning("Only a numeric or integer vector must be called as data",
            call. = FALSE)
  }
}



## The R function UveU calculates the variance matrix V_ud of the random effect u_d.

UveU <- function(thetas) {
  Vuu <- matrix(0, nrow=2, ncol=2)
  sup <- c(thetas[3]*sqrt(thetas[1])*sqrt(thetas[2]))
  Vuu[upper.tri(Vuu)] <- sup
  Vuu <- Vuu + t(Vuu)
  diag(Vuu) <- c(thetas[1], thetas[2])
  return(Vuu)
}



## The R function Thetas calculates the vector of parameters 

Thetas <- function(Vus) {
  Zethas <- vector()
  Zethas <- c(Vus[1,1], Vus[2,2], Vus[1,2]/(sqrt(Vus[1,1]*Vus[2,2])))
  return(Zethas)
}



## The R function FirstDer calculates the matrices of first derivatives
## of V_ud with respect to the three parameters respectively.

FirstDer <- function(thetas){
  Vud1 <- matrix(0, nrow=2, ncol=2)
  Vud1[1,1] <- 1
  Vud1[1,2] <- Vud1[2,1] <- thetas[3]*sqrt(thetas[2])/(2*sqrt(thetas[1]))
  Vud2 <- matrix(0, nrow=2, ncol=2)
  Vud2[2,2] <- 1
  Vud2[1,2] <- Vud2[2,1] <- thetas[3]*sqrt(thetas[1])/(2*sqrt(thetas[2]))
  Vud3 <- matrix(0, nrow=2, ncol=2)
  Vud3[1,2] <- Vud3[2,1] <- sqrt(thetas[1])*sqrt(thetas[2])
  return(list(Vud1, Vud2, Vud3))
}



## The R function REML.BFH calculates the REML estimators of the parameters of the bivariate Fay-Herriot model.

REML.BFH <- function(X, y, D, Ved, Vud, MAXITER = 100) {
  Vud.f  <- Vud
  Xd <- X
  yd <- y
  theta.f <- Thetas(Vud)
  Vuds <- FirstDer(theta.f)
  Bad <- 0
  FLAG <- 0
  p <- ncol(X[[1]])
  # Iteration loop of Fisher-scoring algorithm
  for(ITER in 1:MAXITER){
    Vd.inv <- Vinvyd <- VinvXd <- list()
    VinvVud1 <- XtVinvVud1VinvX <- VinvVud1VinvVud1 <-
      XtVinvVud1VinvVud1VinvX <- list()
    tr.VinvVud1 <- ytVinvX <- ytVinvVud1Vinvy <- ytVinvVud1VinvX <-
      SumXtVinvVud1VinvX <- tr.VinvVud1VinvVud1 <- 0
    VinvVud2 <- XtVinvVud2VinvX <- VinvVud2VinvVud2 <-
      XtVinvVud2VinvVud2VinvX <- list()
    tr.VinvVud2 <- ytVinvVud2Vinvy <- ytVinvVud2VinvX <-
      SumXtVinvVud2VinvX <- tr.VinvVud2VinvVud2 <- 0
    VinvVud3 <- XtVinvVud3VinvX <- VinvVud3VinvVud3 <-
      XtVinvVud3VinvVud3VinvX <- list()
    tr.VinvVud3 <- ytVinvVud3Vinvy <- ytVinvVud3VinvX <-
      SumXtVinvVud3VinvX <- tr.VinvVud3VinvVud3 <- 0
    VinvVud1VinvVud2 <- XtVinvVud1VinvVud2VinvX <- list()
    tr.VinvVud1VinvVud2 <- 0
    VinvVud1VinvVud3 <- XtVinvVud1VinvVud3VinvX <- list()
    tr.VinvVud1VinvVud3 <- 0
    VinvVud2VinvVud3 <- XtVinvVud2VinvVud3VinvX <- list()
    tr.VinvVud2VinvVud3 <- 0
    Q.inv <- matrix(0, nrow=p, ncol=p)
    # Domain loop 1
    for(d in 1:D) {
      Vd <- Vud.f+Ved[[d]]
      if(abs(det(Vd))<1e-09 || abs(det(Vd))>1e+11){
        FLAG <- 1
        Bad <- Bad+1
        break
      }
      Vd.inv[[d]] <- solve(Vd)
      Vinvyd[[d]] <- Vd.inv[[d]]%*%yd[[d]]
      VinvXd[[d]] <- Vd.inv[[d]]%*%Xd[[d]]
      Q.inv <- Q.inv + t(Xd[[d]])%*%VinvXd[[d]]
      # Score S1
      VinvVud1[[d]] <- Vd.inv[[d]]%*%Vuds[[1]]
      tr.VinvVud1 <- tr.VinvVud1 + sum(diag(VinvVud1[[d]]))
      XtVinvVud1VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[1]]%*%VinvXd[[d]]
      ytVinvX <- ytVinvX + t(yd[[d]])%*%VinvXd[[d]]
      ytVinvVud1Vinvy <- ytVinvVud1Vinvy + 
        t(Vinvyd[[d]])%*%Vuds[[1]]%*%Vinvyd[[d]]
      ytVinvVud1VinvX <- ytVinvVud1VinvX + 
        t(Vinvyd[[d]])%*%Vuds[[1]]%*%VinvXd[[d]]
      SumXtVinvVud1VinvX <- SumXtVinvVud1VinvX + XtVinvVud1VinvX[[d]]
      # Score S2
      VinvVud2[[d]] <- Vd.inv[[d]]%*%Vuds[[2]]
      tr.VinvVud2 <-  tr.VinvVud2 + sum(diag(VinvVud2[[d]]))
      XtVinvVud2VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[2]]%*%VinvXd[[d]]
      ytVinvVud2Vinvy <- ytVinvVud2Vinvy + 
        t(Vinvyd[[d]])%*%Vuds[[2]]%*%Vinvyd[[d]]
      ytVinvVud2VinvX <- ytVinvVud2VinvX + 
        t(Vinvyd[[d]])%*%Vuds[[2]]%*%VinvXd[[d]]
      SumXtVinvVud2VinvX <- SumXtVinvVud2VinvX + XtVinvVud2VinvX[[d]]
      # Score S3
      VinvVud3[[d]] <- Vd.inv[[d]]%*%Vuds[[3]]
      tr.VinvVud3 <- tr.VinvVud3 + sum(diag(VinvVud3[[d]]))
      XtVinvVud3VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[3]]%*%VinvXd[[d]]
      ytVinvVud3Vinvy <- ytVinvVud3Vinvy + 
        t(Vinvyd[[d]])%*%Vuds[[3]]%*%Vinvyd[[d]]
      ytVinvVud3VinvX <- ytVinvVud3VinvX + 
        t(Vinvyd[[d]])%*%Vuds[[3]]%*%VinvXd[[d]]
      SumXtVinvVud3VinvX <- SumXtVinvVud3VinvX + XtVinvVud3VinvX[[d]]
      # Fisher information element F11
      VinvVud1VinvVud1[[d]] <- VinvVud1[[d]]%*%VinvVud1[[d]]
      tr.VinvVud1VinvVud1 <- tr.VinvVud1VinvVud1 + 
        sum(diag(VinvVud1VinvVud1[[d]]))
      XtVinvVud1VinvVud1VinvX[[d]] <- 
        t(Xd[[d]])%*%VinvVud1VinvVud1[[d]]%*%VinvXd[[d]]
      # Fisher information element F22
      VinvVud2VinvVud2[[d]] <- VinvVud2[[d]]%*%VinvVud2[[d]]
      tr.VinvVud2VinvVud2 <- tr.VinvVud2VinvVud2 + 
        sum(diag(VinvVud2VinvVud2[[d]]))
      XtVinvVud2VinvVud2VinvX[[d]] <- 
        t(Xd[[d]])%*%VinvVud2VinvVud2[[d]]%*%VinvXd[[d]]
      # Fisher information element F33
      VinvVud3VinvVud3[[d]] <- VinvVud3[[d]]%*%VinvVud3[[d]]
      tr.VinvVud3VinvVud3 <- tr.VinvVud3VinvVud3 + 
        sum(diag(VinvVud3VinvVud3[[d]]))
      XtVinvVud3VinvVud3VinvX[[d]] <- 
        t(Xd[[d]])%*%VinvVud3VinvVud3[[d]]%*%VinvXd[[d]]
      # Fisher information element F12
      VinvVud1VinvVud2[[d]] <- VinvVud1[[d]]%*%VinvVud2[[d]]
      tr.VinvVud1VinvVud2 <- tr.VinvVud1VinvVud2 + 
        sum(diag(VinvVud1VinvVud2[[d]]))
      XtVinvVud1VinvVud2VinvX[[d]] <- 
        t(Xd[[d]])%*%VinvVud1VinvVud2[[d]]%*%VinvXd[[d]]
      # Fisher information element F13
      VinvVud1VinvVud3[[d]] <- VinvVud1[[d]]%*%VinvVud3[[d]]
      tr.VinvVud1VinvVud3 <- tr.VinvVud1VinvVud3 + 
        sum(diag(VinvVud1VinvVud3[[d]]))
      XtVinvVud1VinvVud3VinvX[[d]] <- 
        t(Xd[[d]])%*%VinvVud1VinvVud3[[d]]%*%VinvXd[[d]]
      # Fisher information element F23
      VinvVud2VinvVud3[[d]] <- VinvVud2[[d]]%*%VinvVud3[[d]]
      tr.VinvVud2VinvVud3 <- tr.VinvVud2VinvVud3 + 
        sum(diag(VinvVud2VinvVud3[[d]]))
      XtVinvVud2VinvVud3VinvX[[d]] <- 
        t(Xd[[d]])%*%VinvVud2VinvVud3[[d]]%*%VinvXd[[d]]
      #
      if(FLAG==1){
        FLAG <- 0
        ITER <- MAXITER
        break
      }
    } # End of the domain loop 1
    Q <- solve(Q.inv)
    tr.XtVinvVud1VinvXQ <- tr.XtVinvVud2VinvXQ <- tr.XtVinvVud3VinvXQ <-
      tr.XtVinvVud4VinvXQ <- tr.XtVinvVud5VinvXQ <- tr.XtVinvVud6VinvXQ <- 0
    tr.XtVinvVud1VinvVud1VinvXQ <- tr.XtVinvVud2VinvVud2VinvXQ <-
      tr.XtVinvVud3VinvVud3VinvXQ <- 0
    tr.XtVinvVud1VinvVud2VinvXQ <- tr.XtVinvVud1VinvVud3VinvXQ <-
      tr.XtVinvVud2VinvVud3VinvXQ <- 0
    # Domain loop 2
    for(d in 1:D){
      tr.XtVinvVud1VinvXQ <- tr.XtVinvVud1VinvXQ +
        sum(diag(XtVinvVud1VinvX[[d]]%*%Q))
      tr.XtVinvVud2VinvXQ <- tr.XtVinvVud2VinvXQ +
        sum(diag(XtVinvVud2VinvX[[d]]%*%Q))
      tr.XtVinvVud3VinvXQ <- tr.XtVinvVud3VinvXQ +
        sum(diag(XtVinvVud3VinvX[[d]]%*%Q))
      tr.XtVinvVud1VinvVud1VinvXQ <- tr.XtVinvVud1VinvVud1VinvXQ +
        sum(diag(XtVinvVud1VinvVud1VinvX[[d]]%*%Q))
      tr.XtVinvVud2VinvVud2VinvXQ <- tr.XtVinvVud2VinvVud2VinvXQ +
        sum(diag(XtVinvVud2VinvVud2VinvX[[d]]%*%Q))
      tr.XtVinvVud3VinvVud3VinvXQ <- tr.XtVinvVud3VinvVud3VinvXQ +
        sum(diag(XtVinvVud3VinvVud3VinvX[[d]]%*%Q))
      tr.XtVinvVud1VinvVud2VinvXQ <- tr.XtVinvVud1VinvVud2VinvXQ +
        sum(diag(XtVinvVud1VinvVud2VinvX[[d]]%*%Q))
      tr.XtVinvVud1VinvVud3VinvXQ <- tr.XtVinvVud1VinvVud3VinvXQ +
        sum(diag(XtVinvVud1VinvVud3VinvX[[d]]%*%Q))
      tr.XtVinvVud2VinvVud3VinvXQ <- tr.XtVinvVud2VinvVud3VinvXQ +
        sum(diag(XtVinvVud2VinvVud3VinvX[[d]]%*%Q))
    } # End of the domain loop 2
    tr.PV1 <- tr.VinvVud1 - tr.XtVinvVud1VinvXQ
    tr.PV2 <- tr.VinvVud2 - tr.XtVinvVud2VinvXQ
    tr.PV3 <- tr.VinvVud3 - tr.XtVinvVud3VinvXQ
    SumXtVinvVud1VinvXQ <- SumXtVinvVud1VinvX%*%Q
    SumXtVinvVud2VinvXQ <- SumXtVinvVud2VinvX%*%Q
    SumXtVinvVud3VinvXQ <- SumXtVinvVud3VinvX%*%Q
    tr.PV1PV1  <- tr.VinvVud1VinvVud1 - 2*tr.XtVinvVud1VinvVud1VinvXQ +
      sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud1VinvXQ))
    tr.PV1PV2  <- tr.VinvVud1VinvVud2 - 2*tr.XtVinvVud1VinvVud2VinvXQ +
      sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud2VinvXQ))
    tr.PV1PV3  <- tr.VinvVud1VinvVud3 - 2*tr.XtVinvVud1VinvVud3VinvXQ +
      sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud3VinvXQ))
    tr.PV2PV2  <- tr.VinvVud2VinvVud2 - 2*tr.XtVinvVud2VinvVud2VinvXQ +
      sum(diag(SumXtVinvVud2VinvXQ%*%SumXtVinvVud2VinvXQ))
    tr.PV2PV3  <- tr.VinvVud2VinvVud3 - 2*tr.XtVinvVud2VinvVud3VinvXQ +
      sum(diag(SumXtVinvVud2VinvXQ%*%SumXtVinvVud3VinvXQ))
    tr.PV3PV3  <- tr.VinvVud3VinvVud3 - 2*tr.XtVinvVud3VinvVud3VinvXQ +
      sum(diag(SumXtVinvVud3VinvXQ%*%SumXtVinvVud3VinvXQ))
    ytVinvXQ <- ytVinvX%*%Q
    ytPV1Py  <- ytVinvVud1Vinvy - 2*ytVinvVud1VinvX%*%t(ytVinvXQ) +
      ytVinvXQ%*%SumXtVinvVud1VinvX%*%t(ytVinvXQ)
    ytPV2Py  <- ytVinvVud2Vinvy - 2*ytVinvVud2VinvX%*%t(ytVinvXQ) +
      ytVinvXQ%*%SumXtVinvVud2VinvX%*%t(ytVinvXQ)
    ytPV3Py  <- ytVinvVud3Vinvy - 2*ytVinvVud3VinvX%*%t(ytVinvXQ) +
      ytVinvXQ%*%SumXtVinvVud3VinvX%*%t(ytVinvXQ)
    # Scores
    S1 <- -0.5*tr.PV1 + 0.5*ytPV1Py
    S2 <- -0.5*tr.PV2 + 0.5*ytPV2Py
    S3 <- -0.5*tr.PV3 + 0.5*ytPV3Py
    # Fisher information matrix
    F11  <- 0.5*tr.PV1PV1
    F12  <- 0.5*tr.PV1PV2
    F13  <- 0.5*tr.PV1PV3
    F22  <- 0.5*tr.PV2PV2
    F23  <- 0.5*tr.PV2PV3
    F33  <- 0.5*tr.PV3PV3
    #
    if(ITER>1){
      F.sig.prev <- Fsig
      Q.prev <- Q
      theta.f.prev <- theta.f
    }
    Ssig <- rbind(S1, S2, S3)
    Fsig <- matrix(c(F11, F12, F13, F12, F22, F23, F13, F23, F33),
                   nrow=3, byrow=TRUE)
    #
    Fsig.inv <- solve(Fsig)
    dif <- Fsig.inv%*%Ssig
    theta.rml <- theta.f
    theta.f <- theta.rml + dif
    Vud.f <- UveU(theta.f)
    Vuds <- FirstDer(theta.f)
    # Stop because of warning
    if( theta.f[1]<=0 || theta.f[2]<=0 || abs(theta.f[3])>1){
      if(ITER==1){
        return(list(theta.rml, ITER, ITER))
        cat("WARNING: Stop at the first iteration. Out of parametric space")
      }
      return(list(theta.f.prev, F.sig.prev, ITER, Bad, Q.prev))
      cat("WARNING: Out of parametric space")
    }
    # Stop because of warning
    if(det(Vud.f)<=0){
      return(list(theta.f.prev, F.sig.prev, ITER, Bad, Q.prev))
      cat("WARNING: Singularity of Vud matrix")
    }
    # Stopping rule for iterative algorithm
    if(identical(as.numeric(abs(dif)<1e-06), rep(1,3))){
      break
    }
  }
  # End of Iteration loop (Fisher-scoring algorithm)
  return(list(theta.f, Fsig, ITER, Bad, Q))
}



## By using the REML method the R function BETA.U.BFH calculates
## the mode predictors of the random effects and the estimators
## of the regression parameters of a bivariate Fay-Herriot model.

BETA.U.BFH <- function(X, y, D, Ved, theta.hat) {
  Vd.inv <- Vd.hat <- Xd <- yd <- u <- list()
  p <- ncol(X[[1]]); Xd <- X; yd <- y
  Vudbeta <- UveU(theta.hat)
  Q.inv <- matrix(0, nrow=p, ncol=p)
  XVy <- rep(0,p)
  for(d in 1:D){
    Vd.hat[[d]] <- Vudbeta + Ved[[d]]
    Vd.inv[[d]] <- solve(Vd.hat[[d]])
    Q.inv <- Q.inv + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]
    XVy <- XVy + t(Xd[[d]])%*%Vd.inv[[d]]%*%yd[[d]]
  }
  Q <- solve(Q.inv); betax <- Q%*%XVy
  for(d in 1:D){
    u[[d]] <- Vudbeta%*%Vd.inv[[d]]%*%(yd[[d]]-X[[d]]%*%betax)
  }
  return(list(betax, u))
}



## The R function pvalue calculates the standard error, the z-value
## and the p-value for testing if the regression parameters are equal to zero.

pvalue <- function(beta.fitted, fit) {
  if(!inherits(fit,"try-error")){
    stderror <- sqrt(diag(fit[[5]]))
    z <- abs(as.vector(beta.fitted))/stderror
    p <- pnorm(z, lower.tail=F)
    return(data.frame(stderror, z, "p-value"=2*p))
  }
  else{
    warning("Only a converging algorithm is allowed", call. = FALSE)
  }
}



## The R function CI calculates asymptotic confidence intervals
## for the regression and the variance component parameters.

CI <- function(beta.fitted, fit, conf.level = 0.95) {
  if(!inherits(fit,"try-error")){
    alpha <- 1-conf.level; k <- 1-alpha/2; z <- qnorm(k)
    Finv <- solve(fit[[2]]); Q <- fit[[5]]
    lower.beta <- beta.fitted - z*sqrt(diag(Q))
    lower.theta <- fit[[1]] - z*sqrt(diag(Finv))
    upper.beta <- beta.fitted + z*sqrt(diag(Q))
    upper.theta <- fit[[1]] + z*sqrt(diag(Finv))
    return(list(data.frame(Inf.beta=lower.beta, Sup.beta=upper.beta),
                data.frame(Inf.theta=lower.theta, Sup.theta=upper.theta)))
  }
  else{
    warning("Only a converging algorithm is allowed", call. = FALSE)
  }
}



## The R function MSE.BFH calculates the analytic estimator of the MSE of the BLUP.

MSE.BFH <- function(X, D, Ved, fit){
  p <- ncol(X[[1]])
  Vud <- UveU(fit[[1]])
  Vu <- bdiag(lapply(1:D, function(d) Vud ))
  Ve <- bdiag(Ved); V <- Vu+Ve; V.inv <- solve(V);  Z <- diag(2*D)
  # G1
  Te <- Vu - tcrossprod(Vu,Z)%*%V.inv%*%crossprod(Z,Vu)
  G1 <- Z%*%tcrossprod(Te,Z)
  # G2
  equis.matrix <- sapply(1:p,
                         function(j) sapply(1:D, function(i) X[[i]][,j]))
  Q <- fit[[5]]
  M <- equis.matrix - G1%*%crossprod(solve(Ve), equis.matrix)
  G2 <- M%*%tcrossprod(Q,M)
  # MSEs
  mse.k1 <- diag(G1+G2)[c(T,F)]
  mse.k2 <- diag(G1+G2)[c(F,T)]
  return(list(mse.k1, mse.k2))
}





###
