



## This R code is a set of routines that correspond to chapter 9 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## First, we run the code of chapter 8

source('SAE8_NEReblup.R')


## We do some preliminary matrix calculations

gammad <- nd*sigmau2/(nd*sigmau2+sigmae2)
udom <- sort(unique(dom))
Q <- QQ <- QQQ <- matrix(0,3,3); Vd.inv <- unod <- Xd <- list()
for (d in 1:D) {
  unod[[d]] <- matrix(1,nd[d])
  condition <- dom==udom[[d]]
  # Matrix Xd
  Xd[[d]] <- cbind(1, work[condition], nowork[condition])
  # Identity matrix
  Id <- diag(1,nd[d])
  # Inverse of Vd matrix
  Vd.inv[[d]] <- (Id-(gammad[d]/nd[d])*unod[[d]]%*%t(unod[[d]]))/sigmae2
  Q <- Q + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]
  QQ <- QQ + t(Xd[[d]])%*%Vd.inv[[d]]%*%unod[[d]]%*%t(unod[[d]])%*%
    Vd.inv[[d]]%*%Xd[[d]]
  QQQ <- QQQ + t(Xd[[d]])%*%Vd.inv[[d]]%*%Vd.inv[[d]]%*%Xd[[d]]
}
Q <- solve(Q)



## For calculating the Fisher information matrix Fuu, we run the code

Fuu1 <- Fuu2 <- Fuu3 <- 0
for (d in 1:D) {
  oneVone <- t(unod[[d]])%*%Vd.inv[[d]]%*%unod[[d]]
  Fuu1 <- Fuu1 + oneVone^2
  Fuu2 <- Fuu2 +
    oneVone*t(unod[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]%*%Q%*%t(Xd[[d]])%*%
    Vd.inv[[d]]%*%unod[[d]]
  Fuu3 <- Fuu3 +
    t(unod[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]%*%Q%*%QQ%*%Q%*%t(Xd[[d]])%*%
    Vd.inv[[d]]%*%unod[[d]]
}
Fuu <- 0.5*Fuu1 - Fuu2 + 0.5*Fuu3; Fuu



## For calculating the Fisher information matrix Fue, we run the code

Fue1 <- Fue2 <- Fue3 <- 0
for (d in 1:D) {
  Fue1 <- Fue1 + t(unod[[d]])%*%Vd.inv[[d]]%*%Vd.inv[[d]]%*%unod[[d]]
  Fue2 <- Fue2 + t(unod[[d]])%*%Vd.inv[[d]]%*%Vd.inv[[d]]%*%Xd[[d]]%*%Q%*%
    t(Xd[[d]])%*%Vd.inv[[d]]%*% unod[[d]]
  Fue3 <- Fue3 + t(unod[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]%*%Q%*%QQQ%*%Q%*%
    t(Xd[[d]])%*%Vd.inv[[d]]%*%unod[[d]]
}
Fue <- 0.5*Fue1 - Fue2 + 0.5*Fue3; Fue



## For calculating the Fisher information matrix Fee, we run the code

Fee1 <- Fee2 <- Fee3 <- 0
for (d in 1:D) {
  Fee1 <- Fee1 + sum(diag(Vd.inv[[d]]%*%Vd.inv[[d]]))
  Fee2 <- Fee2 + sum(diag(Vd.inv[[d]]%*%Vd.inv[[d]]%*%Xd[[d]]%*%Q%*%
                            t(Xd[[d]])%*%Vd.inv[[d]]))
  Fee3 <- Fee3 + sum(diag(Vd.inv[[d]]%*%Xd[[d]]%*%Q%*%QQQ%*%Q%*%
                            t(Xd[[d]])%*%Vd.inv[[d]]))
}
Fee <- 0.5*Fee1 - Fee2 + 0.5*Fee3; Fee



## We calculate the REML Fisher information matrix.

F <- matrix(c(Fuu, Fue, Fue, Fee), 2, 2); F
F.inv <- solve(F); F.inv



## We calculate g1

g1 <- (1-gammad)*sigmau2; g1



## We calculate g2

g2 <- vector()
for (d in 1:D) {
  xd.bar <- c(1, Xmean$work[d], Xmean$nowork[d])
  xd.dird <- c(1, dirwork[d], dirnowork[d])
  ad <- xd.bar - gammad[d]*xd.dird
  g2[d] <- t(ad)%*%Q%*%ad
}
g2



## We calculate g3

g31 <- g32 <- g33 <- g3 <- vector()
for (d in 1:D) {
  gfix <- nd[d]/((nd[d]*sigmau2+sigmae2)^3)
  g31[d] <- sigmau2^2*F.inv[2,2]*gfix
  g32[d] <- -2*sigmau2*sigmae2*F.inv[1,2]*gfix
  g33[d] <- sigmae2^2*F.inv[1,1]*gfix
}
g3 <- g31 + g32 + g33; g3



## We calcuate g4 and the mse

g4 <- sigmae2*(aux$TOT-nd)/(aux$TOT^2); g4
mseg <- g1 + g2 + 2*g3 + g4; mseg



## We calculate the relative difference in percentage

R.Bg <- round(100*(ebmse-mseg)/mseg,2); R.Bg



## Saving the results

output <- data.frame(mseB=round(ebmse,0), mseg=round(mseg,0), R.Bg)
output





###
