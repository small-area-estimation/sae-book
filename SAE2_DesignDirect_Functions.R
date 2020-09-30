



## This R code is a set of functions that correspond to chapter 2 of the book 
## A COURSE ON SMALL AREA ESTIMATION AND MIXED MODEL available at
## www.springer.com

## AUTHORS: Morales D., Esteban M.D., Perez A., Hozba T. 



## dir1 Direct function. Horvitz-Thompson version

dir1 <- function(data, w, domain, Nd) {
  if(is.vector(data)){
    last <- length(domain) + 1
    Nd.hat <- aggregate(w, by=domain, sum)[,last]
    nd <- aggregate(rep(1, length(data)), by=domain, sum)[,last]
    tot <- aggregate(w*data, by=domain, sum)
    names(tot) <- c(names(domain), "tot")
    var.tot <- aggregate(w*(w-1)*data^2, by=domain, sum)[,last]
    if(missing(Nd)){
      return(cbind(tot, var.tot, Nd.hat, nd))
    }
    else{
      mean <- tot[,last]/Nd
      var.mean <- var.tot/Nd^2
      return(cbind(tot, var.tot, mean, var.mean, Nd.hat, Nd, nd))
    }
  }
  else{
    warning("Only a numeric or integer vector must be called as data",
            call. = FALSE)
  }
}



## dir2 Direct function. Hajek version

dir2 <- function(data, w, domain, Nd) {
  if(is.vector(data)){
    last <- length(domain) + 1
    Nd.hat <- aggregate(w, by=domain, sum)[,last]
    nd <- aggregate(rep(1, length(data)), by=domain, sum)[,last]
    Sum <- aggregate(w*data, by=domain, sum)
    mean <- Sum[,last]/Nd.hat
    dom <- as.numeric(Reduce("paste0", domain))
    if(length(domain)==1){
      domain.unique <- sort(unique(dom))
    }
    else{
      domain.unique <- as.numeric(Reduce("paste0", Sum[,1:length(domain)]))
    }
    difference <- list()
    for(d in 1:length(mean)){
      condition <- dom==domain.unique[d]
      difference[[d]] <- w[condition]*(w[condition]-1)*(data[condition]-mean[d])^2
    }
    var.mean <- unlist(lapply(difference, sum))/Nd.hat^2
    if(missing(Nd)){
      return(data.frame(dom=Sum[,-last], mean, var.mean, Nd.hat, nd))
    }
    else{
      tot <- mean*Nd
      var.tot <- var.mean*Nd^2
      return(data.frame(dom=Sum[,-last], tot, var.tot, mean, var.mean, Nd.hat, Nd, nd))
    }
  }
  else{
    warning("Only a numeric or integer vector must be called as data",
            call. = FALSE)
  }
}





###
