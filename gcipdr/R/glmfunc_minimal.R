## glmfunc_minimal.R contains modules needed to perform NORTA transformation.
## Copyright (C) 2018 Federico Bonofiglio

## This file is part of gcipdr.

    ## gcipdr is free software: you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation, either version 3 of the License, or
    ## (at your option) any later version.

    ## gcipdr is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.

    ## You should have received a copy of the GNU General Public License
    ## along with gcipdr.  If not, see <https://www.gnu.org/licenses/>.


                 
 ## marginal IPD modelling            

###ancillary space stabilization (in case of convergence error)

stableX <- function(n=n, H=H,
                    expr=function(x)rgmom(n, mx, sdx), type="mean"){

H <- H

Xstar <- matrix(numeric(n*H), ncol=H)


Xstar <- apply(Xstar, 2, expr ) # *repropagation phase*

    Xmode <- apply(

 apply(Xstar, 2, function(x)sort(x)),
  1, function(x)switch(type,
                       mean=mean(x),
                       mode=Mode(x))
   )

    Xmode
   }

##

vecmode <- function(vec){

    type <- typeof(vec)

    type

   }



dvs <- function(frame){

 type <- apply(frame, 2, vecmode)

res <- rbind(colnames(frame), type)
 res    
   }





rgmom <- function(n, mean, sd){

    var <- (sd^2)
    alpha <- (mean^2)/var
    beta <- mean/var

  res <- rgamma(n , alpha, beta)
     res
   }

#

qgmom <- function(x, mean, sd){

    var <- (sd^2)
    alpha <- (mean^2)/var
    beta <- mean/var

  qgamma(x , alpha, beta)
     
   }

#

dgmom <- function(x, mean, sd){

        var <- (sd^2)
    alpha <- (mean^2)/var
    beta <- mean/var

    dgamma(x , alpha, beta)

  }



#

rnbmom <- function(n, mean, sd){

 var <- sd^2
p <- mean/var
    r <- mean*p/(1-p)

    res <- rnbinom(n, size=r, prob=p)
    res

   }



#


qbmom <- function(x, mean, n){

   alpha <- mean*n
    beta <- n-alpha

    qbeta(x, alpha, beta)

  }


dbmom <- function(x, mean, n){

   alpha <- mean*n
    beta <- n-alpha

    dbeta(x, alpha, beta)

  }


##  Johnson system densities (from : http://reference.wolfram.com/language/ref/JohnsonDistribution.html). Replaced: used original definitions from Johnson paper 49


dJohnsonSB <- function(x, gamma, delta, lambda, csi){
  

y <- (x - csi)/lambda  # conversion  0 < y < 1

# errterm <-  (x - csi)
  #  errfac <- (-x + csi + lambda)   # y/(1-y) # conversion: should be always positive (if z > 1, something wrong: FORTRAN errors, McLeod comunication) !!!!!
       # safeguard against FORTRAN errors ?
if ( any(y < 0) | any(y > 1)  ) {
     k <- length(y[y < 0])
     j <- length(y[y > 1])
 lb <- min(y[y > 0], na.rm = T)
     ub <- max(y[y < 1], na.rm = T)
        y[y < 0] <- seq(10e-10, lb, length = k +1)[1:k]
     y[y > 1] <- seq(ub, 0.99999, length = j +1)[2:(j+1)]
        #   y <- sort(y)
    }
    
    comply <- 1 - y  # complement y        
kernel <- gamma + delta * log(y/comply)    # z variate    #  log(errterm/errfac)

 dens <- delta*exp(-0.5*kernel^2)    # drop lambda if u use z

    normalizer <- sqrt(2*pi)*y*comply*lambda    # appearantly Johnson forgot a lambda in the normalizer term         # *errterm*errfac                           
    
    dens/normalizer
    
    }

#
dJohnsonSL <- function(x, gamma, delta, lambda, csi){

    y <- (x - csi)/lambda   # y lognormal,  y > 0
               #errterm <- (x - csi)    # should be always positive !!!!
if (any(y < 0) ) {
     k <- length(y[y < 0])
     lb <- min(y[y > 0], na.rm = T)
        y[y < 0] <- seq(10e-10, lb, length = k +1)[1:k]
   y <- sort(y)
   }
    
    kernel <- gamma + delta* log(y)   # z variate          # log(errterm/lambda)

    dens <- delta*exp(-0.5*kernel^2)

    normalizer <- sqrt(2*pi)*y           #  (x-csi)
    
  dens/normalizer

 }

#
dJohnsonSU <- function(x, gamma, delta, lambda, csi){


## errterm <- (x - csi)/lambda
    
##     kernel <- gamma + delta*asinh(errterm)

##     dens <- delta*exp(-0.5*kernel^2)

##     a <- (x-csi)^2
##     b <- lambda^2
    
##   normalizer <- sqrt(2*pi)*sqrt(a+b)

    trans <- function(y) log(x + sqrt(1+y^2) )
    
  y <- (x - csi)/lambda  
    kernel <- gamma + delta*trans(y)   # z variate
    dens <- delta*exp(-0.5*kernel^2)
    normalizer <- sqrt(2*pi)*sqrt(1+y^2)
    
    dens/normalizer
    
}
#

dJohnsonZ <- function(x, itype, gamma, delta, lambda, csi){

 z <- zJohnsonDistribution(x, itype, gamma, delta, lambda, csi)

    if ( ( itype == 1 | itype ==3 ) & any(z < 0)) {

 k <- length(z[z < 0])
     lb <- min(z[z > 0], na.rm = T)
        z[z < 0] <- seq(10e-10, lb, length = k +1)[1:k]
        z <- sort(z)

    }
            
    percentiles <- pnorm(z)

    density <- c(0, diff(percentiles))

    density

   }
#

dJohnson <- function(x, itype, gamma, delta, lambda, csi, mapping = c("direct", "indirect")){

    mapping <- match.arg(mapping)

    # browser()

    out <- switch( mapping,

     indirect = dJohnsonZ(x, itype, gamma, delta, lambda, csi),

     direct = {

   if (itype == 1)
  dJohnsonSL(x, gamma, delta, lambda, csi)
       else {
   if (itype == 2)
   dJohnsonSU(x, gamma, delta, lambda, csi)
    else {
    if (itype == 3)
   dJohnsonSB(x, gamma, delta, lambda, csi)
       else{
           if (itype == 4)
         dnorm(x, lambda, csi)
           else
               stop("argument 'itype' takes values between 1 and 4")
       }
    }
       }

     }

           )
                                                              
    out

   }


#

reldist <- function(mean, true, range){  # relative distance


    (mean-true)/range

   }



