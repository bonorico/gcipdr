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


##  Johnson system densities: used original definitions from Johnson paper 49
## and from https://reference.wolfram.com/language/ref/JohnsonDistribution.html

dJohnsonSB <- function(x, gamma, delta, lambda, csi){
  

y <- (x - csi)/lambda  # conversion  0 < y < 1

       # safeguard against FORTRAN errors ?
if ( any(y < 0) | any(y > 1)  ) {
     k <- length(y[y < 0])
     j <- length(y[y > 1])
 lb <- min(y[y > 0], na.rm = T)
     ub <- max(y[y < 1], na.rm = T)
        y[y < 0] <- seq(10e-10, lb, length = k +1)[1:k]
     y[y > 1] <- seq(ub, 0.99999, length = j +1)[2:(j+1)]
        
    }
    
    comply <- 1 - y  # complement y        
   kernel <- gamma + delta * log(y/comply)    # z variate   

   dens <- delta*exp(-0.5*kernel^2)    

    normalizer <- sqrt(2*pi)*y*comply*lambda # last lambda term missing in Johnson but not in wolfram
    
    
    dens/normalizer
    
    }


#
dJohnsonSL <- function(x, gamma, delta, lambda, csi){

    y <- (x - csi)/lambda   # y lognormal,  y > 0

if (any(y < 0) ) {
     k <- length(y[y < 0])
     lb <- min(y[y > 0], na.rm = T)
        y[y < 0] <- seq(10e-10, lb, length = k +1)[1:k]
   y <- sort(y)
   }
    
    kernel <- gamma + delta* log(y)   # z variate         

    dens <- delta*exp(-0.5*kernel^2)

    normalizer <- sqrt(2*pi)*y*lambda  # last lambda term missing in Johnson but not in wolfram        
    
  dens/normalizer

 }

#
dJohnsonSU <- function(x, gamma, delta, lambda, csi){

   ## trans <- function(y) asinh(y) 

  ##   y <- (x - csi)/lambda
    
  ##   kernel <- gamma + delta*trans(y)

  ##   dens <- delta*exp(-0.5*kernel^2)

  ##   a <- (x-csi)^2
  ##   b <- lambda^2
    
  ## normalizer <- sqrt(2*pi)*sqrt(a+b)

    
     trans <- function(y) log(y + sqrt(1+y^2) )
    
  y <- (x - csi)/lambda  
    kernel <- gamma + delta*trans(y)   # z variate
    dens <- delta*exp(-0.5*kernel^2)
    normalizer <- sqrt(2*pi)*sqrt(1+y^2)*lambda  # last lambda term missing in Johnson
    
    dens/normalizer
    
  }
#


### deprecated 

dJohnsonZ <- function(x, itype, gamma, delta, lambda, csi){

 z <- zJohnsonDistribution(x, itype, gamma, delta, lambda, csi)

  density <- switch(itype,

      1 ==  delta*exp(-0.5*z^2)/(sqrt(2*pi)*((x-csi)/lambda)),               
    2 == delta*exp(-0.5*z^2)/(sqrt(2*pi)*sqrt(1+((x-csi)/lambda)^2)*lambda),
    3 == delta*exp(-0.5*z^2)/(((x-csi)/lambda)*(1-((x-csi)/lambda))*lambda),
     4 == dnorm(x, lambda, csi)
                    )
            
   return( density )

   }
#


### DEPRECATED:

## dJohnsonZ <- function(x, itype, gamma, delta, lambda, csi){

##  z <- zJohnsonDistribution(x, itype, gamma, delta, lambda, csi)

##     if ( ( itype == 1 | itype ==3 ) & any(z < 0)) {

##  k <- length(z[z < 0])
##      lb <- min(z[z > 0], na.rm = T)
##         z[z < 0] <- seq(10e-10, lb, length = k +1)[1:k]
##         z <- sort(z)

##     }
            
##     percentiles <- pnorm(z)

##     density <- c(0, diff(percentiles))

##     density

##    }
## #

## value 'indirect' of argument mapping is deprecated


dJohnson <- function(x, itype, gamma, delta, lambda, csi, mapping = c("direct", "indirect")){

    mapping <- match.arg(mapping)

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



