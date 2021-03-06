## order_cppWrapper.R contains modules needed to perform permutation-based imposition of incomplete IPD correlation structure (this is *not* NORTA approach).
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







return.MC.extremal.correlation <- function( mc.extremal.corr,
                                           rho, extremal.corr)
{

     #ideally a list with several MC statistics for extremal corr

    if (!is.null(mc.extremal.corr) & is.list(mc.extremal.corr))
    {
        if(abs(rho) <= abs(mc.extremal.corr[[1]]) & sign(mc.extremal.corr[[1]]) == sign(rho) )
            out <- mc.extremal.corr[[1]]
        else
        {
            if(abs(rho) <= abs(mc.extremal.corr[[2]]) & sign(mc.extremal.corr[[2]]) == sign(rho) )
                out <- mc.extremal.corr[[2]]
            else
            {
                if(abs(rho) <= abs(mc.extremal.corr[[3]]) & sign(mc.extremal.corr[[3]]) == sign(rho) )
                    out <- mc.extremal.corr[[3]]
                else
                {
                    warning("|MC extremal correlation| < |observed correlation| OR inverted signs: NAs are expected")
                    out <- extremal.corr
                }
            }
        }
    }
    else
    {
        out <- extremal.corr

    }
    return(out)         

}




 

#' @title Imposition of incomplete IPD correlation structure.
#'
#' @description `rsearch2()` implements a permutation-based re-ordering of first-degree-correlated IPD marginals only w.r.t. some fixed reference variable.
#'
#'
#' @details This function also wraps the C++ module 'findzerocorr' (see source code). This function is *not* a NORTA approach and is primarily meant for experimental and comparative purposes w.r.t. the 'complete-correlation approach' implemented in the NORTA transformation.
#' 


### C++ while loop
# if not null mc.extremal.corr must be a (named) list having the extremal correlation elements in this exact following order: mc.average, mc.median, mc.extremal ...

rsearch2 <- function(xell, x, refstat,
                     mc.extremal.corr = NULL,
                     dz=0.001, random=TRUE,
                     rf.type=c("joint.rank", "joint.moment", "rank.corr", "moment.corr", "tau.corr"),
                     ties=c("min","average", "max") )
{

    N <- length(xell)
    
    if(N!=length(x))
        stop("unequal variables lengths xell x")

    rf.type <- match.arg(rf.type)
    ties <- match.arg(ties)    

    if(  grep("cor", rf.type)  )
        rho <- refstat  
    else
    {
        if( rf.type=="joint.rank" )
            rho <- ( refstat-(
                mean(
                    rank(
                    x,ties=ties),
                    na.rm=T
                )*mean(
                      rank(xell,ties=ties),
                      na.rm=T
                  )
            )
            )/(
                sd(
                    rank(x,
                         ties=ties),
                    na.rm=T
                )*sd(
                      rank(xell,ties),
                      na.rm=T
                  )
            )
        else
            rho <- (refstat-(
                mean(x,
                     na.rm=T
                     )*mean(
                           xell,
                           na.rm=T
                       )
            )
            )/(
                sd(x,
                   na.rm=T
                   )*sd(
                         xell,
                         na.rm=T
                     )
            )

    }

    f1 <- function(xell, x) return(rjm(xell, x, ties))
    f2 <- function(xell, x) return(sum(xell * x))
    f3 <- function(xell, x)  return(spcor(xell, x, ties))
    f4 <- function(xell, x) return(cor(xell, x))
    f5 <- function(xell, x)  return(cor(xell, x, method="kendall"))

    
    lookup <- list("joint.rank"=f1,
                   "join.moment"=f2,
                   "rank.corr"= f3,
                   "moment.corr"=f4,
                   "tau.corr"=f5
                   )
    call_func <- lookup[[rf.type]]

    if ( is.na(call_func(xell, x)) )  # check for existence of correlation between xell and x. if not exists, exit immediatly with NA values
        return( cbind( xell, NA ) )

    zero <- findzerocorr(xell, x, dz, rf.type)  # while loop in C++
    i <- zero$iter 

    z <- zero$vec 
    
    if(i >= N)
    {
        warning(paste(
            "no zero correlation found, minimum was",
            round(zero$cor,4)
        ),
        call.=FALSE
        )
    }
    
    u2 <- sort(x)
    if(rho < 0 )
        u2 <- rev(u2)
    u <- matrix(c( sort(xell), u2 ), ncol=2)
    extremal.corr <- call_func(u[,1], u[,2])    # empirical maximal correlation on observed vector. BEWARE: if u_i is binary and sd = 0, cor may be NA !!!

    if(abs(extremal.corr)> 0.94)
        extremal.corr <- sign(extremal.corr)*1
    if(abs(rho) > abs(extremal.corr)  | sign(extremal.corr)!=sign(rho) )
    {
        warning("|extremal correlation| < |observed correlation|: something is wrong, NAs expected")   # this seems more likely with moment corr when binary variables are employed .....
        extremal.corr <- return.MC.extremal.correlation( mc.extremal.corr, rho, extremal.corr)
        
    }
    
    alpha <- 1-(rho/extremal.corr)  # main relationship
    bi <- rbinom(N, 1, alpha )
    if(!random)
    {
        pos <- sample(1:N, N*alpha )
        bi <- 1:N
        bi[pos] <- NA  # be careful xell & x have no NA when input (they should not because simulated)!!!.
        bi <- ifelse(is.na(bi), 1, 0 )
    }

    res <- I(bi)*z + I(1-bi)*u   # mix
    print(
        noquote(
            paste(
                "zero correlation permutation found after",
                i,
                "iterations"
            )
        )
    )
    return(res)

}





##



