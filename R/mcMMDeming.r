###############################################################################
##
## mcWDeming.r
##
## Function for computing weighted deming regression for two methods with  proportional errors.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
## Copyright (C) 2020 Giorgio Pioda for the MM-Deming
## Copyright (C) 2024 Giorgio Pioda for C version of MM-Deming C algorithm
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Calculate Weighted Deming Regression
#'
#' Calculate weighted deming regression with iterative algorithm suggested by Linnet.
#' This algorithm is avalaible only for positive values. But even in this case there is no guarantee that
#' the algorithm always converges.
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param error.ratio ratio between squared measurement errors of reference- and test method, necessary for Deming regression (Default is 1).
#' @param iter.max maximal number of iterations.
#' @param threshold threshold value.
#' @return a list with elements
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{xw}{average of reference method values.}
#'  \item{iter}{number of iterations.}
#' @references  Linnet K.
#'              Evaluation of Regression Procedures for Methods Comparison Studies.
#'              CLIN. CHEM. 39/3, 424-432 (1993).
#'
#'              Linnet K.
#'              Estimation of the Linear Relationship between the Measurements of two Methods with Proportional Errors.
#'              STATISTICS IN MEDICINE, Vol. 9, 1463-1473 (1990).
mc.mmdemingConstCV <- function(X, Y, error.ratio, iter.max=120, threshold=0.000001)
{
  # Check validity of parameters

  stopifnot(is.numeric(X))
  stopifnot(is.numeric(Y))
  stopifnot(length(X) == length(Y))
  stopifnot(is.numeric(error.ratio))
  stopifnot(error.ratio > 0)
  stopifnot(is.numeric(iter.max))
  stopifnot(round(iter.max) == iter.max)
  stopifnot(iter.max > 0)
  stopifnot(is.numeric(threshold))
  stopifnot(threshold >= 0)

  # This algorithm often doesn't converge if there are negative
  # measurements in data set

#  if (min(X)<0 | min(Y)<0)
#  {
#    return(paste("There are negative values in data set."))
#  }
#  else
#  {
    # 1. Calculate  initial values
    #    (point estimates of unweighted Deming regression)

    #n <- length(X)
    nX <- length(X)

    if (nX >= 100){
      start.n<-20
    } else if(nX >= 76){
      start.n<-50
    } else if (nX >= 46){
      start.n<-100
    } else if (nX >= 36) {
      start.n<-250
    } else {
      #message("there is no convergence warranty below 36 samples")
      start.n<-250
    }

    #mX <- mean(X)
    #mY <- mean(Y)
    #u <- sum((X-mX)^2)
    #q <- sum((Y-mY)^2)
    #p <- sum((X-mX)*(Y-mY))

    ## initial values, b1 as the straight bisector of the two slopes obtained by the robust covariance matrix
    ## and b0 as the derived intercept using its centers. The scale parameter is then obtained with the
    ## mad() of the euclidean calculated residuals.

    #b1 <- ((error.ratio*q-u)+sqrt((u-error.ratio*q)^2+4*error.ratio*p^2))/(2*error.ratio*p)
    #b0 <- mean(Y)-b1*mean(X)
    
    seIntercept <- seSlope <- intercept <- slope <-  0
    mode <- ltsScale <- 1
    xw <- 0
    maxit <- iter.max
    W <- rep(1, nX)

    cov.sest<-rrcov::CovSest(cbind(X,Y),nsamp=start.n,bdp=0.5)
    b1<-mean(c(cov.sest$cov[2,1]/cov.sest$cov[1,1],1/(cov.sest$cov[2,1]/cov.sest$cov[2,2])))
    b0<-cov.sest$center[2]-b1*cov.sest$center[1]

    d <- Y-(b0+b1*X)
    XHAT <- X+(error.ratio*b1*d/(1+error.ratio*b1^2))
    YHAT <- Y-(d/(1+error.ratio*b1^2))
    euclid.d<-sqrt((X-XHAT)^2+(Y-YHAT)^2)
    scale.lts <- mad(euclid.d)
    stdeuclid.d <- euclid.d/scale.lts
    
    k<-4.685
    W<-ifelse(abs(stdeuclid.d) <= k, (1-(stdeuclid.d/k)^2)^2, 0)
    
    
    # If the sum of the weights is too low, use Rocke covariance method
    
    if(sum(W) > (0.4 * nX)) {
      
      model.Deming <- .C("calc_MMDem", 
                         x = as.numeric(X), y = as.numeric(Y), 
                         n = as.integer(nX), 
                         error_ratio = as.numeric(error.ratio), 
                         intercept = as.numeric(b0), 
                         slope = as.numeric(b1), 
                         seIntercept = as.numeric(seIntercept), 
                         seSlope = as.numeric(seSlope), 
                         mode = as.integer(mode), 
                         maxit = as.integer(iter.max), 
                         threshold = as.numeric(threshold), 
                         W = as.numeric(W),
                         xw = as.numeric(xw),
                         ltsScale = as.numeric(scale.lts),
                         PACKAGE="mcrPioda")
      
      if (model.Deming$maxit >= maxit) {
        warning(paste("no convergence after", maxit, "iterations"))
      }
      
      list(b1 = model.Deming$slope, b0 = model.Deming$intercept,
           iter = model.Deming$maxit,
           xw = model.Deming$xw,  weight = model.Deming$W)
      
    } else {
      
      warning(paste("Default CovSest weighting factors not meaningful,
                    Rocke method used instead","b1:",b1,"b0:",b0,"B1:",B1,"sum W:",sum(W)))
      
      cov.sest<-rrcov::CovSest(cbind(X,Y),method="rocke")
      b1<-mean(c(cov.sest$cov[2,1]/cov.sest$cov[1,1],1/(cov.sest$cov[2,1]/cov.sest$cov[2,2])))
      b0<-cov.sest$center[2]-b1*cov.sest$center[1]
      
      d <- Y-(b0+b1*X)
      XHAT <- X+(error.ratio*b1*d/(1+error.ratio*b1^2))
      YHAT <- Y-(d/(1+error.ratio*b1^2))
      euclid.d<-sqrt((X-XHAT)^2+(Y-YHAT)^2)
      scale.lts <- mad(euclid.d)
      stdeuclid.d<-euclid.d/scale.lts
      
      k<-4.685
      W<-ifelse(abs(stdeuclid.d) <= k, (1-(stdeuclid.d/k)^2)^2, 0)
      
    }
    
    if(sum(W) > (0.4 * nX)) {
      
      model.Deming <- .C("calc_MMDem", 
                         x = as.numeric(X), y = as.numeric(Y), 
                         n = as.integer(nX), 
                         error_ratio = as.numeric(error.ratio), 
                         intercept = as.numeric(intercept), 
                         slope = as.numeric(slope), 
                         seIntercept = as.numeric(seIntercept), 
                         seSlope = as.numeric(seSlope), 
                         mode = as.integer(mode), 
                         maxit = as.integer(iter.max), 
                         threshold = as.numeric(threshold), 
                         W = as.numeric(W), 
                         xw = as.numeric(xw),
                         ltsScale = as.numeric(scale.lts),
                         PACKAGE="mcrPioda")
      if (model.Deming$maxit >= maxit) {
        
        warning(paste("no convergence after", maxit, "iterations"))
      }
      
      list(b1 = model.Deming$slope, b0 = model.Deming$intercept,
           iter = model.Deming$maxit,
           xw = model.Deming$xw,  weight = model.Deming$W)
      
    } else {
      
      warning(paste("Not even Rocke weights meaningiful for convergency, convergency loop
                    not executed, Rocke startvalues passed","b1:",b1,"b0:",b0,"sum W:",sum(W)))
      
      ## Do nothing, b0 and b1 are set with Rocke starting values, just return the list
      
      list(b1 = b1, b0 = b0,
           iter = 1,
           xw = xw,  weight = W)
      
    } 
}    

