###############################################################################
##
## mcWDeming.r
##
## Function for computing weighted Deming regression for two methods with proportional errors.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
## Copyright (C) 2020 Giorgio Pioda for the M-Deming adaptation
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
#' Calculate weighted Deming regression with iterative algorithm suggested by Linnet.
#' This algorithm is available only for positive values. But even in this case there is no guarantee that
#' the algorithm always converges.
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param error.ratio ratio between squared measurement errors of reference- and test method, necessary for Deming regression (Default is 1).
#' @param iter.max maximal number of iterations.
#' @param threshold threshold value
#' @param kM Huber's k for the M weighting, default kM = 1.345
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


mc.mdemingConstCV <- function(X, Y, error.ratio, iter.max = 30, threshold = 0.000001, kM = 1.345)
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
  stopifnot(kM > 0)
  
  # This algorithm often doesn't converge if there are negative
  # measurements in data set
  
  if (min(X) < 0 | min(Y) < 0){
    return(paste("There are negative values in data set."))
  }else{
    ######################################################
    ###  call C-function
    intercept <- slope <- seIntercept <- seSlope <- 0
    xw <- 0
    maxit <- iter.max
    nX <- length(X)
    # kM <- 0.95106 rational, diagonal geometric alternative (1.345 * (sqrt(2)/2) )
    # kM <- 1.345
    
    ### mode = 0 - Deming regression
    ### mode = 1 - WDeming regression
    mode <- 1
    W <- rep(1, nX)
    
    model.Deming <- .C("calc_MDem", 
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
                       kM = as.numeric(kM),
                       PACKAGE="mcrPioda")
    if (model.Deming$maxit >= maxit) {
      warning(paste("no convergence after", maxit, "iterations","B1:",
                    model.Deming$slope,"B0:",model.Deming$intercept,"sum(W)",sum(model.Deming$W)))
    }
    
    list(b1 = model.Deming$slope, b0 = model.Deming$intercept,
         se.b0 = model.Deming$seIntercept, se.b1 = model.Deming$seSlope,
         xw = model.Deming$xw,  weight = model.Deming$W)
  }
}

