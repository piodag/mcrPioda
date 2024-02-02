###############################################################################
##
## mcDeming.R
##
## Function for computing Deming regression based method comparison.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
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

#' Calculate Unweighted Deming Regression and Estimate Standard Errors
#' 
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param error.ratio ratio of measurement error of reference method to measurement error of test method.
#' @return a list with elements
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{se.b0}{respective standard error of intercept.}
#'  \item{se.b1}{respective standard error of slope.}
#'  \item{xw}{average of reference method values.}
#' @references  Linnet K.
#'              Evaluation of Regression Procedures for Methods Comparison Studies.
#'              CLIN. CHEM. 39/3, 424-432 (1993).
#' 
#'              Linnet K.
#'              Estimation of the Linear Relationship between the Measurements of two Methods with Proportional Errors.
#'              STATISTICS IN MEDICINE, Vol. 9, 1463-1473 (1990).




mc.deming <- function(X, Y, error.ratio)
{
    # Check validity of parameters
	nX <- length(X)
	nY <- length(Y)
	
    stopifnot(!is.na(X))
    stopifnot(!is.na(Y))
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
	stopifnot(nX > 0)
    stopifnot(nX == nY)
    stopifnot(!is.na(error.ratio))
    stopifnot(is.numeric(error.ratio))
    stopifnot(error.ratio > 0)
    stopifnot(length(error.ratio) > 0)
    
	
	######################################################
	###  call C-function
	intercept <- slope <- seIntercept <- seSlope <- 0
	xw <- maxit <- threshold <- 0
	
	### mode = 0 - Deming regression
	### mode = 1 - WDeming regression
	mode <- 0
	W <- rep(1, nX)
	
	model.Deming <- .C("calc_Deming", 
						x = as.numeric(X), y = as.numeric(Y), 
						n = as.integer(nX), 
						error_ratio = as.numeric(error.ratio), 
						intercept = as.numeric(intercept), 
						slope = as.numeric(slope), 
						seIntercept = as.numeric(seIntercept), 
						seSlope = as.numeric(seSlope), 
						mode = as.integer(mode), 
						maxit = as.integer(maxit), 
						threshold = as.numeric(threshold), 
						W = as.numeric(W), 
						xw = as.numeric(xw), 
						PACKAGE="mcrPioda")

	list(b0 = model.Deming$intercept, b1 = model.Deming$slope, 
		se.b0 = model.Deming$seIntercept, se.b1 = model.Deming$seSlope, 
		xw = model.Deming$xw,  weight = model.Deming$W)
}
