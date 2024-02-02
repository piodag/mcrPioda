###############################################################################
##
## mcPaBa.R
##
## Functions for computing Passing Bablok regression.
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

#' Calculate Matrix of All Pair-wise Slope Angles
#' 
#' This is a very slow R version. It should not be called except for debugging purposes. 
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param posCor should the algorithm assume positive correlation, i.e. symmetry around slope 1? 
#' @return Upper triangular matrix of slopes for all point combinations. Slopes in radian.
mc.calcAngleMat.R <- function(X, Y, posCor=TRUE)
{
	## Check validity of parameters
	stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X)==length(Y))
	
	## Set up
	nData <- length(X)
	angM <- matrix(NA,nrow=nData,ncol=nData)
	
	## Calculate angles
	for(j in 1:(nData-1)) 
    {
		for(k in (j+1):nData) 
        {
			## Calculate x and y difference 
			dx <- calcDiff(X[k],X[j])
			dy <- calcDiff(Y[k],Y[j])	
			if(dx!=0) {
			## x and y != 0 
				angM[j,k] <- atan(dy/dx)
			}
			else if(dy!=0) 
            {
			## x==0, y!=0
			## only positive infinity for pos correlated
				if(posCor) 
                    angM[j,k] <- pi/2
				else 
                    angM[j,k] <- -pi/2
			}
			## else dx == dy == 0 => leave angM as NA, slope undefined 
		}
	}
	return(angM)
}

#' Calculate Matrix of All Pair-wise Slope Angles
#' 
#' This version is implemented in C for computational efficiency.
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param posCor should algorithm assume positive correlation, i.e. symmetry around slope 1? 
#' @return Upper triangular matrix of slopes for all point combinations. Slopes in radian.
mc.calcAngleMat <- function(X,Y,posCor=TRUE) 
{
    ## Check validity of parameters
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X)==length(Y))
    stopifnot(!is.na(posCor))
    ## Call C function
	ans <- .Call("calcAngleMat",X,Y,posCor)
	return(ans)
}


#' Passing-Bablok Regression
#'
#' @param angM upper triangular matrix of slopes for all point combinations (optional). Slopes in radian.
#' @param X measurement values of reference method
#' @param Y measurement values of test method
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level
#' @param posCor should algorithm assume positive correlation, i.e. symmetry around slope 1?
#' @param calcCI should confidence intervals be computed?
#' @param slope.measure angular measure of pairwise slopes  (see \code{\link{mcreg}} for details).\cr   
#'          \code{"radian"} - for data sets with even sample numbers median slope is calculated as average of two central slope angles.\cr
#'          \code{"tangent"} - for data sets with even sample numbers median slope is calculated as average of two central slopes (tan(angle)).\cr
#' @return Matrix of estimates and confidence intervals for intercept and slope. No standard errors provided by this algorithm.  


mc.paba <- function(angM = NULL, X, Y, alpha = 0.05, posCor = TRUE, 
					calcCI = TRUE, slope.measure = c("radian", "tangent")) 
{
	## Check validity of parameters
	slope.measure <- match.arg(slope.measure)
    stopifnot(length(X) == length(Y))
	stopifnot(is.numeric(X))
	stopifnot(is.numeric(Y))
    stopifnot(is.logical(posCor))
	if(slope.measure == "radian"){
		slope.measure <- 0
	}else if(slope.measure == "tangent"){
		slope.measure <- 1
	}else{
		stop(" slope.measure must be one of 'radian' or 'tangent'. \n")
	}
	parallel <- options()$parallel
	if(!is.null(parallel) && parallel == TRUE){
		Ncpu <- detectCores()/2
	}else{
		Ncpu <- 1
	}
	
	######################################################
	###  call C-function
	nX <- length(X)
	intercept <- slope <- 0
	seSlope <- seIntercept <- c(0,0)
	
	model.paba <- .C("calc_PaBa", 
					X = as.numeric(X), Y = as.numeric(Y), N = as.integer(nX), 
					intercept = as.numeric(intercept), slope = as.numeric(slope), 
					seIntercept = as.numeric(seIntercept), seSlope = as.numeric(seSlope), 
					pQuantile = as.numeric(qnorm(1-alpha/2)), pCor = as.integer(posCor), 
					tangent = as.integer(slope.measure), Ncpu = as.integer(Ncpu),
					PACKAGE="mcrPioda")
	
    ##
	## Confidence intervals for slope
	##
	if(calcCI){
		if(is.na(model.paba$seIntercept[1])){
			model.paba$seIntercept[1] <- -Inf
		}
		if(is.na(model.paba$seIntercept[2])){
			model.paba$seIntercept[2] <- Inf
		}
		if(is.na(model.paba$seSlope[1])){
			model.paba$seSlope[1] <- -Inf
		}
		if(is.na(model.paba$seSlope[2])){
			model.paba$seSlope[2] <- Inf
		}
	}else{
		## No theoretical CIs computed, use resampling to get CIs
		model.paba$seSlope[1:2] <- as.numeric(NA)
		model.paba$seIntercept[1:2] <- as.numeric(NA)
	}
	
	######################################################
	## Prepare result matrix
	
	rmat <- matrix(nrow=2,ncol=4)
	rownames(rmat) <- c("Intercept","Slope")
	colnames(rmat) <- c("EST","SE","LCI","UCI")
	rmat[,1] <- c(model.paba$intercept, model.paba$slope)
	rmat[,2] <- NA
	rmat["Intercept","LCI"] <- model.paba$seIntercept[1]
	rmat["Intercept","UCI"] <- model.paba$seIntercept[2]
	rmat["Slope","LCI"] <- model.paba$seSlope[1]
	rmat["Slope","UCI"] <- model.paba$seSlope[2]
	return(rmat)
}
	
