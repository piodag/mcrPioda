//
//C implementation of Deming regression
//

#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "utils.functions.h"




void calc_Deming(const double *X, const double *Y, int *nX,
				double *error_ratio, double *intercept, double *slope,
				double *se_intercept, double *se_slope, int *Wmode,
				int *itermax, double *threshold, double *W, double *xw)
{
	double meanX, meanY, u, q, p, r, x2;
	int i, n;

	n = *nX;
	double lambda = *error_ratio;
	meanX = meanY = u = q = p = r = x2 = 0;

	// calculation mean of X and Y
	get_mean(X, nX, &meanX);
	get_mean(Y, nX, &meanY);
	*xw = meanX;

	// calculation of u, q, p, r
	for(i = 0; i < n; i++){
		x2 += X[i]*X[i];
	    u += (X[i] - meanX) * (X[i] - meanX);
	    q += (Y[i] - meanY) * (Y[i] - meanY);
	    p += (X[i] - meanX) * (Y[i] - meanY);
	}
	r = p / sqrt(u*q);


	// Estimated points
	// [ Ref. K.Linnet. Estimation of the linear relationship between
    //        the measurements of two methods with  Proportional errors.
    //        STATISTICS IN MEDICINE, VOL. 9, 1463-1473 (1990)].

	slope[0] = ((lambda*q - u) + sqrt(pow((u - lambda*q), 2) + 4*lambda*pow(p,2))) / (2*lambda*p);
	intercept[0] = meanY - slope[0]*meanX;


	// if mode = 1, then calculate Deming with weights

	if(Wmode[0] < 1){
		// Standard error
		// [Ref. Strike, P. W. (1991) Statistical Methods in Laboratory Medicine.
		//       Butterworth-Heinemann, Oxford ].
		double r_2 = pow(r, 2);
		se_slope[0] = sqrt(pow(slope[0],2) * (diff(1, r_2) / r_2) / (n - 2));
		se_intercept[0] = sqrt(pow(se_slope[0],2) * (x2/n));
	}else{
		//
		// Iterative Algorithm

		int j = 0;
		double d = 0, XHAT = 0, YHAT = 0;
		double B0 = 0, B1 = 0, U = 0, Q = 0, P = 0;
		i = 0;

		//do loop at least once
		 while(i < itermax[0]){
			double XW = 0, YW = 0, sumW = 0;
			// Calculation of weights
			for(j = 0; j < n; j++){
					d = Y[j] - (intercept[0] + slope[0]*X[j]);
					XHAT = X[j] + (lambda*slope[0]*d / (1 + lambda*pow(slope[0],2)));
					YHAT = Y[j] - (d/(1 + lambda*pow(slope[0],2)));
					W[j] = 1/pow((XHAT + lambda*YHAT)/(1 + lambda), 2);
					sumW += W[j];
					XW += W[j] * X[j];
					YW += W[j] * Y[j];
			}
			XW = XW/sumW;
			YW = YW/sumW;
			*xw = XW;

			//Calculation of regression coefficients
			U = 0, Q = 0, P = 0;
			for(j = 0; j < n; j++){
				    U += W[j]*((X[j] - XW) * (X[j] - XW));
				    Q += W[j]*((Y[j] - YW) * (Y[j] - YW));
				    P += W[j]*((X[j] - XW) * (Y[j] - YW));
			}
			// Estimated points
			B1 = (lambda*Q - U + sqrt(pow((U - lambda*Q), 2) + 4*lambda*pow(P,2))) / (2*lambda*P);
			B0 = YW - B1*XW;

			if(fabs(slope[0] - B1) < threshold[0] && fabs(intercept[0] - B0) < threshold[0]){
				// set new values
				slope[0] = B1;
				intercept[0] = B0;
				break;
			}

			// set new values
			slope[0] = B1;
			intercept[0] = B0;
			i++;
		}

		itermax[0] = i + 1;

		// set standard error to 0
		se_slope[0] = 0;
		se_intercept[0] = 0;
	}
}





void calc_Linreg(const double *X, const double *Y, int *nX,
				double *intercept, double *slope,
				double *se_intercept, double *se_slope,
				double *W, double *xw)
{
	double meanX, meanY, SXXW, SXYW, SXX;
	double SX = 0, SY = 0, SXY = 0;
	int i, n;

	n = *nX;
	meanX = meanY = SXXW = SXYW = SXX = 0;

	// calculation mean of X and Y
	get_Wmean(X, W, nX, &meanX, xw);
	get_Wmean(Y, W, nX, &meanY, xw);


	// calculation of SXXW and SXYW
	for(i = 0; i < n; i++){
		SXX += W[i]*pow(X[i], 2);
		SX += W[i]*X[i];
		SY += W[i]*Y[i];
		SXY += W[i]*Y[i]*X[i];
	}
	SXXW = SXX - pow(SX, 2) / xw[0];
	SXYW = SXY - SX*SY / xw[0];

	// Point estimates
	slope[0] = SXYW/SXXW;
	intercept[0] = SY/ xw[0] - slope[0]*(SX/ xw[0]);

	// Standard errors
    double MSEW = 0;
    for(i = 0; i < n; i++){
    	MSEW += W[i]*pow(Y[i] - (intercept[0] + slope[0]*X[i]),2);
    }
    MSEW = 1.0 / (n-2)*MSEW;

    se_slope[0] = sqrt(MSEW/SXXW);
	se_intercept[0] = se_slope[0]*sqrt(SXX/xw[0]);
	xw[0] = SX/ xw[0];
}






void calc_PaBa(const double *X, const double *Y, const int *nX,
				double *intercept, double *slope,
				double *se_intercept, double *se_slope,
				double *pQuantile, int *pCor, int *tangent, int* Ncpu)
{
	double mean = 0, sigma = 0;
	double tempL = 0, tempR = 0;

	int nNeg = 0, nNeg2 = 0, nPos = 0, nPos2 = 0, nAllItems = 0, Offset = 0;
	int nValIndex2 = 0, LowestIdx = 0, LCLundef = 0, UCLundef = 0;
	int i;
	int N = *nX;

	//omp_set_num_threads(*Ncpu);


	if(N*(N-1)/2 < 2e+8){
		/* allocate storage: all helpers, initially zero */
		double *ans = Calloc(N*(N-1)/2, double);                                  /* (n(n-1)/2) x 1  */

		/*  calculate all angles */
		calc_AngleMat_opt(X, Y, nX, pCor, ans, &nAllItems,
								&nNeg, &nNeg2, &nPos, &nPos2,
								&mean, &sigma);
		if(*pCor == 1){                                                       /* compute Bin-index of the offsetted median */
			Offset = nNeg + nNeg2;
		} else{
			Offset = (-1) * (nPos + nPos2);
		}

		/* determine slope of the regression line */
		nValIndex2 = nAllItems + Offset;
		int half = (int)(nValIndex2 + 1)/2;                                   /* integer part, i.e. floor() */


		if(nValIndex2 % 2 == 0) {                                             /* nValIndex is an even number */
			tempL = quickselect(ans, nAllItems, half - 1);
			tempR = quickselect(ans, nAllItems, half);

			if(*tangent == 1){
				*slope = (tan(tempL) + tan(tempR)) / 2;
			} else {
				*slope = tan((tempL + tempR) / 2);
			}
		} else {                                                              /* nValIndex2 is an odd number */
			*slope = tan(quickselect(ans, nAllItems, half-1));
		}


		/* calculate confidence interval for the slope, the CI for the intercept will be calculated on R-side */
		/******************/
		/* Lower CI-Bound */
		double dConf = (*pQuantile) * sqrt( ((double)*nX) *(((double)*nX)-1.) * (2. * ((double)*nX) + 5) / 18.);
		dConf = round(dConf);                                                /* round to the nearest integer */
		int nInd     = (int) (nAllItems - dConf + Offset);                   /* as in the exact algo */
		int nItems   = (int)((nInd + 1)/2);                                  /* !!! added + 1 as in the exact algo (nInd+1L)%/%2) */

		if(*pCor == 1){
			LowestIdx = 2*(nNeg - nNeg2) + 1;
		} else {
			LowestIdx =  2*(nPos - nPos2) + 1;
		}

		if(nItems >= LowestIdx){                                             /* otherwise the -INF init-value coming from R remains unchanged */
			if( nInd % 2 == 0 ) {                                            /* nInd is an even number */
				tempL = quickselect(ans, nAllItems, nItems - 1);
				tempR = quickselect(ans, nAllItems, nItems);
				if(*tangent == 1){
					se_slope[0] = (tan(tempL) + tan(tempR)) / 2;
				} else {
					se_slope[0] = tan((tempL + tempR) / 2);
				}
			} else {
				se_slope[0] = tan(quickselect(ans, nAllItems, nItems-1));     /* nInd is an odd number */
			}
		}else{
			LCLundef = 1;
		}
		/******************/
		/* Upper CI-Bound */
		nInd   = (int) (nAllItems + dConf + Offset);
		nItems = (int) ((nInd + 1)/2);


		if(nItems <= nAllItems){                                            /* otherwise the INF init-value coming from R remains unchanged */
			if( nInd % 2 == 0 ){                                            /* nInd is an even number */
				tempL = quickselect(ans, nAllItems, nItems - 1);
				tempR = quickselect(ans, nAllItems, nItems);
				if(*tangent == 1){
					se_slope[1] = (tan(tempL) + tan(tempR)) / 2;
				} else {
					se_slope[1] = tan((tempL + tempR) / 2);
				}
		    } else {
		    	se_slope[1] = tan(quickselect(ans, nAllItems, nItems-1));      /* nInd is an odd number */
		    }
		}else{
			UCLundef = 1;
		}
		if(LCLundef || UCLundef){                                            /* signal error(s) computing CI-bounds for slope */
		    se_slope[0] = se_slope[1] = NA_REAL;
		}
		Free(ans);

		/******************/
		/* calculate intercept and confidence interval */
		/******************/
		double *mcres_intercept = Calloc(N, double);
		for(i=0; i < N; i++){
			mcres_intercept[i] = diff(Y[i], slope[0]*X[i]);
		}
		half = (int)(N + 1)/2;
		if( N % 2 == 0 ){
			tempL = quickselect(mcres_intercept, N, half - 1);
			tempR = quickselect(mcres_intercept, N, half);
			*intercept = (tempL + tempR) / 2;
		}else{
			*intercept = quickselect(mcres_intercept, N, half-1);
		}

		/******************/
		/* Lower CI-Bound */
		if(LCLundef == 1){
			se_intercept[0] = NA_REAL;
		}else{
			for(i=0; i < N; i++){
				mcres_intercept[i] = diff(Y[i], se_slope[1]*X[i]);
			}
			if( N % 2 == 0 ){
				tempL = quickselect(mcres_intercept, N, half - 1);
				tempR = quickselect(mcres_intercept, N, half);
				se_intercept[0] = (tempL + tempR) / 2;
			}else{
				se_intercept[0] = quickselect(mcres_intercept, N, half-1);
			}
		}
		/******************/
		/* Upper CI-Bound */
		if(UCLundef == 1){
			se_intercept[1] = NA_REAL;
		}else{
			for(i=0; i < N; i++){
				mcres_intercept[i] = diff(Y[i], se_slope[0]*X[i]);
			}
			if( N % 2 == 0 ){
				tempL = quickselect(mcres_intercept, N, half - 1);
				tempR = quickselect(mcres_intercept, N, half);
				se_intercept[1] = (tempL + tempR) / 2;
			}else{
				se_intercept[1] = quickselect(mcres_intercept, N, half-1);
			}
		}
		Free(mcres_intercept);
	} else {
		*slope = 0;
		*intercept = 0;
	}
}







