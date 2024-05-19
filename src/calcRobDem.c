//
// C implementation of M and MM-Deming regression
//

#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
//#include <gsl.h>
#include "utils.functions.h"


void calc_MDem(const double *X, const double *Y, int *nX,
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
		
		double mad = 0;
		double k = 1.345;
		double euclid[n];
		double work[n];
		
		// first estimation of euclid distances, for estimating a starting MAD value
		
		for(j = 0; j < n; j++){
		  d = Y[j] - (intercept[0] + slope[0]*X[j]);
		  XHAT = X[j] + (lambda*slope[0]*d / (1 + lambda*pow(slope[0],2)));
		  YHAT = Y[j] - (d/(1 + lambda*pow(slope[0],2)));
		  euclid[j] = sqrt(pow((X[j]-XHAT),2)+pow(Y[j]-YHAT,2));
		}
		
		// MAD from GLS library
		
		mad = gsl_stats_mad(euclid,1,n,work);
		
		//do loop at least once
		 while(i < itermax[0]){
			double XW = 0, YW = 0, sumW = 0;
		   
			// Calculation of weights
			for(j = 0; j < n; j++){
					d = Y[j] - (intercept[0] + slope[0]*X[j]);
					XHAT = X[j] + (lambda*slope[0]*d / (1 + lambda*pow(slope[0],2)));
					YHAT = Y[j] - (d/(1 + lambda*pow(slope[0],2)));
					euclid[j] = sqrt(pow((X[j]-XHAT),2)+pow(Y[j]-YHAT,2));
					
					if( euclid[j] / mad <= k){
					  W[j] = 1;
					}else{
					  W[j] = k / (euclid[j]/mad);
					}
					
					sumW += W[j];
					XW += W[j] * X[j];
					YW += W[j] * Y[j];
					
			}
			XW = XW/sumW;
			YW = YW/sumW;
			*xw = XW;
			
			// Update MAD
			mad = gsl_stats_mad(euclid,1,n,work);

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




void calc_MMDem(const double *X, const double *Y, int *nX,
                 double *error_ratio, double *intercept, double *slope,
                 double *se_intercept, double *se_slope, int *Wmode,
                 int *itermax, double *threshold, double *W,
                 double *xw, double *ltsScale)
{
  // double meanX, meanY, u, q, p, r, x2;
  int i, n;
  
  n = *nX;
  double Scale = *ltsScale; 
  double lambda = *error_ratio;
  
  // meanX = meanY = u = q = p = r = x2 = 0;
  
  // calculation mean of X and Y
  // get_mean(X, nX, &meanX);
  // get_mean(Y, nX, &meanY);
  // *xw = meanX;
  
  // calculation of u, q, p, r
  // for(i = 0; i < n; i++){
  //  x2 += X[i]*X[i];
  //  u += (X[i] - meanX) * (X[i] - meanX);
  //  q += (Y[i] - meanY) * (Y[i] - meanY);
  //  p += (X[i] - meanX) * (Y[i] - meanY);
  //}
  // r = p / sqrt(u*q);
  
  
  // Estimated points
  // [ Ref. K.Linnet. Estimation of the linear relationship between
  //        the measurements of two methods with  Proportional errors.
  //        STATISTICS IN MEDICINE, VOL. 9, 1463-1473 (1990)].
  
  // slope[0] = ((lambda*q - u) + sqrt(pow((u - lambda*q), 2) + 4*lambda*pow(p,2))) / (2*lambda*p);
  // intercept[0] = meanY - slope[0]*meanX;
  
  
  // if mode = 1, then calculate Deming with weights
  
  // if(Wmode[0] < 1){
    // Standard error
    // [Ref. Strike, P. W. (1991) Statistical Methods in Laboratory Medicine.
    //       Butterworth-Heinemann, Oxford ].
  //  double r_2 = pow(r, 2);
  //  se_slope[0] = sqrt(pow(slope[0],2) * (diff(1, r_2) / r_2) / (n - 2));
  //  se_intercept[0] = sqrt(pow(se_slope[0],2) * (x2/n));
  // }else{
    //
    
    // Iterative Algorithm
    
    int j = 0;
    double d = 0, XHAT = 0, YHAT = 0;
    double B0 = 0, B1 = 0, U = 0, Q = 0, P = 0;
    i = 0;
    
    //double mad = 0;
    double k = 4.685;
    double euclid[n];
    double steuclid[n];
    
    //do loop at least once
    while(i < itermax[0]){
      double XW = 0, YW = 0, sumW = 0;
      
      // Calculation of weights
      
      for(j = 0; j < n; j++){
        d = Y[j] - (intercept[0] + slope[0]*X[j]);
        XHAT = X[j] + (lambda*slope[0]*d / (1 + lambda*pow(slope[0],2)));
        YHAT = Y[j] - (d/(1 + lambda*pow(slope[0],2)));
        euclid[j] = sqrt(pow((X[j]-XHAT),2)+pow(Y[j]-YHAT,2));
        steuclid[j] = euclid[j]/Scale;

        
        if( steuclid[j] <= k) {
          
          W[j] = pow(1-pow(steuclid[j]/k,2),2);
          
        } else {
          
          W[j] = 0;
        }
        
        sumW += W[j];
        XW += W[j] * X[j];
        YW += W[j] * Y[j];
        
      }
      
      XW = XW/sumW;
      YW = YW/sumW;
      *xw = XW;
      
      // Update MAD
      //mad = gsl_stats_mad(euclid,1,n,work);
      
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
 //  }
}
