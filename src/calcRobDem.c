//
// C implementation of M and MM-Deming regression
// Copyright: G. Pioda, 2024, gfwp@ticino.com
// Inspired from mcr package

#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
#include "utils.functions.h"

void calc_MDem(const double *X, const double *Y, int *nX,
               double *error_ratio, double *intercept, double *slope,
               double *se_intercept, double *se_slope, int *Wmode,
               int *itermax, double *threshold, double *W,
               double *xw, double *kM)
{
  double meanX, meanY, u, q, p;
  int i, n;
  
  n = *nX;
  double lambda = *error_ratio;
  meanX = meanY = u = q = p = 0;
  
  // calculation mean of X and Y
  get_mean(X, nX, &meanX);
  get_mean(Y, nX, &meanY);
  *xw = meanX;
  
  // calculation of u, q, p, r
  for(i = 0; i < n; i++){
    u += (X[i] - meanX) * (X[i] - meanX);
    q += (Y[i] - meanY) * (Y[i] - meanY);
    p += (X[i] - meanX) * (Y[i] - meanY);
  }
  
  // Estimated points
  // [ Ref. K.Linnet. Estimation of the linear relationship between
  //        the measurements of two methods with  Proportional errors.
  //        STATISTICS IN MEDICINE, VOL. 9, 1463-1473 (1990)].
  
  slope[0] = ((lambda*q - u) + sqrt(pow((u - lambda*q), 2) + 4*lambda*pow(p,2))) / (2*lambda*p);
  intercept[0] = meanY - slope[0]*meanX;
  
  int j = 0;
  double d = 0, XHAT = 0, YHAT = 0;
  double B0 = 0, B1 = 0, U = 0, Q = 0, P = 0;
  i = 0;
  
  double mad = 0;
  double k = *kM;
  double euclid[n];
  double work[n];
  
  //do loop at least once
  while(i < itermax[0]){
    double XW = 0, YW = 0, sumW = 0;
    
    // Calculation of weights
    for(j = 0; j < n; j++){
      d = Y[j] - (intercept[0] + slope[0]*X[j]);
      XHAT = X[j] + (lambda*slope[0]*d / (1 + lambda*pow(slope[0],2)));
      YHAT = Y[j] - (d/(1 + lambda*pow(slope[0],2)));
      euclid[j] = sqrt(pow((X[j]-XHAT),2)+pow(Y[j]-YHAT,2));
    }
     
  // Stop loop, calculate MAD
   mad = gsl_stats_mad(euclid,1,n,work);
        
    for(j = 0; j < n; j++){  
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
    
    // This part is set to break resonance in the convergence. It is crucial to 
    // average both slope AND intercept. Makes convergence much faster.
    
    if(i % 2){
      slope[0] = (slope[0] + B1) / 2;
      intercept[0] = (intercept[0] + B0) / 2;
      i++;
      
    }else{
      slope[0] = B1;
      intercept[0] = B0;
      i++;
    }

  }
  
  itermax[0] = i + 1;
  
}


void calc_MMDem(const double *X, const double *Y, int *nX,
               double *error_ratio, double *intercept, double *slope,
               double *se_intercept, double *se_slope, int *Wmode,
               int *itermax, double *threshold, double *W,
               double *xw, double *kM, double *tauMM)
{
  int i, n;
  
  n = *nX;
  double lambda = *error_ratio;
  
  // Estimated points
  // [ Ref. K.Linnet. Estimation of the linear relationship between
  //        the measurements of two methods with  Proportional errors.
  //        STATISTICS IN MEDICINE, VOL. 9, 1463-1473 (1990)].
  
  int j = 0;
  double d = 0, XHAT = 0, YHAT = 0;
  double B0 = 0, B1 = 0, U = 0, Q = 0, P = 0;
  i = 0;
  
  double mad = 0;
  double tau = *tauMM;
  double euclid[n];
  double work[n];
  
  //do loop at least once
  while(i < itermax[0]){
    
    double XW = 0, YW = 0, sumW = 0;
    
    for(j = 0; j < n; j++){
      d = Y[j] - (intercept[0] + slope[0]*X[j]);
      XHAT = X[j] + (lambda*slope[0]*d / (1 + lambda*pow(slope[0],2)));
      YHAT = Y[j] - (d/(1 + lambda*pow(slope[0],2)));
      euclid[j] = sqrt(pow((X[j]-XHAT),2)+pow(Y[j]-YHAT,2));
    }
    
    // Stop loop, calculate MAD
    
     mad = gsl_stats_mad(euclid,1,n,work);
    
    // Calculation of weights
    
    for(j = 0; j < n; j++){
      
     if( euclid[j]/mad <= tau) {
       W[j] = pow(1-pow((euclid[j]/mad)/tau,2),2);
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
    
    // This part is set to break resonance in the convergence. It is crucial to 
    // average both slope AND intercept. Makes convergence much faster.
    
    if(i % 2){
      slope[0] = (slope[0] + B1) / 2;
      intercept[0] = (intercept[0] + B0) / 2;
      i++;
      
    }else{
      slope[0] = B1;
      intercept[0] = B0;
      i++;
    }
    
  }
  
  itermax[0] = i + 1;
  
}


