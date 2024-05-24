//
// C implementation of M and MM-Deming regression
// Copyright: G. Pioda, 2024, gfwp@ticino.com
// Inspired from mcr package


void calc_MDem(const double *X, const double *Y, int *nX,
               double *error_ratio, double *intercept, double *slope,
               double *se_intercept, double *se_slope, int *Wmode,
               int *itermax, double *threshold, double *W,
               double *xw, double *kM);
void calc_MMDem(const double *X, const double *Y, int *nX,
               double *error_ratio, double *intercept, double *slope,
               double *se_intercept, double *se_slope, int *Wmode,
               int *itermax, double *threshold, double *W,
               double *xw, double *kM, double *tauMM);
