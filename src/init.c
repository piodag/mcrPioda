#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/


/* .C calls */
extern void calc_Deming(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calc_Linreg(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calc_PaBa(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PaBaLargeData(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void calc_DemReg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP calcAngleMat(SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(ktau)(void *, void *, void *, void *);


static const R_CMethodDef CEntries[] = {
    {"calc_Deming",   (DL_FUNC) &calc_Deming,   13},
    {"calc_Linreg",   (DL_FUNC) &calc_Linreg,    9},
    {"calc_PaBa",     (DL_FUNC) &calc_PaBa,     11},
    {"PaBaLargeData", (DL_FUNC) &PaBaLargeData, 11},
    
    {"calc_Deming",   (DL_FUNC) &calc_Deming,   13},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"calcAngleMat", (DL_FUNC) &calcAngleMat, 3},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"ktau", (DL_FUNC) &F77_NAME(ktau), 4},
    {NULL, NULL, 0}
};

void R_init_mcr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
