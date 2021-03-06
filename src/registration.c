#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>


extern SEXP dupRowNumMat(SEXP, SEXP);
extern SEXP getGlmBias(SEXP, SEXP, SEXP, SEXP);
extern void initQRdecomp(int*, int*);
extern void finalQRdecomp(void);

static R_CallMethodDef callMethods[]  = {
  {"dupRowNumMat", (DL_FUNC) &dupRowNumMat, 2},
  {"getGlmBias", (DL_FUNC) &getGlmBias, 4},
  {NULL, NULL, 0}
};
static R_CMethodDef cMethods[] = {
  {"initQRdecomp", (DL_FUNC) &initQRdecomp, 2},
  {"finalQRdecomp", (DL_FUNC) &finalQRdecomp, 0},
  {NULL, NULL, 0}
};

void R_init_BinQuasi(DllInfo *info)
{
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
