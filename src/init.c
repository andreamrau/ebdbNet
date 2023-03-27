#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern SEXP RunWrapGen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP,SEXP, SEXP, SEXP);
  
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
/*  Convert package to use registration */
/*  Allow registration by declaring R_init_splines visible */
static const R_CMethodDef CEntries[] = {
  {"RunWrapGen", (DL_FUNC) &RunWrapGen, 28},
  {NULL, NULL, 0}
};

void 
/* attribute_visible */
  R_init_ebdbNet(DllInfo *dll)
    {
      R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
      R_useDynamicSymbols(dll, FALSE);
      R_forceSymbols(dll, TRUE);
    }




