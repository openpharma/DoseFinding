#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void critfunc(double *, int *, int *, double *, int *, double *, double *, double *, double *, double *, double *, int *, int *, double *);
extern void getcomp(int *, int *, int *, int *, int *);
extern void sample(int *, int *, double *, double *, int *, const double *, double *, int *, int *, double *, double *, double *, int *, double *, int *, double *, double *);

static const R_CMethodDef CEntries[] = {
    {"critfunc", (DL_FUNC) &critfunc, 14},
    {"getcomp",  (DL_FUNC) &getcomp,   5},
    {"sample",   (DL_FUNC) &sample,   17},
    {NULL, NULL, 0}
};

void R_init_DoseFinding(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
