#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void critfunc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getcomp(void *, void *, void *, void *, void *);
extern void sample(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

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
