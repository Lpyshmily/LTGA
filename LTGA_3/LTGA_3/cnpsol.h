#ifndef CNPOSOL_H
#define CNPOSOL_H

#include"f2c.h"


#ifdef DLL_API
#else
#define DLL_API extern "C" _declspec(dllimport)
#endif

//test
DLL_API int daxpy_(integer* n, doublereal *alpha, doublereal *x,integer *incx, doublereal *y, integer *incy);

//lssol
DLL_API int lsopti_(char *string, integer *ivalue, ftnlen string_len);
DLL_API int lsoptn_(char *string, ftnlen string_len);
DLL_API int lsoptr_(char *string, doublereal *rvalue, ftnlen string_len);
DLL_API int lsfile_(integer *ioptns, integer *inform__);
DLL_API int lssol_(integer *mm, integer *n, integer *nclin, integer *
	lda, integer *ldr, doublereal *a, doublereal *bl, doublereal *bu, 
	doublereal *cvec, integer *istate, integer *kx, doublereal *x, 
	doublereal *r__, doublereal *b, integer *inform__, integer *iter, 
	doublereal *obj, doublereal *clamda, integer *iw, integer *leniw, 
	doublereal *w, integer *lenw);

//nlsol
DLL_API int nloptn_(char *string, ftnlen string_len);
DLL_API int nlopti_(char *string, integer *ivalue, ftnlen string_len);
DLL_API int nloptr_(char *string, doublereal *rvalue, ftnlen string_len);
DLL_API int nlfile_(integer *ioptns, integer *inform__);
DLL_API int nlssol_(integer *m, integer *n, integer *nclin, integer *
	ncnln, integer *lda, integer *ldcju, integer *ldfju, integer *ldr, 
	doublereal *a, doublereal *bl, doublereal *bu, U_fp funcon, U_fp 
	funobj, integer *inform__, integer *iter, integer *istate, doublereal 
	*c__, doublereal *cjacu, doublereal *y, doublereal *f, doublereal *
	fjacu, doublereal *clamda, doublereal *objf, doublereal *r__, 
	doublereal *x, integer *iw, integer *leniw, doublereal *w, integer *
	lenw);

//npsol
DLL_API int npoptn_(char *string, ftnlen string_len);
DLL_API int npopti_(char *string, integer *ivalue, ftnlen string_len);
DLL_API int npoptr_(char *string, doublereal *rvalue, ftnlen string_len);
DLL_API int npfile_(integer *ioptns, integer *inform__);
DLL_API int npsol_(integer *n, integer *nclin, integer *ncnln, 
	integer *lda, integer *ldju, integer *ldr, doublereal *a, doublereal *
	bl, doublereal *bu, U_fp funcon, U_fp funobj, integer *inform__, 
	integer *iter, integer *istate, doublereal *c__, doublereal *cjacu, 
	doublereal *clamda, doublereal *objf, doublereal *gradu, doublereal *
	r__, doublereal *x, integer *iw, integer *leniw, doublereal *w, 
	integer *lenw);

#endif