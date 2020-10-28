/* Header file for cminpack, by Frederic Devernay.
   The documentation for all functions can be found in the file
   minpack-documentation.txt from the distribution, or in the source
   code of each function. */

#ifndef __CMINPACK_H__
#define __CMINPACK_H__

/* The default floating-point type is "double" for C/C++ and "float" for CUDA,
   but you can change this by defining one of the following symbols when
   compiling the library, and before including cminpack.h when using it:
   __cminpack_double__ for double
   __cminpack_float__ for float
   __cminpack_half__ for half from the OpenEXR library (in this case, you must
                     compile cminpack with a C++ compiler)
*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Cmake will define cminpack_EXPORTS on Windows when it
configures to build a shared library. If you are going to use
another build system on windows or create the visual studio
projects by hand you need to define cminpack_EXPORTS when
building a DLL on windows.
*/


/* Declarations for minpack */

/* Function types: */
/* The first argument can be used to store extra function parameters, thus */
/* avoiding the use of global variables. */
/* the iflag parameter is input-only (with respect to the FORTRAN */
/*  version), the output iflag value is the return value of the function. */
/* If iflag=0, the function shoulkd just print the current values (see */
/* the nprint parameters below). */
  
/* for hybrd1 and hybrd: */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/* return a negative value to terminate hybrd1/hybrd */
int fcnnn(int n, const double *x, double *fvec, int iflag, const double* para);

/* for hybrj1 and hybrj */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate hybrj1/hybrj */
int fcndernn(int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, const double* para);

/* for lmdif1 and lmdif */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/*         if iflag = 1 the result is used to compute the residuals. */
/*         if iflag = 2 the result is used to compute the Jacobian by finite differences. */
/*         Jacobian computation requires exactly n function calls with iflag = 2. */
/* return a negative value to terminate lmdif1/lmdif */
int fcnmn(int m, int n, const double *x, double *fvec, int iflag, const double* para);

/* for lmder1 and lmder */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate lmder1/lmder */
int fcndermn(int m, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, const double* para);

/* for lmstr1 and lmstr */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. */
/*         if iflag = i calculate the (i-1)-st row of the */
/*         jacobian at x and return this vector in fjrow. */
/* return a negative value to terminate lmstr1/lmstr */
int fcnderstrmn(int m, int n, const double *x, double *fvec, double *fjrow, int iflag, const double* para);






/* MINPACK functions: */
/* the info parameter was removed from most functions: the return */
/* value of the function is used instead. */
/* The argument 'p' can be used to store extra function parameters, thus */
/* avoiding the use of global variables. You can also think of it as a */
/* 'this' pointer a la C++. */

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation) */

int hybrd1( int(*fcnnn)(int n, const double *x, double *fvec, int iflag, const double* para), int n, double *x, double *fvec,
		   const double* para, double *wa, double xtol, int nprint, int maxfevno);

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation, more general). */

int hybrd(int(* fcnnn)(int n, const double *x, double *fvec, int iflag, const double* para), int n, double *x, double *fvec, 
		  const double* para, double xtol, int maxfev, int ml, int mu, double epsfcn, double *diag, int mode,
	      double factor, int nprint, int *nfev,
	      double *fjac, int ldfjac, double *r, int lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4);
  
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian) */

int hybrj1( int(* fcndernn)(int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, const double* para),
		   int n, double *x, double *fvec, const double* para, double *fjac, double *wa, double xtol, int nprint, int maxfevno);
          
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian,
   more general) */

int hybrj( int(* fcndernn)(int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, const double* para),
		  int n, double *x, double *fvec, const double* para, double *fjac, int ldfjac, double xtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *nfev, int *njev, double *r,
	      int lr, double *qtf, double *wa1, double *wa2,
	      double *wa3, double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation) */

int lmdif1( int(* fcnmn)(int m, int n, const double *x, double *fvec, int iflag, const double* para),
	       int m, int n, double *x, double *fvec, const double* para, int *iwa, double *wa, double xtol, int nprint, int maxfevno);

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation, more
   general) */

int lmdif(int(* fcnmn)(int m, int n, const double *x, double *fvec, int iflag, const double* para),
	      int m, int n, double *x, double *fvec, const double* para, double ftol,
	      double xtol, double gtol, int maxfev, double epsfcn,
	      double *diag, int mode, double factor, int nprint,
	      int *nfev, double *fjac, int ldfjac, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian) */

int lmder1(int (* fcndermn)(int m, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, const double* para),
	       int m, int n, double *x, double *fvec, const double* para, double *fjac, int *ipvt,
	       double *wa, double xtol, int nprint, int maxfevno);

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, more general) */

int lmder(int (* fcndermn)(int m, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, const double* para),
	      int m, int n, double *x, double *fvec, const double* para, double *fjac,
	      int ldfjac, double ftol, double xtol, double gtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage) */

int lmstr1( int (* fcnderstrmn)(int m, int n, const double *x, double *fvec, double *fjrow, int iflag, const double* para),
		   int m, int n, double *x, double *fvec, const double* para, double *fjac, int *ipvt, double *wa, double xtol, int nprint, int maxfevno);

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage, more general) */

int lmstr( int (* fcnderstrmn)(int m, int n, const double *x, double *fvec, double *fjrow, int iflag, const double* para),
		  int m, int n, double *x, double *fvec, const double* para, double *fjac,
	      int ldfjac, double ftol, double xtol, double gtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );
 

void chkder( int m, int n, const double *x, double *fvec, double *fjac,
	       int ldfjac, double *xp, double *fvecp, int mode,
	       double *err  );


double dpmpar( int i );


double enorm( int n, const double *x );

/* compute a forward-difference approximation to the m by n jacobian
   matrix associated with a specified problem of m functions in n
   variables. */

int fdjac2( int(* fcnmn)(int m, int n, const double *x, double *fvec, int iflag, const double* para),
	     int m, int n, double *x, const double *fvec, const double* para, double *fjac,
	     int ldfjac, double epsfcn, double *wa);

/* compute a forward-difference approximation to the n by n jacobian
   matrix associated with a specified problem of n functions in n
   variables. if the jacobian has a banded form, then function
   evaluations are saved by only approximating the nonzero terms. */

int fdjac1(int(* fcnnn)(int n, const double *x, double *fvec, int iflag, const double* para),
		   int n, double *x, const double *fvec, const double* para, double *fjac, int ldfjac,
	     int ml, int mu, double epsfcn, double *wa1,
	     double *wa2);

/* compute inverse(JtJ) after a run of lmdif or lmder. The covariance matrix is obtained
   by scaling the result by enorm(y)**2/(m-n). If JtJ is singular and k = rank(J), the
   pseudo-inverse is computed, and the result has to be scaled by enorm(y)**2/(m-k). */

void covar(int n, double *r, int ldr, const int *ipvt, double tol, double *wa);

/* covar1 estimates the variance-covariance matrix:
   C = sigma**2 (JtJ)**+
   where (JtJ)**+ is the inverse of JtJ or the pseudo-inverse of JtJ (in case J does not have full rank),
   and sigma**2 = fsumsq / (m - k)
   where fsumsq is the residual sum of squares and k is the rank of J.
   The function returns 0 if J has full rank, else the rank of J.
*/

int covar1(int m, int n, double fsumsq, double *r, int ldr, const int *ipvt, double tol, double *wa);

/* internal MINPACK subroutines */

void dogleg(int n, const double *r, int lr, 
             const double *diag, const double *qtb, double delta, double *x, 
             double *wa1, double *wa2);

void qrfac(int m, int n, double *a, int
            lda, int pivot, int *ipvt, int lipvt, double *rdiag,
            double *acnorm, double *wa);

void qrsolv(int n, double *r, int ldr, 
             const int *ipvt, const double *diag, const double *qtb, double *x, 
             double *sdiag, double *wa);

void qform(int m, int n, double *q, int
            ldq, double *wa);

void r1updt(int m, int n, double *s, int
             ls, const double *u, double *v, double *w, int *sing);

void r1mpyq(int m, int n, double *a, int
             lda, const double *v, const double *w);

void lmpar(int n, double *r, int ldr, 
            const int *ipvt, const double *diag, const double *qtb, double delta, 
            double *par, double *x, double *sdiag, double *wa1, 
            double *wa2);

void rwupdt(int n, double *r, int ldr, 
             const double *w, double *b, double *alpha, double *cos, 
             double *sin);
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __CMINPACK_H__ */
