/*
#######################################################################
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.
*/

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif
#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<R_ext/Lapack.h>
#include<R_ext/Applic.h>


void rank1vec(double *grad, int *nPar, double *alpha, double *A){
  // calculates alpha*grad*grad'+A
  char uplo='U';
  int inc=1;
  F77_CALL(dsyr)(&uplo, nPar, alpha, grad, &inc,
		 A, nPar FCONE);
}

// calculate design matrix
void calcMat(double *grad, int *nPar, double *design, 
             int *nD, double *A, int *incr){
  // nD - number of doses (length of design)
  // nPar - number of parameters = ncol(A) = nrow(A)
  double gradsub[4]={0.0};
  int i,j=0;
  for(i=0;i<*nD;i++){
    for(j=0;j<*nPar;j++){
      gradsub[j] = grad[*incr+*nPar*i+j];
    }
    rank1vec(gradsub, nPar, &design[i], A);
  }
  // complete symmetric matrix from upper triang. part
  for(i=0;i<*nPar;i++){
    for(j=0;j<i;j++){
      A[*nPar*j+i] = A[*nPar*i+j];
    }
  }
}

void calcDetGinv(double *A, int *nPar, double *work, 
		 double *s, double *VT, double *U, 
		 double *tol, int *type, double *resD){
  // calculate svd decomposition and from there
  // the g-inverse and the determinant
  // s - needs to be of size nPar
  // U,VT - dimension nPar*nPar
  // work - min dimension: 5*nPar use 30 here
  // 
  int i,j,nD,nonzero=*nPar;
  char jobu='A';
  char jobvt='A';
  int info,lwork=30;
  
  // Calculate singular value decomposition
  F77_CALL(dgesvd)(&jobu, &jobvt, nPar, nPar, A, nPar,
		   s, U, nPar, VT, nPar, work, &lwork,
		   &info FCONE FCONE);
  
  if((*type == 1) || (*type == 3)){ // calculate g-inverse
    for (i = 1; i < *nPar; i++){
      if (s[i] < *tol*s[0]){
	nonzero = i;
	break;
      }
    }
    for (i = 0; i < *nPar; i++){
      for (j = 0; j < nonzero; j++){
	U[j**nPar + i] = U[j**nPar+i] * 1.0/s[j];
      }
    }
    for (i = 0; i < *nPar; i++){ // g-inverse stored in upper tri part of A
      for (j = i; j < *nPar; j++){
	A[j**nPar+i] = 0.0;
	for (nD=0; nD < nonzero; nD++){
	  A[j**nPar+i] += VT[i**nPar+nD] * U[nD**nPar+j];
	}
      }
    }
  }
  if((*type == 2) || (*type == 3)){ // calculate determinant
    *resD = 1.0;
    for (i = 0; i < *nPar; i++){
      *resD *= s[i];
    }
  }
}

void calcQuadform(double *beta, double *Q, int *nPar, double *out, int *incrbeta){
  // calculates quadratic form beta'Qbeta
  // Q = (Q11,Q12,...,Q22,Q23,... (only upper triangular part of sym. matrix)
  int i,j;
  for(i=0;i<*nPar;i++){
    for(j=i;j<*nPar;j++){
      if(i==j){
        *out += Q[*nPar*j+i]*beta[*incrbeta+i]*beta[*incrbeta+i];
      } else {
        *out += 2*Q[*nPar*j+i]*beta[*incrbeta+i]*beta[*incrbeta+j];
      }
    }
  }
}

void getAllocs(double *w2, double *n2, double *nold, int *nD){
  int i;
  double n1=0;
  for(i=0;i<*nD;i++){
    n1 += nold[i];
  }
  for(i=0;i<*nD;i++){
    w2[i] = (*n2*w2[i] + nold[i])/(*n2 + n1);
  }
}

void setzero(double *x, int nPar){
  int i;
  for(i=0;i<nPar;i++){
    x[i] = 0.0;
  }
}

void critfunc(double *grad, int *nPar, int *nD, double *probs, int *M,
              double *design, double *n2, double *nold,
              double *A, double *tol, double *MEDgrad, int *type,
	      int *stand, double *res){
  // grad - contains gradient vectors (4 cells reserved for each model)
  // nPar - number of parameters (dim A)
  // nD - number of dose-levels
  // design - design
  // type - 1: MED, 2: Dopt, 3: MED&Dopt
  int m,incgrad=0,incb=0;
  double resM=0,resD=0,fracp=0;
  // variables for SVD decomposition, initialize to max possible dimension
  double work[30]={0.0};
  double s[4]={0.0}; 
  double VT[16]={0.0};
  double U[16]={0.0};
  *res = 0.0;

  // calculate weight vector
  getAllocs(design, n2, nold, nD);
  for(m=0;m<*M;m++){
    if(m > 0){
      incgrad+=*nD*nPar[m-1];
      incb+=nPar[m-1];
    }
    setzero(A, 16);resM = 0.0;
    // calulate matrix 
    calcMat(grad, &nPar[m], design, nD, A, &incgrad);
    // calculate det and/or  MP-Inverse 
    calcDetGinv(A, &nPar[m], work, s, VT, U,
		tol, type, &resD);
    if(*type == 1){  // calculate quadratic form (for MED designs)
      calcQuadform(MEDgrad, A, &nPar[m], &resM, &incb);
      *res += probs[m]*log(resM);
    }
    if(*type == 2){
      if(*stand == 1){
	fracp = (double) nPar[m];
	*res += probs[m]*(-log(resD)/fracp);	
      } else {
	*res += probs[m]*(-log(resD));
      }
    }
    if(*type == 3){  // calculate quadratic form (for MED designs)
      calcQuadform(MEDgrad, A, &nPar[m], &resM, &incb);
      if(*stand == 1){
	fracp = (double) nPar[m];
	*res += probs[m]*(-0.5*log(resD)/fracp+0.5*log(resM));
      } else {
	*res += probs[m]*(-0.5*log(resD)+0.5*log(resM));
      }
    }
  }
}
