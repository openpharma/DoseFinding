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


/*
To do:
- rwm - re-introduce random walk metropolis (for efficiency in linear models, 
        optimal tuning parameter can be pre-calculated!, maybe also for linear
        parameters of nonlinear models, or leave it optional)
- dbeta seems to be really slow in R (pre-calculate normalizing constant?)
 */

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include<float.h>
#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<R_ext/Lapack.h>
#include<R_ext/Applic.h>

/* structure to store basic information on problem */
struct modpars{
  double *doses;
  int *modelId;
  int *nPar;
  double *work;
  double *drEst;
  double *clinvCov;
  int *dim;
  double *prior;
  int *prnr;
  int *noint;
};

void R_CheckUserInterrupt(void);

#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<R_ext/Lapack.h>
#include<R_ext/BLAS.h>
#include<R_ext/Applic.h>


/* calculates A*x with A upper triangular */
void trmatvec(double *A, int *dim, double *x){
  char *uplo="U", *trans="N", *diag="N";
  int incx=1;
  F77_CALL(dtrmv)(uplo, trans, diag, dim, A, dim,
		  x, &incx FCONE FCONE FCONE);
}

/* calculates A*x for general A */
void matvec(double *A, int *nrow, int *ncol, 
	    double *x, double *y){
  char *trans="N";
  double alpha = 1.0, beta = 0.0;
  int incx=1;
  F77_CALL(dgemv)(trans, nrow, ncol, &alpha, A, nrow, 
		  x, &incx, &beta, y, &incx FCONE);
}

void crsprod(double *A, double *B, int *nrow, int *ncol, double *C){
  /* calculate A'B */
  /* Nrow - nrows of A*/
  /* ncol - ncols of A*/
  char *transa="T", *transb="N";
  double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(transa, transb, ncol, nrow, nrow,
  		  &alpha, A, nrow, B, nrow, &beta,
  		  C, ncol FCONE FCONE);
}

/* model functions */
void linear(double *doses, const int dim,
	    const double e0, const double delta,
	    double *resp){
  int i;
  for(i=0;i<dim;i++){
    resp[i] = e0 + delta*doses[i];
  }
}

void linlog(double *doses, const int dim,
	    const double e0, const double delta, const double off,
	    double *resp){
  int i;
  for(i=0;i<dim;i++){
    resp[i] = e0 + delta*log(doses[i]+off);
  }
}

void quadratic(double *doses, const int dim,
	       const double e0, const double beta1, const double beta2,
	       double *resp){
  int i;
  for(i=0;i<dim;i++){
    resp[i] = e0 + beta1*doses[i]+beta2*doses[i]*doses[i];
  }
}

void logistic(double *doses, const int dim,
	      const double e0, const double Emax, const double ED50, 
	      const double delta, double *resp){
  int i;
  double tmp1=0.0,tmp2=0.0;
  tmp1 = 1/delta;
  for(i=0;i<dim;i++){
    tmp2 = (ED50-doses[i])*tmp1;
    tmp2 = 1/(1+exp(tmp2));
    resp[i] = e0 + Emax*tmp2;
  }
}

void sigEmax(double *doses, const int dim,
	     const double e0, const double Emax, const double ED50, 
	     const double h, double *resp){
  int i;
  double tmp1=0.0,tmp2=0.0;
  tmp1 = pow(ED50, h);
  for(i=0;i<dim;i++){
    tmp2 = pow(doses[i], h);
    tmp2 = tmp2/(tmp1 + tmp2);
    resp[i] = e0 + Emax*tmp2;
  }
}

void emax(double *doses, const int dim,
	  const double e0, const double Emax, const double ED50, 
	  double *resp){
  sigEmax(doses, dim, e0, Emax, ED50, 1.0, resp);
}

void exponential(double *doses, const int dim,
		 const double e0, const double e1, const double delta,
		 double *resp){
  int i;
  double tmp=0.0;
  tmp = 1/delta;
  for(i=0;i<dim;i++){
    resp[i] = e0 + e1*(exp(doses[i]*tmp)-1);
  }
}

void betaMod(double *doses, const int dim,
	     const double e0, const double Emax, const double delta1, 
	     const double delta2, const double scal,
	     double *resp){
  int i;
  double tmp=0.0, B=0.0;
  B = delta1+delta2;
  B = pow(B, B)/(pow(delta1, delta1)*pow(delta2, delta2));
  for(i=0;i<dim;i++){
    tmp = doses[i]/scal;
    tmp = B*pow(tmp, delta1)*pow(1-tmp, delta2);
    resp[i] = e0 + Emax*tmp;
  }
}

void linInt(double *doses, const int dim,
	    double *par, double *resp){
  int i;
  for(i=0;i<dim;i++)
    resp[i] = par[i];
}


void getResp(double *par, double *doses, int *modelId, double *work,
	     int *dim){
  /* calculate 'approximate' likelihood */

  switch(*modelId){
  case 1:
    linear(doses, *dim, par[0], par[1], work);
    break;
  case 2:
    linlog(doses, *dim, par[0], par[1], par[2], work);
    break;
  case 3:
    quadratic(doses, *dim, par[0], par[1], par[2], 
	      work);
    break;
  case 4:
    linInt(doses, *dim, par, work);
    break;
  case 5:
    emax(doses, *dim, par[0], par[1], par[2], work);
    break;
  case 6:
    logistic(doses, *dim, par[0], par[1], par[2], 
	     par[3], work);
    break;
  case 7:
    exponential(doses, *dim, par[0], par[1], par[2], 
		work);
    break;
  case 8:
    sigEmax(doses, *dim, par[0], par[1], par[2], 
	    par[3], work);
    break;
  case 9:
    betaMod(doses, *dim, par[0], par[1], par[2],
	    par[3], par[4], work);
    break;
    //default: /* needed to remove this due to CRAN policy */
    //printf("invalid model selected"); 
  }
}


void loglik(double *par, double *doses, int *modelId, double *work,
	    double *drEst, double *clinvCov, int *dim, double *out){
  /* calculate 'approximate' likelihood */
  int i;

  getResp(par, doses, modelId, work, dim);
  
  for(i=0;i<*dim;i++){
    work[i] -= drEst[i];
  }
  
  trmatvec(clinvCov, dim, work);
  
  *out = 0.0;
  for(i=0;i<*dim;i++){
    *out -= work[i]*work[i];
  }
  *out *= 0.5;
}

double lg2(double x){
  return (x > 0) ? (log(x)) : (0.0);
}


void logprior(double *par, int *npar, double *prior, int *prnr, 
	      int *noint, double *out){
  /*
    prnr - number for prior (1- normal. 2-t, 3-log-normal, 4-beta)
    prior - prior parameters
    noint - equals 1 if there is no intercept in the model
  */
  *out =0.0;
  int i,count=0,i2=0;
  double p1=0.0,p2=0.0,p3=0.0,p4=0.0;
  for(i=0;i<(*npar-*noint);i++){
    i2 = i+*noint;
    p1 = prior[count];p2 = prior[count+1];
    if(prnr[i] == 1){ // normal-distribution
      *out += dnorm(par[i2],p1,p2, 1);
      count += 2;
    }
    if(prnr[i] == 2){ // t-distribution
      p3 = prior[count+2];
      *out += dt((par[i2]-p1)/p2, p3, 1)-log(p2);
      count += 3;
    }
    if(prnr[i] == 3){ // log-normal-distribution
      *out += dlnorm(par[i2], p1, p2, 1);
      count += 2;
    }
    if(prnr[i] == 4){ // scaled-beta-distribution
      p3 = prior[count+2];p4 = prior[count+3];
      *out += dbeta((par[i2]-p1)/(p2-p1), p3, p4, 1)-log(p2-p1);
      count += 4;
    }
  }
}

/* function to evaluate the pseudo-log-likelihood and log-prior */
double logPost(double *par, struct modpars *mp){
  double out=0.0,out2=0.0;
  logprior(par, mp->nPar, mp->prior, mp->prnr, mp->noint, &out);
  
  if(isfinite(out)){ /* only evaluate likelihood if prior > 0 */
    loglik(par, mp->doses, mp->modelId, mp->work, mp->drEst, 
	   mp->clinvCov, mp->dim, &out2);
    out += out2;
  }
  return out;
}


double logPost1d(double *actpar, int *ind, double *par, struct modpars *mp){
  double out=0.0;
  par[*ind] = *actpar;
  out = logPost(par, mp);
  return out;
}

/* get parameter bounds for non-linear parameters (from info on prior density) */
void getBnds(int *npar, double *prior, int *prnr, 
	     double *lower, double *upper, int *noint){
  /*
    prnr - number for prior (1-normal, 2-t-distribution, 3-log-normal, 4-beta)
    prior - prior parameters
  */
  int i,count=0,i2=0;
  for(i=0;i<*npar-*noint;i++){
    i2 = i+*noint;
    lower[i2] = -DBL_MAX;upper[i2] = DBL_MAX;
    if(prnr[i] == 1) // normal-distribution
      count += 2;
    if(prnr[i] == 2)// t-distribution
      count += 3;
    if(prnr[i] == 3){ // log-normal-distribution
      lower[i2] = 0.0;
      count += 2;
    }
    if(prnr[i] == 4){ // scaled-beta-distribution
      lower[i2] = prior[count];upper[i2] = prior[count+1];
      count += 4;
    }
  }
}



/* slice sampler */
/* stepping out procedure */
void getIntStep(double *par, int *ind, double *L, double *R,
		const double z, const double w, const double lower,
		const double upper, struct modpars *mp){
  double r,temp;
  r = unif_rand();
  temp = par[*ind];
  *L = temp - r*w;
  if(*L < lower)
    *L = lower;
  *R = temp + (1-r)*w;
  if(*R > upper)
    *R = upper;
  while(logPost1d(L, ind, par, mp) > z){
    *L -= w;
    if(*L < lower){
      *L = lower;
      break;
    }
  }
  while(logPost1d(R, ind, par, mp) > z){
    *R += w;
    if(*R > upper){
      *R = upper;
      break;
    }
  }
  par[*ind] = temp;
}

void slice1step(double *par, int *ind,
		const double w, double *lpostx, 
		const double lower, const double upper,
		struct modpars *mp){
  /*
    x - current value of chain (will contain output)
    ind - current dimension
    lpostx - logposterior evaluated at x
    w - tuning parameter of slice sampler
  */
  double z,tmp,xOld,xNew,L,R;
  z = *lpostx - exp_rand();
  xOld = par[*ind];
  /* get enclosing interval */
  getIntStep(par, ind, &L, &R, z, w, 
	     lower, upper, mp);
  while(1){
    xNew = unif_rand()*(R-L) + L;
    tmp = logPost1d(&xNew, ind, par, mp);
    if(tmp >= z - DBL_EPSILON){
      break;
    } else { // shrink interval
      if(xNew > xOld){
	R = xNew;
      } else {
	L = xNew;
      }
    }
  }
  par[*ind] = xNew;
  *lpostx = tmp;
}

/* whole sampler */
void sample(int *nSim, int *thin,
	    double *out, double *par, int *noint,
	    const double *w, 
	    double *doses, int *modelId, int *nPar, double *work, 
	    double *drEst, double *clinvCov, int *dim, double *prior, 
	    int *prnr, double *lower, double *upper){

  int i=0,count=0,j=0,d=0,actSimI=0;
  double lds,actSimD=0;
  /* initialize structural information */
  struct modpars mp = {doses, modelId, nPar, work, drEst, 
		       clinvCov, dim, prior, prnr, noint};
  actSimD = ((double) *nSim) / ((double) *thin);
  actSimI = (int) actSimD;

  /* calculate lower and upper bounds for parameters */
  getBnds(nPar, prior, prnr, lower, upper, noint);

  /* initialize R random number generator */
  GetRNGstate();
  
  lds = logPost(par, &mp); /* starting likelihood value */
  /* actual MCMC loop */
  for(i=0;i< *nSim;i++){
    for(d=*noint;d < *nPar;d++){
      slice1step(par, &d, w[d], &lds, 
		 lower[d], upper[d], &mp);
    }
    /* store information when desired */
    if(!(i%(*thin))){
      for(j = 0; j < *nPar; j++){
	out[count+j*(actSimI)] = par[j];
      }
      count++;
    }
  }
  PutRNGstate();
}
