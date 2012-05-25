#include <stdlib.h>
#include <stdio.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>     /* for Calloc/Free */
#include <R_ext/Applic.h> /* for dgemm */

#define TRUE 1
#define FALSE 0

static inline double soft_threshold(double beta, double gamma)
{
   return sign(beta) * fmax(fabs(beta) - gamma, 0);
}

static inline double hard_threshold(double beta, double gamma)
{
   if(beta > gamma)
      return(gamma);
   else if(beta < -gamma)
      return(-gamma);
   return(beta);
}

/* From R src/main/array.c */
static void tcrossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "T";
    double one = 1.0, zero = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncx, &one,
			x, &nrx, y, &nry, &zero, z, &nrx);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < nrx*nry; i++) z[i] = 0;
    }
}

/* From R src/main/array.c */
static void matprod(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0;
    long double sum;
    Rboolean have_na = FALSE;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	/* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
	 * The test is only O(n) here
	 */
	for (i = 0; i < nrx*ncx; i++)
	    if (isnan(x[i])) {have_na = TRUE; break;}
	if (!have_na)
	    for (i = 0; i < nry*ncy; i++)
		if (isnan(y[i])) {have_na = TRUE; break;}
	if (have_na) {
	    for (i = 0; i < nrx; i++)
		for (k = 0; k < ncy; k++) {
		    sum = 0.0;
		    for (j = 0; j < ncx; j++)
			sum += x[i + j * nrx] * y[j + k * nry];
		    z[i + k * nrx] = sum;
		}
	} else
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}

/* From R src/main/array.c */
static void crossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z)
{
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
			x, &nrx, y, &nry, &zero, z, &ncx);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < ncx*ncy; i++) z[i] = 0;
    }
}

void spg_core(double *xx, double *xy, double *x, double *y,
   int *N_p, int *p_p, int *K_p, double *beta, double *C, int *CE_p,
   double *L_p, double *gamma_p, double *lambda_p,
   double *tol_p, double *mu_p, int *maxiter_p, int *verbose_p,
   int *niter, int *status)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
   int CE = *CE_p;
   int verbose = *verbose_p;
   double L = *L_p,
	  gamma = *gamma_p,
	  lambda = *lambda_p,
	  tol = *tol_p,
	  mu = *mu_p;
   int maxiter = *maxiter_p;
   double oneOnL = 1.0 / L;
   double theta = 1, theta_new, s1, s2, s3;

   double *A = malloc(sizeof(double) * CE * p);
   double *M = malloc(sizeof(double) * CE * p);
   double *W = calloc(p * K, sizeof(double));
   double *Wmu = calloc(p * K, sizeof(double));
   double *tmpNK= malloc(sizeof(double) * N * K);
   double *tmppK = malloc(sizeof(double) * p * K);
   double *tmppK2 = malloc(sizeof(double) * p * K);
   double *tmpCEp = malloc(sizeof(double) * CE * p);
   double *grad = malloc(sizeof(double) * p * K);
   double *XWy = malloc(sizeof(double) * N * p);
   double *beta_new = malloc(sizeof(double) * p * K);
   double *obj = calloc(maxiter, sizeof(double));
   int iter, mod = 0;

   if(verbose)
      Rprintf("spg_core_new\n");

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      for(i = p * K - 1 ; i >= 0 ; --i)
         Wmu[i] = W[i] / mu;

      tcrossprod(C, CE, K, Wmu, p, K, M); 
      for(i = CE * p - 1 ; i >= 0 ; --i)
         A[i] = hard_threshold(M[i], 1);

      if(p < 2 * N && p < 1e4)
      {
         crossprod(xx, p, p, W, p, K, tmppK);
         crossprod(A, CE, p, C, CE, K, tmppK2);
         for(i = p * K - 1 ; i >= 0 ; --i)
            grad[i] = tmppK[i] - xy[i] + tmppK2[i];
      }
      else
      {
	 matprod(x, N, p, W, p, K, tmpNK);
	 for(i = N * K - 1 ; i >= 0 ; --i)
	    tmpNK[i] -= y[i];

	 crossprod(x, N, p, tmpNK, N, K, tmppK); 
	 crossprod(A, CE, p, C, CE, K, tmppK2);

	 for(i = p * K - 1 ; i >= 0 ; --i)
	    grad[i] = tmppK[i] + tmppK2[i];
      }

      for(i = p * K - 1 ; i >= 0 ; --i)
	 beta_new[i] = soft_threshold(W[i] - oneOnL * grad[i],
	       lambda * oneOnL);

      theta_new = (sqrt(pow(theta, 4) + 4 * theta * theta) 
	    - theta * theta) * 0.5;

      for(i = p * K - 1 ; i >= 0 ; --i)
	 W[i] = beta_new[i] + (1 - theta) / theta 
	       * theta_new * (beta_new[i] - beta[i]);

      matprod(x, N, p, beta_new, p, K, tmpNK);
      tcrossprod(C, CE, K, beta_new, p, K, tmpCEp);

      s1 = 0;
      for(i = N * K - 1 ; i >= 0 ; --i)
	 s1 += pow(y[i] - tmpNK[i], 2);
      s1 *= 0.5;

      s2 = 0;
      for(i = CE * p - 1 ; i >= 0 ; --i)
	 s2 += fabs(tmpCEp[i]);

      s3 = 0;
      for(i = K * p - 1 ; i >= 0 ; --i)
	 s3 += fabs(beta_new[i]);
      s3 *= lambda;

      obj[iter] = s1 + s2 + s3;

      if(isinf(obj[iter]) || isnan(obj[iter]))
      {
	 Rprintf("SPG diverged; terminating (%.5f %.5f %.5f)\n", s1, s2, s3);
	 break;
      }

      if(verbose)
	 Rprintf("SPG iter=%d loss %.10f\n", iter, obj[iter]);

      if(iter > 10)
      {
	 if(fabs(obj[iter] - obj[iter-1]) / fabs(obj[iter-1]) < tol)
	 {
	    iter++;
	    break;
	 }
      }

      theta = theta_new;
      for(i = p * K - 1 ; i >= 0 ; --i)
	 beta[i] = beta_new[i];
   }
   *niter = iter;

   if(iter >= maxiter)
   {
      Rprintf("SPG failed to converge after %d iterations (lambda: %.6f, \
gamma: %.6f)\n", lambda, gamma, maxiter);
      *status = FALSE;
   }
   else 
   {
      *status = TRUE;
      if(verbose)
	 Rprintf("SPG converged after %d iterations\n", iter);
   }

   free(A);
   free(M);
   free(W);
   free(Wmu);
   free(tmpNK);
   free(tmppK);
   free(tmppK2);
   free(beta_new);
   free(tmpCEp);
   free(obj);
}

