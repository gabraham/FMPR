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

void tcrossprod(double *x, int nrx, int ncx,
   double *y, int nry, int ncy, double *z)
{
   char *transa = "N",
        *transb = "T";
   double one = 1.0,
          zero = 0.0;
   F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncx, &one,
      x, &nrx, y, &nry, &zero, z, &nrx);
}

void matprod(double *x, int nrx, int ncx,
   double *y, int nry, int ncy, double *z)
{
   char *transa = "N",
        *transb = "N";
   double one = 1.0,
          zero = 0.0;
   F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
      x, &nrx, y, &nry, &zero, z, &nrx);
}

void crossprod(double *x, int nrx, int ncx,
   double *y, int nry, int ncy, double *z)
{
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
      x, &nrx, y, &nry, &zero, z, &ncx);
}

void spg_core(double *xx, double *xy, double *x, double *y,
   int *N_p, int *p_p, int *K_p, double *beta, double *C, int *nE_p,
   double *L_p, double *gamma_p, double *lambda_p,
   double *tol_p, double *mu_p, int *maxiter_p, int *verbose_p,
   int *niter, int *status)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i;
   int nE = *nE_p;
   int verbose = *verbose_p;
   double L = *L_p,
	  gamma = *gamma_p,
	  lambda = *lambda_p,
	  tol = *tol_p,
	  mu = *mu_p;
   int maxiter = *maxiter_p;
   double oneOnL = 1.0 / L;
   double theta = 1, theta_new, s1, s2, s3;

   double *A = malloc(sizeof(double) * nE * p);
   double *M = malloc(sizeof(double) * nE * p);
   double *W = calloc(p * K, sizeof(double));
   double *Wmu = calloc(p * K, sizeof(double));
   double *tmpNK= malloc(sizeof(double) * N * K);
   double *tmppK = malloc(sizeof(double) * p * K);
   double *tmppK2 = malloc(sizeof(double) * p * K);
   double *tmpnEp = malloc(sizeof(double) * nE * p);
   double *grad = malloc(sizeof(double) * p * K);
   double *beta_new = malloc(sizeof(double) * p * K);
   double *obj = calloc(maxiter, sizeof(double));
   int iter;

   if(verbose > 1)
      Rprintf("spg_core\n");

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      for(i = p * K - 1 ; i >= 0 ; --i)
         Wmu[i] = W[i] / mu;

      tcrossprod(C, nE, K, Wmu, p, K, M); 
      for(i = nE * p - 1 ; i >= 0 ; --i)
         A[i] = hard_threshold(M[i], 1);

      if(p < 2 * N && p < 1e4)
      {
         matprod(xx, p, p, W, p, K, tmppK);
         crossprod(A, nE, p, C, nE, K, tmppK2);
         for(i = p * K - 1 ; i >= 0 ; --i)
            grad[i] = tmppK[i] - xy[i] + tmppK2[i];
      }
      else
      {
	 matprod(x, N, p, W, p, K, tmpNK);
	 for(i = N * K - 1 ; i >= 0 ; --i)
	    tmpNK[i] -= y[i];

	 crossprod(x, N, p, tmpNK, N, K, tmppK); 
	 crossprod(A, nE, p, C, nE, K, tmppK2);

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
      tcrossprod(C, nE, K, beta_new, p, K, tmpnEp);

      /* squared loss */
      s1 = 0;
      for(i = N * K - 1 ; i >= 0 ; --i)
	 s1 += pow(y[i] - tmpNK[i], 2);
      s1 *= 0.5;

      /* l1 fusion loss */
      s2 = 0;
      for(i = nE * p - 1 ; i >= 0 ; --i)
	 s2 += fabs(tmpnEp[i]);

      /* l1 lasso loss */
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

      if(verbose > 1)
	 Rprintf("SPG iter=%d loss %.10f (%.5f %.5f %.5f)\n", iter, obj[iter],
	       s1, s2, s3);

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
      if(verbose > 1)
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
   free(tmpnEp);
   free(obj);
}

/* L2 fusion penalty */
void spg_l2_core(double *xx, double *xy, double *x, double *y,
   int *N_p, int *p_p, int *K_p, double *beta, double *C, int *nE_p,
   double *L_p, double *gamma_p, double *lambda_p,
   double *tol_p, double *mu_p, int *maxiter_p, int *verbose_p,
   int *niter, int *status)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i;
   int nE = *nE_p;
   int verbose = *verbose_p;
   double L = *L_p,
	  gamma = *gamma_p,
	  lambda = *lambda_p,
	  tol = *tol_p;
   int iter, maxiter = *maxiter_p;
   double oneOnL = 1.0 / L;
   double theta = 1, theta_new, s1, s2, s3;
   int diverged = FALSE;

   double *W = calloc(p * K, sizeof(double));
   double *xW= malloc(sizeof(double) * N * K);
   double *xxW = malloc(sizeof(double) * p * K);
   double *WCC = malloc(sizeof(double) * p * K);
   double *CB = malloc(sizeof(double) * nE * p);
   double *grad = malloc(sizeof(double) * p * K);
   double *beta_new = malloc(sizeof(double) * p * K);
   double *obj = calloc(maxiter, sizeof(double));
   double *CC = calloc(K * K, sizeof(double));

   /* CC: K by K matrix */
   crossprod(C, nE, K, C, nE, K, CC);

   if(verbose > 1)
      Rprintf("spg_l2_core\n");

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      if(p < 2 * N && p < 1e4)
      {
         matprod(xx, p, p, W, p, K, xxW);
         matprod(W, p, K, CC, K, K, WCC);
         for(i = p * K - 1 ; i >= 0 ; --i)
            grad[i] = xxW[i] - xy[i] + gamma * WCC[i];
      }
      else
      {
	 matprod(x, N, p, W, p, K, xW);
	 for(i = N * K - 1 ; i >= 0 ; --i)
	    xW[i] -= y[i];

	 crossprod(x, N, p, xW, N, K, xxW); 
         matprod(W, p, K, CC, K, K, WCC);
	 for(i = p * K - 1 ; i >= 0 ; --i)
	    grad[i] = xxW[i] + gamma * WCC[i];
      }

      for(i = p * K - 1 ; i >= 0 ; --i)
	 beta_new[i] = soft_threshold(
	    W[i] - oneOnL * grad[i], lambda * oneOnL);

      theta_new = (sqrt(pow(theta, 4) + 4 * theta * theta) 
	    - theta * theta) * 0.5;

      for(i = p * K - 1 ; i >= 0 ; --i)
	 W[i] = beta_new[i] + (1 - theta) / theta 
	       * theta_new * (beta_new[i] - beta[i]);

      matprod(x, N, p, beta_new, p, K, xW);
      tcrossprod(C, nE, K, beta_new, p, K, CB);

      /* squared loss */
      s1 = 0;
      for(i = N * K - 1 ; i >= 0 ; --i)
	 s1 += pow(y[i] - xW[i], 2);
      s1 *= 0.5;

      /* l2 fussion penalty */
      s2 = 0;
      for(i = nE * p - 1 ; i >= 0 ; --i)
	 s2 += pow(CB[i], 2);
      s2 *= gamma;

      /* l1 lasso penalty */
      s3 = 0;
      for(i = K * p - 1 ; i >= 0 ; --i)
	 s3 += fabs(beta_new[i]);
      s3 *= lambda;

      obj[iter] = s1 + s2 + s3;

      if(isinf(obj[iter]) || isnan(obj[iter]) || obj[iter] > 10 * obj[0])
      {
	 Rprintf("SPG diverged; terminating (%.5f %.5f %.5f)\n", s1, s2, s3);
	 diverged = TRUE;
	 break;
      }

      if(verbose > 1)
	 Rprintf("SPG iter=%d loss %.10f (%.5f %.5f %.5f)\n", iter, obj[iter],
	       s1, s2, s3);

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

   if(diverged)
   {
      *status = FALSE;
   }
   else if(iter >= maxiter)
   {
      Rprintf("SPG failed to converge after %d iterations (lambda: %.6f, \
gamma: %.6f)\n", lambda, gamma, maxiter);
      *status = FALSE;
   }
   else 
   {
      *status = TRUE;
      if(verbose > 1)
	 Rprintf("SPG converged after %d iterations\n", iter);
   }

   free(W);
   free(xW);
   free(xxW);
   free(WCC);
   free(beta_new);
   free(CB);
   free(obj);
   free(CC);
}

