#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <R.h>

#define TRUE 1
#define FALSE 0

inline double sign(double x)
{
   if(x < 0) 
      return -1;
   else if(x > 0)
      return 1;
   else
      return 0;
}

inline double soft_threshold(double beta, double gamma)
{
   return sign(beta) * fmax(fabs(beta) - gamma, 0);
}

inline double hard_threshold(double beta, double gamma)
{
   if(beta > gamma)
      return(gamma);
   else if(beta < -gamma)
      return(-gamma);
   return(beta);
}

/*
 * (m by n)^T times (m by p) = (n by p)
 */
void crossprod(double *x, double *y, double *z, int *m_p, int *n_p, int *p_p)
{
   int i, j, k;
   int m = *m_p,
       n = *n_p,
       p = *p_p;

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 k = 0;
	 z[i + j * n] = x[k + m * i] * y[k + m * j];
	 for(k = 1 ; k < m ; k++)
	    z[i + j * n] += x[k + m * i] * y[k + m * j];
      }
   }
}

/*
 * (m by n) times (n by p) = (m by p)
 */
void prod(double *x, double *y, double *z, int *m_p, int *n_p, int *p_p)
{
   int i, j, k;
   int m = *m_p,
       n = *n_p,
       p = *p_p;

   for(i = 0 ; i < m ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 k = 0;
	 z[i + j * m] = x[k * m + i] * y[k + n * j];
	 for(k = 1 ; k < n ; k++)
	    z[i + j * m] += x[k * m + i] * y[k + n * j];
      }
   }
}

/*
 * (m by n) times (p by n)^T = (m by p)
 */
void tcrossprod(double *x, double *y, double *z, int *m_p, int *n_p, int *p_p)
{
   int i, j, k;
   int m = *m_p,
       n = *n_p,
       p = *p_p;

   for(i = 0 ; i < m ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 k = 0;
	 z[i + j * m] = x[k * m + i] * y[k * p + j];
	 for(k = 1 ; k < n ; k++)
	    z[i + j * m] += x[k * m + i] * y[k * p + j];
      }
   }
}

void spg_core(double *xx, double *xy, double *x, double *y,
   int *N_p, int *p_p, int *K_p, double *beta, double *C, int *CE_p,
   double *L_p, double *gamma_p, double *lambda_p,
   double *tol_p, double *mu_p, int *maxiter_p, int *verbose_p,
   int *niter)
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

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      for(i = p * K - 1 ; i >= 0 ; --i)
         Wmu[i] = W[i] / mu;

      tcrossprod(C, Wmu, M, &CE, &K, &p); 
      for(i = CE * p - 1 ; i >= 0 ; --i)
         A[i] = hard_threshold(M[i], 1);

      if(p < 2 * N && p < 1e4)
      {
         crossprod(xx, W, tmppK, &p, &p, &K);
         crossprod(A, C, tmppK2, &CE, &p, &K);
         for(i = p * K - 1 ; i >= 0 ; --i)
            grad[i] = tmppK[i] - xy[i] + tmppK2[i];
      }
      else
      {
	 prod(x, W, tmpNK, &N, &p, &K);
	 for(i = N * K - 1 ; i >= 0 ; --i)
	    tmpNK[i] -= y[i];

	 crossprod(x, tmpNK, tmppK, &N, &p, &K); 
	 crossprod(A, C, tmppK2, &CE, &p, &K);

	 for(i = p * K - 1 ; i >= 0 ; --i)
	    grad[i] = tmppK[i] + tmppK2[i];
      }

      for(i = p * K - 1 ; i >= 0 ; --i)
	 beta_new[i] = soft_threshold(W[i] - oneOnL * grad[i],
	       lambda * oneOnL);

      theta_new = (sqrt(pow(theta, 4) + 4 * theta * theta) 
	    - theta * theta) / 2.0;

      for(i = p * K - 1 ; i >= 0 ; --i)
	 W[i] = beta_new[i] + (1 - theta) / theta 
	       * theta_new * (beta_new[i] - beta[i]);

      prod(x, beta_new, tmpNK, &N, &p, &K);
      tcrossprod(C, beta_new, tmpCEp, &CE, &K, &p);

      s1 = 0;
      for(i = N * K - 1 ; i >= 0 ; --i)
	 s1 += pow(y[i] - tmpNK[i], 2);
      s1 /= 2.0;

      s2 = 0;
      for(i = CE * p - 1 ; i >= 0 ; --i)
	 s2 += fabs(tmpCEp[i]);

      s3 = 0;
      for(i = K * p - 1 ; i >= 0 ; --i)
	 s3 += fabs(beta_new[i]);
      s3 *= lambda;

      obj[iter] = s1 + s2 + s3;

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
   }
   else if(verbose)
   {
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

