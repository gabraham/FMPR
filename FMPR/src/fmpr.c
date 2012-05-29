#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>

#define TRUE 1
#define FALSE 0
#define ZERO_THRESH 1e-15
#define ZERO_VAR_THRESH 1e-6

static inline double sign(double x)
{
   if(x < 0) 
      return -1;
   else if(x > 0)
      return 1;
   else
      return 0;
}

static inline double soft_threshold(double beta, double gamma)
{
   return sign(beta) * fmax(fabs(beta) - gamma, 0);
}

void maxlambda1(double *x, double *y, double *lambda,
      int *N_p, int *p_p, int *K_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k, m;
   double d1, d2, s, xij;
   double const eps = 1e-16;

   for(k = 0 ; k < K ; k++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 d1 = 0;
	 d2 = 0;
	 for(i = N - 1 ; i >= 0 ; --i)
	 {
	    xij = x[i + j * N];
	    d1 += -xij * y[i + k * N];
	    d2 += xij * xij;
	 }
	 s = 0;
	 if(d2 > ZERO_VAR_THRESH)
	    s = fabs(d1 / d2);
	 lambda[k] = (lambda[k] > s ? lambda[k] : s);
      }
      /* small fudge factor to ensure penalty is high enough
       * to truly make everything zero, otherwise numerical error might cause
       * trouble */
      lambda[k] += eps;
   }
}

void lasso(double *x, double *y, double *b,
      int *N_p, int *p_p, double *lambda_p,
      int *maxiter_p, double *eps_p, int *verbose_p)
{
   int N = *N_p,
       p = *p_p;
   int i, j ;
   int iter;
   double d1;
   double loss = INFINITY, oldloss = INFINITY, lossnull = 0, lossnullF;
   int *active = malloc(sizeof(int) * p);
   int *oldactive = malloc(sizeof(int) * p);
   int *ignore = calloc(p, sizeof(int));
   double *d2 = calloc(p, sizeof(double));
   double *LP = calloc(N, sizeof(double));
   double delta, bj;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive, allconverged = 1, numconverged = 0;
   int verbose = *verbose_p;
   double meany = 0;
   double lambda = *lambda_p;
   double s;
   double oneOnN = 1.0 / N;

   for(j = p - 1 ; j >= 0 ; --j)
   {
      for(i = N - 1 ; i >= 0 ; --i)
	 d2[j] +=  x[i + j * N] * x[i + j * N];
      ignore[j] = d2[j] <= ZERO_VAR_THRESH;
   }

   for(j = p - 1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = !ignore[j];

   /* null loss mean((y - mean(y))^2)*/
   for(i = N - 1 ; i >= 0 ; --i)
      meany += y[i];
   meany *= oneOnN;

   for(i = N - 1 ; i >= 0 ; --i)
      lossnull += pow(y[i] - meany, 2);
   lossnull *= oneOnN;
   lossnullF = lossnull * eps;
	 
   if(verbose)
      Rprintf("null loss: %.5f\n", lossnull);

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      numactive = 0;
      numconverged = 0;

      for(j = 0 ; j < p ; j++)
      {
	 oldloss = loss;

	 if(!active[j])
	 {
	    numconverged++;
	    continue;
	 }

	 d1 = 0;
	 for(i = N - 1 ; i >= 0 ; --i)
	    d1 += x[i + j * N] * (LP[i] - y[i]);

	 // different implementation of the soft-thresholding
	 bj = b[j];
	 s = bj - d1 / d2[j];
	 if(fabs(s) <= lambda)
	 {
	    b[j] = 0;
	    delta = -bj;
	 }
	 else
	 {
	    b[j] = s - lambda * sign(s);
	    if(fabs(b[j]) < ZERO_THRESH)
	       b[j] = 0;
	    delta = b[j] - bj;
	 }

	 loss = 0;
	 for(i = N - 1 ; i >= 0 ; --i)
	 {
	    LP[i] += x[i + j * N] * delta;
	    loss += pow(LP[i] - y[i], 2);
	 }
	 loss *= oneOnN;
	 numconverged += fabs(loss - oldloss) < lossnullF;

	 active[j] = b[j] != 0;
	 numactive += active[j];
      }

      if(verbose)
      {
	 Rprintf("%d iter loss %.10f\n", iter, loss);
	 Rprintf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == p)
      {
	 if(verbose)
	    Rprintf("all converged at iter %d\n", iter);
	 if(allconverged == 1)
	 {
	    for(j = p - 1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = !ignore[j];
	    }
	    allconverged = 2;
	 }
	 else
	 {
	    for(j = p - 1 ; j >= 0 ; --j)
	       if(active[j] != oldactive[j])
		  break;
	    if(j < 0)
	    {
	       if(verbose)
		  Rprintf("terminating at iteration %d with %d active vars\n",
		     iter, numactive);
	       break;
	    }

	    if(verbose)
	       Rprintf("active set has changed\n");

	    allconverged = 1;
	    for(j = p - 1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = !ignore[j];
	    }
	 }
      }
   }

   if(iter >= maxiter && verbose)
      Rprintf("failed to converge after %d iterations\n", maxiter);

   free(active);
   free(oldactive);
   free(d2);
   free(LP);
   free(ignore);
}

void fmpr_weighted_warm(double *x, double *y, double *b,
      double *LP, int *N_p, int *p_p, int *K_p,
      double *lambda_p, double *gamma_p,
      double *G, int *maxiter_p, double *eps_p, int *verbose_p,
      int *status, int *iter_p, int *numactive_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
   int iter;
   double d1, d2;
   double delta, bjk;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive,
       allconverged = 1,
       numconverged = 0;
   int verbose = *verbose_p;
   int pK = p * K, pK1 = p * K - 1;
   double s, lambda = *lambda_p, gamma = *gamma_p;
   int g, q, iNk, kqK, jpk;

   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   int *ignore = calloc(p * K, sizeof(int));
   double *Err = calloc(N * K, sizeof(double));
   double *meany = calloc(K, sizeof(double));
   double *loss = calloc(K, sizeof(double));
   double *oldloss = calloc(K, sizeof(double));
   double *lossnull = calloc(K, sizeof(double));
   double *lossnullF = calloc(K, sizeof(double));
   double *d2_0 = calloc(p * K, sizeof(double));
   double *gammaG = calloc(K * K, sizeof(double));
   double *signG = calloc(K * K, sizeof(double));
   double losstotal = 0, oldlosstotal = 0;
   //double *oneOnLambda2PlusOne = malloc(sizeof(double) * K);
   double oneOnN = 1.0 / N,
	  oneOnKp = 1.0 / (K * p);

   /* null loss mean((y - mean(y))^2)*/
   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
         meany[k] += y[i];
      meany[k] *= oneOnN;
   }

   for(k = K * K - 1 ; k >= 0 ; --k)
      signG[k] = (double)sign(G[k]);

   /* setup second derivatives of loss and which variables to ignore due to
    * small variance */
   for(k = 0 ; k < K ; k++)
   {
      for(j = p - 1 ; j >= 0 ; --j) 
      {
	 jpk = j + p * k;
	 for(i = N - 1 ; i >= 0 ; --i)
	    d2_0[jpk] +=  x[i + j * N] * x[i + j * N];
	 ignore[jpk] = d2_0[jpk] <= ZERO_VAR_THRESH;;
      }
   }

   for(j = pK1 ; j >= 0 ; --j)
      oldactive[j] = active[j] = !ignore[j];

   /* ensure that fusion weights for tasks with themselves are zero */
   for(k = 0 ; k < K ; k++)
   {
      for(q = 0 ; q < K ; q++)
      {
	 kqK = k + q * K;
	 /* must take absolute value of G since f(r_ml) is monotonic in
	  * abs(r_ml) but G is signed
	  */
	 if(q != k)
	    gammaG[kqK] = gamma * fabs(G[kqK]);
      }
   }

   for(k = 0 ; k < K ; k++)
   {
      loss[k] = 0;
      for(i = N - 1 ; i >= 0 ; --i)
      {
	 iNk = i + N * k;
	 lossnull[k] += pow(y[iNk] - meany[k], 2);
	 Err[iNk] = LP[iNk] - y[iNk];
	 loss[k] += Err[iNk] * Err[iNk];
      }
      loss[k] *= oneOnN;
      lossnull[k] *= oneOnN;
      lossnullF[k] = lossnull[k] * eps;
   }
   
   for(iter = 0 ; iter < maxiter ; iter++)
   {
      numactive = 0;
      numconverged = 0;
      losstotal = 0;
      oldlosstotal = 0;

      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p ; j++)
	 {
	    jpk = j + p * k;

	    if(!active[jpk])
	       numconverged++;
	    else
	    {
	       d1 = 0;
	       d2 = d2_0[jpk];

	       for(i = N - 1 ; i >= 0 ; --i)
	          d1 += x[i + j * N] * Err[i + N * k];

	       // different implementation of the soft-thresholding
	       bjk = b[jpk];
	       
	       /* Apply inter-task penalty */
	       for(q = 0 ; q < K ; q++)
	       {
	          kqK = k + q * K;
	          d1 += gammaG[kqK] * (bjk - signG[kqK] * b[j + p * q]);
	          d2 += gammaG[kqK];
	       }

	       /* Apply intra-task ridge regression */
	       s = (bjk - d1 / d2); //* oneOnLambda2PlusOne[k];

	       /* Now apply intra-task lasso */
	       if(fabs(s) <= lambda)
	       {
	          b[jpk] = 0;
	          delta = -bjk;
	       }
	       else
	       {
	          b[jpk] = s - lambda * sign(s);
		  if(fabs(b[jpk]) < ZERO_THRESH) // close enough to zero
		     b[jpk] = 0;
	          delta = b[jpk] - bjk;
	       }

	       /* update loss and errors based on new estimates */
	       oldloss[k] = loss[k];
	       loss[k] = 0;
	       for(i = N - 1 ; i >= 0 ; --i)
	       {
		  iNk = i + N * k;
	          LP[iNk] += x[i + N * j] * delta;
		  Err[iNk] = LP[iNk] - y[iNk];
	          loss[k] += Err[iNk] * Err[iNk];
	       }
	       loss[k] *= oneOnN;
	       //loss[k] += lambda * fabs(b[jpk]) * oneOnKp;
	       numconverged += fabs(loss[k] - oldloss[k]) < lossnullF[k];
	    }

	    active[jpk] = b[jpk] != 0;
	    numactive += active[jpk];
	 }

	 losstotal += loss[k];
	 oldlosstotal += oldloss[k];
      }

      if(verbose)
      {
	 Rprintf("%d iter loss %.10f\n", iter, losstotal);
	 Rprintf("%d converged at iter %d\n", numconverged, iter);
      }

      /* active-set convergence */
      if(numconverged == pK)
      {
         if(verbose)
            Rprintf("all (%d) converged at iter %d\n", numconverged, iter);
         if(allconverged == 1)
         {
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = !ignore[j];
            }
            allconverged = 2;
         }
         else
         {
            for(j = pK1 ; j >= 0 ; --j)
               if(active[j] != oldactive[j])
        	  break;
            if(j < 0)
            {
               if(verbose)
        	  Rprintf("terminating at iteration %d with %d active vars\n",
        	     iter, numactive);
	       *status = TRUE;
               break;
            }

	    if(verbose)
	       Rprintf("activeset changed at iter %d\n", iter);

            allconverged = 1;
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = !ignore[j];
            }
         }
      }
   }

   if(iter >= maxiter)
   {
      if(verbose)
	 Rprintf("fmpr_threshold_warm failed to converge after %d iterations\n",
	    maxiter);
      *status = FALSE;
   }
   *iter_p = iter;
   *numactive_p = numactive;

   free(active);
   free(oldactive);
   free(meany);
   free(loss);
   free(oldloss);
   free(lossnull);
   free(lossnullF);
   free(Err);
   free(d2_0);
   free(ignore);
   free(gammaG);
   free(signG);
}

/*
 * Assumes y is centred so there's no intercept
 * 
 * Weighted rather than thresholded correlation matrix
 *
 */
void fmpr_weighted(double *x, double *y, double *b,
      int *N_p, int *p_p, int *K_p,
      double *lambda_p, double *gamma_p,
      double *G, int *maxiter_p, double *eps_p, int *verbose_p,
      int *status, int *iter_p, int *numactive_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
   int iter;
   double d1, d2;
   double delta, bjk;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive,
       allconverged = 1,
       numconverged = 0;
   int verbose = *verbose_p;
   int pK = p * K, pK1 = p * K - 1;
   double s, lambda = *lambda_p, gamma = *gamma_p;
   int q, iNk, kqK, jpk;
   double g;

   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   int *ignore = calloc(p * K, sizeof(int));
   double *LP = calloc(N * K, sizeof(double));
   double *Err = calloc(N * K, sizeof(double));
   double *meany = calloc(K, sizeof(double));
   double *loss = malloc(sizeof(double) * K);
   double *oldloss = malloc(sizeof(double) * K);
   double *lossnull = calloc(K, sizeof(double));
   double *lossnullF = calloc(K, sizeof(double));
   double *d2_0 = calloc(p * K, sizeof(double));
   double losstotal = 0, oldlosstotal = 0;
   double *gammaG = calloc(K * K, sizeof(double));
   double *signG = calloc(K * K, sizeof(double));
   double oneOnN = 1.0 / N;

   /* null loss mean((y - mean(y))^2)*/
   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
         meany[k] += y[i];
      meany[k] *= oneOnN;
   }

   for(k = K * K - 1 ; k >= 0 ; --k)
      signG[k] = (double)sign(G[k]);

   for(k = 0 ; k < K ; k++)
   {
      for(j = p - 1 ; j >= 0 ; --j) 
      {
	 jpk = j + p * k;
	 for(i = N - 1 ; i >= 0 ; --i)
	    d2_0[jpk] +=  x[i + j * N] * x[i + j * N];
	 ignore[jpk] = d2_0[jpk] <= ZERO_VAR_THRESH;;
      }
   }

   for(j = pK1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = !ignore[j];

   for(k = 0 ; k < K ; k++)
   {
      loss[k] = 0;
      for(i = N - 1 ; i >= 0 ; --i)
      {
	 iNk = i + N * k;
	 lossnull[k] += pow(y[iNk] - meany[k], 2);
	 Err[iNk] = LP[iNk] - y[iNk];
	 loss[k] += Err[iNk] * Err[iNk];
      }
      loss[k] *= oneOnN;
      lossnull[k] *= oneOnN;
      lossnullF[k] = lossnull[k] * eps;
   }

   /* ensure that fusion weights for tasks with themselves are zero */
   for(k = 0 ; k < K ; k++)
   {
      for(q = 0 ; q < K ; q++)
      {
	 kqK = k + q * K;
         gammaG[kqK] = 0;
	 /* must take absolute value of G since f(r_ml) is monotonic in
	  * abs(r_ml) but G is signed
	  */
	 if(q != k)
	    gammaG[kqK] = gamma * fabs(G[kqK]);
      }
   }

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      numactive = 0;
      numconverged = 0;
      losstotal = 0;
      oldlosstotal = 0;

      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p ; j++)
	 {
	    jpk = j + p * k;

	    d1 = 0;
	    d2 = d2_0[jpk];

	    if(!active[jpk])
	       numconverged++;
	    else
	    {
	       for(i = N - 1 ; i >= 0 ; --i)
	          d1 += x[i + j * N] * Err[i + N * k];

	       bjk = b[jpk];

	       // Apply inter-task ridge regression
	       for(q = 0 ; q < K ; q++)
	       {
	          kqK = k + q * K;
	          d1 += gammaG[kqK] * (bjk - signG[kqK] * b[j + p * q]);
	          d2 += gammaG[kqK];
	       }

	       s = bjk - d1 / d2;

	       /* Now apply intra-task lasso */
	       if(fabs(s) <= lambda)
	       {
	          b[jpk] = 0;
	          delta = -bjk;
	       }
	       else
	       {
	          b[jpk] = s - lambda * sign(s);

		  // values close enough to zero are considered zero to avoid
		  // silly situations where tiny numbers produced from
		  // numerical error are considered truly non-zero
		  if(fabs(b[jpk]) < ZERO_THRESH)
		     b[jpk] = 0;
	          delta = b[jpk] - bjk;
	       }

	       oldloss[k] = loss[k];
	       loss[k] = 0;
	       for(i = N - 1 ; i >= 0 ; --i)
	       {
		  iNk = i + N * k;
	          LP[iNk] += x[i + N * j] * delta;
		  Err[iNk] = LP[iNk] - y[iNk];
	          loss[k] += Err[iNk] * Err[iNk];
	       }
	       loss[k] *= oneOnN;
	       // TODO: add the fusion loss to the total loss
	       //loss[k] += lambda * fabs(b[jpk]) / K;
	       numconverged += fabs(loss[k] - oldloss[k]) < lossnullF[k];
	    }

	    active[jpk] = b[jpk] != 0;
	    numactive += active[jpk];
	 }

	 losstotal += loss[k];
	 oldlosstotal += oldloss[k];
      }

      if(verbose)
      {
	 Rprintf("%d iter loss %.10f\n", iter, losstotal);
	 Rprintf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == pK)
      {
         if(verbose)
            Rprintf("all (%d) converged at iter %d\n", numconverged, iter);
         if(allconverged == 1)
         {
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = !ignore[j];
            }
            allconverged = 2;
         }
         else
         {
            for(j = pK1 ; j >= 0 ; --j)
               if(active[j] != oldactive[j])
        	  break;
            if(j < 0)
            {
               if(verbose)
        	  Rprintf("terminating at iteration %d with %d active vars\n",
        	     iter, numactive);
	       *status = TRUE;
               break;
            }

	    if(verbose)
	       Rprintf("activeset changed at iter %d\n", iter);

            allconverged = 1;
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = !ignore[j];
            }
         }
      }
   }

   if(iter >= maxiter)
   {
      if(verbose)
	 Rprintf("failed to converge after %d iterations\n", maxiter);
      *status = FALSE;
   }
   *iter_p = iter;
   *numactive_p = numactive;

   free(active);
   free(oldactive);
   free(LP);
   free(meany);
   free(loss);
   free(oldloss);
   free(lossnull);
   free(lossnullF);
   free(Err);
   free(d2_0);
   free(ignore);
   free(gammaG);
   free(signG);
}

