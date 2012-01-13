#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

void maxlambda1(double *x, double *y, double *lambda,
      int *N_p, int *p_p, int *K_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k, m;
   double d1, d2, s, xij;
   double eps = 1e-16;

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
	 s = fabs(d1 / d2);
	 lambda[k] = (lambda[k] > s ? lambda[k] : s);
      }
      /* small fudge factor to ensure penalty is high enough
       * to truly make everything zero */
      lambda[k] += eps;
   }
}

void groupridge_simple(double *x, double *y, double *B,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3,
      int *grp, int *maxiter_p, double *eps_p, int *verbose_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k, q;
   int iter;
   double d1, d2;
   double loss = INFINITY, oldloss = INFINITY;
   double *v = calloc(p, sizeof(double));
   double *LP = calloc(N * K, sizeof(double));
   double s, s2, delta, Bjk;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   double lossk;
   int kN, ikN, itercount = 0;
   double oneOnN = 1.0 / N;
   int numactive;
   int verbose = *verbose_p;

   /* second derivative */
   for(j = p - 1 ; j >= 0 ; --j)
      for(i = N - 1 ; i >= 0 ; --i)
	 v[j] +=  x[i + j * N] * x[i + j * N];

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      if(itercount == 2000)
      {
	 itercount = 0;
	 if(verbose)
	    printf("iter %d\n", iter);
      }
      itercount++;

      loss = 0;
      numactive = 0;
      for(k = 0 ; k < K ; k++)
      {
	 kN = k * N;
	 for(j = 0 ; j < p ; j++)
	 {
	    Bjk = B[j + k * p];
	    d1 = 0;
	    d2 = v[j];
	    for(i = N - 1 ; i >= 0 ; --i)
	       d1 += x[i + j * N] * (LP[i + kN] - y[i + kN]);

	    s = Bjk - d1 / d2;

	    /* soft thresholding, if beta is zero, skip the other penalties */
	    if(fabs(s) <= lambda1[k])
	    {
	       delta = -B[j + k * p];
	       B[j + k * p] = 0;
	       lossk = 0;
	       for(i = N - 1 ; i >= 0 ; --i)
	       {
		  ikN = i + k * N;
	          LP[ikN] += x[i + j * N] * delta;
		  lossk += pow(LP[ikN] - y[ikN], 2);
	       }
	       loss += lossk * oneOnN;
	       continue;
	    }

	    /* lasso penalty when beta!=0, no 2nd derivative */
	    /*d1 -= lambda1[k] * sign(Bjk);*/

	    /* ridge penalty intra-class */
	    /*d1 += lambda2[k] * Bjk;
	    d2 += lambda2[k];*/

	    /* ridge penalty inter-class */
	    /*for(q = 0 ; q < K ; q++)
	    {
	       if(grp[k] == grp[q] && k != q)
	       {
		  d1 += lambda3[k] * (Bjk - B[j + q * p]);
		  d2 += lambda3[k] * N;
	       }
	    }*/
	    
	    /*s = d1 / d2;
	    s2 = Bjk - s;
	    delta = s2 - Bjk;
	    B[j + k * p] = s2;*/

	    B[j + k * p] = s - lambda1[k] * sign(s);
	    delta = B[j + k * p] - Bjk;

	    lossk = 0;
	    for(i = N - 1 ; i >= 0 ; --i)
	    {
	       ikN = i + k * N;
	       LP[ikN] += x[i + j * N] * delta;
	       lossk += pow(LP[ikN] - y[ikN], 2);
	    }
	    loss += lossk * oneOnN;

	    numactive += B[j + k * p] != 0;
	 }
	 printf("iter %d loss %.7f and %d active vars\n",
	       iter, loss, numactive);
      }

      //if(fabs(oldloss - loss) < eps)
      if(fabs(oldloss - loss) / fabs(loss) < eps)
      {
	 if(verbose)
	 {
	    printf(
	       "terminating at iteration %d with loss %.7f and %d active vars\n",
		  iter, loss, numactive);
	    fflush(stdout);
	 }
	 break;
      }
      /*else if(iter > 1 && numactive >= N)
      {
	 printf(
	    "saturated model, terminating at iteration %d \
with loss %.7f and %d active vars\n",
	       iter, loss, numactive);
	 break;
      }*/
      oldloss = loss;
   }
   
   if(iter >= maxiter && verbose)
      printf("failed to converge after %d iterations\n", maxiter);

   free(v);
   free(LP);
}

/*
 * lambda1: lasso penalty within each task
 * lambda2: ridge penalty within each task
 * lambda3: ridge penalty across tasks
 */
void groupridge(double *x, double *y, double *B,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3,
      int *grp, int *maxiter_p, double *eps_p, int *verbose_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k, q;
   int iter;
   double d1, d2;
   double loss = INFINITY, oldloss = INFINITY;
   double *v = calloc(p, sizeof(double));
   double *LP = calloc(N * K, sizeof(double));
   double s, s2, delta, Bjk;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   double lossk;
   int ikN, itercount = 0;
   double oneOnN = 1.0 / N;
   int numactive;
   int verbose = *verbose_p;

   for(j = p - 1 ; j >= 0 ; --j)
      for(i = N - 1 ; i >= 0 ; --i)
	 v[j] +=  x[i + j * N] * x[i + j * N];

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      if(itercount == 2000)
      {
	 itercount = 0;
	 if(verbose)
	    printf("iter %d\n", iter);
      }
      itercount++;

      loss = 0;
      numactive = 0;
      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p ; j++)
	 {
	    d1 = 0;
	    d2 = v[j];
	    for(i = N - 1 ; i >= 0 ; --i)
	    {
	       q = k * N;
	       d1 += x[i + j * N] * (LP[i + q] - y[i + q]);
	    }

	    s = B[j + k * p] - d1 / d2;

	    /* beta is zero, skip the other penalties */
	    if(fabs(s) <= lambda1[k])
	    {
	       delta = - B[j + k * p];
	       B[j + k * p] = 0;
	       lossk = 0;
	       for(i = 0 ; i < N ; i++)
	       {
		  ikN = i + k * N;
	          LP[ikN] += x[i + j * N] * delta;
		  lossk += pow(LP[ikN] - y[ikN], 2);
	       }
	       loss += lossk * oneOnN;
	       continue;
	    }

	    Bjk = B[j + k * p];
	    /* lasso penalty when beta!=0, no 2nd derivative */
	    d1 += lambda1[k] * sign(Bjk);

	    /* ridge penalty intra-class */
	    d1 += lambda2[k] * Bjk;
	    d2 += lambda2[k];

	    /* ridge penalty inter-class */
	    for(q = 0 ; q < K ; q++)
	    {
	       if(grp[k] == grp[q] && k != q)
	       {
		  d1 += lambda3[k] * (Bjk - B[j + q * p]);
		  d2 += lambda3[k] * N;
	       }
	    }
	    
	    s = d1 / d2;
	    s2 = Bjk - s;
	    delta = s2 - Bjk;
	    B[j + k * p] = s2;
	    lossk = 0;
	    for(i = N - 1 ; i >= 0 ; --i)
	    {
	       ikN = i + k * N;
	       LP[ikN] += x[i + j * N] * delta;
	       lossk += pow(LP[ikN] - y[ikN], 2);
	    }
	    loss += lossk * oneOnN;

	    numactive += B[j + k * p] != 0;
	 }
      }

      if(fabs(oldloss - loss) < eps)
      {
	 if(verbose)
	    printf(
	       "terminating at iteration %d with loss %.7f and %d active vars\n",
		  iter, loss, numactive);
	 break;
      }
      /*else if(iter > 1 && numactive >= N)
      {
	 printf(
	    "saturated model, terminating at iteration %d \
with loss %.7f and %d active vars\n",
	       iter, loss, numactive);
	 break;
      }*/
      oldloss = loss;
   }
   
   if(iter >= maxiter && verbose)
      printf("failed to converge after %d iterations\n", maxiter);

   free(v);
   free(LP);
}

void lasso2(double *x, double *y, double *b,
      int *N_p, int *p_p, double *lambda1_p,
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
   double *d2 = calloc(p, sizeof(double));
   double *LP = calloc(N, sizeof(double));
   double delta, bj;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive, allconverged = 1, numconverged = 0;
   int verbose = *verbose_p;
   double meany = 0;
   double lambda1 = *lambda1_p;

   for(j = p - 1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   for(j = p - 1 ; j >= 0 ; --j)
   {
      for(i = N - 1 ; i >= 0 ; --i)
	 d2[j] +=  x[i + j * N] * x[i + j * N];
   }

   /* null loss mean((y - mean(y))^2)*/
   for(i = N - 1 ; i >= 0 ; --i)
      meany += y[i];
   meany /= N;

   for(i = N - 1 ; i >= 0 ; --i)
      lossnull += pow(y[i] - meany, 2);
   lossnull /= N;
   lossnullF = lossnull * 1e-8;
	 
   if(verbose)
      printf("null loss: %.5f\n", lossnull);

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

	 bj = b[j];
	 b[j] = soft_threshold(bj - d1 / d2[j], lambda1);
	 delta = b[j] - bj;

	 loss = 0;
	 for(i = N - 1 ; i >= 0 ; --i)
	 {
	    LP[i] += x[i + j * N] * delta;
	    loss += pow(LP[i] - y[i], 2);
	 }
	 loss /= N;
	 numconverged += fabs(loss - oldloss) < lossnullF;

	 active[j] = b[j] != 0;
	 numactive += active[j];
      }

      if(verbose)
      {
	 printf("%d iter loss %.10f\n", iter, loss);
	 printf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == p)
      {
	 if(verbose)
	    printf("all converged at iter %d\n", iter);
	 if(allconverged == 1)
	 {
	    for(j = p - 1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = TRUE;
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
		  printf("terminating at iteration %d with %d active vars\n",
		     iter, numactive);
	       break;
	    }

	    allconverged = 1;
	    for(j = p - 1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = TRUE;
	    }
	 }
      }
   }

   if(iter >= maxiter && verbose)
      printf("failed to converge after %d iterations\n", maxiter);

   free(active);
   free(oldactive);
   free(d2);
   free(LP);
}

void groupridge2(double *x, double *y, double *b,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3,
      int *grp, int *maxiter_p, double *eps_p, int *verbose_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
   int iter;
   double d1;
   double loss = INFINITY, oldloss = INFINITY, lossnull = 0, lossnullF;
   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *d2 = calloc(p, sizeof(double));
   double *LP = calloc(N, sizeof(double));
   double delta, bj;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive, allconverged = 1, numconverged = 0;
   int verbose = *verbose_p;
   double *meany = calloc(K, sizeof(double));
   int pK = p * K, pK1 = p * K - 1;

   for(j = pK1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   for(j = p - 1 ; j >= 0 ; --j)
      for(i = N - 1 ; i >= 0 ; --i)
	 d2[j] +=  x[i + j * N] * x[i + j * N];

   /* null loss mean((y - mean(y))^2)*/
   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
         meany[k] += y[i];
      meany[k] /= N;
   }

   for(k = 0 ; k < K ; k++)
      for(i = N - 1 ; i >= 0 ; --i)
	 lossnull += pow(y[i + N * k] - meany[k], 2);
   lossnull /= (N * K);
   lossnullF = lossnull * 1e-8;
	 
   if(verbose)
      printf("null loss: %.5f\n", lossnull);

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      numactive = 0;
      numconverged = 0;

      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p ; j++)
	 {
	    oldloss = loss;

	    if(!active[j + p * k])
	    {
	       numconverged++;
	       continue;
	    }

	    d1 = 0;
	    for(i = N - 1 ; i >= 0 ; --i)
	       d1 += x[i + j * N] * (LP[i + N * k] - y[i + N * k]);

	    bj = b[j + p * k];
	    b[j + p * k] = soft_threshold(bj - d1 / d2[j], lambda1[k]);
	    delta = b[j + p * k] - bj;

	    loss = 0;
	    for(i = N - 1 ; i >= 0 ; --i)
	    {
	       LP[i + N * k] += x[i + j * N] * delta;
	       loss += pow(LP[i + N * k] - y[i + N * k], 2);
	    }
	    loss /= N;
	    numconverged += fabs(loss - oldloss) < lossnullF;

	    active[j + p * k] = b[j + p * k] != 0;
	    numactive += active[j + p * k];
	 }
      }

      if(verbose)
      {
	 printf("%d iter loss %.10f\n", iter, loss);
	 printf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == pK)
      {
	 if(verbose)
	    printf("all converged at iter %d\n", iter);
	 if(allconverged == 1)
	 {
	    for(j = pK1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = TRUE;
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
		  printf("terminating at iteration %d with %d active vars\n",
		     iter, numactive);
	       break;
	    }

	    allconverged = 1;
	    for(j = pK1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = TRUE;
	    }
	 }
      }
   }

   if(iter >= maxiter && verbose)
      printf("failed to converge after %d iterations\n", maxiter);

   free(active);
   free(oldactive);
   free(d2);
   free(LP);
   free(meany);
}

void lasso3(double *x, double *y, double *b,
      int *N_p, int *p_p, double *lambda1_p,
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
   double *d2 = calloc(p, sizeof(double));
   double *LP = calloc(N, sizeof(double));
   double delta, bj;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive, allconverged = 1, numconverged = 0;
   int verbose = *verbose_p;
   double meany = 0;
   double lambda1 = *lambda1_p;
   double s;

   for(j = p - 1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   for(j = p - 1 ; j >= 0 ; --j)
   {
      for(i = N - 1 ; i >= 0 ; --i)
	 d2[j] +=  x[i + j * N] * x[i + j * N];
   }

   /* null loss mean((y - mean(y))^2)*/
   for(i = N - 1 ; i >= 0 ; --i)
      meany += y[i];
   meany /= N;

   for(i = N - 1 ; i >= 0 ; --i)
      lossnull += pow(y[i] - meany, 2);
   lossnull /= N;
   lossnullF = lossnull * 1e-8;
	 
   if(verbose)
      printf("null loss: %.5f\n", lossnull);

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
	 if(fabs(s) <= lambda1)
	 {
	    b[j] = 0;
	    delta = -bj;
	 }
	 else
	 {
	    b[j] = s - lambda1 * sign(s);
	    if(fabs(b[j]) < 1e-15)
	       b[j] = 0;
	    delta = b[j] - bj;
	 }

	 loss = 0;
	 for(i = N - 1 ; i >= 0 ; --i)
	 {
	    LP[i] += x[i + j * N] * delta;
	    loss += pow(LP[i] - y[i], 2);
	 }
	 loss /= N;
	 numconverged += fabs(loss - oldloss) < lossnullF;

	 active[j] = b[j] != 0;
	 numactive += active[j];
      }

      if(verbose)
      {
	 printf("%d iter loss %.10f\n", iter, loss);
	 printf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == p)
      {
	 if(verbose)
	    printf("all converged at iter %d\n", iter);
	 if(allconverged == 1)
	 {
	    for(j = p - 1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = TRUE;
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
		  printf("terminating at iteration %d with %d active vars\n",
		     iter, numactive);
	       break;
	    }

	    if(verbose)
	       printf("active set has changed\n");

	    allconverged = 1;
	    for(j = p - 1; j >= 0 ; --j)
	    {
	       oldactive[j] = active[j];
	       active[j] = TRUE;
	    }
	 }
      }
   }

   if(iter >= maxiter && verbose)
      printf("failed to converge after %d iterations\n", maxiter);

   free(active);
   free(oldactive);
   free(d2);
   free(LP);
}

/*
 *  Assumes that columns of x matrix are
 *  standardised to zero-mean and unit-norm
 *  and y is scaled so there's no intercept
 */
void groupridge3(double *x, double *y, double *b,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3_p,
      int *grp, int *maxiter_p, double *eps_p, int *verbose_p,
      int *status)
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
   int numactive, allconverged = 1, numconverged = 0;
   int verbose = *verbose_p;
   int pK = p * K, pK1 = p * K - 1;
   double s, lambda3 = *lambda3_p;
   int q;

   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *LP = calloc(N * K, sizeof(double));
   double *meany = calloc(K, sizeof(double));
   double *loss = malloc(sizeof(double) * K),
	  *oldloss = malloc(sizeof(double) * K),
	  *lossnull = calloc(K, sizeof(double)),
	  *lossnullF = calloc(K, sizeof(double));
   double losstotal = 0, oldlosstotal = 0;

   if(verbose)
   {
      printf("N: %d p: %d K: %d\n", N, p, K);
      printf("l3: %.6f\n", lambda3);
      for(k = 0 ; k < K ; k++)
      {
	 printf("lambda1[%d]: %.6f lambda2[%d]: %.6f grp[%d]: %d\n",
	    k, lambda1[k], k, lambda2[k], k, grp[k]);
      }
   }

   for(j = pK1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   /* null loss mean((y - mean(y))^2)*/
   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
         meany[k] += y[i];
      meany[k] /= N;
   }

   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
	 lossnull[k] += pow(y[i + N * k] - meany[k], 2);
      lossnull[k] /= N;
      lossnullF[k] = lossnull[k] * 1e-8;
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
	    d1 = 0;
	    d2 = 1; // assumes standardised inputs
	    if(!active[j + p * k])
	       numconverged++;
	    else
	    {
	       for(i = N - 1 ; i >= 0 ; --i)
	          d1 += x[i + j * N] * (LP[i + N * k] - y[i + N * k]);

	       // different implementation of the soft-thresholding
	       bjk = b[j + p * k];

	       /* Apply inter-task ridge regression */
	       for(q = 0 ; q < K ; q++)
	       {
	          if(grp[k] == grp[q] && k != q)
	          {
	             d1 += lambda3 * (bjk - b[j + p * q]);
	             d2 += lambda3 * N;
	          }
	       }

	       /* Apply intra-task ridge regression */
	       s = (bjk - d1 / d2) / (1 + lambda2[k]);

	       /* Now apply intra-task lasso */
	       if(fabs(s) <= lambda1[k])
	       {
	          b[j + p * k] = 0;
	          delta = -bjk;
	       }
	       else
	       {
	          b[j + p * k] = s - lambda1[k] * sign(s);
		  if(fabs(b[j + p * k]) < 1e-15)
		     b[j + p * k] = 0;
	          delta = b[j + p * k] - bjk;
	       }

	       oldloss[k] = loss[k];
	       loss[k] = 0;
	       for(i = N - 1 ; i >= 0 ; --i)
	       {
	          LP[i + N * k] += x[i + N * j] * delta;
	          loss[k] += pow(LP[i + N * k] - y[i + N * k], 2);
	       }
	       loss[k] /= N;
	       numconverged += fabs(loss[k] - oldloss[k]) < lossnullF[k];
	    }

	    active[j + p * k] = b[j + p * k] != 0;
	    numactive += active[j + p * k];
	 }

	 losstotal += loss[k];
	 oldlosstotal += oldloss[k];
      }

      if(verbose)
      {
	 printf("%d iter loss %.10f\n", iter, losstotal);
	 printf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == pK)
      {
         if(verbose)
            printf("all converged at iter %d\n", iter);
         if(allconverged == 1)
         {
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = TRUE;
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
        	  printf("terminating at iteration %d with %d active vars\n",
        	     iter, numactive);
	       *status = TRUE;
               break;
            }

	    if(verbose)
	       printf("activeset changed at iter %d\n", iter);

            allconverged = 1;
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = TRUE;
            }
         }
      }
   }

   if(iter >= maxiter && verbose)
   {
      printf("failed to converge after %d iterations\n", maxiter);
      *status = FALSE;
   }

   free(active);
   free(oldactive);
   free(LP);
   free(meany);
   free(loss);
   free(oldloss);
   free(lossnull);
   free(lossnullF);
}

/*
 * Assumes that columns of x matrix are
 * standardised to zero-mean and unit-norm
 * and y is scaled so there's no intercept
 * 
 * Uses a different formulation of the groups, based on 
 * the symmetric K * K matrix G (diagonal is zero), G_{ij} \in {-1, 0, 1} to
 * denote positive, zero, and negative correlation, resp.,
 * rather than the group K-vector grp.
 *
 */
void groupridge4(double *x, double *y, double *b,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3_p,
      int *G, int *maxiter_p, double *eps_p, int *verbose_p,
      int *status)
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
   double s, lambda3 = *lambda3_p;
   int q, iNk;

   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *LP = calloc(N * K, sizeof(double));
   double *Err = calloc(N * K, sizeof(double));
   double *meany = calloc(K, sizeof(double));
   double *loss = malloc(sizeof(double) * K);
   double *oldloss = malloc(sizeof(double) * K);
   double *lossnull = calloc(K, sizeof(double));
   double *lossnullF = calloc(K, sizeof(double));
   double losstotal = 0, oldlosstotal = 0;
   int g;

   for(j = pK1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   /* null loss mean((y - mean(y))^2)*/
   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
         meany[k] += y[i];
      meany[k] /= N;
   }

   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
	 lossnull[k] += pow(y[i + N * k] - meany[k], 2);
      lossnull[k] /= N;
      lossnullF[k] = lossnull[k] * eps;
      loss[k] = INFINITY;
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
	    d1 = 0;
	    d2 = 1; // assumes standardised inputs
	    if(!active[j + p * k])
	       numconverged++;
	    else
	    {
	       for(i = N - 1 ; i >= 0 ; --i)
	          d1 += x[i + j * N] * Err[i + N * k];

	       // different implementation of the soft-thresholding
	       bjk = b[j + p * k];

	       /* Apply inter-task ridge regression */
	       // TODO: what about negative correlation?
	       for(q = 0 ; q < K ; q++)
	       {
		  g = G[k + q * K];
	          if(k != q && g != 0)
	          {
	             d1 += lambda3 * (bjk - sign(g) * b[j + p * q]);
	             d2 += lambda3 * N;
	          }
	       }

	       /* Apply intra-task ridge regression */
	       s = (bjk - d1 / d2) / (1 + lambda2[k]);

	       /* Now apply intra-task lasso */
	       if(fabs(s) <= lambda1[k])
	       {
	          b[j + p * k] = 0;
	          delta = -bjk;
	       }
	       else
	       {
	          b[j + p * k] = s - lambda1[k] * sign(s);
		  if(fabs(b[j + p * k]) < 1e-15) // close enough to zero
		     b[j + p * k] = 0;
	          delta = b[j + p * k] - bjk;
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
	       loss[k] /= N;
	       numconverged += fabs(loss[k] - oldloss[k]) < lossnullF[k];
	    }

	    active[j + p * k] = b[j + p * k] != 0;
	    numactive += active[j + p * k];
	 }

	 losstotal += loss[k];
	 oldlosstotal += oldloss[k];
      }

      if(verbose)
      {
	 printf("%d iter loss %.10f\n", iter, losstotal);
	 printf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == pK)
      {
         if(verbose)
            printf("all converged at iter %d\n", iter);
         if(allconverged == 1)
         {
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = TRUE;
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
        	  printf("terminating at iteration %d with %d active vars\n",
        	     iter, numactive);
	       *status = TRUE;
               break;
            }

	    if(verbose)
	       printf("activeset changed at iter %d\n", iter);

            allconverged = 1;
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = TRUE;
            }
         }
      }
   }

   if(iter >= maxiter)
   {
      if(verbose)
	 printf("failed to converge after %d iterations\n", maxiter);
      *status = FALSE;
   }

   free(active);
   free(oldactive);
   free(LP);
   free(meany);
   free(loss);
   free(oldloss);
   free(lossnull);
   free(lossnullF);
   free(Err);
}

/*
 * Assumes that columns of x matrix are
 * standardised to zero-mean and unit-norm
 * and y is scaled so there's no intercept
 * 
 * Uses a different formulation of the groups, based on 
 * the symmetric K * K matrix G (diagonal is zero), G_{ij} \in {-1, 0, 1} to
 * denote positive, zero, and negative correlation, resp.,
 * rather than the group K-vector grp.
 *
 * With warm restarts
 *
 */
void groupridge5(double *x, double *y, double *b, double *LP,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3_p,
      int *G, int *maxiter_p, double *eps_p, int *verbose_p,
      int *status)
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
   double s, lambda3 = *lambda3_p;
   int q, iNk;

   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *Err = calloc(N * K, sizeof(double));
   double *meany = calloc(K, sizeof(double));
   double *loss = malloc(sizeof(double) * K);
   double *oldloss = malloc(sizeof(double) * K);
   double *lossnull = calloc(K, sizeof(double));
   double *lossnullF = calloc(K, sizeof(double));
   double losstotal = 0, oldlosstotal = 0;
   int g;

   for(j = pK1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   /* null loss mean((y - mean(y))^2)*/
   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
         meany[k] += y[i];
      meany[k] /= N;
   }

   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
	 lossnull[k] += pow(y[i + N * k] - meany[k], 2);
      lossnull[k] /= N;
      lossnullF[k] = lossnull[k] * eps;
      //loss[k] = INFINITY;
   }

   for(k = K - 1 ; k >= 0 ; --k)
   {
      loss[k] = 0;
      for(i = N - 1 ; i >= 0 ; --i)
      {
         iNk = i + N * k;
         Err[iNk] = LP[iNk] - y[iNk];
         loss[k] += Err[iNk] * Err[iNk];
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
	    d1 = 0;
	    d2 = 1; // assumes standardised inputs
	    if(!active[j + p * k])
	       numconverged++;
	    else
	    {
	       for(i = N - 1 ; i >= 0 ; --i)
	          d1 += x[i + j * N] * Err[i + N * k];

	       // different implementation of the soft-thresholding
	       bjk = b[j + p * k];

	       /* Apply inter-task ridge regression */
	       // TODO: what about negative correlation?
	       for(q = 0 ; q < K ; q++)
	       {
		  g = G[k + q * K];
	          if(k != q && g != 0)
	          {
	             d1 += lambda3 * (bjk - sign(g) * b[j + p * q]);
	             d2 += lambda3 * N;
	          }
	       }

	       /* Apply intra-task ridge regression */
	       s = (bjk - d1 / d2) / (1 + lambda2[k]);

	       /* Now apply intra-task lasso */
	       if(fabs(s) <= lambda1[k])
	       {
	          b[j + p * k] = 0;
	          delta = -bjk;
	       }
	       else
	       {
	          b[j + p * k] = s - lambda1[k] * sign(s);
		  if(fabs(b[j + p * k]) < 1e-15) // close enough to zero
		     b[j + p * k] = 0;
	          delta = b[j + p * k] - bjk;
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
	       loss[k] /= N;
	       numconverged += fabs(loss[k] - oldloss[k]) < lossnullF[k];
	    }

	    active[j + p * k] = b[j + p * k] != 0;
	    numactive += active[j + p * k];
	 }

	 losstotal += loss[k];
	 oldlosstotal += oldloss[k];
      }

      if(verbose)
      {
	 printf("%d iter loss %.10f\n", iter, losstotal);
	 printf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == pK)
      {
         if(verbose)
            printf("all converged at iter %d\n", iter);
         if(allconverged == 1)
         {
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = TRUE;
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
        	  printf("terminating at iteration %d with %d active vars\n",
        	     iter, numactive);
	       *status = TRUE;
               break;
            }

	    if(verbose)
	       printf("activeset changed at iter %d\n", iter);

            allconverged = 1;
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = TRUE;
            }
         }
      }
   }

   if(iter >= maxiter)
   {
      if(verbose)
	 printf("failed to converge after %d iterations\n", maxiter);
      *status = FALSE;
   }

   free(active);
   free(oldactive);
   free(meany);
   free(loss);
   free(oldloss);
   free(lossnull);
   free(lossnullF);
   free(Err);
}

