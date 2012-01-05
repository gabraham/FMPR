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

void maxlambda1(double *x, double *y, double *lambda,
      int *N_p, int *p_p, int *K_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k, m;
   double d1, d2, s, xij;

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
	 /*s = fabs((d1 / N) / (d2 / (N - 1)));*/
	 s = fabs(d1 / d2);
	 lambda[k] = (lambda[k] > s ? lambda[k] : s);
      }
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

      if(fabs(oldloss - loss) < eps)
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

/*
 * With active set convergence
 *
 * lambda1: lasso penalty within each task
 * lambda2: ridge penalty within each task
 * lambda3: ridge penalty across tasks
 */
void groupridge2(double *x, double *y, double *B,
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
   double loss = INFINITY, oldloss = INFINITY, lossk;
   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *v = calloc(p, sizeof(double));
   double *LP = calloc(N * K, sizeof(double));
   double s, s2, delta, Bjk;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive, allconverged = 1, numconverged = 0;
   int pK1 = p * K - 1;
   int verbose = *verbose_p;
   int *converged = calloc(p * K, sizeof(int));

   for(j = pK1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   for(j = p - 1 ; j >= 0 ; --j)
      for(i = N - 1 ; i >= 0 ; --i)
	 v[j] +=  x[i + j * N] * x[i + j * N];
   
   for(iter = 0 ; iter < maxiter ; iter++)
   {
      numactive = 0;
      oldloss = loss;
      loss = 0;
      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p ; j++)
	 {
	    if(!active[j + k * p])
	       continue;

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
	    }
	    else
	    {
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
	    }

	    lossk = 0;
	    for(i = N - 1 ; i >= 0 ; --i)
	    {
	       LP[i + k * N] += x[i + j * N] * delta;
	       lossk += pow(LP[i + k * N] - y[i + k * N], 2);
	    }
	    loss += lossk / N; 
	    converged[j + k * p] = fabs(loss - oldloss) < 1e-6;
	    numconverged += converged[j + k * p];

	    active[j + k * p] = B[j + k * p] != 0;
	    numactive += active[j + k * p];
	 }
      }

      if(numconverged == p * K)
      {
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
   free(v);
   free(LP);
   free(converged);
}

