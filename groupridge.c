#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

/*
 * lambda1: lasso penalty within each task
 * lambda2: ridge penalty within each task
 * lambda3: ridge penalty across tasks
 */
void groupridge(double *x, double *y, double *B,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3,
      int *grp, int *maxiter_p, double *eps_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k, q;
   int iter;
   double d1, d2;
   double loss = INFINITY, oldloss = INFINITY;
   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *v = malloc(sizeof(double) * p);
   double *LP = calloc(N * K, sizeof(double));
   double s, s2, delta, Bjk;
   double eps = *eps_p;
   int maxiter = *maxiter_p;

   printf("eps: %.6f\n", eps);

   for(i = p * K - 1 ; i >= 0 ; --i)
      active[i] = oldactive[i] = 1.0;

   for(j = p - 1 ; j >= 0 ; --j)
   {
      v[j] = 0;
      for(i = N - 1 ; i >= 0 ; --i)
	 v[j] +=  x[i + j * N] * x[i + j * N];
   }

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p ; j++)
	 {
	    d1 = 0;
	    d2 = v[j];
	    for(i = N - 1 ; i >= 0 ; --i)
	       d1 += x[i + j * N] * (LP[i + k * N] - y[i + k * N]);

	    s = B[j + k * p] - d1 / d2;

	    /* beta is zero, skip the other penalties */
	    if(fabs(s) <= lambda1[k])
	    {
	       delta = - B[j + k * p];
	       B[j + k * p] = 0;
	       for(i = 0 ; i < N ; i++)
	          LP[i + k * N] += x[i + j * N] * delta;
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
	    for(i = N - 1 ; i >= 0 ; --i)
	       LP[i + k * N] += x[i + j * N] * delta;
	 }
      }

      loss = 0;
      for(i = N * K - 1 ; i >= 0 ; --i)
	 loss += pow(LP[i] - y[i], 2);
      loss /= (N * K);
      printf("iter: %d oldloss: %.5f loss: %.5f\n", iter, oldloss,
	    loss);
      if(fabs(oldloss - loss) < eps)
      {
	 printf("converged after %d iterations\n", iter);
	 break;
      }
      oldloss = loss;
   }

   printf("done\n");
   printf("final loss: %.7f\n", loss);

   free(active);
   free(oldactive);
   free(v);
   free(LP);
}

