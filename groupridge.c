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
   int i, j, k;
   double d1, d2, s;

   for(k = 0 ; k < K ; k++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 d1 = 0;
	 d2 = 0;
	 for(i = 0 ; i < N ; i++)
	 {
	    d1 += -x[i + j * N] * y[i + k * N];
	    d2 += pow(x[i + j * N], 2);
	 }
	 s = fabs((d1 / N) / (d2 / (N - 1)));
	 lambda[k] = (lambda[k] > s ? lambda[k] : s);
      }
   }
}

void groupridge(double *x, double *y, double *B, int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3, int *g)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
   int iter;
   double d1, d2;
   double loss = INFINITY, oldloss = 1e9;
   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *var = malloc(sizeof(double) * p);
   double *LP = calloc(N * K, sizeof(double));
   double s, s2, delta;
   int const maxiter = 5e2L;

   for(i = p * K - 1 ; i >= 0 ; --i)
      active[i] = oldactive[i] = 1.0;

   for(j = p - 1 ; j >= 0 ; --j)
   {
      var[j] = 0;
      for(i = N - 1 ; i >= 0 ; --i)
	 var[j] +=  x[i + j * N] * x[i + j * N];
      var[j] = var[j] / (N - 1);
   }

   for(iter = 0 ; iter < maxiter ; iter++)
   {
      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p ; j++)
	 {
	    d1 = 0;
	    d2 = var[j];
	    for(i = 0 ; i < N ; i++)
	       d1 += x[i + j * N] * (LP[i + k * N] - y[i + k * N]);
	    d1 /= N;

	    s = B[j + k * p] - d1 / d2;
	    s2 = sign(s) * fmax(fabs(s) - lambda1[k], 0);
	    if(s2 == 0)
	    {
	       delta = s2 - B[j + k * p];
	       B[j + k * p] = s2;
	       for(i = 0 ; i < N ; i++)
	          LP[i + k * N] += x[i + j * N] * delta;
	       continue;
	    }

	    /* from here on, beta isn't zero */

	    printf("s2: %.5f\n", s2); 

	    s = (d1 + lambda1[k] * sign(B[j + k * p])) / d2;
	    s2 = B[j + k * p] - s;
	    delta = s2 - B[j + k * p];
	    B[j + k * p] = s2;
	    for(i = 0 ; i < N ; i++)
	       LP[i + k * N] += x[i + j * N] * delta;

	 }
      }
      printf("iter: %d loss: %.5f\n", iter, loss);

      loss = 0;
      for(i = N * K - 1 ; i >= 0 ; --i)
	 loss += pow(LP[i] - y[i], 2);
      loss /= N;
      if(fabs(oldloss - loss) < 1e-5)
	 break;
      oldloss = loss;
   }

   printf("done\n");

   free(active);
   free(oldactive);
   free(var);
   free(LP);
}

