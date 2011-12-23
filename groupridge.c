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

void groupridge(double *x, double *y, double *B, int *N_p, int *p_p, int *K_p,
      double *lambda1_p, double *lambda2_p, double *lambda3_p, int *g)
{
   double lambda1 = *lambda1_p,
	  lambda2 = *lambda2_p,
	  lambda3 = *lambda3_p;
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
   int iter;
   double d1, d2;
   double loss = 0, olsloss = 1e9;
   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   double *var = malloc(sizeof(double) * p);
   double *LP = calloc(N * K, sizeof(double));
   double s, s2, delta;

   for(i = p * K - 1 ; i >= 0 ; --i)
      active[i] = oldactive[i] = 1.0;

   for(j = p - 1 ; j >= 0 ; --j)
   {
      var[j] = 0;
      for(i = N - 1 ; i >= 0 ; --i)
	 var[j] +=  x[i + j * N] * x[i + j * N];
      var[j] = var[j] / (N - 1);
   }

   for(iter = 0 ; iter < 500 ; iter++)
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

	    /*if(lambda1 > 0)
	    {
	       s = B[j + k * p] - d1 / d2;
	       s2 = sign(s) * fmax(fabs(s) - lambda1, 0);

	       if(s2 == 0)
	       {
		  B[j + k * p] = active[j + k * p] = 0;
		  delta = B[j + k * p];
		  for(i = 0 ; i < N ; i++)
		     LP[i + k * N] = x[i * j * N] * delta;
		  break;
	       }

	       d1 += lambda1 * sign(B[j + k * p]);
	    }*/

	    delta = d1 / d2;
	    B[j + k * p] -= delta;
	    for(i = 0 ; i < N ; i++)
	       LP[i + k * N] -= x[i + j * N] * delta;
	 }
      }
      printf("iter %d\r", iter);
   }

   printf("done\n");

   free(active);
   free(oldactive);
   free(var);
   free(LP);
}

