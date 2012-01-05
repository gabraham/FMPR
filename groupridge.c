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
//void groupridge2(double *x, double *y, double *B,
//      int *N_p, int *p_p, int *K_p,
//      double *lambda1, double *lambda2, double *lambda3,
//      int *grp, int *maxiter_p, double *eps_p, int *verbose_p)
//{
//   int N = *N_p,
//       p = *p_p,
//       K = *K_p;
//   int i, j, k, q;
//   int iter;
//   double d1, d2;
//   double loss = INFINITY, oldloss = INFINITY, lossk, lossnull = 0, lossnullF;
//   int *active = malloc(sizeof(int) * p * K);
//   int *oldactive = malloc(sizeof(int) * p * K);
//   double *v = calloc(p, sizeof(double));
//   double *LP = calloc(N * K, sizeof(double));
//   double s, s2, delta, Bjk;
//   double eps = *eps_p;
//   int maxiter = *maxiter_p;
//   int numactive, allconverged = 1, numconverged = 0;
//   int pK = p * K, pK1 = p * K - 1, ikN;
//   int verbose = *verbose_p;
//   double *meany = calloc(K, sizeof(double));
//
//   for(j = pK1 ; j >= 0 ; --j)
//      active[j] = oldactive[j] = TRUE;
//
//   for(j = p - 1 ; j >= 0 ; --j)
//      for(i = N - 1 ; i >= 0 ; --i)
//	 v[j] +=  x[i + j * N] * x[i + j * N];
//
//   /* null loss mean((y - mean(y))^2)*/
//   for(k = K - 1 ; k >= 0 ; --k)
//   {
//      for(i = N - 1 ; i >= 0 ; --i)
//	 meany[k] += y[i + k * N];
//      meany[k] /= N;
//   }
//
//   for(k = K - 1 ; k >= 0 ; --k)
//      for(i = N - 1 ; i >= 0 ; --i)
//	 lossnull += pow(y[i + k * N] - meany[k], 2);
//   lossnull /= (N * K);
//   lossnullF = lossnull * 1e-6;
//	 
//   if(verbose)
//      printf("null loss: %.5f\n", lossnull);
//
//   for(iter = 0 ; iter < maxiter ; iter++)
//   {
//      numactive = 0;
//      loss = 0;
//      numconverged = 0;
//
//      for(k = 0 ; k < K ; k++)
//      {
//	 for(j = 0 ; j < p ; j++)
//	 {
//	    oldloss = loss;
//	    if(!active[j + k * p])
//	    {
//	       numconverged++;
//	       continue;
//	    }
//
//	    d1 = 0;
//	    d2 = v[j];
//	    for(i = N - 1 ; i >= 0 ; --i)
//	    {
//	       q = k * N;
//	       d1 += x[i + j * N] * (LP[i + q] - y[i + q]);
//	    }
//
//	    //s = B[j + k * p] - d1 / d2;
//
//	    /* beta is zero, skip the other penalties */
//	    //if(fabs(s) <= lambda1[k])
//	    //{
//	    //   delta = - B[j + k * p];
//	    //   B[j + k * p] = 0;
//	    //}
//	    //else
//	    //{
//	    //   Bjk = B[j + k * p];
//	    //   /* lasso penalty when beta!=0, no 2nd derivative */
//	    //   //d1 += lambda1[k] * sign(Bjk);
//
//	    //   /* ridge penalty intra-class */
//	    //   //d1 += lambda2[k] * Bjk;
//	    //   //d2 += lambda2[k];
//
//	    //   /* ridge penalty inter-class */
//	    //   //for(q = 0 ; q < K ; q++)
//	    //   //{
//	    //   //   if(grp[k] == grp[q] && k != q)
//	    //   //   {
//	    //   //      d1 += lambda3[k] * (Bjk - B[j + q * p]);
//	    //   //      d2 += lambda3[k] * N;
//	    //   //   }
//	    //   //}
//	    //   //
//	    //   //s = d1 / d2;
//	    //   s2 = Bjk - s;
//	    //   delta = s2 - Bjk;
//	    //   B[j + k * p] = s2;
//	    //}
//
//	    Bjk = B[j + k * p];
//	    B[j + k * p] = soft_threshold(Bjk - d1 / d2, lambda1[k]);
//	    delta = B[j + k * p] - Bjk;
//
//	    lossk = 0;
//	    for(i = N - 1 ; i >= 0 ; --i)
//	    {
//	       ikN = i + k * N;
//	       LP[ikN] += x[i + j * N] * delta;
//	       lossk += pow(LP[ikN] - y[ikN], 2);
//	    }
//	    loss += lossk / N; 
//	    numconverged += fabs(loss - oldloss) < lossnullF;
//
//	    active[j + k * p] = B[j + k * p] != 0;
//	    numactive += active[j + k * p];
//	 }
//      }
//
//      if(verbose)
//      {
//	 printf("%d iter loss %.10f\n", iter, loss);
//	 printf("%d converged at iter %d\n", numconverged, iter);
//      }
//
//      if(numconverged == pK)
//      {
//	 if(verbose)
//	    printf("all converged at iter %d\n", iter);
//	 if(allconverged == 1)
//	 {
//	    for(j = pK1; j >= 0 ; --j)
//	    {
//	       oldactive[j] = active[j];
//	       active[j] = TRUE;
//	    }
//	    allconverged = 2;
//	 }
//	 else
//	 {
//	    for(j = pK1 ; j >= 0 ; --j)
//	       if(active[j] != oldactive[j])
//		  break;
//	    if(j < 0)
//	    {
//	       if(verbose)
//		  printf("terminating at iteration %d with %d active vars\n",
//		     iter, numactive);
//	       break;
//	    }
//
//	    allconverged = 1;
//	    for(j = pK1; j >= 0 ; --j)
//	    {
//	       oldactive[j] = active[j];
//	       active[j] = TRUE;
//	    }
//	 }
//      }
//   }
//
//   if(iter >= maxiter && verbose)
//      printf("failed to converge after %d iterations\n", maxiter);
//
//   free(active);
//   free(oldactive);
//   free(v);
//   free(LP);
//   free(meany);
//}

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
 */
void groupridge3(double *x, double *y, double *b,
      int *N_p, int *p_p, int *K_p,
      double *lambda1, double *lambda2, double *lambda3_p,
      int *grp, int *maxiter_p, double *eps_p, int *verbose_p)
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
   //double *d2 = calloc(p, sizeof(double));
   double *LP = calloc(N * K, sizeof(double));
   double *meany = calloc(K, sizeof(double));
   double *loss = malloc(sizeof(double) * K),
	  *oldloss = malloc(sizeof(double) * K),
	  *lossnull = calloc(K, sizeof(double)),
	  *lossnullF = calloc(K, sizeof(double));
   double losstotal = 0, oldlosstotal = 0;

   printf("N: %d p: %d K: %d\n", N, p, K);

   for(j = pK1 ; j >= 0 ; --j)
      active[j] = oldactive[j] = TRUE;

   //for(j = p - 1 ; j >= 0 ; --j)
   // for(i = N - 1 ; i >= 0 ; --i)
   // d2[j] +=  x[i + j * N] * x[i + j * N];

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
	 
   //   if(verbose)
   //      printf("null loss: %.5f lossnullF: %.10f\n",
   //	    lossnull, lossnullF);

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
	    {
	       numconverged++;
	       //printf("j: %d k: %d lossdiff: (NA)\n", j, k);
	    }
	    else
	    {
	       for(i = N - 1 ; i >= 0 ; --i)
	          d1 += x[i + j * N] * (LP[i + N * k] - y[i + N * k]);

	       // different implementation of the soft-thresholding
	       bjk = b[j + p * k];
	       //b[j + p * k] = soft_threshold(bjk - d1 / d2, lambda1[k]);
	       //delta = b[j + p * k] - bjk;

	       /* Apply inter-task ridge regression */
	       for(q = 0 ; q < K ; q++)
	       {
	          if(grp[k] == grp[q] && k != q)
	          {
	             d1 += lambda3 * (bjk - b[j + p * q]);
	             // TODO: fix the second derivative
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
	       //printf("j: %d k: %d loss: %.10f oldloss: %.10f diff: %.10f d1: %.10f\n", j, k,
	       //	     loss[k], oldloss[k], fabs(loss[k] - oldloss[k]), d1);
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
      printf("failed to converge after %d iterations\n", maxiter);

   free(active);
   free(oldactive);
   //free(d2);
   free(LP);
   free(meany);
   free(loss);
   free(oldloss);
   free(lossnull);
   free(lossnullF);
}

