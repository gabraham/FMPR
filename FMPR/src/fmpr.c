#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>     /* for Calloc/Free */
#include <R_ext/Applic.h> /* for dgemm */

#define TRUE 1
#define FALSE 0
#define ZERO_THRESH 1e-10
#define ZERO_VAR_THRESH 1e-6

static void crossprod(double *x, int nrx, int ncx,
   double *y, int nry, int ncy, double *z)
{
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
      x, &nrx, y, &nry, &zero, z, &ncx);
}

static double dotprod(double *x, double *y, int n)
{
    int one = 1.0;
    return F77_CALL(ddot)(&n, x, &one, y, &one);
}

/* really, it's just equivalent to max(abs(crossprod(X, Y))) and probably
 * slower */
void maxlambda1(double *x, double *y, double *lambda,
      int *N_p, int *p_p, int *K_p)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
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
	 
   if(verbose > 1)
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

	 // soft-thresholding
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

      if(verbose > 1)
      {
	 Rprintf("%d iter loss %.10f\n", iter, loss);
	 Rprintf("%d converged at iter %d\n", numconverged, iter);
      }

      if(numconverged == p)
      {
	 if(verbose > 1)
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
	       if(verbose > 1)
		  Rprintf("terminating at iteration %d with %d active vars\n",
		     iter, numactive);
	       break;
	    }

	    if(verbose > 1)
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

   if(iter >= maxiter && verbose > 1)
      Rprintf("failed to converge after %d iterations\n", maxiter);

   free(active);
   free(oldactive);
   free(d2);
   free(LP);
   free(ignore);
}

/* Uses the C from gennetwork
 */
void fmpr(double *X, double *Y, double *B,
      double *LP, int *N_p, int *p_p, int *K_p,
      double *lambda_p, double *lambda2_p, double *gamma_p,
      double *C, int *pairs, int *edges, int *maxiter_p,
      double *eps_p, int *verbose_p, int *status, int *iter_p,
      int *numactive_p, int *divbyN)
{
   int N = *N_p,
       p = *p_p,
       K = *K_p;
   int i, j, k;
   int iter;
   double d1;
   double delta, Bjk;
   double eps = *eps_p;
   int maxiter = *maxiter_p;
   int numactive = 0,
       allconverged = 1,
       numconverged = 0;
   int verbose = *verbose_p;
   int pK = p * K, pK1 = p * K - 1;
   double s, lambda = *lambda_p, gamma = *gamma_p;
   double oneOnN = (*divbyN) ? 1.0 / N : 1.0;
   int iNk, jpk, iNj;
   int nE = K * (K - 1) / 2, e;
   double df1, df2;
   double oneOn2N = 1.0 / (2.0 * N);
   double sv;
   int v1, v2, l, kK1;

   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   int *ignore = calloc(p * K, sizeof(int));
   double *Err = calloc(N * K, sizeof(double));
   double *meanY = calloc(K, sizeof(double));
   double *loss = calloc(K, sizeof(double));
   double *oldloss = calloc(K, sizeof(double));
   double *lossnull = calloc(K, sizeof(double));
   double *lossnullF = calloc(K, sizeof(double));
   double *d2 = calloc(p * K, sizeof(double));
   double losstotal = 0, oldlosstotal = 0;
   double *oneOnLambda2PlusOne = malloc(sizeof(double) * K);
   double *CC = calloc(K * K, sizeof(double));
   double *diagCC = calloc(K, sizeof(double));
   double *BCC = calloc(p * K, sizeof(double));
   double *colsumsB = calloc(K, sizeof(double));

   int dofusion = K > 1 && gamma > 0;

   /* per task sum of losses */
   //double *sumsC = calloc(K, sizeof(double));

   if(nE > 1)
   {
      crossprod(C, nE, K, C, nE, K, CC);
      for(k = 0 ; k < K ; k++)
	 diagCC[k] = CC[k * K + k];
   }

   /* null loss mean((y - mean(y))^2)*/
   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
         meanY[k] += Y[i];
      meanY[k] *= oneOn2N;

      oneOnLambda2PlusOne[k] = 1.0 / (1.0 + *lambda2_p);
   }

   /* setup second derivatives of loss and which variables to ignore due to
    * small variance */
   for(k = 0 ; k < K ; k++)
   {
      for(j = p - 1 ; j >= 0 ; --j) 
      {
	 jpk = j + p * k;
	 for(i = N - 1 ; i >= 0 ; --i)
	    d2[jpk] +=  X[i + j * N] * X[i + j * N] * oneOnN;
	 //ignore[jpk] = (d2[jpk] <= 1.0 / N);
	 ignore[jpk] = (d2[jpk] <= ZERO_VAR_THRESH);
      }
   }

   for(j = pK1 ; j >= 0 ; --j)
      oldactive[j] = active[j] = !ignore[j];

   for(k = 0 ; k < K ; k++)
   {
      loss[k] = 0;
      for(i = N - 1 ; i >= 0 ; --i)
      {
	 iNk = i + N * k;
	 lossnull[k] += pow(Y[iNk] - meanY[k], 2);
	 Err[iNk] = LP[iNk] - Y[iNk];
	 loss[k] += Err[iNk] * Err[iNk];
      }
      loss[k] *= oneOn2N;
      lossnull[k] *= oneOn2N;
      lossnullF[k] = lossnull[k] * eps;
   }

   for(iter = 1 ; iter <= maxiter ; iter++)
   {
      numactive = 0;
      numconverged = 0;
      losstotal = 0;
      oldlosstotal = 0;

      for(k = 0 ; k < K ; k++)
      {
	 kK1 = (K - 1) * k;

	 for(j = 0 ; j < p ; j++)
	 {
	    jpk = j + p * k;

	    if(!active[jpk])
	    {
	       numconverged++;
#ifdef DEBUG
	       Rprintf("skipping inactive k=%d j=%d\n", k, j+1);
#endif
	    }
	    else
	    {
	       /* 1st derivative of loss function wrt B_jk */
	       //crossprod(X + N * j, N, 1, Err + N * k, N, 1, &d1);
	       d1 = dotprod(X + N * j, Err + N * k, N);
	       d1 *= oneOnN; 
	       Bjk = B[jpk];
	       
	       /* 1st/2nd derivative of fusion loss */
	       df1 = 0;
	       df2 = 0;
	       
	       if(dofusion)
	       {
		  // For most tasks, C is are zero for each pair,
	       	  // and can be skipped, as for each task k there are
	       	  // only K-1 other tasks it shares an edge with
	       	  for(l = 0 ; l < K - 1 ; l++)
	       	  {
	       	     e = edges[l + kK1];
	       	     v1 = pairs[e];
	       	     v2 = pairs[e + nE];
	       	     sv = B[j + v1 * p] * C[e + v1 * nE]
	       	        + B[j + v2 * p] * C[e + v2 * nE];
	       	     df1 += sv * C[e + k * nE];
	       	  }

	       	  df1 *= gamma;

	       	  /* 2nd derivative of fusion loss */
	       	  df2 = gamma * diagCC[k];
	       }

	       //s = Bjk - (d1 + pd1) / L;
	       s = Bjk - (d1 + df1) / (d2[jpk] + df2);

	       /* lasso soft-thresholding */
	       //if(fabs(s) <= lambda)
	       //{
	       //   B[jpk] = 0;
	       //   delta = -Bjk;
	       //}
	       //else
	       //{
	       //   B[jpk] = (s - lambda * sign(s)) * oneOnLambda2PlusOne[k];
	       //   //if(fabs(b[jpk]) < ZERO_THRESH) // numerically close enough to zero
	       //   // b[jpk] = 0;
	       //   delta = B[jpk] - Bjk;
	       //}
	       //B[jpk] = sign(s) * fmax(fabs(s) - lambda / (d2[jpk] + df2), 0);
	       B[jpk] = sign(s) * fmax(fabs(s) - lambda, 0);

	       /* close enough to zero */
	       if(fabs(B[jpk]) < ZERO_THRESH)
		  B[jpk] = 0;
	       delta = B[jpk] - Bjk;

	       #ifdef DEBUG
	       Rprintf("[k=%d j=%d] d1=%.6f df1=%.6f d2=%.6f df2=%.6f s=%.6f \
delta=%.6f beta_old=%.6f beta_new=%.6f active:%d\n",
		     k, j+1, d1, df1, d2[jpk], df2, s, delta, Bjk, B[jpk], active[jpk]);
	       #endif
 
	       //if(B[jpk] * Bjk < 0)
	       //{
	       //   Rprintf("sign mismatch, s: %.6f old: %.6f new: %.6f lambda: %.6f\n",
	       //      s, Bjk, B[jpk], lambda);
	       //   return;
	       //}


	       /* Update loss and errors based on new estimates. */
	       oldloss[k] = loss[k];
	       loss[k] = 0;
	       iNk = N * k + N - 1;
	       iNj = N * j + N - 1;
	       for(i = N - 1 ; i >= 0 ; --i)
	       {
		  LP[iNk] += X[iNj] * delta;
		  Err[iNk] = LP[iNk] - Y[iNk];
		  //loss[k] += Err[iNk] * Err[iNk];
		  iNk--;
		  iNj--;
	       }
	       loss[k] *= oneOn2N;

	       /* update l1 loss */
	       colsumsB[k] += fabs(B[jpk]) - fabs(Bjk);
	       loss[k] += lambda * colsumsB[k];

	       /* update fusion loss */
	       //tmp = delta * (2 * B[jpk] + delta);
	       //for(e = 0 ; e < nE ; e++)
	       //   sumsC[k] += C[e + k * nE] * C[e + k * nE] * tmp;
	       //loss[k] += 0.5 * gamma * sumsC[k];
	       
	       //numconverged += fabs(loss[k] - oldloss[k]) / fabs(loss[k]) < eps;
	       //if(isnan(delta))
	       //{
	       //   Rprintf("delta is bad: %.6f\n", delta);
	       //   return;
	       //}
	       if((delta == 0 && Bjk == 0) || fabs(delta) / fabs(Bjk) < eps)
		  numconverged++;
	       //else
		  //Rprintf("delta: %.6f s: %.6f\n", delta,  fabs(delta) / fabs(Bjk));
	       
	       active[jpk] = (B[jpk] != 0);
	    }

	    numactive += active[jpk];
	 }

#ifdef DEBUG
	 Rprintf("[end of task k=%d] numconverged=%d numactive=%d\n",
	    k, numconverged, numactive);
	 Rprintf("----------------------\n");
#endif

	 losstotal += loss[k];
	 oldlosstotal += oldloss[k];
      }

      if(verbose > 1)
      {
	 Rprintf("%d iter loss %.10f\n", iter, losstotal);
	 Rprintf("%d converged at iter %d\n", numconverged, iter);
	 Rprintf("%d numactive at iter %d\n", numactive, iter);
      }

      /* active-set convergence */
      if(numconverged == pK)
      {
         if(verbose > 1)
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
               if(verbose > 1)
        	  Rprintf("terminating at iteration %d with %d active vars\n",
        	     iter, numactive);
	       *status = TRUE;
               break;
            }

	    if(verbose > 1)
	       Rprintf("activeset changed at iter %d\n", iter);

            allconverged = 1;
            for(j = pK1; j >= 0 ; --j)
            {
               oldactive[j] = active[j];
               active[j] = !ignore[j];
            }
         }
      }
      else if(verbose > 1)
	 Rprintf("active set not yet converged at iter %d (numconverged=%d)\n",
	    iter, numconverged);
   }

   if(iter >= maxiter)
   {
      if(verbose > 1)
	 Rprintf("fmpr_threshold_warm failed to converge after %d iterations\n",
	    maxiter);
      *status = FALSE;
   }
   *iter_p = iter;
   *numactive_p = numactive;

   free(active);
   free(oldactive);
   free(meanY);
   free(loss);
   free(oldloss);
   free(lossnull);
   free(lossnullF);
   free(Err);
   free(d2);
   free(ignore);
   free(oneOnLambda2PlusOne);
   free(CC);
   free(diagCC);
   free(colsumsB);
   //free(sumsC);
   free(BCC);
}

