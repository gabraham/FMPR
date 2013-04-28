#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>     /* for Calloc/Free */
#include <R_ext/Applic.h> /* for dgemm */

#define TRUE 1
#define FALSE 0
#define ZERO_THRESH 1e-10
#define ZERO_VAR_THRESH 1e-6

void timestamp()
{
   time_t t = time(NULL);
   char *s = asctime(localtime(&t));
   s[strlen(s) - 1] = '\0';
   printf("[%s]", s);
}

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
       * to truly make everything zero, otherwise numerical
       * error might cause trouble */
      lambda[k] += eps;
   }
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
       allconverged = 0;
   int verbose = *verbose_p;
   int pK = p * K, pK1 = p * K - 1;
   double s, lambda = *lambda_p, gamma = *gamma_p;
   double oneOnN = (*divbyN) ? 1.0 / N : 1.0;
   int iNk, jpk, iNj;
   int nE = K * (K - 1) / 2, e, m, jep;
   double df1, df2;
   double oneOn2N = 1.0 / (2.0 * N);
   double sv;
   int v1, v2, l, kK1;
   int dofusion = K > 1 && gamma > 0;
   double tmp;
   double loss = 0, l1loss = 0, floss = 0, lossold = 0;

   int *active = malloc(sizeof(int) * p * K);
   int *oldactive = malloc(sizeof(int) * p * K);
   int *ignore = calloc(p * K, sizeof(int));
   double *Err = calloc(N * K, sizeof(double));
   double *d2 = calloc(p * K, sizeof(double));
   double *oneOnLambda2PlusOne = calloc(K, sizeof(double));
   /*double *BCT = NULL; */
   double *CC = NULL;
   double *diagCC = NULL;

   if(dofusion)
   {
      CC = calloc(K * K, sizeof(double));
      diagCC = calloc(K, sizeof(double));
      /*BCT = calloc(p * nE, sizeof(double)); */

      crossprod(C, nE, K, C, nE, K, CC);
      for(k = 0 ; k < K ; k++)
	 diagCC[k] = CC[k * K + k];
      free(CC);
   }

   for(k = 0 ; k < K ; k++)
      oneOnLambda2PlusOne[k] = 1.0 / (1.0 + *lambda2_p);

   /* setup second derivatives of loss and which
    * variables to ignore due to small variance */
   //for(k = 0 ; k < K ; k++)
   //{
   //   for(j = p - 1 ; j >= 0 ; --j) 
   //   {
   //      jpk = j + p * k;
   //      for(i = N - 1 ; i >= 0 ; --i)
   //         d2[jpk] +=  X[i + j * N] * X[i + j * N] * oneOnN;
   //      ignore[jpk] = (d2[jpk] <= ZERO_VAR_THRESH);
   //   }
   //}

   for(j = pK1 ; j >= 0 ; --j)
   {
      oldactive[j] = active[j] = (B[j] != 0);
   }

   for(k = 0 ; k < K ; k++)
   {
      for(i = N - 1 ; i >= 0 ; --i)
      {
	 iNk = i + N * k;
	 Err[iNk] = LP[iNk] - Y[iNk];
      }
   }

   for(iter = 1 ; iter <= maxiter ; iter++)
   {
      numactive = 0;

      for(k = 0 ; k < K ; k++)
      {
	 kK1 = (K - 1) * k;

	 for(j = 0 ; j < p ; j++)
	 {
	    jpk = j + p * k;

	    if(active[jpk])
	    {
	       /* 1st derivative of loss function wrt B_jk */
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
	       	  for(l = K - 2 ; l >= 0 ; --l)
	       	  {
	       	     e = edges[l + kK1];
	       	     v1 = pairs[e];
	       	     v2 = pairs[e + nE];
	       	     sv = B[j + v1 * p] * C[e + v1 * nE]
	       	        + B[j + v2 * p] * C[e + v2 * nE];
	       	     df1 += sv * C[e + k * nE];
	       	  }

	       	  df1 *= gamma;
	       	  df2 = gamma * diagCC[k];
	       }

	       /*s = Bjk - (d1 + df1) / (d2[jpk] + df2);*/
	       s = Bjk - (d1 + df1) / (1 + df2);
	       B[jpk] = sign(s) * fmax(fabs(s) - lambda, 0) 
		  * oneOnLambda2PlusOne[k];

	       /* close enough to zero */
	       if(fabs(B[jpk]) < ZERO_THRESH)
		  B[jpk] = 0;

	       delta = B[jpk] - Bjk;

	       if(delta != 0)
	       {
		  iNk = N * k + N - 1;
	       	  iNj = N * j + N - 1;
	       	  for(i = N - 1 ; i >= 0 ; --i)
	       	  {
	       	     LP[iNk] += X[iNj] * delta;
	       	     Err[iNk] = LP[iNk] - Y[iNk];
	       	     iNk--;
	       	     iNj--;
	       	  }
	       }
	    	       
	       active[jpk] = (B[jpk] != 0);
	    }

	    numactive += active[jpk];
	 }
      }

      loss = 0;
      for(i = 0 ; i < N ; i++)
      {
	 for(k = 0 ; k < K ; k++)
	 {
	    m = i + k * N;
	    tmp = (LP[m] - Y[m]) * (LP[m] - Y[m]);
	    loss += tmp;
	 }
      }
      loss *= oneOn2N;

      l1loss = 0;
      for(j = 0 ; j < p ; j++)
	 for(k = 0 ; k < K ; k++)
	    l1loss += fabs(B[j + k * p]);
      loss += lambda * l1loss;
      
      if(dofusion)
      {
	 floss = 0;
	 for(j = p - 1 ; j >= 0 ; --j)
	 {
	    for(e = nE - 1 ; e >= 0 ; --e)
	    {
	       v1 = pairs[e];
	       v2 = pairs[e + nE];
	       sv = B[j + v1 * p] * C[e + v1 * nE] 
	          + B[j + v2 * p] * C[e + v2 * nE];
	       floss += sv * sv;
	    }
	 }

	 loss += gamma * 0.5 * floss;
      }

      if(fabs(loss - lossold) / fabs(lossold) < eps)
      {
	 if(verbose > 1)
	 {
	    timestamp();
	    Rprintf(" converged, lossold: %.6f loss: %.6f\n", lossold, loss);
	 }
	 allconverged++;
      }
      else
	 allconverged = 0;


       if(allconverged == 1)
      {
	 /* reset active-set to contain all (non monomorphic) coordinates, in
	  * order to check whether non-active coordinates become active again 
	  * or vice-versa */
         for(j = pK1 ; j >= 0 ; --j)
         {
            oldactive[j] = active[j];
            active[j] = !ignore[j];
         }
	 if(verbose > 1)
	 {
	    timestamp();
	    Rprintf(" resetting activeset at iter %d, loss: %.6f floss: %.6f\n",
	       iter, loss, floss);
	 }
	 //mult = 2;
      }
      else if(allconverged == 2)
      {
         for(j = pK1 ; j >= 0 ; --j)
            if(active[j] != oldactive[j])
               break;

         if(j < 0)
         {
            if(verbose > 1)
	    {
	       timestamp();
               Rprintf(" terminating at iter %d with %d active vars\n",
		  iter, numactive);
	    }
            *status = TRUE;
            break;
         }

         if(verbose > 1)
	 {
	    timestamp();
            Rprintf(" active set changed, %d active vars", numactive);
	 }

	 /* keep iterating over existing active set, keep track
	  * of the current active set */
         for(j = pK1 ; j >= 0 ; --j)
	 {
            oldactive[j] = active[j];
	    active[j] = !ignore[j]; 
	 }

	 ///* double the size of the active set */
	 //for(k = 0 ; k < K ; k++)
	 //{
	 //   mx = fminl(mult * numactiveK[k], p);
	 //   printf("%d ", mx);
	 //   for(j = mx - 1 ; j >= 0 ; --j)
	 //   {
	 //      idx = g->grad_array[j + k * p].index + k * p;
	 //      g->active[idx] = !g->ignore[idx];
	 //   }
	 //}

	 //printf("\n");

         allconverged = 0;
	 //mult *= 2;
      }     

      iter++;
      lossold = loss;
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
   free(Err);
   free(d2);
   free(ignore);
   free(oneOnLambda2PlusOne);
   /*if(CC)
      free(CC);*/
   if(diagCC)
      free(diagCC);
   /*if(BCT)
      free(BCT);*/
}

