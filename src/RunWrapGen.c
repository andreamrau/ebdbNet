/* $Id$ */
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
/* Load global subprograms */
void MatrixInv(double**, int, double**, double*);
void MatrixMult(double**, int, int, double**, int, double**);
void MatrixSum(double**, double**, double**, int*, int*);
void MatrixTrans(double**, double**, int*, int*);

/******************************************************************/
/* This is the wrapper function for to RunCode.R. Its inputs are  */
/* R = number of replicates, P = number of genes, T = number of   */
/* times, K = number of hidden states (0 means no x's estimated), */
/* x0 = initial values of x, yorig = gene expression, a0-sigma0=  */
/* initial values of alpha-sigma, conv1-conv3 = convergence       */
/* criterai, DEst and DvarEst will return the value of the        */
/* posterior mean and variance of D.                              */
/******************************************************************/

void RunWrapGen(int *R, int *P, int *T, int *K, int *M, double *xx,
    double *yy, double *uu, double *alpha, double *beta, double *gamma,
    double *delta, double *v, double *mu, double *sigma,
    double *conv1, double *conv2, double *conv3,
    double *APost, double *BPost, double *CPost, double *DPost,
    double *CvarPost, double *DvarPost, int *alliterations, int *maxiterations,
    int *subiterations, int *verboseInd)
{
    /* Load subprograms */
    void RunCode(int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*,
    	double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
        double*, double*, double*, int*, int*, int*, int*);

    if(*verboseInd == 1) {
        Rprintf("Running EBDBN algorithm ...\n");
       Rprintf("\n");
        if(*K > 0) {Rprintf("Iterations:\n");}
    }

    RunCode(R, P, T, K, M, xx, yy, uu, alpha, beta, gamma, delta, v, mu,
        sigma, conv1, conv2, conv3, APost, BPost, CPost, DPost, CvarPost, DvarPost, alliterations, maxiterations,
        subiterations, verboseInd);
}

/* Include other programs */
/******************************************************************/
/*  This is the counterpart of Overall_Posterior_Mean.R. It has   */
/*  inputs alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, & M */
/*  and returns the value of the overall posterior mean of A, B,  */
/*  C, and D (in the case where x's are estimated), or the        */
/*  overall posterior mean and variance of D (in the case where   */
/*  x's are not estimated).                                       */
/******************************************************************/

void PostMeanOverall(double *alpha, double *beta, double *gamma, double *delta,
	double *v, double ***x, double ***y, double ***u, int *K, int *P, int *T,
	int *R, int *M, double *AA, double *BB, double *CC, double *DD, double *CCvar, double *DDvar)
{
	int i, *all, j, jj, index, *KK, index2, m, mm;
    double **MNLNinv, **LNMNinv, **JGFGinv, **FGJGinv, ***HNLS,
            ***SNMH, ***EGFQ, ***QGJE, **matk1, **matm1, **Dvartemp,
            **D, **C, **B, **A, ***Dvar, ***Cvar;

	/* Load subprograms */
	void SimplifyNoX(double*, double*, double***, double***, int*, int*, int*, int*, int*, double**,
        double**);
	void SimplifyX(double*, double*, double*, double*, double*, double***, double***, double***, int*,
        int*, int*, int*, int*, int*, double**, double**, double**, double**, double***,
        double***, double***, double***);

    /* Allocate memory */
	all = (int*) calloc(1, sizeof(int));
	KK = (int*) calloc(1, sizeof(int));
	*all = 1;
	*KK = 1;
	if(*K > 0)
	{
	    *KK = *K;
	}
	A = (double**) calloc(*KK, sizeof(double*));
	B = (double**) calloc(*KK, sizeof(double*));
    C = (double**) calloc(*P, sizeof(double*));
	Cvar = (double***) calloc(*P, sizeof(double**));
	D = (double**) calloc(*P, sizeof(double*));
	Dvar = (double***) calloc(*P, sizeof(double**));
    for(j=0; j<*KK; j++)
	{
	    *(A+j) = (double*) calloc(*KK, sizeof(double));
	    *(B+j) = (double*) calloc(*M, sizeof(double));
	}
	for(i=0; i<*P; i++)
	{
        *(C+i) = (double*) calloc(*KK, sizeof(double));
	    *(D+i) = (double*) calloc(*M, sizeof(double));
	    *(Dvar+i) = (double**) calloc(*M, sizeof(double*));
	    *(Cvar+i) = (double**) calloc(*KK, sizeof(double*));
	    for(m=0; m<*M; m++)
	    {
	        *(*(Dvar+i)+m) = (double*) calloc(*M, sizeof(double));
	    }
        for(j=0; j<*KK; j++)
        {
            *(*(Cvar+i)+j) = (double*) calloc(*KK, sizeof(double));
        }
	}
    Dvartemp = (double**) calloc(*M, sizeof(double*));
    for(m=0; m<*M; m++)
    {
        *(Dvartemp+m) = (double*) calloc(*M, sizeof(double));
    }

    /* Begin function -- no x's */
	if(*K==0)
	{
        SimplifyNoX(delta, v, y, u, P, T, M, R, all, D, Dvartemp);
        for(i=0; i<*P; i++)
        {
            for(m=0; m<*M; m++)
            {
                for(mm=0; mm<*M; mm++)
                {
                   *(*(*(Dvar+i)+m)+mm) = (1.0/(*(v+i))) * (*(*(Dvartemp+m)+mm));
                }
            }
        }
	}

    /* Free memory */
    for(m=0; m<*M; m++)
    {
        free(*(Dvartemp+m));
    }
    free(Dvartemp);

	/* Begin function -- x's */
    /* Allocate MNLNinf, LNMNinf, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QFJE */
    MNLNinv = (double**) calloc(*KK, sizeof(double*));
    LNMNinv = (double**) calloc(*M, sizeof(double*));
    JGFGinv = (double**) calloc(*KK, sizeof(double*));
    FGJGinv = (double**) calloc(*M, sizeof(double*));
    HNLS = (double***) calloc(*KK, sizeof(double**));
    SNMH = (double***) calloc(*KK, sizeof(double**));
    EGFQ = (double***) calloc(*P, sizeof(double**));
    QGJE = (double***) calloc(*P, sizeof(double**));
    for(j=0; j<*KK; j++)
    {
        *(MNLNinv+j) = (double*) calloc(*KK, sizeof(double));
        *(JGFGinv+j) = (double*) calloc(*KK, sizeof(double));
        *(HNLS+j) = (double**) calloc(*KK, sizeof(double*));
        *(SNMH+j) = (double**) calloc(*M, sizeof(double*));
        for(jj=0; jj<*KK; jj++)
        {
            *(*(HNLS+j)+jj) = (double*) calloc(1, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(SNMH+j)+m) = (double*) calloc(1, sizeof(double));
        }
    }
    for(m=0; m<*M; m++)
    {
        *(LNMNinv+m) = (double*) calloc(*M, sizeof(double));
        *(FGJGinv+m) = (double*) calloc(*M, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(EGFQ+i) = (double**) calloc(*KK, sizeof(double*));
        *(QGJE+i) = (double**) calloc(*M, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(EGFQ+i)+j) = (double*) calloc(1, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(QGJE+i)+m) = (double*) calloc(1, sizeof(double));
        }
    }
    matk1 = (double**) calloc(*KK, sizeof(double*));
    matm1 = (double**) calloc(*M, sizeof(double*));
    for(j=0; j<*KK; j++)
    {
        *(matk1+j) = (double*) calloc(1, sizeof(double));
    }
    for(m=0; m<*M; m++)
    {
        *(matm1+m) = (double*) calloc(1, sizeof(double));
    }

	if(*K > 0)
	{
        SimplifyX(alpha, beta, gamma, delta, v, x, y, u, K, P, T, M, all, R,
            MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QGJE);
        for(j=0; j<*K; j++)
        {
            MatrixMult(MNLNinv, *K, *K, *(HNLS+j), 1, matk1);
            MatrixMult(LNMNinv, *M, *M, *(SNMH+j), 1, matm1);
            for(jj=0; jj<*K; jj++)
            {
                *(*(A+j)+jj) = *(*(matk1+jj));
            }
            for(m=0; m<*M; m++)
            {
                *(*(B+j)+m) = *(*(matm1+m));
            }
        }

        for(i=0; i<*P; i++)
        {
            MatrixMult(FGJGinv, *M, *M, *(QGJE+i), 1, matm1);
            MatrixMult(JGFGinv, *K, *K, *(EGFQ+i), 1, matk1);
            for(m=0; m<*M; m++)
            {
                *(*(D+i)+m) = *(*(matm1+m));
                for(mm=0; mm<*M; mm++)
                {
                    *(*(*(Dvar+i)+m)+mm) = (1.0/(*(v+i))) * (*(*(FGJGinv+m)+mm));
                }
            }
            for(j=0; j<*K; j++)
            {
               *(*(C+i)+j) =  *(*(matk1+j));
                for(jj=0; jj<*K; jj++)
                {
                    *(*(*(Cvar+i)+j)+jj) = (1.0/(*(v+i))) * (*(*(JGFGinv+j)+jj));
                }
            }
        }
	}

	/* Set CCvar = Cvar, and DDvar = Dvar */
	index = 0;
	index2 = 0;
	for(i=0; i<*P; i++)
	{
	    for(m=0; m<*M; m++)
	    {
	        for(mm=0; mm<*M; mm++)
	        {
	           *(DDvar+index) = *(*(*(Dvar+i)+m)+mm);
	           index++;
	        }
	    }
	    for(j=0; j<*K; j++)
	    {
            for(jj=0; jj<*K; jj++)
            {
                *(CCvar+index2) = *(*(*(Cvar+i)+j)+jj);
                index2++;
            }
	    }
	}

	/* Set CC = C, and DD = D */
	index = 0;
	index2 = 0;
	for(i=0; i<*P; i++)
	{
	    for(m=0; m<*M; m++)
	    {
	        *(DD + index) = *(*(D+i)+m);
	        index++;
	    }
	    for(j=0; j<*K; j++)
	    {
	        *(CC + index2) = *(*(C+i)+j);
		index2++;
	    }
	}

    /* Set AA = A, BB = B */
	index = 0;
    index2 = 0;
    for(j=0; j<*K; j++)
    {
        for(jj=0; jj<*K; jj++)
        {
            *(AA + index) = *(*(A+j)+jj);
            index++;
        }
        for(m=0; m<*M; m++)
        {
            *(BB + index2) = *(*(B+j)+m);
            index2++;
        }
    }

	/* Free memory */
    for(j=0; j<*KK; j++)
    {
        for(jj=0; jj<*KK; jj++)
        {
            free(*(*(HNLS+j)+jj));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(SNMH+j)+m));
        }
        free(*(HNLS+j));
        free(*(SNMH+j));
        free(*(MNLNinv+j));
        free(*(JGFGinv+j));
        free(*(matk1+j));
        free(*(A+j));
	    free(*(B+j));
    }
    for(i=0; i<*P; i++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(EGFQ+i)+j));
            free(*(*(Cvar+i)+j));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(QGJE+i)+m));
            free(*(*(Dvar+i)+m));
        }
        free(*(EGFQ+i));
        free(*(QGJE+i));
	    free(*(Cvar+i));
	    free(*(Dvar+i));
	    free(*(D+i));
        free(*(C+i));
    }
    for(m=0; m<*M; m++)
    {
        free(*(LNMNinv+m));
        free(*(FGJGinv+m));
        free(*(matm1+m));
    }
    free(HNLS);
    free(SNMH);
    free(MNLNinv);
    free(JGFGinv);
    free(EGFQ);
    free(QGJE);
    free(LNMNinv);
    free(FGJGinv);
    free(matk1);
    free(matm1);
	free(Cvar);
	free(Dvar);
	free(D);
    free(C);
	free(KK);
	free(all);
    free(A);
    free(B);

}

/******************************************************************/
/* This is counterpart to RunCode.R.  It takes as inputs          */
/* R = number of replicates, P = number of genes, T = number of   */
/* times, K = number of hidden states (0 means no x's estimated), */
/* M = dimension of input variable,                               */
/* xx = initial values of x, yy = gene expression, a-sigma=       */
/* initial values of alpha-sigma, conv1-conv3 = convergence       */
/* criteria, DPost and DvarPost will return the value of the      */
/* posterior mean and variance of D.                              */
/* ****************************************************************/
void RunCode(int *R, int *P, int *T, int *K, int *M, double *xx,
    double *yy, double *uu, double *alpha, double *beta, double *gamma,
    double *delta, double *v, double *mu, double *sigma,
    double *conv1, double *conv2, double *conv3,
    double *APost, double *BPost, double *CPost, double *DPost, double *CvarPost,
    double *DvarPost, int *alliterations, int *maxiterations, int *subiterations,
    int *verboseInd)
{
    int i, j, r, t, index, *KK, m;
    double ***x, ***y, ***u;

    /* Load subprograms */
    void FullAlgorithm(double***, double***, double***, double*, double*, double*, double*, double*,
        double*, double*, int*, int*, int*, int*, int*, double*, double*, double*, int*, int*, int*, int*);
    void PostMeanOverall(double*, double*, double*, double*, double*, double***, double***, double***,
        int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*);

    /* Allocate y, read in data */
    y = (double***) calloc(*R, sizeof(double**));
    KK = (int*) calloc(1, sizeof(int));
    index = 0;
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }
    for(r=0; r<*R; r++)
    {
        *(y+r) = (double**) calloc(*P, sizeof(double*));
        for(i=0; i<*P; i++)
        {
            *(*(y+r)+i) = (double*) calloc(*T, sizeof(double));
            for(t=0; t<*T; t++)
            {
                *(*(*(y+r)+i)+t) = *(yy+index);
                index++;
            }
        }
    }

    /* Allocate x's, read in data */
    x = (double***) calloc(*R, sizeof(double**));
    for(r=0; r<*R; r++)
    {
        *(x+r) = (double**) calloc(*KK, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(x+r)+j) = (double*) calloc(*T, sizeof(double));
        }
    }
    index = 0;
    if(*K > 0)
    {
        for(r=0; r<*R; r++)
        {
            for(j=0; j<*K; j++)
            {
                for(t=0; t<*T; t++)
                {
                    *(*(*(x+r)+j)+t) = *(xx+index);
                    index++;
                }
            }
        }
    }

    /* Allocate u's, read in data */
    u = (double***) calloc(*R, sizeof(double**));
    index = 0;
    for(r=0; r<*R; r++)
    {
        *(u+r) = (double**) calloc(*M, sizeof(double*));
        for(m=0; m<*M; m++)
        {
            *(*(u+r)+m) = (double*) calloc(*T, sizeof(double));
            for(t=0; t<*T; t++)
            {
                *(*(*(u+r)+m)+t) = *(uu+index);
                index++;
            }
        }
    }

    /*****************************/
    /*  Run Full algorithm       */
    /*****************************/

    FullAlgorithm(y, x, u, alpha, beta, gamma, delta, v, mu, sigma,
        K, P, R, T, M, conv1, conv2, conv3, alliterations, maxiterations, subiterations,
        verboseInd);

    /****************************************/
    /*  Find posterior mean & variance      */
    /****************************************/

    if(*verboseInd == 1)
    {
        Rprintf("EBDBN Algorithm complete! \n");
    }

    PostMeanOverall(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M, APost, BPost, CPost, DPost, CvarPost, DvarPost);

    /* Read in final estimate of x's */
    if(*K > 0)
    {
        index = 0;
        for(r=0; r<*R; r++)
        {
            for(j=0; j<*KK; j++)
            {
                for(t=0; t<*T; t++)
                {
                    *(xx+index) = *(*(*(x+r)+j)+t);
                    index++;
                }
            }
        }
    }

    /* Release memory of x's */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(x+r)+j));
        }
        free(*(x+r));
    }
    free(x);

    /* Release memory of y's */
    for(r=0; r<*R; r++)
    {
        for(i=0; i<*P; i++)
        {
            free(*(*(y+r)+i));
        }
        free(*(y+r));
    }
    free(y);

    /* Release memory of u's */
    for(r=0; r<*R; r++)
    {
        for(m=0; m<*M; m++)
        {
            free(*(*(u+r)+m));
        }
        free(*(u+r));
    }
    free(u);
    free(KK);
}

/******************************************************************/
/* FullAlgorithm implements the EBDBN algorithm.  It has as       */
/* arguments y = gene expression, x = hidden states, u = inputs,  */
/* alpha - sigma = hyperparameter initial values, K = hidden state*/
/* dimension, P = # of genes, R = # of replicates, T = time pts,  */
/* M = dimension of inputs, conv1 - conv3 = convergence criteria, */
/* alliterations = count of total iterations, maxiterations =     */
/* maximum number of allowed overall iterations, subiterations =  */
/* maximum number of allowed iterations in EM-type algorithm      */
/* ****************************************************************/
void FullAlgorithm(double ***y, double ***x, double ***u, double *alpha,
    double *beta, double *gamma, double *delta, double *v,
    double *mu, double *sigma, int *K, int *P, int *R, int *T, int *M,
    double *conv1, double *conv2, double *conv3, int *alliterations, int *maxiterations,
    int *subiterations, int *verboseInd)
{
    int iter=1, j, i, r, *KK, m;
    double overalldiff, ***A, ***B, ***C, ***D, ***Dvar,
        *alpha0, *beta0, *gamma0, *delta0, *v0,
        *alphaold, *betaold, *gammaold, *deltaold, *vold;
    double alphadiff, betadiff, gammadiff, deltadiff, sumnum, sumden;

    /* Load subprograms */
    void EmTypeConv(double*, double*, double*, double*, double*, double***, double***, double***,
        int*, int*, int*, int*, int*, double*, double*, int*);
    void PostMeanR(double*, double*, double*, double*, double*, double***, double***, double***, int*,
        int*, int*, int*, int*, double***, double***, double***, double***, double***);
    void Kalman(double***, double***, double***, double***, double***, double*, double*, double*, int*,
        int*, int*, int*, int*, double***, double***);

    /* Allocate memory */
    KK = (int*) calloc(1, sizeof(int));
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }
    A = (double***) calloc(*R, sizeof(double**));
    B = (double***) calloc(*R, sizeof(double**));
    C = (double***) calloc(*R, sizeof(double**));
    D = (double***) calloc(*R, sizeof(double**));
    Dvar = (double***) calloc(*R, sizeof(double**));
    for(r=0; r<*R; r++)
    {
        *(A+r) = (double**) calloc(*KK, sizeof(double*));
        *(B+r) = (double**) calloc(*KK, sizeof(double*));
        *(C+r) = (double**) calloc(*P, sizeof(double*));
        *(D+r) = (double**) calloc(*P, sizeof(double*));
        *(Dvar+r) = (double**) calloc(*M, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(A+r)+j) = (double*) calloc(*KK, sizeof(double));
            *(*(B+r)+j) = (double*) calloc(*M, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(C+r)+i) = (double*) calloc(*KK, sizeof(double));
            *(*(D+r)+i) = (double*) calloc(*M, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(Dvar+r)+m) = (double*) calloc(*M, sizeof(double));
        }
    }
    alpha0 = (double*) calloc(*KK, sizeof(double));
    beta0 = (double*) calloc(*M, sizeof(double));
    gamma0 = (double*) calloc(*KK, sizeof(double));
    delta0 = (double*) calloc(*M, sizeof(double));
    v0 = (double*) calloc(*P, sizeof(double));
    alphaold = (double*) calloc(*KK, sizeof(double));
    betaold = (double*) calloc(*M, sizeof(double));
    gammaold = (double*) calloc(*KK, sizeof(double));
    deltaold = (double*) calloc(*M, sizeof(double));
    vold = (double*) calloc(*P, sizeof(double));

    /* If no x's */
    if(*K == 0)
    {
        EmTypeConv(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M, conv1, conv2, subiterations);
    }

    /* If x estimates */
    if(*K > 0)
    {
        /* Three convergence criteria:
           Convergence for first sub-loop of EM (pre-v) is conv1
           Convergence for second sub-loop of EM (post-v) is conv2
           Overall algorithm convergence is conv3 */
        /*********************************************************/
        /*  Initial hyperparameter estimates, based on initial x */
        /*  and initial hyperparameters                          */
        /*  New values are stored back in alpha, beta, etc.      */
        /*********************************************************/
        for(j=0; j<*K; j++)
        {
            *(alpha0+j) = *(alpha+j);
            *(gamma0+j) = *(gamma+j);
        }
        for(m=0; m<*M; m++)
        {
            *(beta0+m) = *(beta+m);
            *(delta0+m) = *(delta+m);
        }
        for(i=0; i<*P; i++)
        {
            *(v0+i) = *(v+i);
        }
        EmTypeConv(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M,
            conv1, conv2, subiterations);

        for(j=0; j<*K; j++)
        {
            *(alphaold+j) = *(alpha+j);
            *(gammaold+j) = *(gamma+j);
        }
        for(m=0; m<*M; m++)
        {
            *(betaold+m) = *(beta+m);
            *(deltaold+m) = *(delta+m);
        }
        for(i=0; i<*P; i++)
        {
            *(vold+i) = *(v+i);
        }
        overalldiff = 100.0;

        /**************************************************************************************/
        /* Stop if we reach convergence criterion OR if we get up to 100 overall iterations   */
        /* If we are over 100 iterations, we will start over with a new set of initial values */
        /**************************************************************************************/
        while(overalldiff > *conv3)
        {
            if(*alliterations > *maxiterations) break;
            PostMeanR(alpha, beta, gamma, delta, v, x, y, u, K, P, T,
                R, M, A, B, C, D, Dvar);

            /******************************************************/
            /*  Kalman filter and smoother estimates of x         */
            /******************************************************/
            Kalman(y, A, B, C, D, v, mu, sigma, K, P, T, R, M, x, u);

            /*  Normally we would re-estimate mu and sigma here */

            /******************************************************/
            /*  Running EM algorithm based on new x values        */
            /*  Start with original initial values                */
            /******************************************************/

            for(j=0; j<*K; j++)
            {
                *(alpha+j) = *(alpha0+j);
                *(gamma+j) = *(gamma0+j);
            }
            for(m=0; m<*M; m++)
            {
                *(beta+m) = *(beta0+m);
                *(delta+m) = *(delta0+m);
            }
            for(i=0; i<*P; i++)
            {
                *(v+i) = *(v0+i);
            }
            EmTypeConv(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M,
                conv1, conv2, subiterations);

            /******************************************************/
            /*  Check convergence of hyperparameters              */
            /******************************************************/
            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(alphaold+j)) * (*(alphaold+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(alpha+j) - *(alphaold+j))*(*(alpha+j) - *(alphaold+j)))/sumden;
            }
            alphadiff = sqrt(sumnum);                                       /* alpha.diff */

            sumnum = 0;
            sumden = 0;
            for(m=0; m<*M; m++)
            {
                sumden += (*(betaold+m)) * (*(betaold+m));
            }
            for(m=0; m<*M; m++)
            {
                sumnum += ((*(beta+m) - *(betaold+m))*(*(beta+m) - *(betaold+m)))/sumden;
            }
            betadiff = sqrt(sumnum);                                       /* beta.diff */

            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(gammaold+j)) * (*(gammaold+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(gamma+j) - *(gammaold+j))*(*(gamma+j) - *(gammaold+j)))/sumden;
            }
            gammadiff = sqrt(sumnum);                                       /* gamma.diff */

            sumnum = 0;
            sumden = 0;
            for(m=0; m<*M; m++)
            {
                sumden += (*(deltaold+m)) * (*(deltaold+m));
            }
            for(m=0; m<*M; m++)
            {
                sumnum += ((*(delta+m) - *(deltaold+m))*(*(delta+m) - *(deltaold+m)))/sumden;
            }
            deltadiff = sqrt(sumnum);                                       /* delta.diff */

            /* Find maximum convergence criterion */
            overalldiff = alphadiff;
            if(betadiff>overalldiff) {overalldiff = betadiff;}
            if(gammadiff>overalldiff) {overalldiff = gammadiff;}
            if(deltadiff>overalldiff) {overalldiff = deltadiff;}

            if(*verboseInd == 1) {Rprintf("Max difference = %f, ", overalldiff);}
            /* Set new values = to old values of hyperparameters */
            for(j=0; j<*K; j++)
            {
                *(alphaold+j) = *(alpha+j);
                *(gammaold+j) = *(gamma+j);
            }
            for(m=0; m<*M; m++)
            {
                *(betaold+m) = *(beta+m);
                *(deltaold+m) = *(delta+m);
            }

            /* Update the number of iterations, and loop */
           if(*verboseInd == 1) {Rprintf("** Iteration %d complete! ** \n", iter);}
           *alliterations = iter;
            iter++;
        }
    }

    /* Free memory */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(A+r)+j));
            free(*(*(B+r)+j));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(C+r)+i));
            free(*(*(D+r)+i));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(Dvar+r)+m));
        }
        free(*(A+r));
        free(*(B+r));
        free(*(C+r));
        free(*(D+r));
        free(*(Dvar+r));
    }
    free(A);
    free(B);
    free(C);
    free(D);
    free(Dvar);
    free(alpha0);
    free(beta0);
    free(gamma0);
    free(delta0);
    free(v0);
    free(alphaold);
    free(betaold);
    free(gammaold);
    free(deltaold);
    free(vold);
    free(KK);
}


/******************************************************************/
/*  This is the counterpart of Kalman_Filter.R.  It takes as      */
/*  inputs y, the current values of A, B, C, D, v, mu, and sigma  */
/*  (where sigma is a matrix), K, P, T, and R.  It returns the    */
/*  new values of x into the original triple pointer x.           */
/******************************************************************/
void Kalman(double ***y, double ***A, double ***B,
    double ***C, double ***D, double *v, double *mu, double *sigma,
    int *K, int *P, int *T, int *R, int *M, double ***x, double ***u)
{
    int r, t, i, j, jj, index, m;
    double **yr, **Ar, **Br, **Cr, **Dr, **sigmamat, **xminus,
        **filter, **Pminus, **Pk, **smoother, **Ps, **ur;

    /* Load subprograms */
    void KalmanFilter(double**, double**, double**, double**, double**, double**, double*, double*,
        double**, int*, int*, int*, int*, double**, double**, double**, double**);
    void KalmanSmoother(double**, double**, double**, double**,
        double**, int*, int*, double**, double**);

    /* Allocate memory */
    ur = (double**) calloc(*M, sizeof(double*));
    yr = (double**) calloc(*P, sizeof(double*));
    Ar = (double**) calloc(*K, sizeof(double*));
    Br = (double**) calloc(*K, sizeof(double*));
    Cr = (double**) calloc(*P, sizeof(double*));
    Dr = (double**) calloc(*P, sizeof(double*));
    sigmamat = (double**) calloc(*K, sizeof(double*));
    xminus = (double**) calloc(*K, sizeof(double*));
    filter = (double**) calloc(*K, sizeof(double*));
    Pminus = (double**) calloc(*K, sizeof(double*));
    Pk = (double**) calloc(*K, sizeof(double*));
    smoother = (double**) calloc(*K, sizeof(double*));
    Ps = (double**) calloc(*K, sizeof(double*));
    for(j=0; j<*K; j++)
    {
        *(Ar+j) = (double*) calloc(*K, sizeof(double));
        *(Br+j) = (double*) calloc(*M, sizeof(double));
        *(sigmamat+j) = (double*) calloc(*K, sizeof(double));
        *(xminus+j) = (double*) calloc(*T, sizeof(double));
        *(filter+j) = (double*) calloc(*T, sizeof(double));
        *(Pminus+j) = (double*) calloc(*K, sizeof(double));
        *(Pk+j) = (double*) calloc(*K, sizeof(double));
        *(smoother+j) = (double*) calloc(*T, sizeof(double));
        *(Ps+j) = (double*) calloc(*K, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(yr+i) = (double*) calloc(*T, sizeof(double));
        *(Cr+i) = (double*) calloc(*K, sizeof(double));
        *(Dr+i) = (double*) calloc(*M, sizeof(double));
    }
    for(m=0; m<*M; m++)
    {
        *(ur+m) = (double*) calloc(*T, sizeof(double));
    }

    /* Change sigma into matrix form */
    index = 0;
    for(j=0; j<*K; j++)
    {
        for(jj=0; jj<*K; jj++)
        {
            *(*(sigmamat+j)+jj) = *(sigma+index);
            index++;
        }
    }

    /* Do Kalman filter for each replicate r independently */
    /* Set xminus, filter, Pminus, Pk, smoother, and PS equal to zero at start */
    for(r=0; r<*R; r++)
    {
        /* Read in correct values for y, A, B, C, and D */
        for(i=0; i<*P; i++)
        {
            for(t=0; t<*T; t++)
            {
                *(*(yr+i)+t) = *(*(*(y+r)+i)+t);
            }
            for(j=0; j<*K; j++)
            {
                *(*(Cr+i)+j) = *(*(*(C+r)+i)+j);
            }
            for(m=0; m<*M; m++)
            {
                *(*(Dr+i)+m) = *(*(*(D+r)+i)+m);
            }
        }
        for(j=0; j<*K; j++)
        {
            for(jj=0; jj<*K; jj++)
            {
                *(*(Ar+j)+jj) = *(*(*(A+r)+j)+jj);
            }
            for(m=0; m<*M; m++)
            {
                *(*(Br+j)+m) = *(*(*(B+r)+j)+m);
            }
        }
        for(m=0; m<*M; m++)
        {
            for(t=0; t<*T; t++)
            {
                *(*(ur+m)+t) = *(*(*(u+r)+m)+t);
            }
        }
        KalmanFilter(yr, ur, Ar, Br, Cr, Dr, v, mu, sigmamat, K, P,
            T, M, xminus, filter, Pminus, Pk);                                      /* Kalman filter */
        KalmanSmoother(Ar, xminus, filter, Pminus, Pk, K, T,
            smoother, Ps);                                                         /* Kalman smoother */
        for(j=0; j<*K; j++)
        {
            for(t=0; t<*T; t++)
            {
                *(*(*(x+r)+j)+t) = *(*(smoother+j)+t);                         /* Set x's to smoothed value */
           }
        }
    }
     /* Release memory */
    for(j=0; j<*K; j++)
    {
        free(*(Ar+j));
        free(*(Br+j));
        free(*(sigmamat+j));
        free(*(xminus+j));
        free(*(filter+j));
        free(*(Pminus+j));
        free(*(Pk+j));
        free(*(smoother+j));
        free(*(Ps+j));
    }
    for(i=0; i<*P; i++)
    {
        free(*(yr+i));
        free(*(Cr+i));
        free(*(Dr+i));
    }
    for(m=0; m<*M; m++)
    {
        free(*(ur+m));
    }
    free(ur);
    free(yr);
    free(Ar);
    free(Br);
    free(Cr);
    free(Dr);
    free(sigmamat);
    free(xminus);
    free(filter);
    free(Pminus);
    free(Pk);
    free(smoother);
    free(Ps);
}

/* Kalman Filter function */
void KalmanFilter(double **yr, double **ur, double **Ar, double **Br, double **Cr,
    double **Dr, double *v, double *mu, double **sigmamat, int *K, int *P,
    int *T, int *M, double **xminus, double **filter, double **Pminus, double **Pk)
{
    int t, j, jj, i, m;
    double **gain, **matp1, **matk1, **matk1b, **matp1b, **matkk,
        **xtemp, **ytemp, **xtemp2, **ytemp2, **Art, **matkk2, **utemp, **utemp2;

    /* Load subprograms */
    void KalmanGain(double**, double**, double*, int*, int*, double**);

    /* Allocate memory */
    gain = (double**) calloc(*K, sizeof(double*));
    matp1 = (double**) calloc(*P, sizeof(double*));
    matp1b = (double**) calloc(*P, sizeof(double*));
    matk1 = (double**) calloc(*K, sizeof(double*));
    matk1b = (double**) calloc(*K, sizeof(double*));
    matkk = (double**) calloc(*K, sizeof(double*));
    xtemp = (double**) calloc(*K, sizeof(double*));
    ytemp = (double**) calloc(*P, sizeof(double*));
    xtemp2 = (double**) calloc(*K, sizeof(double*));
    ytemp2 = (double**) calloc(*P, sizeof(double*));
    Art = (double**) calloc(*K, sizeof(double*));
    matkk2 = (double**) calloc(*K, sizeof(double*));
    utemp = (double**) calloc(*M, sizeof(double*));
    utemp2 = (double**) calloc(*M, sizeof(double*));
    for(j=0; j<*K; j++)
    {
        *(gain+j) = (double*) calloc(*P, sizeof(double));
        *(matk1+j) = (double*) calloc(1, sizeof(double));
        *(matk1b+j) = (double*) calloc(1, sizeof(double));
        *(matkk+j) = (double*) calloc(*K, sizeof(double));
        *(xtemp+j) = (double*) calloc(1, sizeof(double));
        *(xtemp2+j) = (double*) calloc(1, sizeof(double));
        *(Art+j) = (double*) calloc(*K, sizeof(double));
        *(matkk2+j) = (double*) calloc(*K, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(matp1+i) = (double*) calloc(1, sizeof(double));
        *(matp1b+i) = (double*) calloc(1, sizeof(double));
        *(ytemp+i) = (double*) calloc(1, sizeof(double));
        *(ytemp2+i) = (double*) calloc(1, sizeof(double));
    }
    for(m=0; m<*M; m++)
    {
        *(utemp+m) = (double*) calloc(1, sizeof(double));
        *(utemp2+m) = (double*) calloc(1, sizeof(double));
    }

    /* Begin Kalman filter */
    for(t=0; t<*T; t++)
    {
        if(t==0)                                                            /* First time point */
        {
            for(j=0; j<*K; j++)
            {
                *(*(xminus+j)) = *(mu+j);                                   /* x.minus */
                for(jj=0; jj<*K; jj++)
                {
                    *(*(Pminus+j)+jj) = *(*(sigmamat+j)+jj);                /* P.minus */
                }
            }
        }
        if(t > 0)                                                           /* All other time points */
       {
            for(j=0; j<*K; j++)
            {
                *(*(xtemp+j)) = *(*(filter+j)+(t-1));
            }
            for(m=0; m<*M; m++)
            {
                *(*(utemp+m)) = *(*(ur+m)+t);
            }
            MatrixMult(Ar, *K, *K, xtemp, 1, matk1);
            MatrixMult(Br, *K, *M, utemp, 1, matk1b);
            for(j=0; j<*K; j++)
            {
                *(*(xminus+j)+t) = *(*(matk1+j)) + *(*(matk1b+j));      /* x.minus */
            }
            MatrixMult(Ar, *K, *K, Pk, *K, matkk);
            MatrixTrans(Ar, Art, K, K);
            MatrixMult(matkk, *K, *K, Art, *K, Pminus);
            for(j=0; j<*K; j++)
            {
                *(*(Pminus+j)+j) += 1;                                  /* P.minus */
           }
        }
        /* Set y */
        for(i=0; i<*P; i++)
        {
            *(*(ytemp+i)) = *(*(yr+i)+t);
            *(*(ytemp2+i)) = 0;
        }

        /* Set gain = 0 to set up for KalmanGain */
        for(j=0; j<*K; j++)
        {
            for(i=0; i<*P; i++)
            {
                *(*(gain+j)+i) = 0;
            }
        }
        KalmanGain(Pminus, Cr, v, K, P, gain);                         /* Kalman Gain */
        for(j=0; j<*K; j++)
        {
            *(*(xtemp2+j)) = *(*(xminus+j)+t);
        }
        for(m=0; m<*M; m++)
        {
            *(*(utemp2+m)) = *(*(ur+m)+t);
        }
        MatrixMult(Cr, *P, *K, xtemp2, 1, matp1);
        MatrixMult(Dr, *P, *M, utemp2, 1, matp1b);
        for(i=0; i<*P; i++)
        {
            *(*(ytemp2+i)) = *(*(ytemp+i)) - *(*(matp1+i)) - *(*(matp1b+i));
        }
        MatrixMult(gain, *K, *P, ytemp2, 1, matk1);
        for(j=0; j<*K; j++)
        {
            *(*(filter+j)+t) = *(*(xminus+j)+t) + *(*(matk1+j));        /* Filter (x.k) */
        }
        MatrixMult(gain, *K, *P, Cr, *K, matkk);
        for(j=0; j<*K; j++)
        {
            for(jj=0; jj<*K; jj++)
            {
                if(j != jj)
                {
                    *(*(matkk2+j)+jj) = 0-(*(*(matkk+j)+jj));
                }
                if(j == jj)
                {
                    *(*(matkk2+j)+jj) = 1-(*(*(matkk+j)+jj));
                }
            }
        }
        MatrixMult(matkk2, *K, *K, Pminus, *K, Pk);                  /* Pk */
    }

    /* Release memory */
    for(j=0; j<*K; j++)
    {
        free(*(gain+j));
        free(*(matk1+j));
        free(*(matk1b+j));
        free(*(matkk+j));
        free(*(xtemp+j));
        free(*(xtemp2+j));
        free(*(Art+j));
        free(*(matkk2+j));
    }
    for(i=0; i<*P; i++)
    {
        free(*(matp1+i));
        free(*(matp1b+i));
        free(*(ytemp+i));
        free(*(ytemp2+i));
    }
    for(m=0; m<*M; m++)
    {
        free(*(utemp+m));
        free(*(utemp2+m));
    }
    free(utemp);
    free(utemp2);
    free(Art);
    free(xtemp);
    free(xtemp2);
    free(ytemp);
    free(ytemp2);
    free(matkk);
    free(matk1);
    free(matk1b);
    free(matp1);
    free(matp1b);
    free(gain);
    free(matkk2);
}

/* Kalman Gain calculator */
void KalmanGain(double **Pminus, double **Cr, double *v, int *K, int *P, double **gain)
{
    double **matpp, **matpk, **matkp, **Crt, **inv, *det;
    int j, i;

    /* Allocate memory */
    matpp = (double**) calloc(*P, sizeof(double*));
    matpk = (double**) calloc(*P, sizeof(double*));
    matkp = (double**) calloc(*K, sizeof(double*));
    Crt = (double**) calloc(*K, sizeof(double*));
    inv = (double**) calloc(*P, sizeof(double*));
    det = (double*) calloc(1, sizeof(double));
    for(j=0; j<*K; j++)
    {
        *(matkp+j) = (double*) calloc(*P, sizeof(double));
        *(Crt+j) = (double*) calloc(*P, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(matpp+i) = (double*) calloc(*P, sizeof(double));
        *(matpk+i) = (double*) calloc(*P, sizeof(double));
        *(inv+i) = (double*) calloc(*P, sizeof(double));
    }

    MatrixMult(Cr, *P, *K, Pminus, *K, matpk);
    MatrixTrans(Cr, Crt, P, K);
    MatrixMult(matpk, *P, *K, Crt, *P, matpp);

    for(i=0; i<*P; i++)
    {
        *(*(matpp+i)+i) += 1/(*(v+i));
    }
    MatrixInv(matpp, *P, inv, det);
    MatrixMult(Pminus, *K, *K, Crt, *P, matkp);
    MatrixMult(matkp, *K, *P, inv, *P, gain);

    /* Free memory */
    for(i=0; i<*P; i++)
    {
        free(*(matpp+i));
        free(*(matpk+i));
        free(*(inv+i));
    }
    for(j=0; j<*K; j++)
    {
        free(*(matkp+j));
        free(*(Crt+j));
    }
    free(matkp);
    free(matpp);
    free(matpk);
    free(Crt);
    free(inv);
    free(det);
}


/* Kalman Smoother function */
void KalmanSmoother(double **Ar, double **xminus, double **filter,
    double **Pminus, double **Pk, int *K, int *T, double **smoother,
    double **Ps)
{
    int t, j, jj;
    double **smooth, **xtemp, **matk1, **Art, **matkk;

    /* Load subprograms */
    void KalmanSmooth(double**, double**, double**, int*, double**);

    /* Allocate memory */
    smooth = (double**) calloc(*K, sizeof(double*));
    xtemp = (double**) calloc(*K, sizeof(double*));
    matk1 = (double**) calloc(*K, sizeof(double*));
    Art = (double**) calloc(*K, sizeof(double*));
    matkk = (double**) calloc(*K, sizeof(double*));
    for(j=0; j<*K; j++)
    {
        *(smooth+j) = (double*) calloc(*K, sizeof(double));
        *(xtemp+j) = (double*) calloc(1, sizeof(double));
        *(matk1+j) = (double*) calloc(1, sizeof(double));
        *(Art+j) = (double*) calloc(*K, sizeof(double));
        *(matkk+j) = (double*) calloc(*K, sizeof(double));
    }

    for(t=((*T)-1); t>=0; t--)
    {
        if(t==((*T)-1))                                                    /* Last time point */
        {
            for(j=0; j<*K; j++)
            {
                *(*(smoother+j)+t) = *(*(filter+j)+t);                     /* Smoother */
                for(jj=0; jj<*K; jj++)
                {
                    *(*(Ps+j)+jj) = *(*(Pk+j)+jj);
                }
            }
        }

        if(t<((*T)-1))                                                      /* Remaining time points */
        {
            KalmanSmooth(Pminus, Pk, Ar, K, smooth);                        /* Smooth (A.s) */
            for(j=0; j<*K; j++)
            {
                *(*(xtemp+j)) = *(*(smoother+j)+(t+1)) - *(*(xminus+j)+(t+1));
            }
            MatrixMult(smooth, *K, *K, xtemp, 1, matk1);
            for(j=0; j<*K; j++)
            {
                *(*(smoother+j)+t) = *(*(filter+j)+t) + *(*(matk1+j));      /* Smoother */
            }
            MatrixTrans(Ar, Art, K, K);
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    *(*(matkk+j)+jj) = *(*(Ps+j)+jj) - *(*(Pminus+j)+jj);
                }
            }
            MatrixMult(smooth, *K, *K, matkk, *K, matkk);
            MatrixMult(matkk, *K, *K, Art, *K, matkk);
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    *(*(Ps+j)+jj) = *(*(Pk+j)+jj) + *(*(matkk+j)+jj);          /* Pa */
                }
            }
        }
    }

   /* Release memory */
    for(j=0; j<*K; j++)
    {
       free(*(smooth+j));
       free(*(xtemp+j));
       free(*(matk1+j));
       free(*(Art+j));
       free(*(matkk+j));
    }
    free(smooth);
    free(xtemp);
    free(matk1);
    free(Art);
    free(matkk);
}

/* Kalman Smooth calculator */
void KalmanSmooth(double **Pminus, double **Pk, double **Ar, int *K, double **smooth)
{
    int j;
    double **Art, **matkk, **inv, *det;

    /* Allocate memory */
    Art = (double**) calloc(*K, sizeof(double*));
    matkk = (double**) calloc(*K, sizeof(double*));
    inv = (double**) calloc(*K, sizeof(double*));
    det = (double*) calloc(1, sizeof(double));
    for(j=0; j<*K; j++)
    {
        *(Art+j) = (double*) calloc(*K, sizeof(double));
        *(matkk+j) = (double*) calloc(*K, sizeof(double));
        *(inv+j) = (double*) calloc(*K, sizeof(double));
    }

    MatrixTrans(Ar, Art, K, K);
    MatrixInv(Pminus, *K, inv, det);
    MatrixMult(Pk, *K, *K, Art, *K, matkk);
    MatrixMult(matkk, *K, *K, inv, *K, smooth);

    /* Release memory */
   for(j=0; j<*K; j++)
    {
        free(*(matkk+j));
        free(*(Art+j));
        free(*(inv+j));
    }
    free(Art);
    free(matkk);
    free(inv);
    free(det);
}


/******************************************************************/
/* EmTypeConv runs the EM-type alogrithm.  Inputs are alpha - v = */
/* hyperparameters, x = hidden states, y = gene expression, u =   */
/* inputs, K = dimension of hidden state, P = # of genes, T = time*/
/* pts, R = # of replicates, M = input dimension, conv1 - conv2 = */
/* convergence criteria, subiterations = maximum allowed          */
/* iterations for EM-type loop.                                   */
/* ****************************************************************/
void EmTypeConv(double *alpha, double *beta, double *gamma,
    double *delta, double *v, double ***x, double ***y, double ***u,
    int *K, int *P, int *T, int *R, int *M, double *conv1, double *conv2,
    int *subiterations)
{
    double ***A, ***B, ***C, ***D, ***Dvar;
    int r, i, j, *KK, m;

    /* Load subprograms */
    void HyperMax(double*, double*, double*, double*, double*, double***, double***, double***, int*,
        int*, int*, int*, int*, double*, int*);
    void PostMeanR(double*, double*, double*, double*, double*, double***, double***, double***, int*,
        int*, int*, int*, int*, double***, double***, double***, double***, double***);
    void VarMaxR(double***, double***, double***, double***, double***, int*, int*, int*, int*, int*,
        double*);

    /* Allocate space for vnew, A, B, C, D, and Dvar */
    KK = (int*) calloc(1, sizeof(int));
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }
    A = (double***) calloc(*R, sizeof(double**));
    B = (double***) calloc(*R, sizeof(double**));
    C = (double***) calloc(*R, sizeof(double**));
    D = (double***) calloc(*R, sizeof(double**));
    Dvar = (double***) calloc(*R, sizeof(double**));
    for(r=0; r<*R; r++)
    {
        *(A+r) = (double**) calloc(*KK, sizeof(double*));
        *(B+r) = (double**) calloc(*KK, sizeof(double*));
        *(C+r) = (double**) calloc(*P, sizeof(double*));
        *(D+r) = (double**) calloc(*P, sizeof(double*));
        *(Dvar+r) = (double**) calloc(*M, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(A+r)+j) = (double*) calloc(*KK, sizeof(double));
            *(*(B+r)+j) = (double*) calloc(*M, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(C+r)+i) = (double*) calloc(*KK, sizeof(double));
            *(*(D+r)+i) = (double*) calloc(*M, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(Dvar+r)+m) = (double*) calloc(*M, sizeof(double));
        }
    }

    /* Initial estimation of hyperparameter estimates
       New estimates are put back into alpha, beta, gamma, delta
       Use convergence criterion conv1 */
    HyperMax(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M, conv1, subiterations);

    /* Estimation of gene precisions, v
       New estimate of v is put back into v */
    PostMeanR(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M, A, B, C,
        D, Dvar);
    VarMaxR(x, y, u, C, D, P, R, T, K, M, v);

    /* Final estimation of hyperparameter estimates
       New esimates are put back into alpha, beta, gamma, delta
       Use convergence criterion conv2 */
    HyperMax(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M, conv2, subiterations);

    /* Release memory */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*K; j++)
        {
            free(*(*(A+r)+j));
            free(*(*(B+r)+j));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(C+r)+i));
            free(*(*(D+r)+i));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(Dvar+r)+m));
        }
        free(*(A+r));
        free(*(B+r));
        free(*(C+r));
        free(*(D+r));
        free(*(Dvar+r));
    }
    free(A);
    free(B);
    free(C);
    free(D);
    free(Dvar);
    free(KK);
}

/******************************************************************/
/*  Function to conduct EM algorithm to estimate alpha - delta.   */
/*  Corresponds to the hyper.max function in R.  Inputs are       */
/*  initial values of alpha - delta, current value of v and x,    */
/*  y, u, K, P, T, R, M, and a convergence criterion. Updated     */
/*  values of alpha - delta are returned in the original pointers.*/
/******************************************************************/
void HyperMax(double *alpha, double *beta, double *gamma, double *delta,
    double *v, double ***x, double ***y, double ***u, int *K, int *P, int *T,
    int *R, int *M, double *conv, int *subiterations)
{
    int r, j, i, jj, t, iter, *all, *Rchoice, *KK, m, mm;
    double convergence, K2, P2;
    double **MNLNinv, **LNMNinv, **JGFGinv, **FGJGinv,
        ***HNLS, ***SNMH, ***EGFQ, ***QGJE, **matk1, **matm1, **matp1,
        *tempa, *tempb, *tempc, *tempd,
        **alphanewmat, **betanewmat, **gammanewmat, **deltanewmat,
        *alphanew, *betanew, *gammanew, *deltanew,
        alphadiff, betadiff, gammadiff, deltadiff, sumnum, sumden,
        ***A, ***B, ***C, ***D, ***Dvar, *temp, *temp2;

    /* Load subprograms */
    double VecMedian(double*, int*);
    void PostMeanR(double*, double*, double*, double*, double*, double***, double***, double***, int*,
        int*, int*, int*, int*, double***, double***, double***, double***, double***);
    void SimplifyX(double*, double*, double*, double*, double*, double***, double***, double***, int*,
        int*, int*, int*, int*, int*, double**, double**, double**, double**, double***,
        double***, double***, double***);

    all = (int*) calloc(1, sizeof(int));
    Rchoice = (int*) calloc(1, sizeof(int));
    KK = (int*) calloc(1, sizeof(int));
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }
    *all = 0;
    iter = 1;
    convergence = *conv + 1.0;
    K2 = (double) (*K)/2.0;
    P2 = (double) (*P)/2.0;

    /* Allocate MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, and QGJE and set to 0 */
    MNLNinv = (double**) calloc(*KK, sizeof(double*));
    LNMNinv = (double**) calloc(*M, sizeof(double*));
    JGFGinv = (double**) calloc(*KK, sizeof(double*));
    FGJGinv = (double**) calloc(*M, sizeof(double*));
    HNLS = (double***) calloc(*KK, sizeof(double**));
    SNMH = (double***) calloc(*KK, sizeof(double**));
    EGFQ = (double***) calloc(*P, sizeof(double**));
    QGJE = (double***) calloc(*P, sizeof(double**));

    matk1 = (double**) calloc(*KK, sizeof(double*));
    matm1 = (double**) calloc(*M, sizeof(double*));
    matp1 = (double**) calloc(*P, sizeof(double*));
    tempa = (double*) calloc(*KK, sizeof(double));
    tempb = (double*) calloc(*M, sizeof(double));
    tempc = (double*) calloc(*KK, sizeof(double));
    tempd = (double*) calloc(*M, sizeof(double));

    alphanewmat = (double**) calloc(*KK, sizeof(double*));
    betanewmat = (double**) calloc(*M, sizeof(double*));
    gammanewmat = (double**) calloc(*KK, sizeof(double*));
    deltanewmat = (double**) calloc(*M, sizeof(double*));
    alphanew = (double*) calloc(*KK, sizeof(double));
    betanew = (double*) calloc(*M, sizeof(double));
    gammanew = (double*) calloc(*KK, sizeof(double));
    deltanew = (double*) calloc(*M, sizeof(double));

    for(j=0; j<*KK; j++)
    {
        *(MNLNinv+j) = (double*) calloc(*KK, sizeof(double));
        *(JGFGinv+j) = (double*) calloc(*KK, sizeof(double));
        *(HNLS+j) = (double**) calloc(*KK, sizeof(double*));
        *(SNMH+j) = (double**) calloc(*M, sizeof(double*));
        for(jj=0; jj<*KK; jj++)
        {
            *(*(HNLS+j)+jj) = (double*) calloc(1, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(SNMH+j)+m) = (double*) calloc(1, sizeof(double));
        }
        *(matk1+j) = (double*) calloc(1, sizeof(double));
        *(alphanewmat+j) = (double*) calloc(*R, sizeof(double));
        *(gammanewmat+j) = (double*) calloc(*R, sizeof(double));
    }
    for(m=0; m<*M; m++)
    {
        *(LNMNinv+m) = (double*) calloc(*M, sizeof(double));
        *(FGJGinv+m) = (double*) calloc(*M, sizeof(double));
        *(betanewmat+m) = (double*) calloc(*R, sizeof(double));
        *(deltanewmat+m) = (double*) calloc(*R, sizeof(double));
        *(matm1+m) = (double*) calloc(1, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(EGFQ+i) = (double**) calloc(*KK, sizeof(double*));
        *(QGJE+i) = (double**) calloc(*M, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(EGFQ+i)+j) = (double*) calloc(1, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(QGJE+i)+m) = (double*) calloc(1, sizeof(double));
        }
        *(matp1+i) = (double*) calloc(1, sizeof(double));
    }

    /* If estimating x's */
    if(*K > 0)
    {
        while(convergence > *conv)
       {
            /* Rprintf("Sub convergence = %f, \n", convergence); */
            if(iter > *subiterations) break;

            for(r=0; r<*R; r++)
            {
                *Rchoice = r;

                /* Set MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QGJE to 0 */
                for(j=0; j<*K; j++)
                {
                    for(jj=0; jj<*K; jj++)
                    {
                        *(*(MNLNinv+j)+jj)=0;
                        *(*(JGFGinv+j)+jj)=0;
                        *(*(*(HNLS+j)+jj))=0;
                    }
                    for(m=0; m<*M; m++)
                    {
                        *(*(*(SNMH+j)+m))=0;
                    }
                    *(*(matk1+j)) = 0;
                }
                for(m=0; m<*M; m++)
                {
                    for(mm=0; mm<*M; mm++)
                    {
                        *(*(LNMNinv+m)+mm)=0;
                        *(*(FGJGinv+m)+mm)=0;
                    }
                }
                for(i=0; i<*P; i++)
                {
                    for(m=0; m<*M; m++)
                    {
                        *(*(*(QGJE+i)+m))=0;
                    }
                    for(j=0; j<*K; j++)
                    {
                        *(*(*(EGFQ+i)+j))=0;
                    }
                }
                SimplifyX(alpha, beta, gamma, delta, v, x, y, u, K, P,             /* Simplify function */
                    T, M, all, Rchoice, MNLNinv, LNMNinv, JGFGinv,
                    FGJGinv, HNLS, SNMH, EGFQ, QGJE);
                for(j=0; j<*K; j++)
                {
                    *(tempa+j) = 0;
                    *(tempc+j) = 0;
                }
                for(m=0; m<*M; m++)
                {
                    *(tempb+m) = 0;
                    *(tempd+m) = 0;
                }

                for(j=0; j<*K; j++)
                {
                    MatrixMult(MNLNinv, *K, *K, *(HNLS+j), 1, matk1);
                    for(jj=0; jj<*K; jj++)
                    {
                        *(tempa+jj) += (*(*(matk1+jj))) * (*(*(matk1+jj)));     /* temp.a */
                    }
                    MatrixMult(LNMNinv, *M, *M, *(SNMH+j), 1, matm1);
                    for(m=0; m<*M; m++)
                    {
                        *(tempb+m) += (*(*(matm1+m))) * (*(*(matm1+m)));        /* temp.b */
                    }
                }

                for(i=0; i<*P; i++)
                {
                    MatrixMult(JGFGinv, *K, *K, *(EGFQ+i), 1, matk1);
                    for(j=0; j<*K; j++)
                    {
                        *(*(matk1+j)) = (*(*(matk1+j))) * (*(*(matk1+j))) * (*(v+i));
                        *(tempc+j) += (*(*(matk1+j)));                          /* temp.c */
                    }
                    MatrixMult(FGJGinv, *M, *M, *(QGJE+i), 1, matm1);
                    for(m=0; m<*M; m++)
                    {
                        *(*(matm1+m)) = (*(*(matm1+m))) * (*(*(matm1+m))) * (*(v+i));
                        *(tempd+m) += (*(*(matm1+m)));                         /* temp.d */
                    }
                }

                for(j=0; j<*K; j++)
                {
                    *(*(alphanewmat+j)+r) = K2 * (1.0/(K2 * (*(*(MNLNinv+j)+j)) + (1.0/2.0)*(*(tempa+j))));
                    *(*(gammanewmat+j)+r) = P2 * (1.0/(P2 * (*(*(JGFGinv+j)+j)) + (1.0/2.0)*(*(tempc+j))));
                }
                for(m=0; m<*M; m++)
                {
                    *(*(betanewmat+m)+r) = K2 * (1/(K2 * (*(*(LNMNinv+m)+m)) + (1.0/2.0)*(*(tempb+m))));
                    *(*(deltanewmat+m)+r) = P2 * (1/(P2 * (*(*(FGJGinv+m)+m)) + (1.0/2.0)*(*(tempd+m))));
                }
            }

            for(j=0; j<*K; j++)
            {
                *(alphanew+j) = VecMedian(*(alphanewmat+j), R);             /* alpha.new */
                *(gammanew+j) = VecMedian(*(gammanewmat+j), R);             /* beta.new */
            }
            for(m=0; m<*M; m++)
            {
                *(betanew+m) = VecMedian(*(betanewmat+m),R);                /* gamma.new */
                *(deltanew+m) = VecMedian(*(deltanewmat+m),R);              /* delta.new */
            }

            /* Check convergence */
            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(alpha+j)) * (*(alpha+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(alpha+j) - *(alphanew+j))*(*(alpha+j) - *(alphanew+j)))/sumden;
            }
            alphadiff = sqrt(sumnum);                                       /* alpha.diff */

            sumnum = 0;
            sumden = 0;
            for(m=0; m<*M; m++)
            {
                sumden += (*(beta+m)) * (*(beta+m));
            }
            for(m=0; m<*M; m++)
            {
                sumnum += ((*(beta+m) - *(betanew+m))*(*(beta+m) - *(betanew+m)))/sumden;
            }
            betadiff = sqrt(sumnum);                                       /* beta.diff */

            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(gamma+j)) * (*(gamma+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(gamma+j) - *(gammanew+j))*(*(gamma+j) - *(gammanew+j)))/sumden;
            }
            gammadiff = sqrt(sumnum);                                       /* gamma.diff */

            sumnum = 0;
            sumden = 0;
            for(m=0; m<*M; m++)
            {
                sumden += (*(delta+m)) * (*(delta+m));
            }
            for(m=0; m<*M; m++)
            {
                sumnum += ((*(delta+m) - *(deltanew+m))*(*(delta+m) - *(deltanew+m)))/sumden;
            }
            deltadiff = sqrt(sumnum);                                       /* delta.diff */

            /* Find maximum convergence criterion */
            convergence = alphadiff;                                        /* Convergence criterion */
            if(betadiff>convergence) {convergence = betadiff;}
            if(gammadiff>convergence) {convergence = gammadiff;}
            if(deltadiff>convergence) {convergence = deltadiff;}

            for(j=0; j<*K; j++)                                             /* Update hyperparameters */
            {
                *(alpha+j) = *(alphanew+j);
                *(gamma+j) = *(gammanew+j);
            }
            for(m=0; m<*M; m++)
            {
                *(beta+m) = *(betanew+m);
                *(delta+m) = *(deltanew+m);
            }
            iter++;                                                            /* Update iteration number */
        }
    }

    /* Free MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, and QGJE */
    for(j=0; j<*KK; j++)
    {
        for(jj=0; jj<*KK; jj++)
        {
            free(*(*(HNLS+j)+jj));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(SNMH+j)+m));
        }
        free(*(HNLS+j));
        free(*(SNMH+j));
        free(*(MNLNinv+j));
        free(*(JGFGinv+j));
        free(*(matk1+j));
        free(*(alphanewmat+j));
        free(*(gammanewmat+j));
    }
    for(m=0; m<*M; m++)
    {
        free(*(LNMNinv+m));
        free(*(FGJGinv+m));
        free(*(matm1+m));
        free(*(betanewmat+m));
        free(*(deltanewmat+m));
    }

    for(i=0; i<*P; i++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(EGFQ+i)+j));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(QGJE+i)+m));
        }
        free(*(EGFQ+i));
        free(*(QGJE+i));
        free(*(matp1+i));
    }
    free(HNLS);
    free(SNMH);
    free(MNLNinv);
    free(JGFGinv);
    free(EGFQ);
    free(QGJE);
    free(LNMNinv);
    free(FGJGinv);
    free(matm1);
    free(matk1);
    free(matp1);
    free(tempa);
    free(tempb);
    free(tempc);
    free(tempd);
    free(alphanewmat);
    free(betanewmat);
    free(gammanewmat);
    free(deltanewmat);
    free(alphanew);
    free(betanew);
    free(gammanew);
    free(deltanew);

    /* If NOT estimating x's */
    A = (double***) calloc(*R, sizeof(double**));
    B = (double***) calloc(*R, sizeof(double**));
    C = (double***) calloc(*R, sizeof(double**));
    D = (double***) calloc(*R, sizeof(double**));
    Dvar = (double***) calloc(*R, sizeof(double**));
    temp = (double*) calloc(*M, sizeof(double));
    temp2 = (double*) calloc(*M, sizeof(double));
    deltanewmat = (double**) calloc(*M, sizeof(double*));
    deltanew = (double*) calloc(*M, sizeof(double));
    matm1 = (double**) calloc(*M, sizeof(double*));

    for(m=0; m<*M; m++)
    {
        *(deltanewmat+m) = (double*) calloc(*R, sizeof(double));
        *(matm1+m) = (double*) calloc(1, sizeof(double));
    }
    for(r=0; r<*R; r++)
    {
        *(A+r) = (double**) calloc(*KK, sizeof(double*));
        *(B+r) = (double**) calloc(*KK, sizeof(double*));
        *(C+r) = (double**) calloc(*P, sizeof(double*));
        *(D+r) = (double**) calloc(*P, sizeof(double*));
        *(Dvar+r) = (double**) calloc(*M, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(A+r)+j) = (double*) calloc(*KK, sizeof(double));
            *(*(B+r)+j) = (double*) calloc(*M, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(C+r)+i) = (double*) calloc(*KK, sizeof(double));
            *(*(D+r)+i) = (double*) calloc(*M, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(Dvar+r)+m) = (double*) calloc(*M, sizeof(double));
        }
    }

    if(*K == 0)
    {
        while(convergence > *conv)
        {
            if(iter > *subiterations) break;
            PostMeanR(alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M,
                A, B, C, D, Dvar);
            for(r=0; r<*R; r++)
            {
                for(m=0; m<*M; m++)
                {
                    *(temp+m) = 0;
                }
                for(i=0; i<*P; i++)
                {
                    for(m=0; m<*M; m++)
                    {
                        *(temp2+m) = 0;
                    }
                    for(t=0; t<*T; t++)
                    {
                        for(m=0; m<*M; m++)
                        {
                            *(temp2+m) += (*(*(*(u+r)+m)+t)) * (*(*(*(y+r)+i)+t));
                        }
                    }

                    /* Set matp1 = temp2 */
                    for(m=0; m<*M; m++)
                    {
                        *(*(matm1+m)) = *(temp2+m);
                    }

                    MatrixMult(*(Dvar+r), *M, *M, matm1, 1, matm1);
                    for(m=0; m<*M; m++)
                    {
                        *(temp+m) += (*(v+i)) * (*(*(matm1+m))) * (*(*(matm1+m)));
                    }
                }

                for(m=0; m<*M; m++)
                {
                    *(*(deltanewmat+m)+r) = P2 * (1/(P2* (*(*(*(Dvar+r)+m)+m)) + (1.0/2.0) * (*(temp+m))));
                }
            }
            for(m=0; m<*M; m++)
            {
                *(deltanew+m) = VecMedian(*(deltanewmat+m),R);                  /* delta.new */
            }

            sumnum = 0;
            sumden = 0;
            for(m=0; m<*M; m++)
            {
                sumden += (*(delta+m)) * (*(delta+m));
            }
            for(m=0; m<*M; m++)
            {
                sumnum += ((*(delta+m) - *(deltanew+m))*(*(delta+m) - *(deltanew+m)))/sumden;
            }
            convergence = sqrt(sumnum);                                         /* Convergence criterion */
            for(m=0; m<*M; m++)                                                 /* Update delta */
            {
                *(delta+m) = *(deltanew+m);
            }
            iter++;                                                             /* Update iteration number */
        }
    }

    /* Release memory */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(A+r)+j));
            free(*(*(B+r)+j));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(C+r)+i));
            free(*(*(D+r)+i));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(Dvar+r)+m));
        }
        free(*(A+r));
        free(*(B+r));
        free(*(C+r));
        free(*(D+r));
        free(*(Dvar+r));
    }
    for(m=0; m<*M; m++)
    {
        free(*(deltanewmat+m));
        free(*(matm1+m));
    }
    free(A);
    free(B);
    free(C);
    free(D);
    free(Dvar);
    free(temp);
    free(temp2);
    free(deltanewmat);
    free(deltanew);
    free(matm1);
    free(KK);
    free(Rchoice);
    free(all);
}

/******************************************************************/
/*  This is the counterpart of Posterior_Mean.R.  It takes as     */
/*  inputs alpha, beta, gamma, delta, v, x, y, u, K, P, T, R, M   */
/*  and returns the value of the posterior mean of A, B, C, and D */
/*  (in the case where x's are estimated), or the posterior mean  */
/*  and variance of D (in the case where x's are not estimated).  */
/******************************************************************/
void PostMeanR(double *alpha, double *beta, double *gamma, double *delta,
	double *v, double ***x, double ***y, double ***u, int *K, int *P, int *T,
	int *R,	int *M, double ***A, double ***B, double ***C, double ***D, double ***Dvar)
{
	int r, i, *all, j, jj, *Rchoice, *KK, m, mm;
	double **Dmean, **Dv;
    double **MNLNinv, **LNMNinv, **JGFGinv, **FGJGinv, ***HNLS,
            ***SNMH, ***EGFQ, ***QGJE, **matk1, **matm1;

	/* Load subprograms */
	void SimplifyNoX(double*, double*, double***, double***, int*, int*, int*, int*, int*, double**,
        double**);
	void SimplifyX(double*, double*, double*, double*, double*, double***, double***, double***, int*,
        int*, int*, int*, int*, int*, double**, double**, double**, double**, double***,
        double***, double***, double***);

	all = (int*) calloc(1, sizeof(int));
	Rchoice = (int*) calloc(1, sizeof(int));
	KK = (int*) calloc(1, sizeof(int));
	*all = 0;
	*KK = 1;
    if(*K > 0)
    {
        *KK = *K;
    }

    /* Allocate memory */
    Dmean = (double**) calloc(*P, sizeof(double*));
    Dv = (double**) calloc(*M, sizeof(double*));
    for(i=0; i<*P; i++)
    {
        *(Dmean + i) = (double*) calloc(*M, sizeof(double));
    }
    for(m=0; m<*M; m++)
    {
        *(Dv + m) = (double*) calloc(*M, sizeof(double));
    }

    /* Allocate MNLNinf, LNMNinf, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QFJE */
    MNLNinv = (double**) calloc(*KK, sizeof(double*));
    LNMNinv = (double**) calloc(*M, sizeof(double*));
    JGFGinv = (double**) calloc(*KK, sizeof(double*));
    FGJGinv = (double**) calloc(*M, sizeof(double*));
    HNLS = (double***) calloc(*KK, sizeof(double**));
    SNMH = (double***) calloc(*KK, sizeof(double**));
    EGFQ = (double***) calloc(*P, sizeof(double**));
    for(j=0; j<*KK; j++)
    {
        *(MNLNinv+j) = (double*) calloc(*KK, sizeof(double));
        *(JGFGinv+j) = (double*) calloc(*KK, sizeof(double));
        *(HNLS+j) = (double**) calloc(*KK, sizeof(double*));
        *(SNMH+j) = (double**) calloc(*M, sizeof(double*));
        for(jj=0; jj<*KK; jj++)
        {
            *(*(HNLS+j)+jj) = (double*) calloc(1, sizeof(double));
        }
        for(m=0; m<*M; m++)
        {
            *(*(SNMH+j)+m) = (double*) calloc(1, sizeof(double));
        }
    }
    for(m=0; m<*M; m++)
    {
        *(LNMNinv+m) = (double*) calloc(*M, sizeof(double));
        *(FGJGinv+m) = (double*) calloc(*M, sizeof(double));
    }

    for(i=0; i<*P; i++)
    {
        *(EGFQ+i) = (double**) calloc(*KK, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(EGFQ+i)+j) = (double*) calloc(1, sizeof(double));
        }
    }
    QGJE = (double***) calloc(*P, sizeof(double**));
    for(i=0; i<*P; i++)
    {
        *(QGJE+i) = (double**) calloc(*M, sizeof(double*));
        for(m=0; m<*M; m++)
        {
            *(*(QGJE+i)+m) = (double*) calloc(1, sizeof(double));
        }
    }

    matk1 = (double**) calloc(*KK, sizeof(double*));
    matm1 = (double**) calloc(*M, sizeof(double*));
    for(j=0; j<*KK; j++)
    {
        *(matk1+j) = (double*) calloc(1, sizeof(double));
    }
    for(m=0; m<*M; m++)
    {
        *(matm1+m) = (double*) calloc(1, sizeof(double));
    }

	/* Begin function -- no x's */
	if(*K==0)
	{
		for(r=0; r<*R; r++)
		{
			*Rchoice = r;
			/* Set Dmean and Dv to 0 */
			for(i=0; i<*P; i++)
			{
			    for(m=0; m<*M; m++)
			    {
			        *(*(Dmean+i)+m) = 0;
			    }
			}
			for(m=0; m<*M; m++)
			{
			    for(mm=0; mm<*M; mm++)
			    {
			        *(*(Dv+m)+mm) = 0;
			    }
			}
			SimplifyNoX(delta, v, y, u, P, T, M, Rchoice, all, Dmean, Dv);

			for(i=0; i<*P; i++)
			{
				for(m=0; m<*M; m++)
				{
					*(*(*(D+r)+i)+m) = *(*(Dmean+i)+m);
				}
			}
			for(m=0; m<*M; m++)
			{
			    for(mm=0; mm<*M; mm++)
			    {
					*(*(*(Dvar+r)+m)+mm) = *(*(Dv+m)+mm);
				}
			}
		}
	}

	/* Begin function -- x's */
	if(*K > 0)
	{
		for(r=0; r<*R; r++)
		{
			*Rchoice = r;

            /* Set MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QGJE to 0 */
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    *(*(MNLNinv+j)+jj)=0;
                    *(*(JGFGinv+j)+jj)=0;
                    *(*(*(HNLS+j)+jj))=0;
                }
                for(m=0; m<*M; m++)
                {
                    *(*(*(SNMH+j)+m))=0;
                }
            }
            for(m=0; m<*M; m++)
            {
                for(mm=0; mm<*M; mm++)
                {
                    *(*(LNMNinv+m)+mm)=0;
                    *(*(FGJGinv+m)+mm)=0;
                }
            }
            for(i=0; i<*P; i++)
            {
                for(m=0; m<*M; m++)
                {
                    *(*(*(QGJE+i)+m))=0;
                }
                for(j=0; j<*K; j++)
                {
                    *(*(*(EGFQ+i)+j))=0;
                }
            }

			SimplifyX(alpha, beta, gamma, delta, v, x, y, u, K, P,
				T, M, all, Rchoice, MNLNinv, LNMNinv, JGFGinv,
                FGJGinv, HNLS, SNMH, EGFQ, QGJE);

			for(j=0; j<*K; j++)
			{
				MatrixMult(MNLNinv, *K, *K, *(HNLS+j), 1, matk1);
				MatrixMult(LNMNinv, *M, *M, *(SNMH+j), 1, matm1);
				for(jj=0; jj<*K; jj++)
				{
					*(*(*(A+r)+j)+jj) = *(*(matk1+jj));
				}
				for(m=0; m<*M; m++)
				{
					*(*(*(B+r)+j)+m) = *(*(matm1+m));
				}
			}
			for(i=0; i<*P; i++)
			{
				MatrixMult(JGFGinv, *K, *K, *(EGFQ+i), 1, matk1);
				MatrixMult(FGJGinv, *M, *M, *(QGJE+i), 1, matm1);
				for(j=0; j<*K; j++)
				{
					*(*(*(C+r)+i)+j) = *(*(matk1+j));
				}
				for(m=0; m<*M; m++)
				{
					*(*(*(D+r)+i)+m) = *(*(matm1+m));
				}
			}
		}
	}

	/* Free memory */
    for(j=0; j<*KK; j++)
    {
        for(jj=0; jj<*KK; jj++)
        {
            free(*(*(HNLS+j)+jj));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(SNMH+j)+m));
        }
        free(*(HNLS+j));
        free(*(SNMH+j));
        free(*(MNLNinv+j));
        free(*(JGFGinv+j));
        free(*(matk1+j));
    }

    for(i=0; i<*P; i++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(EGFQ+i)+j));
        }
        for(m=0; m<*M; m++)
        {
            free(*(*(QGJE+i)+m));
        }
        free(*(EGFQ+i));
        free(*(QGJE+i));
        free(*(Dmean+i));

    }
    for(m=0; m<*M; m++)
    {
        free(*(LNMNinv+m));
        free(*(FGJGinv+m));
        free(*(matm1+m));
        free(*(Dv+m));
    }
    free(HNLS);
    free(SNMH);
    free(MNLNinv);
    free(JGFGinv);
    free(EGFQ);
    free(QGJE);
    free(LNMNinv);
    free(FGJGinv);
    free(matk1);
    free(matm1);
	free(all);
	free(Rchoice);
	free(KK);
    free(Dmean);
    free(Dv);
}


/* *************************************************** */
/* Function to calculate the median variance estimator */
/* This function is the C counterpart of the program   */
/* Variance_Estimator.R                                */
/* *************************************************** */
void VarMaxR(double ***x, double ***y, double ***u, double ***C, double ***D, int *P,
	int *R, int *T, int *K, int *M, double *vEst)
{
	int r, t, i, j, m, *KK;
	double **vTemp, **vTempt, *temp, **tempDy, **tempCx, **tempu,
		*tempyDy, *tempyCx, **tempxold, *tempyCxDy;

	/* Load subprograms */
	double VecMedian(double*, int*);

	/* Allocate temporary and vTemp matrices */
	KK = (int*) calloc(1, sizeof(int));
	*KK = 1;
	if(*K > 0)
	{
	    *KK = *K;
	}
	vTemp = (double**) calloc(*R, sizeof(double*));
	vTempt = (double**) calloc(*P, sizeof(double*));
	temp = (double*) calloc(*P, sizeof(double));
	tempDy = (double**) calloc(*P, sizeof(double*));
	tempCx = (double**) calloc(*P, sizeof(double*));
	tempu= (double**) calloc(*M, sizeof(double*));
	tempyDy = (double*) calloc(*P, sizeof(double));
	tempyCx = (double*) calloc(*P, sizeof(double));
	tempxold = (double**) calloc(*KK, sizeof(double*));
	tempyCxDy = (double*) calloc(*P, sizeof(double));

	for(r=0; r<*R; r++)
	{
		*(vTemp+r) = (double*) calloc(*P, sizeof(double));
	}
	for(i=0; i<*P; i++)
	{
		*(tempDy+i) = (double*) calloc(1, sizeof(double));
		*(tempCx+i) = (double*) calloc(1, sizeof(double));
		*(vTempt+i) = (double*) calloc(*R, sizeof(double));
	}
	for(m=0; m<*M; m++) {
		*(tempu+m) = (double*) calloc(1, sizeof(double));
	}
	for(j=0; j<*K; j++)
	{
		*(tempxold+j) = (double*) calloc(1, sizeof(double));
	}

	/* Calculate the variance estimator if no x's */
	if(*K == 0)
	{
		for(r=0; r<*R; r++)
		{
			for(i=0; i<*P; i++)
			{
				*(temp+i) = 0;
			}
			for(t=0; t<*T; t++)
			{
			    for(m=0; m<*M; m++)
			    {
			        *(*(tempu+m)) =*(*(*(u+r)+m)+t);
			    }
				MatrixMult(*(D+r),*P,*M,tempu,1,tempDy);
				for(i=0; i<*P; i++)
				{
					*(tempyDy+i)=*(*(*(y+r)+i)+t) - *(*(tempDy+i));
					*(temp+i) += (*(tempyDy+i)) * (*(tempyDy+i));
				}
			}

			for(i=0; i<*P; i++)
			{
				*(*(vTemp+r)+i) = (*(temp+i)) / ((*T)-1);
			}
		}
	}

	/* Calculate the variance estimator if x's */
    if(*K > 0)
    {
        for(r=0; r<*R; r++)
        {
            for(i=0; i<*P; i++)
            {
                *(temp+i) = 0;
            }
            for(t=0; t<*T; t++)
            {
                for(m=0; m<*M; m++)
                {
                    *(*(tempu+m))=*(*(*(u+r)+m)+t);
                }
				for(j=0; j<*K; j++)
				{
					*(*(tempxold+j))=*(*(*(x+r)+j)+t);
				}
                MatrixMult(*(D+r),*P,*M,tempu,1,tempDy);
				MatrixMult(*(C+r),*P,*K,tempxold,1,tempCx);
                for(i=0; i<*P; i++)
                {
                    *(tempyCxDy+i)=*(*(*(y+r)+i)+t) - *(*(tempCx+i)) - *(*(tempDy+i));
                    *(temp+i) += (*(tempyCxDy+i)) * (*(tempyCxDy+i));
                }
            }
            for(i=0; i<*P; i++)
            {
                *(*(vTemp+r)+i) = (*(temp+i)) / ((*T)-1);
            }
        }
    }

	/* Take inverse of all elements of vTemp matrix */
	for(r=0; r<*R; r++)
	{
		for(i=0; i<*P; i++)
		{
			*(*(vTemp+r)+i) = 1/(*(*(vTemp+r)+i));
		}
	}

	/* Rewrite the vTemp matrix as its transpose, find the median */
	MatrixTrans(vTemp, vTempt, R, P);
	for(i=0; i<*P; i++)
	{
		*(vEst+i) = VecMedian(*(vTempt+i), R);
	}

	/* Free memory */
	for(r=0; r<*R; r++)
	{
		free(*(vTemp+r));
	}
	for(i=0; i<*P; i++)
	{
		free(*(tempDy+i));
		free(*(tempCx+i));
		free(*(vTempt+i));
	}
	for(m=0; m<*M; m++)
	{
	    free(*(tempu+m));
	}
	for(j=0; j<*KK; j++)
	{
		free(*(tempxold+j));
	}
	free(tempyCx);
	free(tempxold);
	free(tempyCxDy);
	free(vTemp);
	free(vTempt);
	free(temp);
	free(tempDy);
	free(tempCx);
	free(tempu);
	free(tempyDy);
    free(KK);
}


/*********************************************************************/
/*  Function to calculate the median of a vector of a given length.  */
/*  Will return the median as a double.                              */
/*********************************************************************/
double VecMedian(double *vector, int *length)
{
	int index1, index2;
 	double len, tempnum = 1.0, median = 0;
	len = (double) (*length)/tempnum;

	/* Sort the vector using the rsort function in R */
	R_rsort(vector, *length);

	/* Odd case */
	if((len)/2 != floor((len)/2))
	{
		index1 = (int) floor((len)/2);
		median = *(vector + index1);
	}

	/* Even case */
	if((len)/2 == floor((len)/2))
	{
		index1 = (int) (((len)/2)-1);
		index2 = (int) ((len)/2);
		median = ((*(vector + index1)) + (*(vector + index2)))/2;
	}
	return(median);
}


/*************************************************/
/*  This is the C version of Simplifying_code.R  */
/*  in the case where x's are estimated.         */
/*************************************************/
void SimplifyX(double *alpha, double *beta, double *gamma, double *delta,
	double *v, double ***x, double ***y, double ***u, int *K, int *P, int *T, int *M, int *all,
	int *Rchoice, double **MNLNinv, double **LNMNinv, double **JGFGinv,
	double **FGJGinv, double ***HNLS, double ***SNMH, double ***EGFQ,
	double ***QGJE)
{
	int rlower = 0, rupper = 0, i, j, jj, r, t, *one, m, mm;
	double **xXx, **uXx, **uXu, **xXu, **xtemp,
		**txtemp, **utemp, **tutemp, *xsingle, *ysingle, *det,
       **J, **G, ***H, ***S, **MM, **N, **L,
		**F, ***E, ***Q, **Linv, **Minv, **Finv, **Jinv, **Nt, **Gt,
		**matkk, **matmm, **matkm, **matmk;

	/* Allocate memory, set M, N, L, H, and S = 0 */
	/* Allocate memory, set J, G, F, E, and Q = 0 */

	one = (int*) calloc(1, sizeof(int));
    det = (double*) calloc(1, sizeof(double));
    H = (double***) calloc(*K, sizeof(double**));
    S = (double***) calloc(*K, sizeof(double**));
    MM = (double**) calloc(*K, sizeof(double*));
    N = (double**) calloc(*M, sizeof(double*));
    L = (double**) calloc(*M, sizeof(double*));
    J = (double**) calloc(*K, sizeof(double*));
    G = (double**) calloc(*K, sizeof(double*));
    F = (double**) calloc(*M, sizeof(double*));
    Q = (double***) calloc(*P, sizeof(double**));
    E = (double***) calloc(*P, sizeof(double**));

	xXx = (double**) calloc(*K, sizeof(double*));
	uXx = (double**) calloc(*M, sizeof(double*));
	uXu = (double**) calloc(*M, sizeof(double*));
	xXu = (double**) calloc(*K, sizeof(double*));
	xtemp = (double**) calloc(*K, sizeof(double*));
	txtemp = (double**) calloc(1, sizeof(double*));
	utemp = (double**) calloc(*M, sizeof(double*));
	tutemp = (double**) calloc(1, sizeof(double*));
	xsingle = (double*) calloc(1, sizeof(double));
	ysingle = (double*) calloc(1, sizeof(double));

    for(m=0; m<*M; m++)
    {
        *(N+m) = (double*) calloc(*K, sizeof(double));
		*(L+m) = (double*) calloc(*M, sizeof(double));
		*(F+m) = (double*) calloc(*M, sizeof(double));
		*(utemp+m) = (double*) calloc(1, sizeof(double));
		*(uXx+m) = (double*) calloc(*K, sizeof(double));
		*(uXu+m) = (double*) calloc(*M, sizeof(double));
		for(j=0; j<*K; j++)
		{
		    *(*(N+m)+j) = 0;
		}
		for(mm=0; mm<*M; mm++)
		{
            *(*(L+m)+mm) = 0;
            *(*(F+m)+mm) = 0;

		}
    }
	for(i=0; i<*P; i++)
	{
		*(E + i) = (double**) calloc(*K, sizeof(double*));
		/* ERROR: MAY 30, 2012 */
        /* *(Q + i) = (double**) calloc(*P, sizeof(double*)); */
		*(Q + i) = (double**) calloc(*M, sizeof(double*));
		for(j=0; j<*K; j++)
		{
			*(*(E+i)+j) = (double*) calloc(1, sizeof(double));
			*(*(*(E+i)+j)) = 0;
		}
		for(m=0; m<*M; m++)
		{
			*(*(Q+i)+m) = (double*) calloc(1, sizeof(double));
			*(*(*(Q+i)+m)) = 0;
		}
	}
	for(j=0; j<*K; j++)
	{
		*(xXx + j) = (double*) calloc(*K, sizeof(double));
		*(xXu + j) = (double*) calloc(*M, sizeof(double));
		*(xtemp + j) = (double*) calloc(1, sizeof(double));
        *(H + j) = (double**) calloc(*K, sizeof(double*));
        *(S + j) = (double**) calloc(*M, sizeof(double*));
		*(MM + j) = (double*) calloc(*K, sizeof(double));
		*(J + j) = (double*) calloc(*K, sizeof(double));
		*(G + j) = (double*) calloc(*M, sizeof(double));
		for(jj=0; jj<*K; jj++)
		{
			*(*(H+j)+jj) = (double*) calloc(1, sizeof(double));
			*(*(MM+j)+jj) = 0;
			*(*(*(H+j)+jj)) = 0;
			*(*(J+j)+jj) = 0;
		}
		for(m=0; m<*M; m++)
		{
			*(*(S+j)+m) = (double*) calloc(1, sizeof(double));
			*(*(*(S+j)+m)) = 0;
			*(*(G+j)+m) = 0;
		}
	}
	*(txtemp) = (double*) calloc(*K, sizeof(double));
	*(tutemp) = (double*) calloc(*M, sizeof(double));

	/* Set limits of function */
	if(*all == 1)
	{
		rlower = 0;
		rupper = *Rchoice;
	}
	if(*all == 0)
	{
		rlower = *Rchoice;
		rupper = rlower + 1;
	}
	*one = 1;

    /* Begin function part 1 (corresponds to simplify.1 in R code) */
	for(r=rlower; r<rupper; r++)
	{
        for(t=1; t<*T; t++)
		{
			for(m=0; m<*M; m++)
			{
				*(*(utemp+m)) = *(*(*(u+r)+m)+t);
				*(*(tutemp)+m) = *(*(*(u+r)+m)+t);
			}
			for(j=0; j<*K; j++)
            {
                *(*(xtemp+j)) = *(*(*(x+r)+j)+t);           /* x_{t} */
            }
            MatrixMult(xtemp,*K,1,tutemp,*M,xXu);           /* temp.kp */

  			for(j=0; j<*K; j++)                             /* x_{t-1} */
            {
                *(*(xtemp+j)) = *(*(*(x+r)+j)+(t-1));
                *(*(txtemp)+j) = *(*(*(x+r)+j)+(t-1));
            }

			MatrixMult(xtemp,*K,1,txtemp,*K,xXx);		    /* temp.kk */
			MatrixMult(utemp,*M,1,txtemp,*K,uXx);			/* temp.pk */
			MatrixMult(utemp,*M,1,tutemp,*M,uXu);			/* temp.pp */

			for(j=0; j<*K; j++)
			{
				for(jj=0; jj<*K; jj++)
				{
					*(*(MM+j)+jj) += *(*(xXx+j)+jj);			/* M */
					*(*(J+j)+jj) += *(*(xXx+j)+jj);			/* J, need to add t = T */
				}

				for(m=0; m<*M; m++)
				{
					*(*(G+j)+m) += *(*(xXu+j)+m);			/* G, need to add t = 1 */
				}
			}

			for(m=0; m<*M; m++)
			{
				for(j=0; j<*K; j++)
				{
					*(*(N+m)+j) += *(*(uXx+m)+j);			/* N */
				}
				for(mm=0; mm<*M; mm++)
				{
					*(*(L+m)+mm) += *(*(uXu+m)+mm);			/* L */
					*(*(F+m)+mm) += *(*(uXu+m)+mm);			/* F, need to add t = 1 */
				}
			}

			for(j=0; j<*K; j++)
			{
				*xsingle = *(*(*(x+r)+j)+t);
				for(jj=0; jj<*K; jj++)
				{							                /* H */
					*(*(*(H+j)+jj)) += (*(*(xtemp+jj))) *
						(*xsingle);
				}
				for(m=0; m<*M; m++)
				{							                /* S */
					*(*(*(S+j)+m)) += (*(*(utemp+m))) *
						(*xsingle);
				}
			}

			for(i=0; i<*P; i++)
			{
				*ysingle = *(*(*(y+r)+i)+(t-1));
				for(j=0; j<*K; j++)					        /* E, need to add t = T */
				{
					*(*(*(E+i)+j)) += (*(*(xtemp+j))) *
						(*ysingle);
				}
				*ysingle = *(*(*(y+r)+i)+t);
				for(m=0; m<*M; m++)
				{							                /* Q, need to add t = 1 */
					*(*(*(Q+i)+m)) += (*(*(utemp+m))) *
						(*ysingle);
				}
			}
		}

		/* Add time t = T to the J matrix and E list */
  		for(j=0; j<*K; j++)
        {
            *(*(xtemp+j)) = *(*(*(x+r)+j)+(*T-1));
            *(*(txtemp)+j) = *(*(*(x+r)+j)+(*T-1));
        }
        MatrixMult(xtemp,*K,1,txtemp,*K,xXx);
		MatrixSum(J,xXx,J,K,K);
		for(i=0; i<*P; i++)
		{
			for(j=0; j<*K; j++)
			{
				*(*(xtemp+j)) = (*(*(*(x+r)+j)+(*T-1))) * (*(*(*(y+r)+i)+(*T-1)));
			}
			MatrixSum(*(E+i), xtemp, *(E+i), K, one);
		}

		/* Add first time t = 1 to the G and F matrices and the Q list*/
        for(j=0; j<*K; j++)
        {
            *(*(xtemp+j)) = *(*(*(x+r)+j));
        }
        for(m=0; m<*M; m++)
        {
            *(*(utemp+m)) = *(*(*(u+r)+m));
            *(*(tutemp)+m) = *(*(*(u+r)+m));
        }
        MatrixMult(xtemp,*K,1,tutemp,*M,xXu);
        MatrixSum(G,xXu,G,K,M);
        MatrixMult(utemp,*M,1,tutemp,*M,uXu);
        MatrixSum(F,uXu,F,M,M);
        for(i=0; i<*P; i++)
        {
            for(m=0; m<*M; m++)
            {
                *(*(utemp+m)) = (*(*(*(u+r)+m))) * (*(*(*(y+r)+i)));
            }
            MatrixSum(*(Q+i), utemp, *(Q+i), M, one);
        }

        for(j=0; j<*K; j++)
		{
            *(*(MM+j)+j) += *(alpha+j);
            *(*(J+j)+j) += *(gamma+j);
        }
        for(m=0; m<*M; m++)
        {
            *(*(L+m)+m) += *(beta+m);
            *(*(F+m)+m) += *(delta+m);
        }
	}

    /* Free Unnecessary Memory */
    /* Keep M, N, L, S, H, J, G, F, E, and Q for later */
	for(m=0; m<*M; m++)
    {
        free(*(utemp+m));
        free(*(uXx+m));
        free(*(uXu+m));
    }
    for(j=0; j<*K; j++)
    {
        free(*(xtemp+j));
        free(*(xXx+j));
		free(*(xXu+j));
    }
    free(*(txtemp));
    free(*(tutemp));
    free(xsingle);
	free(ysingle);
    free(utemp);
    free(xtemp);
    free(txtemp);
    free(tutemp);
    free(xXx);
	free(xXu);
    free(uXx);
    free(uXu);

	/* Begin function part 2 (corresponds to simplify in R code) */
	/* Allocate vectors Linv, Minv, Finv, and Jinv */
	Linv = (double**) calloc(*M, sizeof(double*));
	Minv = (double**) calloc(*K, sizeof(double*));
	Finv = (double**) calloc(*M, sizeof(double*));
	Jinv = (double**) calloc(*K, sizeof(double*));
	Nt = (double**) calloc(*K, sizeof(double*));
	Gt = (double**) calloc(*M, sizeof(double*));
	matkk = (double**) calloc(*K, sizeof(double*));
	matmm = (double**) calloc(*M, sizeof(double*));
	matkm = (double**) calloc(*K, sizeof(double*));
	matmk =	(double**) calloc(*M, sizeof(double*));
    for(m=0; m<*M; m++)
    {
        *(Linv+m) = (double*) calloc(*M, sizeof(double));
        *(Finv+m) = (double*) calloc(*M, sizeof(double));
		*(Gt+m) = (double*) calloc(*K, sizeof(double));
        *(matmm+m) = (double*) calloc(*M, sizeof(double));
		*(matmk+m) = (double*) calloc(*K, sizeof(double));
    }
	for(j=0; j<*K; j++)
	{
		*(Minv+j) = (double*) calloc(*K, sizeof(double));
		*(Jinv+j) = (double*) calloc(*K, sizeof(double));
		*(Nt+j) = (double*) calloc(*M, sizeof(double));
		*(matkk+j) = (double*) calloc(*K, sizeof(double));
		*(matkm+j) = (double*) calloc(*M, sizeof(double));
	}

	/* Find matrix inverses */
	MatrixInv(L, *M, Linv, det);							/* L.inv */
	MatrixInv(MM, *K, Minv, det);							/* M.inv */
	MatrixInv(F, *M, Finv, det);							/* F.inv */
	MatrixInv(J, *K, Jinv, det);							/* J.inv */

	MatrixTrans(N, Nt, M, K);
	MatrixTrans(G, Gt, K, M);

	MatrixMult(Nt, *K, *M, Linv, *M, matkm);
	MatrixMult(matkm, *K, *M, N, *K, matkk);
	for(j=0; j<*K; j++)
	{
		for(jj=0; jj<*K; jj++)
		{
			*(*(matkk+j)+jj) = *(*(MM+j)+jj) - *(*(matkk+j)+jj);
		}
	}
	MatrixInv(matkk, *K, MNLNinv, det);						/* MNLN.inv */

	MatrixMult(N, *M, *K, Minv, *K, matmk);
	MatrixMult(matmk, *M, *K, Nt, *M, matmm);
	for(m=0; m<*M; m++)
	{
		for(mm=0; mm<*M; mm++)
		{
			*(*(matmm+m)+mm) = *(*(L+m)+mm) - *(*(matmm+m)+mm);
		}
	}
	MatrixInv(matmm, *M, LNMNinv, det);						/* LNMN.inv */

	MatrixMult(G, *K, *M, Finv, *M, matkm);
	MatrixMult(matkm, *K, *M, Gt, *K, matkk);
	for(j=0; j<*K; j++)
	{
		for(jj=0; jj<*K; jj++)
		{
			*(*(matkk+j)+jj) = *(*(J+j)+jj) - *(*(matkk+j)+jj);
		}
	}
	MatrixInv(matkk, *K, JGFGinv, det);						/*  JGFG.inv */

	MatrixMult(Gt, *M, *K, Jinv, *K, matmk);
	MatrixMult(matmk, *M, *K, G, *M, matmm);
	for(m=0; m<*M; m++)
	{
		for(mm=0; mm<*M; mm++)
		{
			*(*(matmm+m)+mm) = *(*(F+m)+mm) - *(*(matmm+m)+mm);
		}
	}
	MatrixInv(matmm, *M, FGJGinv, det);						/* FGJG.inv */

	for(j=0; j<*K; j++)
	{
		MatrixMult(Nt, *K, *M, Linv, *M, matkm);
		MatrixMult(matkm, *K, *M, *(S+j), 1, *(HNLS+j));
		MatrixMult(N, *M, *K, Minv, *K, matmk);
		MatrixMult(matmk, *M, *K, *(H+j), 1, *(SNMH+j));
		for(jj=0; jj<*K; jj++)
		{
			*(*(*(HNLS+j)+jj)) = *(*(*(H+j)+jj)) - *(*(*(HNLS+j)+jj)); 	/* HNLS */
		}
		for(m=0; m<*M; m++)
		{
			*(*(*(SNMH+j)+m)) = *(*(*(S+j)+m)) - *(*(*(SNMH+j)+m));	 	/* SNMH */
		}
	}

	for(i=0; i<*P; i++)
	{
		MatrixMult(G, *K, *M, Finv, *M, matkm);
		MatrixMult(matkm, *K, *M, *(Q+i), 1, *(EGFQ+i));
		MatrixMult(Gt, *M, *K, Jinv, *K, matmk);
		MatrixMult(matmk, *M, *K, *(E+i), 1, *(QGJE+i));
		for(j=0; j<*K; j++)
		{
			*(*(*(EGFQ+i)+j)) = *(*(*(E+i)+j)) - *(*(*(EGFQ+i)+j));
		}
		for(m=0; m<*M; m++)
		{
			*(*(*(QGJE+i)+m)) = *(*(*(Q+i)+m)) - *(*(*(QGJE+i)+m));
		}
	}

	/* Free Memory */
	for(i=0; i<*P; i++)
	{
	    for(j=0; j<*K; j++)
		{
			free(*(*(E+i)+j));
		}
		for(m=0; m<*M; m++)
		{
			free(*(*(Q+i)+m));
		}
		free(*(E+i));
		free(*(Q+i));
	}

    for(m=0; m<*M; m++)
    {
        free(*(N+m));
        free(*(L+m));
        free(*(F+m));
        free(*(Linv+m));
		free(*(Finv+m));
        free(*(Gt+m));
        free(*(matmm+m));
		free(*(matmk+m));
    }

	for(j=0; j<*K; j++)
	{
		free(*(MM+j));
		for(jj=0; jj<*K; jj++)
		{
			free(*(*(H+j)+jj));
		}
		for(m=0; m<*M; m++)
		{
			free(*(*(S+j)+m));
		}
		free(*(H+j));
		free(*(S+j));
		free(*(J+j));
		free(*(G+j));
		free(*(Minv+j));
		free(*(Jinv+j));
		free(*(Nt+j));
		free(*(matkk+j));
		free(*(matkm+j));
	}
	free(matmm);
	free(matmk);
	free(Gt);
	free(Linv);
	free(Finv);
	free(matkk);
	free(matkm);
	free(Nt);
	free(Minv);
	free(Jinv);
	free(E);
    free(Q);
    free(F);
	free(G);
	free(J);
	free(H);
	free(S);
	free(MM);
	free(N);
	free(L);
	free(det);
	free(one);
}

/*********************************************************************************/
/*  Simplify code, corresponds to the simplify.2 and simplify functions in R     */
/*  in the case where there are no x's to be estimated.  The inputs are          */
/*  delta = current value of delta, v = current value of v, y, u, P, T, M,       */
/*  Rchoice (=a given r if all = FALSE, =R if all = TRUE), all = whether 	     */
/*  estimation is done for a single replicate (FALSE) or all replicates (TRUE),  */
/*  DmeanNox = matrix to hold D mean, DvarNox = matrix to hold Dvar.		     */
/*										                                         */
/*  Note: This is the base variance, to get actual variance for each row, need to*/
/*  multiply by v_i^(-1).							                             */
/*  DmeanNox is PxM, and DvarNox is MxM                                          */
/*********************************************************************************/
void SimplifyNoX(double *delta, double *v, double ***y, double ***u, int *P, int *T,
    int *M, int *Rchoice, int *all, double **DmeanNox, double **DvarNox)
{
	int lower = 0, upper = 0, r, t, i, m, mm;
	double **m1, **m2, **utemp, **tutemp, **uXu, **DmeanNoxt, *det;
    if(*all == 1) {lower = 0; upper = *Rchoice;}
    if(*all == 0) {lower = *Rchoice; upper = (*Rchoice)+1;}

	/* Allocate memory and initialize */
	det = (double*) calloc(1, sizeof(double));
	m1 = (double**) calloc(*M, sizeof(double*));
	m2 = (double**) calloc(*M, sizeof(double*));
	uXu = (double**) calloc(*M, sizeof(double*));
	utemp = (double**) calloc(*M, sizeof(double*));
	tutemp = (double**) calloc(1, sizeof(double*));
	DmeanNoxt = (double**) calloc(*M, sizeof(double*));

	for(m=0; m<*M; m++)
	{
	    *(m1+m) = (double*) calloc(*M, sizeof(double));
        *(m2+m) = (double*) calloc(*P, sizeof(double));
        *(utemp+m) = (double*) calloc(1, sizeof(double));
        *(uXu+m) = (double*) calloc(*M, sizeof(double));
        *(DmeanNoxt+m) = (double*) calloc(*P, sizeof(double));
	    for(mm=0; mm<*M; mm++)
	    {
	        *(*(m1+m)+mm) = 0;
	    }
	    for(i=0; i<*P; i++)
	    {
	        *(*(m2+m)+i) = 0;
	    }
	}
	*(tutemp) = (double*) calloc(*M, sizeof(double));

	/* Begin function, Part 1 */
	/* Corresponds to simplify.2 in R code */

	for(r=lower; r<upper; r++)
	{
		for(t=0; t<*T; t++)
		{
			for(m=0; m<*M; m++)
			{
				*(*(utemp+m)) = *(*(*(u+r)+m)+t);
				*(*(tutemp)+m) = *(*(*(u+r)+m)+t);
			}
			MatrixMult(utemp, *M, 1, tutemp, *M, uXu);
			for(m=0; m<*M; m++)
			{
				for(mm=0; mm<*M; mm++)
				{
					*(*(m1+m)+mm) += *(*(uXu+m)+mm);
				}
			}
			for(m=0; m<*M; m++)
			{
				for(i=0; i<*P; i++)
				{
					*(*(m2+m)+i) += (*(*(utemp+m))) *
						(*(*(*(y+r)+i)+t));
				}
			}
		}
	}

	/* Begin function, Part 2
	 Corresponds to simplify in R code */

	for(m=0; m<*M; m++)
	{
		*(*(m1+m)+m) += *(delta + m);
	}

	/* This is the base variance, to get actual variance for each row, need to multiply by v_i^(-1) */
	MatrixInv(m1, *M, DvarNox, det);

	/* This is the mean */
	MatrixMult(DvarNox,*M,*M,m2,*P,DmeanNoxt);
	MatrixTrans(DmeanNoxt, DmeanNox, M, P);

	/* Free memory */
	for(m=0; m<*M; m++)
	{
	    free(*(m1+m));
	    free(*(m2+m));
        free(*(utemp+m));
        free(*(uXu+m));
        free(*(DmeanNoxt+m));
	}
	free(m1);
	free(m2);
	free(DmeanNoxt);
	free(uXu);
	free(utemp);
	free(tutemp);
	free(det);
}


/******************************************************************/
/*  Function to calculate the inverse and determinant of a matrix */
/*  Code thanks to Cherie Ochsenfeld, using LAPACK function       */
/*  "mat" is an n x n matrix, "invmat" is the double pointer to   */
/*  hold the inverse                                              */
/* 								                                  */
/*  Code needs to be compiled with -llapack.                      */
/******************************************************************/
void MatrixInv(double **mat, int n, double **invmat, double *det)
{
	double **u, *w, **v;
	int i, j, k;
	double *fvect, *u2, *vt, *work;
	char jobu = 'A', jobvt = 'A';
	int lwork = -1, info = 0, m = n;

	/* Allocate memory */
	fvect = (double*) calloc(n*n, sizeof(double));
	u2 = (double*) calloc(n*n, sizeof(double));
	w = (double*) calloc(n, sizeof(double));
	vt = (double*) calloc(n*n, sizeof(double));
	work = (double*) calloc(1, sizeof(double));
	u = (double**) calloc(n, sizeof(double*));
	v = (double**) calloc(n, sizeof(double*));

	for(i=0; i<n; i++)
	{
		*(u+i) = (double*) calloc(n, sizeof(double));
		*(v+i) = (double*) calloc(n, sizeof(double));
	}

	/* Convert matrix to vector form for Lapack SVD call */
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(fvect+(j+i*n)) = *(*(mat+j)+i);

	/* Singular value decomposition using LAPACK function */
	F77_CALL(dgesvd)(&jobu,&jobvt,&m,&m,fvect,&m,w,u2,&m,vt,&m,work,&lwork,&info);
	lwork = *work;
	free(work);
	work = (double*) calloc(lwork, sizeof(double));
	F77_CALL(dgesvd)(&jobu,&jobvt,&m,&m,fvect,&m,w,u2,&m,vt,&m,work,&lwork,&info);

	/* Convert U and V to matrix form */
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(*(u+j)+i)=*(u2+(j+i*n));
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(*(v+i)+j)=*(vt+(j+i*n));

	/* Inversion of SVD */
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(*(v+i)+j) = *(*(v+i)+j) *(1/(*(w+j)));

	/* Multiply inverted SVD to get inverted matrix */
	for(i=0; i<n; i++) for(j=0; j<n; j++) for(k=0; k<n; k++) *(*(invmat+i)+j) = *(*(invmat+i)+j) + (*(*(v+i)+k)) * (*(*(u+j)+k));

	/* Calculate the determinate of the matrix */
	*det=0;
	for(i=0; i<n; i++) *det += log(*(w+i));
	/* Free memory */
	for(i=0; i<n; i++)
	{
		free(*(u+i));
		free(*(v+i));
	}
	free(w);
	free(u);
	free(v);
	free(fvect);
	free(u2);
	free(vt);
	free(work);
}

/****************************************************************************/
/*  This is a function to perform matrix multiplication, m1 %*% m2 = sol.   */
/*  This code is thanks to Cherie Ochsenfeld, and uses the BLAS function.   */
/*  m1r and m1c are the number of rows and cols of m1, and m2c is the       */
/*  number of columns of m2. 					`	                        */
/****************************************************************************/
void MatrixMult(double **m1, int m1r, int m1c, double **m2,
	int m2c, double **sol)
{
	int i,j;
	double *A, *B, *C, alph=1.0, bta=0.0;
	char transa='N', transb='N';

	/* Allocate memory */
	A = (double*) calloc((m1r*m1c), sizeof(double));
	B = (double*) calloc((m1c*m2c), sizeof(double));
	C = (double*) calloc((m1r*m2c), sizeof(double));

	/* Turn matrix into a vector */
	for(i=0; i<m1c; i++) for(j=0; j<m1r; j++) *(A+(i*m1r+j))=*(*(m1+j)+i);
	for(i=0; i<m2c; i++) for(j=0; j<m1c; j++) *(B+(i*m1c+j))=*(*(m2+j)+i);

	/* Call BLAS function */
	F77_CALL(dgemm)(&transa, &transb, &m1r, &m2c, &m1c, &alph, A, &m1r, B, &m1c, &bta, C, &m1r);

	/* set the solution as a matrix */
	for(i=0; i<m2c; i++) for(j=0; j<m1r; j++) *(*(sol+j)+i) = *(C+(i*m1r+j));

	/* Free memory */
	free(A);
	free(B);
	free(C);
}

/************************************************************************/
/*  Function to perform element-wise addition of two matrices,          */
/*  m1 + m2 = sum.                                                      */
/************************************************************************/
void MatrixSum(double **m1, double **m2, double **sum, int *row, int *col)
{
	int i, j;
	for(i=0; i<*row; i++)
	{
		for(j=0; j<*col; j++)
		{
			*(*(sum+i)+j) = (*(*(m1+i)+j)) + (*(*(m2+i)+j));
		}
	}
}

/*************************************************************/
/*  Function to transpose a matrix contained in a double     */
/*  pointer, mt = t(m).  "row" and "col" are the number of   */
/*  rows and columns of the original matrix, m.              */
/*************************************************************/
void MatrixTrans(double **m, double **mt, int *row, int *col)
{
	int i,j;
	for(i=0; i<*row; i++)
	{
		for(j=0; j<*col; j++)
		{
			*(*(mt+j)+i) = *(*(m+i)+j);
		}
	}
}

#ifdef win32
            R_FlushConsole();
#endif
