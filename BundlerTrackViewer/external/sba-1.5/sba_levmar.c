/////////////////////////////////////////////////////////////////////////////////
//// 
////  Expert drivers for sparse bundle adjustment based on the
////  Levenberg - Marquardt minimization algorithm
////  Copyright (C) 2004-2008 Manolis Lourakis (lourakis at ics forth gr)
////  Institute of Computer Science, Foundation for Research & Technology - Hellas
////  Heraklion, Crete, Greece.
////
////  This program is free software; you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation; either version 2 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
///////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "compiler.h"
#include "sba.h"
#include "sba_chkjac.h"

#include "matrix.h"

#define SBA_EPSILON       1E-12
#define SBA_EPSILON_SQ    ( (SBA_EPSILON)*(SBA_EPSILON) )

#define SBA_ONE_THIRD     0.3333333334 /* 1.0/3.0 */


#define emalloc(sz)       emalloc_(__FILE__, __LINE__, sz)

#define FABS(x)           (((x)>=0)? (x) : -(x))

#define ROW_MAJOR         0
#define COLUMN_MAJOR      1
#define MAT_STORAGE       COLUMN_MAJOR


// #define TIMINGS

#ifdef TIMINGS
#include <time.h>
#endif

/* contains information necessary for computing a finite difference approximation to a jacobian,
 * e.g. function to differentiate, problem dimensions and pointers to working memory buffers
 */
struct fdj_data_x_ {
    void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata); /* function to differentiate */
    int cnp, pnp, mnp;  /* parameter numbers */
    int *func_rcidxs,
        *func_rcsubs;   /* working memory for func invocations.
                         * Notice that this has to be different
                         * than the working memory used for
                         * evaluating the jacobian!
                         */
    double *hx, *hxx;   /* memory to save results in */
    void *adata;
};

/* auxiliary memory allocation routine with error checking */
inline static void *emalloc_(char *file, int line, size_t sz)
{
    void *ptr;

    ptr=(void *)malloc(sz);
    if(ptr==NULL){
        fprintf(stderr, "SBA: memory allocation request for %u bytes failed in file %s, line %d, exiting", sz, file, line);
        exit(1);
    }

    return ptr;
}

/* auxiliary routine for clearing an array of doubles */
inline static void _dblzero(register double *arr, register int count)
{
    while(--count>=0)
        *arr++=0.0;
}

/* auxiliary routine for computing the mean reprojection error; used for debugging */
static double sba_mean_repr_error(int n, int mnp, double *x, double *hx, struct sba_crsm *idxij, int *rcidxs, int *rcsubs)
{
    register int i, j;
    int nnz, nprojs;
    double *ptr1, *ptr2;
    double err;

    for(i=nprojs=0, err=0.0; i<n; ++i){
        nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
        nprojs+=nnz;
        for(j=0; j<nnz; ++j){ /* point i projecting on camera rcsubs[j] */
			ptr1=x + idxij->val[rcidxs[j]]*mnp;
            ptr2=hx + idxij->val[rcidxs[j]]*mnp;

            err+=sqrt((ptr1[0]-ptr2[0])*(ptr1[0]-ptr2[0]) + (ptr1[1]-ptr2[1])*(ptr1[1]-ptr2[1]));
        }
    }

    return err/((double)(nprojs));
}

/* print the solution in p using sba's text format. If cnp/pnp==0 only points/cameras are printed */
static void sba_print_sol(int n, int m, double *p, int cnp, int pnp, double *x, int mnp, struct sba_crsm *idxij, int *rcidxs, int *rcsubs)
{
    register int i, j, ii;
    int nnz;
    double *ptr;

    if(cnp){
        /* print camera parameters */
        for(j=0; j<m; ++j){
            ptr=p+cnp*j;
            for(ii=0; ii<cnp; ++ii)
                printf("%g ", ptr[ii]);
            printf("\n");
        }
    }

    if(pnp){
        /* 3D & 2D point parameters */
        printf("\n\n\n# X Y Z  nframes  frame0 x0 y0  frame1 x1 y1 ...\n");
        for(i=0; i<n; ++i){
            ptr=p+cnp*m+i*pnp;
            for(ii=0; ii<pnp; ++ii) // print 3D coordinates
                printf("%g ", ptr[ii]);

            nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
            printf("%d ", nnz);

            for(j=0; j<nnz; ++j){ /* point i projecting on camera rcsubs[j] */
                ptr=x + idxij->val[rcidxs[j]]*mnp;

                printf("%d ", rcsubs[j]);
                for(ii=0; ii<mnp; ++ii) // print 2D coordinates
                    printf("%g ", ptr[ii]);
            }
            printf("\n");
        }
    }
}

/* Compute e=x-y for two n-vectors x and y and return the squared L2 norm of e.
 * e can coincide with either x or y. 
 * Uses loop unrolling and blocking to reduce bookkeeping overhead & pipeline
 * stalls and increase instruction-level parallelism; see http://www.abarnett.demon.co.uk/tutorial.html
 */
static double nrmL2xmy(double *const e, const double *const x, const double *const y, const int n)
{
    const int blocksize=8, bpwr=3; /* 8=2^3 */
    register int i;
    int j1, j2, j3, j4, j5, j6, j7;
    int blockn;
    register double sum0=0.0, sum1=0.0, sum2=0.0, sum3=0.0;

    /* n may not be divisible by blocksize, 
     * go as near as we can first, then tidy up.
     */
    blockn = (n>>bpwr)<<bpwr; /* (n / blocksize) * blocksize; */

    /* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
    for(i=blockn-1; i>0; i-=blocksize){
        e[i ]=x[i ]-y[i ]; sum0+=e[i ]*e[i ];
        j1=i-1; e[j1]=x[j1]-y[j1]; sum1+=e[j1]*e[j1];
        j2=i-2; e[j2]=x[j2]-y[j2]; sum2+=e[j2]*e[j2];
        j3=i-3; e[j3]=x[j3]-y[j3]; sum3+=e[j3]*e[j3];
        j4=i-4; e[j4]=x[j4]-y[j4]; sum0+=e[j4]*e[j4];
        j5=i-5; e[j5]=x[j5]-y[j5]; sum1+=e[j5]*e[j5];
        j6=i-6; e[j6]=x[j6]-y[j6]; sum2+=e[j6]*e[j6];
        j7=i-7; e[j7]=x[j7]-y[j7]; sum3+=e[j7]*e[j7];
    }

    /*
     * There may be some left to do.
     * This could be done as a simple for() loop, 
     * but a switch is faster (and more interesting) 
     */

    i=blockn;
    if(i<n){ 
        /* Jump into the case at the place that will allow
         * us to finish off the appropriate number of items. 
         */
        switch(n - i){ 
        case 7 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
        case 6 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
        case 5 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
        case 4 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
        case 3 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
        case 2 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
        case 1 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
        }
    }

    return sum0+sum1+sum2+sum3;
}

/* Compute e=W*(x-y) for two n-vectors x and y and return the squared L2 norm of e.
 * This norm equals the squared C norm of x-y with C=W^T*W, W being block diagonal
 * matrix with nvis mnp*mnp blocks: e^T*e=(x-y)^T*W^T*W*(x-y)=||x-y||_C.
 * Note that n=nvis*mnp; e can coincide with either x or y.
 *
 * Similarly to nrmL2xmy() above, uses loop unrolling and blocking
 */
static double nrmCxmy(double *const e, const double *const x, const double *const y,
                      const double *const W, const int mnp, const int nvis)
{
    const int n=nvis*mnp;
    const int blocksize=8, bpwr=3; /* 8=2^3 */
    register int i, ii, k;
    int j1, j2, j3, j4, j5, j6, j7;
    int blockn, mnpsq;
    register double norm, sum;
    register const double *Wptr, *auxptr;
    register double *eptr;

    /* n may not be divisible by blocksize, 
     * go as near as we can first, then tidy up.
     */
    blockn = (n>>bpwr)<<bpwr; /* (n / blocksize) * blocksize; */

	printf("a\n");
	/* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
    for(i=blockn-1; i>0; i-=blocksize){
        e[i ]=x[i ]-y[i ];
        j1=i-1; e[j1]=x[j1]-y[j1];
        j2=i-2; e[j2]=x[j2]-y[j2];
        j3=i-3; e[j3]=x[j3]-y[j3];
        j4=i-4; e[j4]=x[j4]-y[j4];
        j5=i-5; e[j5]=x[j5]-y[j5];
        j6=i-6; e[j6]=x[j6]-y[j6];
        j7=i-7; e[j7]=x[j7]-y[j7];
    }

	printf("b\n");
    /*
     * There may be some left to do.
     * This could be done as a simple for() loop, 
     * but a switch is faster (and more interesting) 
     */

    i=blockn;
    if(i<n){ 
        /* Jump into the case at the place that will allow
         * us to finish off the appropriate number of items. 
         */
        switch(n - i){ 
        case 7 : e[i]=x[i]-y[i]; ++i;
        case 6 : e[i]=x[i]-y[i]; ++i;
        case 5 : e[i]=x[i]-y[i]; ++i;
        case 4 : e[i]=x[i]-y[i]; ++i;
        case 3 : e[i]=x[i]-y[i]; ++i;
        case 2 : e[i]=x[i]-y[i]; ++i;
        case 1 : e[i]=x[i]-y[i]; ++i;
        }
    }

	printf("c\n");
    /* compute w_x_ij e_ij in e and its L2 norm.
     * Since w_x_ij is upper triangular, the products can be safely saved
     * directly in e_ij, without the need for intermediate storage
     */
    mnpsq=mnp*mnp;
    /* Wptr, eptr point to w_x_ij, e_ij below */
    for(i=0, Wptr=W, eptr=e, norm=0.0; i<nvis; ++i, Wptr+=mnpsq, eptr+=mnp){
        for(ii=0, auxptr=Wptr; ii<mnp; ++ii, auxptr+=mnp){ /* auxptr=Wptr+ii*mnp */
            for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
                sum+=auxptr[k]*eptr[k]; //Wptr[ii*mnp+k]*eptr[k];
            eptr[ii]=sum;
            norm+=sum*sum;
        }
    }

	printf("d\n");
    return norm;
}

/* search for & print image projection components that are infinite; useful for identifying errors */
static void sba_print_inf(double *hx, int nimgs, int mnp, struct sba_crsm *idxij, int *rcidxs, int *rcsubs)
{
    register int i, j, k;
    int nnz;
    double *phxij;

    for(j=0; j<nimgs; ++j){
        nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */
        for(i=0; i<nnz; ++i){
            phxij=hx + idxij->val[rcidxs[i]]*mnp;
            for(k=0; k<mnp; ++k)
                if(!SBA_FINITE(phxij[k]))
                    printf("SBA: component %d of the estimated projection of point %d on camera %d is invalid!\n", k, rcsubs[i], j);
        }
    }
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is approximated with the aid of finite differences and is returned in the order
 * (A_11, B_11, ..., A_1m, B_1m, ..., A_n1, B_n1, ..., A_nm, B_nm),
 * where A_ij=dx_ij/da_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 * Problem-specific information is assumed to be stored in a structure pointed to by "dat".
 *
 * NOTE: The jacobian (for n=4, m=3) in matrix form has the following structure:
 *       A_11  0     0     B_11 0    0    0
 *       0     A_12  0     B_12 0    0    0
 *       0     0     A_13  B_13 0    0    0
 *       A_21  0     0     0    B_21 0    0
 *       0     A_22  0     0    B_22 0    0
 *       0     0     A_23  0    B_23 0    0
 *       A_31  0     0     0    0    B_31 0
 *       0     A_32  0     0    0    B_32 0
 *       0     0     A_33  0    0    B_33 0
 *       A_41  0     0     0    0    0    B_41
 *       0     A_42  0     0    0    0    B_42
 *       0     0     A_43  0    0    0    B_43
 *       To reduce the total number of objective function evaluations, this structure can be
 *       exploited as follows: A certain d is added to the j-th parameters of all cameras and
 *       the objective function is evaluated at the resulting point. This evaluation suffices
 *       to compute the corresponding columns of *all* A_ij through finite differences. A similar
 *       strategy allows the computation of the B_ij. Overall, only cnp+pnp+1 objective function
 *       evaluations are needed to compute the jacobian, much fewer compared to the m*cnp+n*pnp+1
 *       that would be required by the naive strategy of computing one column of the jacobian
 *       per function evaluation. See Nocedal-Wright, ch. 7, pp. 169. Although this approach is
 *       much faster compared to the naive strategy, it is not preferable to analytic jacobians,
 *       since the latter are considerably faster to compute and result in fewer LM iterations.
 */
static void sba_fdjac_x(
                        double *p,                /* I: current parameter estimate, (m*cnp+n*pnp)x1 */
                        struct sba_crsm *idxij,   /* I: sparse matrix containing the location of x_ij in hx */
                        int    *rcidxs,           /* work array for the indexes of nonzero elements of a single sparse matrix row/column */
                        int    *rcsubs,           /* work array for the subscripts of nonzero elements in a single sparse matrix row/column */
                        double *jac,              /* O: array for storing the approximated jacobian */
                        void   *dat)              /* I: points to a "fdj_data_x_" structure */
{
    register int i, j, ii, jj;
    double *pa, *pb, *pqr, *ppt;
    register double *pAB, *phx, *phxx;
    int n, m, nm, nnz, Asz, Bsz, ABsz, idx;

    double *tmpd;
    register double d;

    struct fdj_data_x_ *fdjd;
    void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
    double *hx, *hxx;
    int cnp, pnp, mnp;
    void *adata;


    /* retrieve problem-specific information passed in *dat */
    fdjd=(struct fdj_data_x_ *)dat;
    func=fdjd->func;
    cnp=fdjd->cnp; pnp=fdjd->pnp; mnp=fdjd->mnp;
    hx=fdjd->hx;
    hxx=fdjd->hxx;
    adata=fdjd->adata;

    n=idxij->nr; m=idxij->nc;
    pa=p; pb=p+m*cnp;
    Asz=mnp*cnp; Bsz=mnp*pnp; ABsz=Asz+Bsz;

    nm=(n>=m)? n : m; // max(n, m);
    tmpd=(double *)emalloc(nm*sizeof(double));

    (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hx, adata); // evaluate supplied function on current solution

    if(cnp){ // is motion varying?
        /* compute A_ij */
        for(jj=0; jj<cnp; ++jj){
            for(j=0; j<m; ++j){
                pqr=pa+j*cnp; // j-th camera parameters
                /* determine d=max(SBA_DELTA_SCALE*|pqr[jj]|, SBA_MIN_DELTA), see HZ */
                d=(double)(SBA_DELTA_SCALE)*pqr[jj]; // force evaluation
                d=FABS(d);
                if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;

                tmpd[j]=d;
                pqr[jj]+=d;
            }

            (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);

            for(j=0; j<m; ++j){
                pqr=pa+j*cnp; // j-th camera parameters
                pqr[jj]-=tmpd[j]; /* restore */
                d=1.0/tmpd[j]; /* invert so that divisions can be carried out faster as multiplications */

                nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
                for(i=0; i<nnz; ++i){
                    idx=idxij->val[rcidxs[i]];
                    phx=hx + idx*mnp; // set phx to point to hx_ij
                    phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    pAB=jac + idx*ABsz; // set pAB to point to A_ij

                    for(ii=0; ii<mnp; ++ii)
                        pAB[ii*cnp+jj]=(phxx[ii]-phx[ii])*d;
                }
            }
        }
    }

    if(pnp){ // is structure varying?
        /* compute B_ij */
        for(jj=0; jj<pnp; ++jj){
            for(i=0; i<n; ++i){
                ppt=pb+i*pnp; // i-th point parameters
                /* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
                d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
                d=FABS(d);
                if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;

                tmpd[i]=d;
                ppt[jj]+=d;
            }

            (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);

            for(i=0; i<n; ++i){
                ppt=pb+i*pnp; // i-th point parameters
                ppt[jj]-=tmpd[i]; /* restore */
                d=1.0/tmpd[i]; /* invert so that divisions can be carried out faster as multiplications */

                nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
                for(j=0; j<nnz; ++j){
                    idx=idxij->val[rcidxs[j]];
                    phx=hx + idx*mnp; // set phx to point to hx_ij
                    phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    pAB=jac + idx*ABsz + Asz; // set pAB to point to B_ij

                    for(ii=0; ii<mnp; ++ii)
                        pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
                }
            }
        }
    }

    free(tmpd);
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras and a grid of repetitions, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is approximated with the aid of finite differences and is returned in the order
 * (A_11, B_11, ..., A_1m, B_1m, ..., A_n1, B_n1, ..., A_nm, B_nm),
 * where A_ij=dx_ij/da_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 * Problem-specific information is assumed to be stored in a structure pointed to by "dat".
 *
 * NOTE: The jacobian (for n=1, m=3, r=3) in matrix form has the following structure:
 *       A_11  0     0     B_11 0    0    0		0		0		0		
 *       0     A_12  0     B_12 0    0    0		0		0		0
 *       0     0     A_13  B_13 0    0    0		0		0		0
 *       A_21  0     0     0    0	 0    0		Pt_21	THor_21	TVer_21
 *       0     A_22  0     0    0	 0    0		Pt_22	THor_22	TVer_22
 *       0     0     A_23  0    0	 0    0		Pt_23	THor_23	TVer_23
 *       A_31  0     0     0    0    0	  0		Pt_31	THor_31	TVer_31
 *       0     A_32  0     0    0    0	  0		Pt_32	THor_32	TVer_32
 *       0     0     A_33  0    0    0	  0		Pt_33	THor_33	TVer_33
 *       A_41  0     0     0    0    0    0		Pt_41	THor_41	TVer_41
 *       0     A_42  0     0    0    0    0		Pt_42	THor_42	TVer_42
 *       0     0     A_43  0    0    0    0		Pt_43	THor_43	TVer_43
 *       To reduce the total number of objective function evaluations, this structure can be
 *       exploited as follows: A certain d is added to the j-th parameters of all cameras and
 *       the objective function is evaluated at the resulting point. This evaluation suffices
 *       to compute the corresponding columns of *all* A_ij through finite differences. A similar
 *       strategy allows the computation of the B_ij, P_ij, THor_ij, TVer_ij. Overall, only cnp+pnp+1 objective function
 *       evaluations are needed to compute the jacobian, much fewer compared to the m*cnp+n*pnp+1
 *       that would be required by the naive strategy of computing one column of the jacobian
 *       per function evaluation. See Nocedal-Wright, ch. 7, pp. 169. Although this approach is
 *       much faster compared to the naive strategy, it is not preferable to analytic jacobians,
 *       since the latter are considerably faster to compute and result in fewer LM iterations.
 */
static void sba_fdjac_grid_x(
							 double *p,                /* I: current parameter estimate, (m*cnp+n*pnp+3*pnp*g)x1 */
							 int n,					   /* n: number of points not on the grid */
							 int nvis,				   /* nvis: number of measurements from points not on the grid */
							 int n_grids,			   /* n_grids: number of grids */
							 int *n_grid_rows,		   /* n_grid_rows: number of rows in the grids*/
							 int *n_grid_cols,		   /* n_grid_cols: number of cols in the grids*/
							 struct sba_crsm *idxij,   /* I: sparse matrix containing the location of x_ij in hx */
							 int    *rcidxs,           /* work array for the indexes of nonzero elements of a single sparse matrix row/column */
							 int    *rcsubs,           /* work array for the subscripts of nonzero elements in a single sparse matrix row/column */
							 double *jac,              /* O: array for storing the approximated jacobian */
							 void   *dat)              /* I: points to a "fdj_data_x_" structure */
{
    fprintf(stderr, "Computing the derivatives...\n");
	
	register int i, j, g, ii, jj;
    double *pa, *pb, *pt, *pqr, *ppt;
    register double *pAB, *phx, *phxx;
    int tn, m, tnm, nnz, Asz, Bsz, Tsz, ABsz, ATsz, idx;
	
    double *tmpd;
    register double d;
	
    struct fdj_data_x_ *fdjd;
    void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
    double *hx, *hxx;
    int cnp, pnp, mnp;
    void *adata;
	
	
    /* retrieve problem-specific information passed in *dat */
    fdjd=(struct fdj_data_x_ *)dat;
    func=fdjd->func;
    cnp=fdjd->cnp; pnp=fdjd->pnp; mnp=fdjd->mnp;
    hx=fdjd->hx;
    hxx=fdjd->hxx;
    adata=fdjd->adata;
	
    tn=idxij->nr; m=idxij->nc; //tn:total no points
    pa=p; pb=p+m*cnp; pt = p+m*cnp + n*pnp;
    Asz=mnp*cnp; Bsz=mnp*pnp; Tsz = mnp*pnp; ABsz=Asz+Bsz; ATsz = Asz + 3*Tsz;
	
    tnm=(tn>=m)? tn : m; // max(n, m);
    tmpd=(double *)emalloc(tnm*sizeof(double));
	
    (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hx, adata); // evaluate supplied function on current solution
	
    if(cnp){ // is motion varying?
        /* compute A_ij */
        for(jj=0; jj<cnp; ++jj){
            for(j=0; j<m; ++j){
                pqr=pa+j*cnp; // j-th camera parameters
                /* determine d=max(SBA_DELTA_SCALE*|pqr[jj]|, SBA_MIN_DELTA), see HZ */
                d=(double)(SBA_DELTA_SCALE)*pqr[jj]; // force evaluation
                d=FABS(d);
                if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
				
                tmpd[j]=d;
                pqr[jj]+=d;
            }
			
            (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
			
            for(j=0; j<m; ++j){
                pqr=pa+j*cnp; // j-th camera parameters
                pqr[jj]-=tmpd[j]; /* restore */
                d=1.0/tmpd[j]; /* invert so that divisions can be carried out faster as multiplications */
				
                nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...tn-1 */
                for(i=0; i<nnz; ++i){
                    idx=idxij->val[rcidxs[i]];
                    phx=hx + idx*mnp; // set phx to point to hx_ij
                    phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    if(idx < nvis)
					{
						pAB=jac + idx*ABsz; // set pAB to point to A_ij
						//fprintf(stderr, "cnp: vis:%d\n", idx);
					}
					else
						pAB=jac + nvis*ABsz + (idx-nvis)*ATsz;
					
                    for(ii=0; ii<mnp; ++ii)
					{
                        pAB[ii*cnp+jj]=(phxx[ii]-phx[ii])*d;
						//printf("jac-A:%f\n", pAB[ii*cnp+jj]);
					}
                }
            }
        }
    }
	printf("derivatives for cameras computed\n...");
	
    if(pnp){ // is structure varying?
        /* compute B_ij */
        for(jj=0; jj<pnp; ++jj){
            for(i=0; i<n; ++i){
                ppt=pb+i*pnp; // i-th point parameters
                /* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
                d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
                d=FABS(d);
                if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
				
                tmpd[i]=d;
                ppt[jj]+=d;
            }
			
            (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
			
            for(i=0; i<n; ++i){
                ppt=pb+i*pnp; // i-th point parameters
                ppt[jj]-=tmpd[i]; /* restore */
                d=1.0/tmpd[i]; /* invert so that divisions can be carried out faster as multiplications */
				
                nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
                for(j=0; j<nnz; ++j){
                    idx=idxij->val[rcidxs[j]];
                    phx=hx + idx*mnp; // set phx to point to hx_ij
                    phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    pAB=jac + idx*ABsz + Asz; // set pAB to point to B_ij
					//fprintf(stderr, "pnp: vis:%d\n", idx);
                    for(ii=0; ii<mnp; ++ii)
					{
                        pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
						//printf("jac-B:%f\n", pAB[ii*pnp+jj]);
					}
                }
            }
        }
    }
	
	printf("derivatives for points computed\n...");
	
	if(n_grids > 0)
	{ 
		// is there a grid repetition?
        /* compute Pt_ij*/
		int prev_points = n;
		for(g=0; g<n_grids; g++)
		{
			ppt=pt + g*pnp*3; // pt parameters
			
			for(jj=0; jj<pnp; ++jj)
			{
				/* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
				d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
				d=FABS(d);
				if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
				
				tmpd[0]=d;
				ppt[jj]+=d;
			
				(*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
			
				ppt[jj]-=tmpd[0]; /* restore */
				d=1.0/tmpd[0]; /* invert so that divisions can be carried out faster as multiplications */
			
				for(i=prev_points; i<prev_points+n_grid_cols[g]*n_grid_rows[g]; i++)
				{
					nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero Pt_ij j=0...m-1 */
					for(j=0; j<nnz; ++j)
					{
						idx=idxij->val[rcidxs[j]];
						phx=hx + idx*mnp; // set phx to point to hx_ij
						phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    
						//fprintf(stderr, "pt:nvis:%d\n", idx);
						pAB=jac + nvis*ABsz + (idx-nvis)*ATsz + Asz; // set pAB to point to P_ij
					
						for(ii=0; ii<mnp; ++ii)
						{
							pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
							//printf("jac-G:%f\n", pAB[ii*pnp+jj]);
						}
					}
				}
			}
		
			/* compute THor_ij*/
			ppt=pt + g*pnp*3 + pnp; // pt parameters
			for(jj=0; jj<pnp; ++jj)
			{
				/* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
				d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
				d=FABS(d);
				if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
			
				tmpd[0]=d;
				ppt[jj]+=d;
			
				(*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
			
				ppt[jj]-=tmpd[0]; /* restore */
				d=1.0/tmpd[0]; /* invert so that divisions can be carried out faster as multiplications */
			
				for(i=prev_points; i<prev_points+n_grid_cols[g]*n_grid_rows[g]; i++)
				{
					nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero Pt_ij j=0...m-1 */
					for(j=0; j<nnz; ++j)
					{
						idx=idxij->val[rcidxs[j]];
						phx=hx + idx*mnp; // set phx to point to hx_ij
						phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    
						//fprintf(stderr, "thor:nvis:%d\n", idx);
					
						pAB=jac + nvis*ABsz + (idx-nvis)*ATsz + Asz + Tsz; // set pAB to point to THor_ij
					
						for(ii=0; ii<mnp; ++ii)
						{
							pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
							//printf("jac-G:%f\n", pAB[ii*pnp+jj]);
						}
					}
				}
			}
		
			/* compute TVer_ij*/
			ppt=pt + g*pnp*3 + pnp + pnp; // pt parameters
			for(jj=0; jj<pnp; ++jj)
			{
				/* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
				d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
				d=FABS(d);
				if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
			
				tmpd[0]=d;
				ppt[jj]+=d;
			
				(*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
			
				ppt[jj]-=tmpd[0]; /* restore */
				d=1.0/tmpd[0]; /* invert so that divisions can be carried out faster as multiplications */
			
				for(i=prev_points; i<prev_points+n_grid_cols[g]*n_grid_rows[g]; i++)
				{
					nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero THor_ij j=0...m-1 */
					for(j=0; j<nnz; ++j)
					{
						idx=idxij->val[rcidxs[j]];
						phx=hx + idx*mnp; // set phx to point to hx_ij
						phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    
						//fprintf(stderr, "pt:tver:%d\n", idx);
					
						pAB=jac + nvis*ABsz + (idx-nvis)*ATsz + Asz + Tsz + Tsz; // set pAB to point to TVer_ij
					
						for(ii=0; ii<mnp; ++ii)
						{
							pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
							//printf("jac-G:%f\n", pAB[ii*pnp+jj]);
						}
					}
				}
			}
			prev_points += n_grid_cols[g]*n_grid_rows[g];
			printf("derivatives for grid %d computed\n...", g);
		}
	}
	
    free(tmpd);
}

static void sba_fdjac_str_grid_x(
							 double *p,                /* I: current parameter estimate, (m*cnp+n*pnp+3*pnp*g)x1 */
							 int n,					   /* n: number of points not on the grid */
							 int nvis,				   /* nvis: number of measurements from points not on the grid */
							 int n_grids,			   /* n_grids: number of grids */
							 int *n_grid_rows,		   /* n_grid_rows: number of rows in the grids*/
							 int *n_grid_cols,		   /* n_grid_cols: number of cols in the grids*/
							 struct sba_crsm *idxij,   /* I: sparse matrix containing the location of x_ij in hx */
							 int    *rcidxs,           /* work array for the indexes of nonzero elements of a single sparse matrix row/column */
							 int    *rcsubs,           /* work array for the subscripts of nonzero elements in a single sparse matrix row/column */
							 double *jac,              /* O: array for storing the approximated jacobian */
							 void   *dat)              /* I: points to a "fdj_data_x_" structure */
{
    fprintf(stderr, "Computing the derivatives...\n");
	
	register int i, j, g, ii, jj;
    double *pa, *pb, *pt, *pqr, *ppt;
	double *pAB;
    register double *phx, *phxx;
    int tn, m, tnm, nnz, Bsz, Tsz,  idx;
	
    double *tmpd;
    register double d;
	
    struct fdj_data_x_ *fdjd;
    void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
    double *hx, *hxx;
    int pnp, mnp;
    void *adata;
	
	
    /* retrieve problem-specific information passed in *dat */
    fdjd=(struct fdj_data_x_ *)dat;
    func=fdjd->func;
    pnp=fdjd->pnp; mnp=fdjd->mnp;
    hx=fdjd->hx;
    hxx=fdjd->hxx;
    adata=fdjd->adata;
	
    tn=idxij->nr; m=idxij->nc; //tn:total no points
    pb=p; pt = p+ n*pnp;
    Bsz=mnp*pnp; Tsz = mnp*pnp;
	
    tnm=(tn>=m)? tn : m; // max(n, m);
    tmpd=(double *)emalloc(tnm*sizeof(double));
	
    (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hx, adata); // evaluate supplied function on current solution
	
	printf("compute derivatives for points....\n");
    if(pnp){ // is structure varying?
        /* compute B_ij */
        for(jj=0; jj<pnp; ++jj){
            for(i=0; i<n; ++i){
                ppt=pb+i*pnp; // i-th point parameters
                /* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
                d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
                d=FABS(d);
                if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
				
                tmpd[i]=d;
                ppt[jj]+=d;
            }
			
            (*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
			
            for(i=0; i<n; ++i){
                ppt=pb+i*pnp; // i-th point parameters
                ppt[jj]-=tmpd[i]; /* restore */
                d=1.0/tmpd[i]; /* invert so that divisions can be carried out faster as multiplications */
				
                nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
                for(j=0; j<nnz; ++j){
                    idx=idxij->val[rcidxs[j]];
                    phx=hx + idx*mnp; // set phx to point to hx_ij
                    phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
                    pAB=jac + idx*Bsz; // set pAB to point to B_ij
					//fprintf(stderr, "pnp: vis:%d\n", idx);
                    for(ii=0; ii<mnp; ++ii)
					{
                        pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
						//printf("jac-B:%f\n", pAB[ii*pnp+jj]);
					}
                }
            }
        }
    }
	
	if(n_grids > 0)
	{ 
		// is there a grid repetition?
        /* compute Pt_ij*/
		//printf("*******Pt******\n");
		int prev_points = n;
		for(g=0; g<n_grids; g++)
		{
			ppt=pt + g*pnp*3; // pt parameters
			
			for(jj=0; jj<pnp; ++jj)
			{
				/* determine d=max(SBA_DELTA_SCALE*|ppt[jj]|, SBA_MIN_DELTA), see HZ */
				d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
				d=FABS(d);
				if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
				
				tmpd[0]=d;
				ppt[jj]+=d;
				
				(*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
				
				ppt[jj]-=tmpd[0]; /* restore */
				d=1.0/tmpd[0]; /* invert so that divisions can be carried out faster as multiplications */
				
				for(i=prev_points; i<prev_points+n_grid_cols[g]*n_grid_rows[g]; i++)
				{
					nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); /* find nonzero Pt_ij j=0...m-1 */
					for(j=0; j<nnz; ++j)
					{
						idx=idxij->val[rcidxs[j]];
						phx=hx + idx*mnp; // set phx to point to hx_ij
						phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
						
						//fprintf(stderr, "pt:nvis:%d\n", idx);
						pAB=jac + nvis*Bsz + (idx-nvis)*Tsz*3; // set pAB to point to P_ij
						
						for(ii=0; ii<mnp; ++ii)
						{
							//printf("hx[%d]:%f, hxx[%d]:%f\n", idx*mnp+ii, phx[ii], idx*mnp+ii, phxx[ii]);
							pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
							//printf("jac-G:%f\n", (phxx[ii]-phx[ii])*d);
						}
					}
				}
			}
			
			/* compute THor_ij*/
			//printf("*******TH******\n");
			ppt=pt + g*pnp*3 + pnp; // pt parameters
			for(jj=0; jj<pnp; ++jj)
			{
				d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
				d=FABS(d);
				if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
				
				tmpd[0]=d;
				ppt[jj]+=d;
				
				(*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
				
				ppt[jj]-=tmpd[0]; // restore 
				d=1.0/tmpd[0]; // invert so that divisions can be carried out faster as multiplications 
				
				for(i=prev_points; i<prev_points+n_grid_cols[g]*n_grid_rows[g]; i++)
				{
					nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); // find nonzero Pt_ij j=0...m-1 
					for(j=0; j<nnz; ++j)
					{
						idx=idxij->val[rcidxs[j]];
						//printf("idx:%d\n", idx);
						phx=hx + idx*mnp; // set phx to point to hx_ij
						phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
						
						//fprintf(stderr, "thor:nvis:%d\n", idx);
						
						pAB=jac + nvis*Bsz + (idx-nvis)*Tsz*3 + Tsz; // set pAB to point to THor_ij
						
						for(ii=0; ii<mnp; ++ii)
						{
							//printf("hx[%d]:%f, hxx[%d]:%f\n", idx*mnp+ii, phx[ii], idx*mnp+ii, phxx[ii]);
							pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
							//printf("jac-G:%f\n", (phxx[ii]-phx[ii])*d);
						}
					}
				}
			}
			
			/* compute TVer_ij*/
			//printf("*******TV******\n");
			ppt=pt + g*pnp*3 + pnp + pnp; // pt parameters
			for(jj=0; jj<pnp; ++jj)
			{
				d=(double)(SBA_DELTA_SCALE)*ppt[jj]; // force evaluation
				d=FABS(d);
				if(d<SBA_MIN_DELTA) d=SBA_MIN_DELTA;
				
				tmpd[0]=d;
				ppt[jj]+=d;
				
				(*func)(p, idxij, fdjd->func_rcidxs, fdjd->func_rcsubs, hxx, adata);
				
				ppt[jj]-=tmpd[0]; // restore 
				d=1.0/tmpd[0]; // invert so that divisions can be carried out faster as multiplications 
				
				for(i=prev_points; i<prev_points+n_grid_cols[g]*n_grid_rows[g]; i++)
				{
					nnz=sba_crsm_row_elmidxs(idxij, i, rcidxs, rcsubs); // find nonzero THor_ij j=0...m-1 
					for(j=0; j<nnz; ++j)
					{
						idx=idxij->val[rcidxs[j]];
						//printf("idx:%d\n", idx);
						phx=hx + idx*mnp; // set phx to point to hx_ij
						phxx=hxx + idx*mnp; // set phxx to point to hxx_ij
						//printf("phx:%f %f, phxx:%f %f\n", phx[0], phx[1], phxx[0], phxx[1]);
						//fprintf(stderr, "pt:tver:%d\n", idx);
						
						pAB=jac + nvis*Bsz + (idx-nvis)*Tsz*3 + Tsz + Tsz; // set pAB to point to TVer_ij
						
						for(ii=0; ii<mnp; ++ii)
						{
							//printf("hx[%d]:%f, hxx[%d]:%f\n", idx*mnp+ii, phx[ii], idx*mnp+ii, phxx[ii]);
							pAB[ii*pnp+jj]=(phxx[ii]-phx[ii])*d;
							//printf("jac-G:%f\n", (phxx[ii]-phx[ii])*d);
						}
					}
				}
			}
			prev_points += n_grid_cols[g]*n_grid_rows[g];
			//printf("derivatives for grid %d computed\n...", g);
		}
	}
	
    free(tmpd);
}


typedef int (*PLS)(double *A, double *B, double *x, int m, int iscolmaj);

/* Bundle adjustment on camera and structure parameters 
 * with grid based repetitions using the sparse 
 * Levenberg-Marquardt as described in HZ p. 568
 *
 * Returns the number of iterations (>=0) if successfull, SBA_ERROR if failed
 */

int sba_motstr_levmar_grid_x(const int n,			/*number of non-grid points*/
							 const int noGrids,		/* number of grids*/
							 const int* grid_rows, /*number of rows in each grid*/
							 const int* grid_cols, /*number of cols in each grid*/
							 const int m,   /* number of images */
							 const int mcon,/* number of images (starting from the 1st) whose parameters should not be modified.
											 * All A_ij (see below) with j<mcon are assumed to be zero
											 */
							 char *vmask,  /* visibility mask: vmask[i, j]=1 if point i visible in image j, 0 otherwise. (n+rows1*cols1+rows2*cols2...)xm */
							 double *p,    /* initial parameter vector p0: (a1, ..., am, b1, ..., bn, pt1, thor1, tver1, pt2, thor2, tver2...).
											* aj are the image j parameters, bi are the i-th point parameters, pt is the grid
											* representative point, thor is the grid horizontal translation, tver is the grid
											* vertical translation
											* size m*cnp + n*pnp + pnp*3*noGrids
											*/
							 const int cnp,/* number of parameters for ONE camera; e.g. 6 for Euclidean cameras */
							 const int pnp,/* number of parameters for ONE point; e.g. 3 for Euclidean points */
							 double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
											* x_ij is the projection of the i-th point on the j-th image.
											* NOTE: some of the x_ij might be missing, if point i is not visible in image j;
											* see vmask[i, j], max. size (n+grid_rows*grid_cols)*m*mnp
											*/
							 double *covx, /* measurements covariance matrices: (Sigma_x_11, .. Sigma_x_1m, ..., Sigma_x_(n+grid_rows*grid_cols)1, .. Sigma_x_(n+grid_rows*grid_cols)m),
											* where Sigma_x_ij is the mnp x mnp covariance of x_ij stored row-by-row. Set to NULL if no
											* covariance estimates are available (identity matrices are implicitly used in this case).
											* NOTE: a certain Sigma_x_ij is missing if the corresponding x_ij is also missing;
											* see vmask[i, j], max. size (n+grid_rows*grid_cols)*m*mnp*mnp
											*/
							 const int mnp,/* number of parameters for EACH measurement; usually 2 */
							 void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
							 /* functional relation describing measurements. Given a parameter vector p,
							  * computes a prediction of the measurements \hat{x}. p is (m*cnp + n*pnp + pnp*3)x1,
							  * \hat{x} is ((n+grid_rows*grid_cols)*m*mnp)x1, maximum
							  * rcidxs, rcsubs are max(m, n+grid_rows*grid_cols) x 1, allocated by the caller and can be used
							  * as working memory
							  */
							 void (*fjac)(double *p, int n, int nvis, int n_grid, int *grid_rows, int *grid_cols, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
							 /* function to evaluate the sparse jacobian dX/dp.
							  * The Jacobian is returned in jac as
							  * (dx_11/da_1, ..., dx_1m/da_m, ..., dx_(n+grid_rows*grid_cols)1/da_1, ..., dx_(n+grid_rows*grid_cols)m/da_m,
							  *  dx_11/db_1, ..., dx_1m/db_1, ..., dx_(n+grid_rows*grid_cols)1/db_n, ..., dx_(n+grid_rows*grid_cols)m/db_n,
							  *  dx_11/d_pt, ..., dx_1m/dpt, ..., dx_(n+grid_rows*grid_cols)1/d_pt, ..., dx_(n+grid_rows*grid_cols)m/d_pt
							  *  dx_11/d_thor, ..., dx_1m/dthor, ..., dx_(n+grid_rows*grid_cols)1/dthor, ..., dx_(n+grid_rows*grid_cols)m/dthor,
							  *  dx_11/d_tver, ..., dx_1m/dtver, ..., dx_(n+grid_rows*grid_cols)1/dtver, ..., dx_(n+grid_rows*grid_cols)m/dtver) or (using HZ's notation),
							  * jac=(A_11, B_11, Pt_11, THor_11, TVer_11, ..., A_1m, B_1m, Pt_1m, Thor_1m, Tver_1m, ..., A_n1, B_n1, Pt_n1, Thor_n1, Tver_n1, ..., A_nm, B_nm, Pt_nm, THor_nm, TVer_nm)
							  * Notice that depending on idxij, some of the A_ij and B_ij might be missing.
							  * Note also that A_ij and B_ij are mnp x cnp and mnp x pnp matrices resp. and they
							  * should be stored in jac in row-major order.
							  * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
							  * as working memory
							  *
							  * If NULL, the jacobian is approximated by repetitive func calls and finite
							  * differences. This is computationally inefficient and thus NOT recommended.
							  */
							 void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 
							 const int itmax,   /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
							 const int verbose, /* I: verbosity */
							 const double opts[SBA_OPTSSZ],
							 /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \epsilon4]. Respectively the scale factor for initial \mu,
							  * stopping thresholds for ||J^T e||_inf, ||dp||_2, ||e||_2 and (||e||_2-||e_new||_2)/||e||_2
							  */
							 double info[SBA_INFOSZ],
							 /* O: information regarding the minimization. Set to NULL if don't care
							  * info[0]=||e||_2 at initial p.
							  * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
							  * info[5]= # iterations,
							  * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
							  *                                 2 - stopped by small dp
							  *                                 3 - stopped by itmax
							  *                                 4 - stopped by small relative reduction in ||e||_2
							  *                                 5 - stopped by small ||e||_2
							  *                                 6 - too many attempts to increase damping. Restart with increased mu
							  *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
							  * info[7]= # function evaluations
							  * info[8]= # jacobian evaluations
							  * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
							  */
							 int use_constraints, camera_constraints_t *constraints,  /* Constraints on camera parameters */
							 int use_point_constraints, point_constraints_t *point_constraints
							 )
{
	
	register int i, j, ii, jj, k, l, gridIt;
    int nvis, nnz, retval;
	int *nVisGrids;
	
    /* The following are work arrays that are dynamically allocated by sba_motstr_levmar_x() */
    double *jac;  /* work array for storing the jacobian, max. size (n+rows1*cols1+rows2*cols2...)*m*(mnp*cnp) + m*(n+(rows1*cols1+rows2*cols2...)*3)*(mnp*pnp) */
    
	double *U;    /* work array for storing the U_j in the order U_1, ..., U_m, size m*cnp*cnp */
    double *V;    /* work array for storing the *strictly upper triangles* of V_i in the order V_1, ..., V_n, size n*pnp*pnp.
                   * V also stores the lower triangles of (V*_i)^-1 in the order (V*_1)^-1, ..., (V*_n)^-1.
                   * Note that diagonal elements of V_1 are saved in diagUV
                   */
	/* Notice that the blocks W_ij, Y_ij are zero iff A_ij (equivalently B_ij) is zero. This means
     * that the matrix consisting of blocks W_ij is itself sparse, similarly to the
     * block matrix made up of the A_ij and B_ij (i.e. jac)
     */
    double *W;    /* work array for storing the W_ij in the order W_11, ..., W_1m, ..., W_n1, ..., W_nm,
				   max. size n*m*cnp*pnp */
	double *M;	  /*size m*cnp*3*pnp*noGrids*/
	double *N;	  /*size 3*pnp*3*pnp*noGrids*/
	double *Ninv; /*size 3*pnp*3*pnp*noGrids*/
	double *Yj;   /* work array for storing the Y_ij for a *fixed* j in the order Y_1j, Y_nj,
				   max. size n*cnp*pnp */
	double *YWt;  /* work array for storing \sum_i Y_ij W_ik^T, size cnp*cnp */
	double *Zj;   /* work array for storing the Z_ij for a *fixed* j in the order Z_1j, Z_nj,
				   max. size 3*cnp*pnp*noGrids */
	double *ZMt;  /* work array for storing \sum_i Z_ij M_ik^T, size cnp*cnp*noGrids */
	double *S;    /* work array for storing the block array S_jk, size m*m*cnp*cnp */
	double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_(n+grid_rows*grid_cols)1, ..., e_(n+grid_rows*grid_cols)m,
				   max. size (n+grid_rows*grid_cols)*m*mnp */
    double *eabtrans;  /* work array for storing the ea_j & eb_i & e_pt & e_thor & e_tver in the order ea_1, .. ea_m eb_1, .. eb_n size m*cnp + n*pnp + 3*pnp*noGrids */
	double *dp;   /* work array for storing the parameter vector updates da_1, ..., da_m, db_1, ..., db_n,d_p,d_thor,d_tver size m*cnp + n*pnp + 3*pnp*noGrids*/
	double *Wtda; /* work array for storing \sum_j W_ij^T da_j, size pnp */
    double *Mtda; /* work array for storing \sum_j M_ij^T da_j, size pnp*noGrids */
	double *wght= NULL; /* work array for storing the weights computed from the covariance inverses, max. size (n+rows1*cols1+rows2*cols2...)*m*mnp*mnp */
	double *E;   /* work array for storing the e_j in the order e_1, .. e_m, size m*cnp */
	
    /* Of the above arrays, jac, e, W, Yj, wght are sparse and
     * U, V, eabtrans, E, S, dp are dense. Sparse arrays (except Yj) are indexed
     * through idxij (see below), that is with the same mechanism as the input 
     * measurements vector x
     */
	
    double *pa, *pb, *pt, *ea, *eb, *et, *dpa, *dpb, *dpt; /* pointers into p, jac, eab and dp respectively */
	
    /* submatrices sizes */
    int Asz, Bsz, Ptsz, THorsz, TVersz, ABsz, ATsz, Usz, Vsz,
	Wsz, Msz, Nsz, Ysz, Zsz, esz, easz, ebsz, etranssz,
	YWtsz, ZMtsz, Wtdasz, Mtdasz, Sblsz, covsz;
	
    int Sdim; /* S matrix actual dimension */
	
    register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
	register double *ptr1_0, *ptr1_1, *ptr1_2;
	register double *ptr2_0, *ptr2_1, *ptr2_2;
	register double *ptr3_0, *ptr3_1, *ptr3_2;
	register double *ptr1_psum, *ptr1_THorP, *ptr1_TVerP, *ptr1_PTHor, *ptr1_hsum, *ptr1_TVerTHor, *ptr1_PTVer, *ptr1_THorTVer, *ptr1_vsum; 
	register double *ptr3_pt, *ptr3_thor, *ptr3_tver;
	
    struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also
                            * the location of A_ij, B_ij in jac, etc.
                            * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                            * index k and it is used to efficiently lookup the memory locations where the non-zero
                            * blocks of a sparse matrix/vector are stored
                            */
    int maxCvis, /* max. of projections of a single point  across cameras, <=m */
	maxPvis, /* max. of projections in a single camera across points,  <=(n+rows1*cols1+rows2*cols2...) */
	maxCPvis, /* max. of the above */
	*rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
			   column in a sparse matrix, size max((n+grid_rows*grid_cols), m) */
	*rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
			   sparse matrix, size max((n+grid_rows*grid_cols), m) */
	
    /* The following variables are needed by the LM algorithm */
    register int itno;  /* iteration counter */
    int issolved;
    /* temporary work arrays that are dynamically allocated */
    double *hx,         /* \hat{x}_i, max. size m*(n+grid_rows*grid_cols)*mnp */
	*diagUVN,     /* diagonals of U_j, V_i, N_i size m*cnp + n*pnp + 3*pnp*noGrids */
	*pdp;        /* p + dp, size m*cnp + n*pnp + 3*pnp*/
	
    double *diagU, *diagV, *diagN; /* pointers into diagUVN */
	
    register double mu,  /* damping constant */
	tmp; /* mainly used in matrix & vector multiplications */
    double p_eL2, eabtrans_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
    double p_L2, dp_L2=DBL_MAX, dF, dL;
    double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2],
	eps3_sq=opts[3]*opts[3], eps4_sq=opts[4]*opts[4], eps5=opts[5];
    double init_p_eL2;
    int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;
    int nobs, nvars;
    const int mmcon=m-mcon;
    PLS linsolver=NULL;
    int (*matinv)(double *A, int m)=NULL;
	
    struct fdj_data_x_ fdj_data;
    void *jac_adata;
	
	float gridWeight = 1.0;
	
    /* Initialization */
	
#ifdef TIMINGS
    clock_t start = clock();
#endif
	
    mu=eabtrans_inf=0.0; /* -Wall */
	
    /* block sizes */
    Asz=mnp * cnp; Bsz=mnp * pnp; Ptsz = mnp * pnp; THorsz = mnp * pnp; TVersz = mnp * pnp; 
	ABsz=Asz + Bsz;
	ATsz =Asz + Ptsz + THorsz + TVersz; 
    Usz=cnp * cnp; Vsz=pnp * pnp;
    Wsz=cnp * pnp; Ysz=cnp * pnp;
	Msz=cnp * pnp; Nsz=pnp * pnp; 
	Zsz = cnp * pnp;
	esz=mnp;
    easz=cnp; ebsz=pnp; etranssz = pnp;
    YWtsz=cnp * cnp;
	ZMtsz = cnp * cnp;
    Wtdasz=pnp;
	Mtdasz=pnp;
    Sblsz=cnp * cnp;
    Sdim=mmcon * cnp;
    covsz=mnp * mnp;
	
	int totalNoPoints = n;
	for(i=0; i<noGrids; i++)
	{
		totalNoPoints += grid_cols[i]*grid_rows[i];
	}
	
    /* count total number of visible image points */
    int totalNVis = 0;
	nvis = 0;
	nVisGrids = (int*) malloc(sizeof(int)*noGrids);
	for(i=0; i<noGrids; i++)
	{
		nVisGrids[i] = 0;
	}
	
	for(i=nvis=0, jj=totalNoPoints*m; i<jj; ++i)
	{
		totalNVis += (vmask[i]!=0);
		
        if(i<n*m)
		{
			nvis += (vmask[i]!=0);
			continue;
		}
		int prevPoints = n*m;
		for(j=0; j<noGrids; j++)
		{
			prevPoints += grid_cols[j]*grid_rows[j]*m;
			if(i<prevPoints)
			{
				nVisGrids[j] += (vmask[i]!=0);
				break;
			}
		}
	}
	
	for(j=0; j<noGrids; j++)
	{
		printf("vis for grid %d: %d\n", j, nVisGrids[j]);
	}
	
	nobs=totalNVis*mnp;
    nvars=m*cnp + n*pnp + 3*pnp*noGrids;
	
	if(nobs<nvars){
        fprintf(stderr, "SBA: sba_motstr_levmar_grid_x() cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
        return SBA_ERROR;
    }
	
	fprintf(stderr, "SBA: sba_motstr_levmar_grid_x() has [%d] measurements [%d] unknowns\n", nobs, nvars);
	
    /* allocate & fill up the idxij structure */
    sba_crsm_alloc(&idxij, totalNoPoints, m, totalNVis);
    for(i=k=0; i<totalNoPoints; ++i){
        idxij.rowptr[i]=k;
        ii=i*m;
        for(j=0; j<m; ++j)
            if(vmask[ii+j]){
                idxij.val[k]=k;
                idxij.colidx[k++]=j;
            }
    }
    idxij.rowptr[totalNoPoints]=totalNVis;
	
    /* find the maximum number (for all cameras) of visible image projections coming from a single 3D point */
    for(i=maxCvis=0; i<totalNoPoints; ++i)
        if((k=idxij.rowptr[i+1]-idxij.rowptr[i])>maxCvis) maxCvis=k;
	
    /* find the maximum number (for all points) of visible image projections in any single camera */
    for(j=maxPvis=0; j<m; ++j){
        for(i=ii=0; i<totalNoPoints; ++i)
            if(vmask[i*m+j]) ++ii;
        if(ii>maxPvis) maxPvis=ii;
    }
    maxCPvis=(maxCvis>=maxPvis)? maxCvis : maxPvis;
	
	fprintf(stderr, "maxCvis:%d, maxPvis:%d\n", maxCvis, maxPvis);
	
	jac=(double *)emalloc(((nvis*ABsz)+((totalNVis-nvis)*ATsz))*sizeof(double));
	W=(double *)emalloc((nvis*Wsz)*sizeof(double));
	U=(double *)emalloc(m*Usz*sizeof(double));
	V=(double *)emalloc(n*Vsz*sizeof(double));
    M = (double *)emalloc(noGrids*3*m*Msz*sizeof(double));
	N = (double *)emalloc(noGrids*3*3*Nsz*sizeof(double));
	Ninv = (double *)emalloc(noGrids*3*3*Nsz*sizeof(double));
	e=(double *)emalloc(nobs*sizeof(double));
    eabtrans=(double *)emalloc(nvars*sizeof(double));
    E=(double *)emalloc(m*cnp*sizeof(double));
    Yj=(double *)emalloc(maxPvis*Ysz*sizeof(double));
    YWt=(double *)emalloc(YWtsz*sizeof(double));
    Zj = (double *)emalloc(noGrids*3*Zsz*sizeof(double));
	ZMt = (double *)emalloc(noGrids*ZMtsz*sizeof(double));
	S=(double *)emalloc(m*m*Sblsz*sizeof(double));
    dp=(double *)emalloc(nvars*sizeof(double));
    Wtda=(double *)emalloc(pnp*sizeof(double));
	Mtda=(double *)emalloc(noGrids*pnp*3*sizeof(double));
    rcidxs=(int *)emalloc(maxCPvis*sizeof(int));
    rcsubs=(int *)emalloc(maxCPvis*sizeof(int));
	
	covx = (double*)emalloc(totalNVis*covsz*sizeof(double));
	for(i=0; i<nvis; i++)
	{
		for(j=0; j<mnp; j++ )
		{
			for(k=0; k<mnp; k++)
			{
				if(j==k)
					covx[i*mnp*mnp + j*mnp + k] = 1.0;
				else
					covx[i*mnp*mnp + j*mnp + k] = 0.0;
			}
		}
	}
	for(i=nvis; i<totalNVis; i++)
	{
		for(j=0; j<mnp; j++ )
		{
			for(k=0; k<mnp; k++)
			{
				if(j==k)
					covx[i*mnp*mnp + j*mnp + k] = 1.0 / gridWeight;
				else
					covx[i*mnp*mnp + j*mnp + k] = 0.0;
			}
		}
	}
	
#ifndef SBA_DESTROY_COVS
    if(covx!=NULL) wght=(double *)emalloc(totalNVis*covsz*sizeof(double));
#else
    if(covx!=NULL) wght=covx;
#endif /* SBA_DESTROY_COVS */
	
	
    hx=(double *)emalloc(nobs*sizeof(double));
    diagUVN = (double *)emalloc(nvars*sizeof(double));
	
    pdp=(double *)emalloc(nvars*sizeof(double));
	
#ifdef TIMINGS
    clock_t end = clock();
    printf("[sba_motstr_levmar_x] initialization took %0.3fs\n", 
           (end - start) / (float) CLOCKS_PER_SEC);
#endif
	
    /* set up auxiliary pointers */
    pa=p; pb=p+m*cnp; pt=p+m*cnp+n*pnp;
    ea=eabtrans; eb=eabtrans+m*cnp; et=eabtrans+m*cnp+n*pnp;
    dpa=dp; dpb=dp+m*cnp; dpt=dp+m*cnp+n*pnp;
	
    diagU=diagUVN; diagV=diagUVN + m*cnp; diagN=diagUVN + m*cnp + n*pnp;
	
    /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
    if(!fjac){
        fdj_data.func=func;
        fdj_data.cnp=cnp;
        fdj_data.pnp=pnp;
        fdj_data.mnp=mnp;
        fdj_data.hx=hx;
        fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
        fdj_data.func_rcidxs=(int *)emalloc(2*maxCPvis*sizeof(int));
        fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxCPvis;
        fdj_data.adata=adata;
		
        fjac=sba_fdjac_grid_x;
        jac_adata=(void *)&fdj_data;
    }
    else{
        fdj_data.hxx=NULL;
        jac_adata=adata;
    }
	
	if(itmax==0){ /* verify jacobian */
        sba_motstr_chkjac_grid_x(func, fjac, p, &idxij, rcidxs, rcsubs, n, nvis, noGrids, grid_rows, grid_cols, mcon, cnp, pnp, mnp, adata, jac_adata);
		retval=0;
        goto freemem_and_return;
    }
	
    /* compute the error vectors e_ij in hx */
	//for(i=0; i<nvars; i++)
	//{
	//	printf("p[%d]:%f\n", i, p[i]);
	//}
	
	if(covx!=NULL){
        for(i=0; i<totalNoPoints; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1, ptr2 to point to cov_x_ij, w_x_ij resp. */
                ptr1=covx + (k=idxij.val[rcidxs[j]]*covsz);
                ptr2=wght + k;
                if(!sba_mat_cholinv(ptr1, ptr2, mnp)){ /* compute w_x_ij s.t. w_x_ij^T w_x_ij = cov_x_ij^-1 */
                    fprintf(stderr, "SBA: invalid covariance matrix for x_ij (i=%d, j=%d) in sba_motstr_levmar_x()\n", i, rcsubs[j]);
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }
            }
        }
        sba_mat_cholinv(NULL, NULL, 0); /* cleanup */
    }
	
    /* compute the error vectors e_ij in hx */
    (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
    /* ### compute e=x - f(p) [e=w*(x - f(p)] and its L2 norm */
    if(covx==NULL)
		p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
    else
		p_eL2=nrmCxmy(e, x, hx, wght, mnp, totalNVis); /* e=wght*(x-hx), p_eL2=||e||=||x-hx||_Sigma^-1 */
	
	printf("c\n");
    /* Add in the camera constraints */
    if (use_constraints) {
        for (j = 0; j < m; j++) {
            /* Find number of points projecting into this camera */
            int nnz = 0;
			
            for(i = 0; i < totalNoPoints; i++)
                nnz += vmask[i * m + j];
			
            for (jj = 0; jj < cnp; jj++) {
                if (constraints[j].constrained[jj]) {
                    double diff = constraints[j].constraints[jj] - p[j * cnp + jj];
					
                    p_eL2 += constraints[j].weights[jj] * diff * diff;
                }
            }
        }
    }
	
    /* Add in the point constraints */
    if (use_point_constraints) {
        for (i = 0; i < totalNoPoints; i++) {
            if (point_constraints[i].constrained) {
                for (ii=0; ii < pnp; ii++) {
					double diff;
					if(i<n)
					{
						diff = point_constraints[i].constraints[ii] - p[m * cnp + i * pnp + ii];
					}
                    else
					{
						int noCurPoints = n;
						for(j=0; j<noGrids; j++)
						{
							noCurPoints += grid_rows[j]*grid_cols[j];
							
							if(i < noCurPoints)
							{
								int noPrevPoints = n;
								for(jj=0; jj<j; jj++)
									noPrevPoints += grid_cols[jj]*grid_rows[jj];
								
								int row = (i-noPrevPoints)/grid_cols[j];
								int col = (i-noPrevPoints)%grid_cols[j];
								
								diff = point_constraints[i].constraints[ii] - (p[m*cnp+n*pnp+j*3*3+ii]+
																			   p[m*cnp+n*pnp+j*3*3+ii+3]*col+
																			   p[m*cnp+n*pnp+j*3*3+ii+6]*row);
								break;
							}
						}
					}
					
                    p_eL2 += totalNVis* point_constraints[i].weight * diff * diff;
                }
            }
        }
    }
	
    if(verbose) printf("initial motstr-SBA error %g [%g]\n", p_eL2, p_eL2/totalNVis);
    init_p_eL2=p_eL2;
    if(!SBA_FINITE(p_eL2)) stop=7;
	
    for(itno=0; itno<itmax && !stop; ++itno){
        /* Note that p, e and ||e||_2 have been updated at the previous iteration */
		
#ifdef TIMINGS
        start = clock();
#endif
		
        /* compute derivative submatrices A_ij, B_ij, Pt_ij, THor_ij, TVer_ij */
        (*fjac)(p, n, nvis, noGrids, grid_rows, grid_cols, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;
		
#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing derivs took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif
		
		for(i=nvis; i<totalNVis; ++i)
		{
			ptr1=wght + i*covsz;
		 
			//A
			ptr2=jac  + nvis*ABsz + (i-nvis)*ATsz;
			for(ii=0; ii<mnp; ++ii)
			{
				for(jj=0; jj<cnp; ++jj){
					for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
						sum+=ptr1[ii*mnp+k]*ptr2[k*cnp+jj];
					ptr2[ii*cnp+jj]=sum;
				}
			}
		 
			//P
			ptr2=jac  + nvis*ABsz + (i-nvis)*ATsz + Asz;
			for(ii=0; ii<mnp; ++ii)
			{
				for(jj=0; jj<pnp; ++jj){
					for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
						sum+=ptr1[ii*mnp+k]*ptr2[k*pnp+jj];
					ptr2[ii*pnp+jj]=sum;
				}
			}
		 
			//H
			ptr2=jac  + nvis*ABsz + (i-nvis)*ATsz + Asz + Ptsz;
			for(ii=0; ii<mnp; ++ii)
			{
				for(jj=0; jj<pnp; ++jj){
					for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
						sum+=ptr1[ii*mnp+k]*ptr2[k*pnp+jj];
					ptr2[ii*pnp+jj]=sum;
				}
			}
		 
			//V
			ptr2=jac  + nvis*ABsz + (i-nvis)*ATsz + Asz + Ptsz + THorsz;
			for(ii=0; ii<mnp; ++ii)
			{
				for(jj=0; jj<pnp; ++jj){
					for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
						sum+=ptr1[ii*mnp+k]*ptr2[k*pnp+jj];
					ptr2[ii*pnp+jj]=sum;
				}
			}
		 
		}
		
		
		
		/* compute U_j = \sum_i A_ij^T A_ij  // \Sigma here!
         * U_j is symmetric, therefore its computation can be sped up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that A_ij is mnp x cnp
         *
         * Also compute ea_j = \sum_i A_ij^T e_ij  // \Sigma here!
         * Recall that e_ij is mnp x 1
         */
		
#ifdef TIMINGS
        start = clock();
#endif
		
        _dblzero(U, m*Usz); /* clear all U_j */
        _dblzero(ea, m*easz); /* clear all ea_j */
        for(j=mcon; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            ptr2=ea + j*easz; // set ptr2 to point to ea_j
			
            nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
            for(i=0; i<nnz; ++i){
                /* set ptr3 to point to A_ij, actual row number in rcsubs[i] */
                int r = idxij.val[rcidxs[i]];
				if(r<nvis)
					ptr3=jac + r*ABsz;
				else
					ptr3=jac + nvis*ABsz + (r-nvis)*ATsz;
				
                /* compute the UPPER TRIANGULAR PART of A_ij^T A_ij and add it to U_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=ii; jj<cnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*cnp+ii]*ptr3[k*cnp+jj];
                        ptr1[ii*cnp+jj]+=sum;
                    }
					
                    /* copy the LOWER TRIANGULAR PART of U_j from the upper one */
                    for(jj=0; jj<ii; ++jj)
                        ptr1[ii*cnp+jj]=ptr1[jj*cnp+ii];
                }
				
                ptr4=e + idxij.val[rcidxs[i]]*esz; /* set ptr4 to point to e_ij */
                /* compute A_ij^T e_ij and add it to ea_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*cnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }
			
            /* Add in the camera constraints */
            if (use_constraints) {
                for (jj=0; jj < cnp; jj++) {
                    if (constraints[j].constrained[jj]) {
                        double diff = constraints[j].constraints[jj] - p[j * cnp + jj];
                        /* Add to the U matrix */
                        ptr1[jj*cnp+jj] += /*nnz **/ constraints[j].weights[jj];
                        ptr2[jj] += /*nnz **/ constraints[j].weights[jj] * diff;
                    }
                }
            }
			
			/*printf("***********U_%d***********\n", j);
			 for(ii=0; ii<cnp; ii++)
			 {
			 for(jj=0; jj<cnp; jj++)
			 {
			 printf("%f ", ptr1[ii*cnp+jj]);
			 }
			 printf("\n");
			 }*/
        }
        
#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing U took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif
		
		/* compute M_j = [\sum_i A_ij^T Pt_ij \sum_i A_ij^T THor_ij \sum_i A_ij^T TVer_ij
         * Recall that A_ij is mnp x cnp and Pt_ij, THor_ij, TVer_ij is mnp x pnp
         */
		
#ifdef TIMINGS
        start = clock();
#endif
		
        _dblzero(M, noGrids*m*3*Msz); /* clear all M_j */
        for(gridIt = 0; gridIt<noGrids; gridIt++)
		{
			int noPrevVisPoints = nvis;
			
			for(j=0; j<gridIt; j++)
				noPrevVisPoints += nVisGrids[j];
			
			int noCurVisPoints = noPrevVisPoints + nVisGrids[gridIt];
			
			for(j=mcon; j<m; ++j)
			{
				ptr1_0=M + gridIt*3*m*Msz + j*3*Msz; // set ptr1 to point to \sum_i A_ij^T Pt_ij in M_j
				ptr1_1=M + gridIt*3*m*Msz + j*3*Msz + Msz; // set ptr1 to point to \sum_i A_ij^T THor_ij in M_j
				ptr1_2=M + gridIt*3*m*Msz + j*3*Msz + 2*Msz; // set ptr1 to point to \sum_i A_ij^T TVer_ij in M_j
				
				nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
				for(i=0; i<nnz; ++i)
				{
					/* set ptr3 to point to A_ij, actual row number in rcsubs[i] */
					int r = idxij.val[rcidxs[i]];
					if(r<noPrevVisPoints || r >= noCurVisPoints)
					{
						continue;
					}
					else
					{
						//fprintf(stderr, "computeM:nvis:%d\n", r);
						ptr3=jac + nvis*ABsz + (r-nvis)*ATsz;
						ptr3_pt = ptr3 + Asz;
						ptr3_thor = ptr3_pt + Ptsz;
						ptr3_tver = ptr3_thor + THorsz;
					}
					
					double sum_pt, sum_thor, sum_tver;
					
					for(ii=0; ii<cnp; ++ii)
					{
						for(jj=0; jj<pnp; ++jj)
						{
							sum_pt = sum_thor = sum_tver = 0.0;
							for(k=0; k<mnp; ++k)
							{
								sum_pt += ptr3[k*cnp+ii]*ptr3_pt[k*pnp+jj];
								sum_thor += ptr3[k*cnp+ii]*ptr3_thor[k*pnp+jj];
								sum_tver += ptr3[k*cnp+ii]*ptr3_tver[k*pnp+jj];
							}
							ptr1_0[ii*pnp+jj]+=sum_pt;
							ptr1_1[ii*pnp+jj]+=sum_thor;
							ptr1_2[ii*pnp+jj]+=sum_tver;
						}
					}
				}
			}
			/*printf("***********M0_%d***********\n", j);
			 for(ii=0; ii<cnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1_0[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********M1_%d***********\n", j);
			 for(ii=0; ii<cnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1_1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********M2_%d***********\n", j);
			 for(ii=0; ii<cnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1_2[ii*pnp+jj]);
			 }
			 printf("\n");
			 }*/
        }
        
#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing M took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif
		
		
        /* compute V_i = \sum_j B_ij^T B_ij */ // \Sigma here!
        /* V_i is symmetric, therefore its computation can be sped up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that B_ij is mnp x pnp
         */
        /* Also compute eb_i = \sum_j B_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */
		
#ifdef TIMINGS
        start = clock();
#endif
		
        _dblzero(V, n*Vsz); /* clear all V_i */
        _dblzero(eb, n*ebsz); /* clear all eb_i */
        for(i=0; i<n; ++i)
		{
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
			
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr3 to point to B_ij, actual column number in rcsubs[j] */
                int r = idxij.val[rcidxs[j]];
				if(r>nvis)
				{
					printf("Error in computation of matrix V!");
				}
				
				//printf("computeV:vis: %d\n", r);
				ptr3=jac + r*ABsz + Asz;
				
                /* compute the UPPER TRIANGULAR PART of B_ij^T B_ij and add it to V_i */
                for(ii=0; ii<pnp; ++ii)
				{
                    for(jj=ii; jj<pnp; ++jj)
					{
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*pnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]+=sum;
                    }
                }
				
                ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
                /* compute B_ij^T e_ij and add it to eb_i */
                for(ii=0; ii<pnp; ++ii)
				{
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*pnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }
			
            /* Add in the constraints */
            if (use_point_constraints && point_constraints[i].constrained) {
                for (ii=0; ii < pnp; ii++) {
                    double diff = point_constraints[i].constraints[ii] - p[m * cnp + i * pnp + ii];
                    
                    // printf("  p[%d][%d] diff: %0.3f\n", i, ii, diff);
                    
                    /* Add to the V matrix */
                    ptr1[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
                    ptr2[ii] += totalNVis * point_constraints[i].weight * diff;
                }
            }
			/*printf("***********V_%d***********\n", i);
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=ii; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }*/
        }
        
#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing V took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif
		
        /* compute N = [\sum_i \sum_j Pt_ij^T Pt_ij  \sum_i \sum_j Pt_ij^T THor_ij \sum_i \sum_j Pt_ij^T TVer_ij 
		 \sum_i \sum_j THor_ij^T Pt_ij  \sum_i \sum_j THor_ij^T THor_ij \sum_i \sum_j THor_ij^T TVer_ij
		 \sum_i \sum_j TVer_ij^T Pt_ij  \sum_i \sum_j TVer_ij^T THor_ij \sum_i \sum_j TVer_ij^T TVer_ij]
		 * Recall that Pt_ij, THor_ij, TVer_ij is mnp x pnp
         */
        /* Also compute et_i = Pt_ij^T e_ij + THor_ij^T e_ij + TVer_ij^T e_ij
         * Recall that e_ij is mnp x 1
         */
		
#ifdef TIMINGS
        start = clock();
#endif
		
        _dblzero(N, noGrids*3*3*Nsz); /* clear all N_i */
        _dblzero(et, noGrids*3*etranssz); /* clear all et_i */
		
		for(gridIt=0; gridIt<noGrids; gridIt++)
		{
			ptr1_psum=N+gridIt*Nsz*9; // set ptr1_psum to point to Pt^T*Pt
			ptr1_THorP=ptr1_psum+Nsz; // set ptr1_THorP to point to THor^T*Pt
			ptr1_TVerP=ptr1_THorP+Nsz; // set ptr1_TVerP to point to TVer^T*Pt
			ptr1_PTHor=ptr1_TVerP+Nsz; // set ptr1_PTHor to point to Pt^T*THor
			ptr1_hsum=ptr1_PTHor+Nsz; // set ptr1_hsum to point to THor^T*THor
			ptr1_TVerTHor=ptr1_hsum+Nsz; // set ptr1_TVerTHor to point to TVer^T*THor
			ptr1_PTVer=ptr1_TVerTHor+Nsz; //set ptr1_PTVer to point to Pt^T*TVer
			ptr1_THorTVer=ptr1_PTVer+Nsz; //set ptr1_THorTVer to point to THor^T*TVer
			ptr1_vsum=ptr1_THorTVer+Nsz; //set ptr1_vsum to point to TVer^T*TVer
			
			ptr2_0 = et+gridIt*etranssz*3; // set ptr2_0 to point to et_0
			ptr2_1 = ptr2_0 + etranssz; // set ptr2_1 to point to et_1
			ptr2_2 = ptr2_1 + etranssz; // set ptr2_2 to point to et_2
			
			int curPoints;
			int prevPoints = n;
			int noPrevVisPoints = nvis;
			int noCurVisPoints;
			
			for(i=0; i<gridIt; i++)
			{
				prevPoints += grid_cols[i]*grid_rows[i];
				noPrevVisPoints += nVisGrids[i];
			}
			curPoints = prevPoints;
			curPoints += grid_cols[gridIt]*grid_rows[gridIt];
			noCurVisPoints = noPrevVisPoints;
			noCurVisPoints += nVisGrids[gridIt];
			
			for(i=prevPoints; i<curPoints; ++i)
			{
				nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero Pt_ij,THor_ij, TVer_ij j=0...m-1 */
				
				for(j=0; j<nnz; ++j)
				{
					/* set ptr3 to point to Pt_ij, actual column number in rcsubs[j] */
					int r = idxij.val[rcidxs[j]];
					//printf("Point:%d, r:%d, prevVis:%d\n", i, r, noPrevVisPoints);
					if(r<noPrevVisPoints || r>=noCurVisPoints)
					{
						printf("Error in computation of matrix N!");
					}
					
					if(rcsubs[j]<mcon)
					{ 
						/* Pt_ij is zero */
						continue;
					}
					
					//printf("Computing N: vis:%d\n", r);
					
					ptr3_pt=jac + nvis*ABsz + (r-nvis)*ATsz + Asz;
					ptr3_thor=jac + nvis*ABsz + (r-nvis)*ATsz + Asz + Ptsz;
					ptr3_tver=jac + nvis*ABsz + (r-nvis)*ATsz + Asz + Ptsz + THorsz;
					
					/*printf("********pt%d%d********\n", i, rcsubs[j]);
					 for(ii=0; ii<mnp; ii++)
					 {
					 for(jj=0; jj<pnp; jj++)
					 {
					 printf("%f ", ptr3_pt[ii*pnp+jj]);
					 }
					 printf("\n");
					 }
					 printf("********thor%d%d********\n", i, rcsubs[j]);
					 for(ii=0; ii<mnp; ii++)
					 {
					 for(jj=0; jj<pnp; jj++)
					 {
					 printf("%f ", ptr3_thor[ii*pnp+jj]);
					 }
					 printf("\n");
					 }
					 printf("********tver%d%d********\n", i, rcsubs[j]);
					 for(ii=0; ii<mnp; ii++)
					 {
					 for(jj=0; jj<pnp; jj++)
					 {
					 printf("%f ", ptr3_tver[ii*pnp+jj]);
					 }
					 printf("\n");
					 }*/
					
					double sum_p, sum_ph, sum_pv, sum_h, sum_hv, sum_v;
					
					for(ii=0; ii<pnp; ++ii)
					{
						for(jj=0; jj<pnp; ++jj)
						{
							sum_p = sum_ph = sum_pv = sum_h = sum_hv = sum_v = 0.0;
							for(k=0; k<mnp; ++k)
							{
								sum_p += ptr3_pt[k*pnp+ii] * ptr3_pt[k*pnp+jj];
								sum_ph += ptr3_pt[k*pnp+ii] * ptr3_thor[k*pnp+jj];
								sum_pv += ptr3_pt[k*pnp+ii] * ptr3_tver[k*pnp+jj];
								sum_h += ptr3_thor[k*pnp+ii] * ptr3_thor[k*pnp+jj];
								sum_hv += ptr3_thor[k*pnp+ii] * ptr3_tver[k*pnp+jj];
								sum_v += ptr3_tver[k*pnp+ii] * ptr3_tver[k*pnp+jj];
							}
							ptr1_psum[ii*pnp+jj]+=sum_p;
							ptr1_THorP[ii*pnp+jj]+=sum_ph;
							ptr1_TVerP[ii*pnp+jj]+=sum_pv;
							ptr1_PTHor[ii*pnp+jj]+=sum_ph;
							ptr1_hsum[ii*pnp+jj]+=sum_h;
							ptr1_TVerTHor[ii*pnp+jj]+=sum_hv;
							ptr1_PTVer[ii*pnp+jj]+=sum_pv;
							ptr1_THorTVer[ii*pnp+jj]+=sum_hv;
							ptr1_vsum[ii*pnp+jj]+=sum_v;
							//printf("ptr1_psum[ii*pnp+jj]:%f\n", ptr1_psum[ii*pnp+jj]);
						}
					}
					
					
					double sum_0, sum_1, sum_2;
					
					ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
					/* compute Pt_ij^T e_ij and add it to et_0 */
					for(ii=0; ii<pnp; ++ii)
					{
						sum_0 = sum_1 = sum_2 = 0.0;
						for(jj=0, sum=0.0; jj<mnp; ++jj)
						{
							sum_0+=ptr3_pt[jj*pnp+ii]*ptr4[jj];
							sum_1+=ptr3_thor[jj*pnp+ii]*ptr4[jj];
							sum_2+=ptr3_tver[jj*pnp+ii]*ptr4[jj];
						}
						ptr2_0[ii]+=sum_0;
						ptr2_1[ii]+=sum_1;
						ptr2_2[ii]+=sum_2;
					}
					
					/* Add in the constraints */
					if (use_point_constraints && point_constraints[i].constrained) 
					{
						for (ii=0; ii < pnp; ii++) 
						{
							int r = (i-prevPoints)/grid_cols[gridIt];
							int c = (i-prevPoints)%grid_cols[gridIt];
							
							double diff = point_constraints[i].constraints[ii] - (p[m*cnp+n*pnp+gridIt*3*3+ii]+
																				  p[m*cnp+n*pnp+gridIt*3*3+ii+3]*c+
																				  p[m*cnp+n*pnp+gridIt*3*3+ii+6]*r);
							// printf("  p[%d][%d] diff: %0.3f\n", i, ii, diff);
							
							/* Add to the N matrix */
							ptr1_psum[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_THorP[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_TVerP[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_PTHor[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_hsum[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_TVerTHor[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_PTVer[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_THorTVer[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							ptr1_vsum[ii*pnp+ii] += totalNVis * point_constraints[i].weight;
							
							ptr2_0[ii] += totalNVis * point_constraints[i].weight * diff;
							ptr2_1[ii] += totalNVis * point_constraints[i].weight * diff;
							ptr2_2[ii] += totalNVis * point_constraints[i].weight * diff;
						}
					}
				}
			}
			/*printf("**********et_0*********\n");
			 ptr1 = et+gridIt*etranssz*3;
			 for(ii=0; ii<etranssz; ii++)
			 {
			 printf("%f ", ptr1[ii]);
			 }
			 printf("\n");
			 printf("**********et_1*********\n");
			 ptr1 = et+gridIt*etranssz*3 + etranssz;
			 for(ii=0; ii<etranssz; ii++)
			 {
			 printf("%f ", ptr1[ii]);
			 }
			 printf("\n");
			 printf("**********et_2*********\n");
			 ptr1 = et+gridIt*etranssz*3 + 2*etranssz;
			 for(ii=0; ii<etranssz; ii++)
			 {
			 printf("%f ", ptr1[ii]);
			 }
			 printf("\n");
			 
			 printf("***********N_0***********\n");
			 ptr1=N+gridIt*Nsz*9;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_1***********\n");
			 ptr1=N+gridIt*Nsz*9+1*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_2***********\n");
			 ptr1=N+gridIt*Nsz*9+2*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_3***********\n");
			 ptr1=N+gridIt*Nsz*9+3*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_4***********\n");
			 ptr1=N+gridIt*Nsz*9+4*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_5***********\n");
			 ptr1=N+gridIt*Nsz*9+5*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_6***********\n");
			 ptr1=N+gridIt*Nsz*9+6*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_7***********\n");
			 ptr1=N+gridIt*Nsz*9+7*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }
			 printf("***********N_8***********\n");
			 ptr1=N+gridIt*Nsz*9+8*Nsz;
			 for(ii=0; ii<pnp; ii++)
			 {
			 for(jj=0; jj<pnp; jj++)
			 {
			 printf("%f ", ptr1[ii*pnp+jj]);
			 }
			 printf("\n");
			 }*/
		}
        
#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing N took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif
		
		/* compute W_ij =  A_ij^T B_ij */ // \Sigma here!
        /* Recall that A_ij is mnp x cnp and B_ij is mnp x pnp
         */
        for(i=0; i<n; ++i)
		{
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j)
			{
				int r = idxij.val[rcidxs[j]];
				
				if(r >= nvis)
				{
					printf("Error when computing W!\n");
				}
				
				//printf("compute W: nvis:%d\n", r);
				
                /* set ptr1 to point to W_ij, actual column number in rcsubs[j] */
				ptr1=W + idxij.val[rcidxs[j]]*Wsz;
				
                if(rcsubs[j]<mcon)
				{ 
					/* A_ij is zero */
                    _dblzero(ptr1, Wsz); /* clear W_ij */
                    continue;
                }
				
                /* set ptr2 & ptr3 to point to A_ij & B_ij resp. */
                ptr2=jac  + idxij.val[rcidxs[j]]*ABsz;
                ptr3=ptr2 + Asz;
                /* compute A_ij^T B_ij and store it in W_ij
                 * Recall that storage for A_ij, B_ij does not overlap with that for W_ij,
                 * see the comments related to the initialization of jac above
                 */
                /* assert(ptr2-ptr1>=Wsz); */
                for(ii=0; ii<cnp; ++ii)
                    for(jj=0; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr2[k*cnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]=sum;
                    }
				/*printf("***********W%d%d***********\n", i, rcsubs[j]);
				 for(ii=0; ii<cnp; ii++)
				 {
				 for(jj=0; jj<pnp; jj++)
				 {
				 printf("%f ", ptr1[ii*pnp+jj]);
				 }
				 printf("\n");
				 }*/
            }
        }
		
        /* Compute ||J^T e||_inf and ||p||^2 */
        for(i=0, p_L2=eabtrans_inf=0.0; i<nvars; ++i)
		{
            if(eabtrans_inf < (tmp=FABS(eabtrans[i]))) eabtrans_inf=tmp;
            p_L2+=p[i]*p[i];
        }
        //p_L2=sqrt(p_L2);
		
        /* save diagonal entries so that augmentation can be later canceled.
         * Diagonal entries are in U_j and V_i, N_i
         */
        for(j=mcon; j<m; ++j)
		{
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
            for(i=0; i<cnp; ++i)
                ptr2[i]=ptr1[i*cnp+i];
        }
        for(i=0; i<n; ++i)
		{
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
            for(j=0; j<pnp; ++j)
                ptr2[j]=ptr1[j*pnp+j];
        }
		for(gridIt=0; gridIt<noGrids; gridIt++)
		{
			for(i=0; i<3; ++i)
			{
				ptr1=N + gridIt*9*Nsz + (i*4)*Nsz; // set ptr1 to point to N_i
				ptr2=diagN + gridIt*3*pnp + i*pnp; // set ptr2 to point to diagN_i
				for(j=0; j<pnp; ++j)
				{
					ptr2[j]=ptr1[j*pnp+j];
					//printf("diagN:%f\n", ptr2[j]);
				}
			}
		}
		
        /* check for convergence */
        if((eabtrans_inf <= eps1)){
            dp_L2=0.0; /* no increment for p in this case */
            stop=1;
            break;
        }
		
        /* compute initial damping factor */
        if(itno==0)
		{
            for(i=mcon*cnp, tmp=DBL_MIN; i<nvars; ++i)
			{
				if(diagUVN[i]>tmp) tmp=diagUVN[i]; /* find max diagonal element */
			}
            mu=tau*tmp;
        }
		
		printf("damping factor:%f\n", mu);
        /* determine increment using adaptive damping */
        while(1)
		{
			
#ifdef TIMINGS
            start = clock();
#endif
			
            /* augment U, V, N */
            for(j=mcon; j<m; ++j)
			{
                ptr1=U + j*Usz; // set ptr1 to point to U_j
                for(i=0; i<cnp; ++i)
                    ptr1[i*cnp+i]+=mu;
            }
            for(i=0; i<n; ++i)
			{
                ptr1=V + i*Vsz; // set ptr1 to point to V_i
                for(j=0; j<pnp; ++j)
                    ptr1[j*pnp+j]+=mu;
				
                /* compute V*_i^-1.
                 * Recall that only the upper triangle of the symmetric pnp x pnp matrix V*_i
                 * is stored in ptr1; its (also symmetric) inverse is saved in the lower triangle of ptr1
                 */
                /* inverting V*_i with LDLT seems to result in faster overall execution compared to when using LU or Cholesky */
                //j=sba_symat_invert_LU(ptr1, pnp); matinv=sba_symat_invert_LU;
                //j=sba_symat_invert_Chol(ptr1, pnp); matinv=sba_symat_invert_Chol;
                j=sba_symat_invert_BK(ptr1, pnp); matinv=sba_symat_invert_BK;
                if(!j){
                    fprintf(stderr, "SBA: singular matrix V*_i (i=%d) in sba_motstr_levmar_x(), increasing damping\n", i);
                    goto moredamping; // increasing damping will eventually make V*_i diagonally dominant, thus nonsingular
                    //retval=SBA_ERROR;
                    //goto freemem_and_return;
                }
            }
			
			_dblzero(Ninv, Nsz*3*3*noGrids);
			for(gridIt=0; gridIt<noGrids; gridIt++)
			{
				for(i=0; i<3; ++i)
				{
					ptr1=N + gridIt*9*Nsz + (i*4)*Nsz; // set ptr1 to point to N_i
					for(j=0; j<pnp; ++j)
						ptr1[j*pnp+j]+=mu;
				}
				
				int count=0;
				ptr1=Ninv + gridIt*9*Nsz;
				for(i=0; i<3; i++)
				{
					ptr1_0=N+gridIt*9*Nsz+i*3*Nsz;
					ptr1_1=ptr1_0+Nsz;
					ptr1_2=ptr1_0+Nsz*2;
					for(j=0; j<pnp; j++)
					{
						for(k=0; k<pnp; k++)
						{
							ptr1[count]=ptr1_0[j*pnp+k];
							count+=1;
						}
						for(k=0; k<pnp; k++)
						{
							ptr1[count]=ptr1_1[j*pnp+k];
							count+=1;
						}
						for(k=0; k<pnp; k++)
						{
							ptr1[count]=ptr1_2[j*pnp+k];
							count+=1;
						}
					}
				}
				/* compute N*_i^-1.
				 * Recall that only the upper triangle of the symmetric pnp x pnp matrix N*_i
				 * is stored in ptr1; its (also symmetric) inverse is saved in the lower triangle of ptr1
				 */
				/* inverting N*_i with LDLT seems to result in faster overall execution compared to when using LU or Cholesky */
				//j=sba_symat_invert_LU(ptr1, pnp); matinv=sba_symat_invert_LU;
				//j=sba_symat_invert_Chol(ptr1, pnp); matinv=sba_symat_invert_Chol;
				/*printf("**************N************\n");
				 for(ii=0; ii<3*pnp; ii++)
				 {
				 for(jj=0; jj<3*pnp; jj++)
				 {
				 printf("%f ", ptr1[ii*3*pnp+jj]);
				 }
				 printf("\n");
				 }*/
				j=sba_symat_invert_full_BK(ptr1, 3*pnp); matinv=sba_symat_invert_BK;
				
				if(!j)
				{
					fprintf(stderr, "SBA: singular matrix N in sba_motstr_levmar_x(), increasing damping\n");
					goto moredamping; // increasing damping will eventually make V*_i diagonally dominant, thus nonsingular
					//retval=SBA_ERROR;
					//goto freemem_and_return;
				}
				
				/*printf("**************Ninv************\n");
				 for(ii=0; ii<3*pnp; ii++)
				 {
				 for(jj=0; jj<3*pnp; jj++)
				 {
				 printf("%f ", ptr1[ii*3*pnp+jj]);
				 }
				 printf("\n");
				 }*/
				//change Ninv to block structure
				count = 0;
				double *tmpNInv = malloc(sizeof(double)*3*pnp*3*pnp);
				for(i=0; i<3; i++)
				{
					for(j=0; j<3; j++)
					{
						for(ii=i*3; ii<(i+1)*3; ii++)
						{
							for(jj=j*3; jj<(j+1)*3; jj++)
							{
								tmpNInv[count] = ptr1[ii*3*pnp+jj];
								count +=1;
							}
						}
					}
				}
				memcpy(ptr1, tmpNInv, sizeof(double)*3*pnp*3*pnp);
				/*printf("**************Ninv************\n");
				 for(ii=0; ii<3*pnp; ii++)
				 {
				 for(jj=0; jj<3*pnp; jj++)
				 {
				 printf("%f ", ptr1[ii*3*pnp+jj]);
				 }
				 printf("\n");
				 }
				 free(tmpNInv);*/
			}
			
#ifdef TIMINGS
            end = clock();
            printf("[sba_motstr_levmar_x] computing damping factor took %0.3fs\n", 
                   (end - start) / (float) CLOCKS_PER_SEC);
#endif
			
            _dblzero(E, m*easz); /* clear all e_j */
            /* compute the mmcon x mmcon block matrix S and e_j */
			
            /* Note that S is symmetric, therefore its computation can be
             * sped up by computing only the upper part and then reusing
             * it for the lower one.
             */
			
#ifdef TIMINGS
            start = clock();
#endif
			
            for(j=mcon; j<m; ++j)
			{
                int mmconxUsz=mmcon*Usz;
				
                nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero Y_ij, i=0...n-1 */
				
                /* compute all Y_ij = W_ij (V*_i)^-1 for a *fixed* j.
                 * To save memory, the block matrix consisting of the Y_ij
                 * is not stored. Instead, only a block column of this matrix
                 * is computed & used at each time: For each j, all nonzero
                 * Y_ij are computed in Yj and then used in the calculations
                 * involving S_jk and e_j.
                 * Recall that W_ij is cnp x pnp and (V*_i) is pnp x pnp
                 */
                for(i=0; i<nnz; ++i)
				{
                    if(idxij.val[rcidxs[i]] >= nvis)
						continue;
					
					/* set ptr3 to point to (V*_i)^-1, actual row number in rcsubs[i] */
                    ptr3=V + rcsubs[i]*Vsz;
					
                    /* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
                    ptr1=Yj + i*Ysz;
                    /* set ptr2 to point to W_ij resp. */
                    ptr2=W + idxij.val[rcidxs[i]]*Wsz;
                    /* compute W_ij (V*_i)^-1 and store it in Y_ij.
                     * Recall that only the lower triangle of (V*_i)^-1 is stored
                     */
                    for(ii=0; ii<cnp; ++ii)
					{
                        ptr4=ptr2+ii*pnp;
                        for(jj=0; jj<pnp; ++jj)
						{
                            for(k=0, sum=0.0; k<=jj; ++k)
                                sum+=ptr4[k]*ptr3[jj*pnp+k]; //ptr2[ii*pnp+k]*ptr3[jj*pnp+k];
                            for( ; k<pnp; ++k)
                                sum+=ptr4[k]*ptr3[k*pnp+jj]; //ptr2[ii*pnp+k]*ptr3[k*pnp+jj];
                            ptr1[ii*pnp+jj]=sum;
                        }
                    }
					
					/*printf("***********Y%d%d***********\n", i, j);
					 for(ii=0; ii<cnp; ii++)
					 {
					 for(jj=0; jj<pnp; jj++)
					 {
					 printf("%f ", ptr1[ii*pnp+jj]);
					 }
					 printf("\n");
					 }*/
					
                }
				
				/* compute all Z_ij = M_ij (N*_i)^-1 for a *fixed* j.
				 * To save memory, the block matrix consisting of the Y_ij
				 * is not stored. Instead, only a block column of this matrix
				 * is computed & used at each time: For each j, all nonzero
				 * Y_ij are computed in Yj and then used in the calculations
				 * involving S_jk and e_j.
				 * Recall that W_ij is cnp x pnp and (V*_i) is pnp x pnp
				 */
				for(gridIt=0; gridIt<noGrids; gridIt++)
				{
					for(i=0; i<3; ++i)
					{
						/* set ptr3 to point to (N*_i)^-1, actual row number in rcsubs[i] */
						ptr3_0=Ninv + gridIt*Nsz*9 + i*Nsz;
						ptr3_1=Ninv + gridIt*Nsz*9 + (3+i)*Nsz;
						ptr3_2=Ninv + gridIt*Nsz*9 + (6+i)*Nsz;   
						
						/* set ptr1 to point to Z_ij, actual row number in rcsubs[i] */
						ptr1=Zj + gridIt*Zsz*3 + i*Zsz;
						
						/* set ptr2 to point to M_ij resp. */
						ptr2_0=M + gridIt*Msz*3*m + j*3*Msz;
						ptr2_1=M + gridIt*Msz*3*m + j*3*Msz + Msz;
						ptr2_2=M + gridIt*Msz*3*m + j*3*Msz + 2*Msz;
						
						/*printf("***********M0***********\n", i,j);
						 for(ii=0; ii<cnp; ii++)
						 {
						 for(jj=0; jj<pnp; jj++)
						 {
						 printf("%f ", ptr2_0[ii*pnp+jj]);
						 }
						 printf("\n");
						 }
						 
						 printf("***********M1***********\n", i,j);
						 for(ii=0; ii<cnp; ii++)
						 {
						 for(jj=0; jj<pnp; jj++)
						 {
						 printf("%f ", ptr2_1[ii*pnp+jj]);
						 }
						 printf("\n");
						 }
						 
						 printf("***********M2***********\n", i,j);
						 for(ii=0; ii<cnp; ii++)
						 {
						 for(jj=0; jj<pnp; jj++)
						 {
						 printf("%f ", ptr2_2[ii*pnp+jj]);
						 }
						 printf("\n");
						 }
						 
						 printf("***********N0***********\n", i,j);
						 for(ii=0; ii<pnp; ii++)
						 {
						 for(jj=0; jj<pnp; jj++)
						 {
						 printf("%f ", ptr3_0[ii*pnp+jj]);
						 }
						 printf("\n");
						 }
						 printf("***********N1***********\n", i,j);
						 for(ii=0; ii<pnp; ii++)
						 {
						 for(jj=0; jj<pnp; jj++)
						 {
						 printf("%f ", ptr3_1[ii*pnp+jj]);
						 }
						 printf("\n");
						 }
						 printf("***********N2***********\n", i,j);
						 for(ii=0; ii<pnp; ii++)
						 {
						 for(jj=0; jj<pnp; jj++)
						 {
						 printf("%f ", ptr3_2[ii*pnp+jj]);
						 }
						 printf("\n");
						 }*/
						
						/* compute M_ij (N*_i)^-1 and store it in Z_ij.
						 */
						for(ii=0; ii<cnp; ++ii)
						{
							for(jj=0; jj<pnp; ++jj)
							{
								for(k=0,sum=0.0 ; k<pnp; ++k)
								{
									sum+=ptr2_0[ii*pnp+k]*ptr3_0[k*pnp+jj];
									sum+=ptr2_1[ii*pnp+k]*ptr3_1[k*pnp+jj];
									sum+=ptr2_2[ii*pnp+k]*ptr3_2[k*pnp+jj];
								}
								ptr1[ii*pnp+jj]=sum;
							}
						}
						/*printf("***********Z%d%d***********\n", i,j);
						 for(ii=0; ii<cnp; ii++)
						 {
						 for(jj=0; jj<pnp; jj++)
						 {
						 printf("%f ", ptr1[ii*pnp+jj]);
						 }
						 printf("\n");
						 }*/
					}
				}
				
				/* compute S */
				for(k=0; k<m; ++k)
				{ 
					// j>=mcon
                    /* compute \sum_i Y_ij W_ik^T in YWt. Note that
                     * for an off-diagonal block defined by j, k YWt
                     * (and thus S_jk) is nonzero only if there exists
                     * a point that is visible in both the j-th and
                     * k-th images
                     */
					
                    /* Recall that Y_ij is cnp x pnp and W_ik is 
                     * cnp x pnp */ 
                    _dblzero(YWt, YWtsz); /* clear YWt */
					
                    for(i=0; i<nnz; ++i)
					{
						if(idxij.val[rcidxs[i]] >= nvis)
							continue;
						
                        register double *pYWt;
						
                        /* find the min and max column indices of the elements in row i (actually rcsubs[i])
                         * and make sure that k falls within them. This test handles W_ik's which are
                         * certain to be zero without bothering to call sba_crsm_elmidx()
                         */
                        ii=idxij.colidx[idxij.rowptr[rcsubs[i]]];
                        jj=idxij.colidx[idxij.rowptr[rcsubs[i]+1]-1];
                        if(k<ii || k>jj) continue; /* W_ik == 0 */
						
                        /* set ptr2 to point to W_ik */
                        l=sba_crsm_elmidxp(&idxij, rcsubs[i], k, j, rcidxs[i]);
                        //l=sba_crsm_elmidx(&idxij, rcsubs[i], k);
                        if(l==-1) continue; /* W_ik == 0 */
						
                        ptr2=W + idxij.val[l]*Wsz;
                        /* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
                        ptr1=Yj + i*Ysz;
						
                        for(ii=0; ii<cnp; ++ii)
						{
                            ptr3=ptr1+ii*pnp;
                            pYWt=YWt+ii*cnp;
							
                            ptr4=ptr2;
                            for(jj=0; jj<cnp; ++jj)
							{
                                //ptr4=ptr2+jj*pnp;
                                for(l=0, sum=0.0; l<pnp; ++l)
                                    sum+=ptr3[l]*ptr4[l]; //ptr1[ii*pnp+l]*ptr2[jj*pnp+l];
                                pYWt[jj]+=sum; //YWt[ii*cnp+jj]+=sum;
                                ptr4+=pnp;
                            }
                        }
                    }
					
					/*printf("***********YWt***********\n");
					 for(ii=0; ii<cnp; ii++)
					 {
					 for(jj=0; jj<cnp; jj++)
					 {
					 printf("%f ", YWt[ii*cnp+jj]);
					 }
					 printf("\n");
					 }*/
					
					/* Recall that Z_ij is cnp x pnp and M_ik is 
                     * cnp x pnp */ 
                    _dblzero(ZMt, ZMtsz*noGrids); /* clear ZMt */
					
					for(gridIt=0; gridIt<noGrids; gridIt++)
					{
						for(i=0; i<3; ++i)
						{
							register double *pZMt;
							
							ptr2=M + gridIt*Msz*3*m + k*3*Msz + i*Msz;
							/* set ptr1 to point to Z_ij, actual row number in rcsubs[i] */
							ptr1=Zj + gridIt*Zsz*3 + i*Zsz;
							
							for(ii=0; ii<cnp; ++ii)
							{
								pZMt=ZMt+gridIt*ZMtsz+ii*cnp;
								
								for(jj=0; jj<cnp; ++jj)
								{
									//ptr4=ptr2+jj*pnp;
									for(l=0, sum=0.0; l<pnp; ++l)
										sum+=ptr1[ii*pnp+l]*ptr2[jj*pnp+l];
									pZMt[jj]+=sum;
								}
							}
						}
						
						/*printf("***********ZMt***********\n");
						 for(ii=0; ii<cnp; ii++)
						 {
						 for(jj=0; jj<cnp; jj++)
						 {
						 printf("%f ", ZMt[gridIt*ZMtsz+ii*cnp+jj]);
						 }
						 printf("\n");
						 }*/
					}
					
					
                    /* since the linear system involving S is solved with lapack,
                     * it is preferable to store S in column major (i.e. fortran)
                     * order, so as to avoid unecessary transposing/copying.
                     */
#if MAT_STORAGE==COLUMN_MAJOR
                    ptr2=S + (k-mcon)*mmconxUsz + (j-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#else
                    ptr2=S + (j-mcon)*mmconxUsz + (k-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#endif
					
                    if(j!=k)
					{ 
						/* Kronecker */
                        for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
                            for(jj=0; jj<cnp; ++jj)
							{
								ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
								-YWt[jj*cnp+ii];
								
#else
								-YWt[ii*cnp+jj];
#endif
							}
						for(gridIt=0; gridIt<noGrids; gridIt++)
						{
							ptr2 -= Sdim*cnp;
							for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
								for(jj=0; jj<cnp; ++jj)
								{
									ptr2[jj] -=
#if MAT_STORAGE==COLUMN_MAJOR
									ZMt[gridIt*ZMtsz + jj*cnp+ii];
#else
									ZMt[gridIt*ZMtsz + ii*cnp+jj];
#endif
								}
						}
                    }
                    else
					{
                        ptr1=U + j*Usz; // set ptr1 to point to U_j
						
                        for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
                            for(jj=0; jj<cnp; ++jj)
							{
								ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
								ptr1[jj*cnp+ii] - YWt[jj*cnp+ii];
#else
								ptr1[ii*cnp+jj] - YWt[ii*cnp+jj];
#endif
							}
						for(gridIt=0; gridIt<noGrids; gridIt++)
						{
							ptr2 -= Sdim*cnp;
							for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
								for(jj=0; jj<cnp; ++jj)
								{
									ptr2[jj] -=
#if MAT_STORAGE==COLUMN_MAJOR
									ZMt[gridIt*ZMtsz + jj*cnp+ii];
#else
								    ZMt[gridIt*ZMtsz + ii*cnp+jj];
#endif
								}
						}
                    }
                }
				
                /* compute e_j=ea_j - \sum_i Y_ij eb_i - \sum_i Z_ij et_i */
                /* Recall that Y_ij is cnp x pnp and eb_i is pnp x 1 */
                ptr1=E + j*easz; // set ptr1 to point to e_j
				
                for(i=0; i<nnz; ++i)
				{
                    if(idxij.val[rcidxs[i]] >= nvis)
						break;
					/* set ptr2 to point to Y_ij, actual row number in rcsubs[i] */
                    ptr2=Yj + i*Ysz;
					
                    /* set ptr3 to point to eb_i */
                    ptr3=eb + rcsubs[i]*ebsz;
                    for(ii=0; ii<cnp; ++ii)
					{
                        ptr4=ptr2+ii*pnp;
                        for(jj=0, sum=0.0; jj<pnp; ++jj)
                            sum+=ptr4[jj]*ptr3[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
                        ptr1[ii]+=sum;
                    }
                }
				
				
				/* Recall that Z_ij is cnp x pnp and et_i is pnp x 1 */
				for(gridIt=0; gridIt<noGrids; gridIt++)
				{
					for(i=0; i<3; ++i)
					{
						/* set ptr2 to point to Z_ij, actual row number in rcsubs[i] */
						ptr2=Zj + gridIt*Zsz*3 + i*Zsz;
						
						/* set ptr3 to point to et_i */
						ptr3=et + gridIt*etranssz*3 + i*etranssz;
						for(ii=0; ii<cnp; ++ii)
						{
							ptr4=ptr2+ii*pnp;
							for(jj=0, sum=0.0; jj<pnp; ++jj)
								sum+=ptr4[jj]*ptr3[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
							ptr1[ii]+=sum;
						}
					}
				}
				
                ptr2=ea + j*easz; // set ptr2 to point to ea_j
                for(i=0; i<easz; ++i)
                    ptr1[i]=ptr2[i] - ptr1[i];
            }
			
#ifdef TIMINGS
            end = clock();
            printf("[sba_motstr_levmar_x] computing S and e_j took %0.3fs\n", 
                   (end - start) / (float) CLOCKS_PER_SEC);
#endif
			
            /* solve the linear system S dpa = E to compute the da_j.
             *
             * Note that if MAT_STORAGE==1 S is modified in the following call;
             * this is OK since S is recomputed for each iteration
             */
			
#ifdef TIMINGS
            start = clock();
#endif
			
			//issolved=sba_Axb_LU(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_LU;
            issolved=sba_Axb_Chol(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_Chol;
            //issolved=sba_Axb_BK(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_BK;
            //issolved=sba_Axb_QRnoQ(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_QRnoQ;
            //issolved=sba_Axb_QR(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_QR;
			//issolved=sba_Axb_SVD(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_SVD;
			//issolved=sba_Axb_CG(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, (3*Sdim)/2, 1E-10, SBA_CG_JACOBI, MAT_STORAGE); linsolver=(PLS)sba_Axb_CG;
			
            ++nlss;
			
			_dblzero(dpa, mcon*cnp); /* no change for the first mcon camera params */
			
#ifdef TIMINGS
            end = clock();
            printf("[sba_motstr_levmar_x] solving linear system took %0.3fs\n", 
                   (end - start) / (float) CLOCKS_PER_SEC);
#endif
			
            if(issolved){
				
                /* compute the db_i */
				
#ifdef TIMINGS
                start = clock();
#endif
				
				for(i=0; i<n; ++i)
				{
                    ptr1=dpb + i*ebsz; // set ptr1 to point to db_i
					
                    /* compute \sum_j W_ij^T da_j */
                    /* Recall that W_ij is cnp x pnp and da_j is cnp x 1 */
                    _dblzero(Wtda, Wtdasz); /* clear Wtda */
                    nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
                    for(j=0; j<nnz; ++j)
					{
                        /* set ptr2 to point to W_ij, actual column number in rcsubs[j] */
                        if(rcsubs[j]<mcon) continue; /* W_ij is zero */
						
                        ptr2=W + idxij.val[rcidxs[j]]*Wsz;
						
                        /* set ptr3 to point to da_j */
                        ptr3=dpa + rcsubs[j]*cnp;
						
                        for(ii=0; ii<pnp; ++ii)
						{
                            ptr4=ptr2+ii;
                            for(jj=0, sum=0.0; jj<cnp; ++jj)
                                sum+=ptr4[jj*pnp]*ptr3[jj]; //ptr2[jj*pnp+ii]*ptr3[jj];
                            Wtda[ii]+=sum;
                        }
                    }
					
                    /* compute eb_i - \sum_j W_ij^T da_j = eb_i - Wtda in Wtda */
                    ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
                    for(ii=0; ii<pnp; ++ii)
                        Wtda[ii]=ptr2[ii] - Wtda[ii];
					
                    /* compute the product (V*_i)^-1 Wtda = (V*_i)^-1 (eb_i - \sum_j W_ij^T da_j).
                     * Recall that only the lower triangle of (V*_i)^-1 is stored
                     */
                    ptr2=V + i*Vsz; // set ptr2 to point to (V*_i)^-1
                    for(ii=0; ii<pnp; ++ii)
					{
                        for(jj=0, sum=0.0; jj<=ii; ++jj)
                            sum+=ptr2[ii*pnp+jj]*Wtda[jj];
                        for( ; jj<pnp; ++jj)
                            sum+=ptr2[jj*pnp+ii]*Wtda[jj];
                        ptr1[ii]=sum;
                    }
                }
				
#ifdef TIMINGS
                end = clock();
                printf("[sba_motstr_levmar_x] computing db_i's took %0.3fs\n", 
                       (end - start) / (float) CLOCKS_PER_SEC);
#endif
				
				/* compute the dt_i */
				
#ifdef TIMINGS
                start = clock();
#endif
				/************************************************/
				/* compute \sum_j M_ij^T da_j */
				/* Recall that M_ij is cnp x pnp and da_j is cnp x 1 */
				_dblzero(Mtda, Mtdasz*3*noGrids); /* clear Mtda */
				
				for(gridIt=0; gridIt<noGrids; gridIt++)
				{
					for(i=0; i<3; i++)
					{
						for(j=0; j<m; ++j)
						{
							/* set ptr2 to point to M_ij */
							ptr2=M + gridIt*Msz*m*3 + (j*3+i)*Msz;
							
							/* set ptr3 to point to da_j */
							ptr3=dpa + j*cnp;
							
							for(ii=0; ii<pnp; ++ii)
							{
								ptr4=ptr2+ii;
								for(jj=0, sum=0.0; jj<cnp; ++jj)
								{
									sum+=ptr4[jj*pnp]*ptr3[jj]; //ptr2[jj*pnp+ii]*ptr3[jj];
								}
								Mtda[gridIt*Mtdasz*3 + i*pnp+ii]+=sum;
							}
						}
						
						/* compute et_i - \sum_j M_ij^T da_j = et_i - Mtda in Mtda */
						ptr2=et + gridIt*etranssz*3 + i*etranssz; // set ptr2 to point to et_i
						for(ii=0; ii<pnp; ++ii)
						{
							//printf("Mtda[%d] b:%f\n", gridIt*Mtdasz*3 + i*pnp+ii, Mtda[gridIt*Mtdasz*3 + i*pnp+ii]);
							//printf("ptr2[%d]:%f\n]", etranssz*3 + i*etranssz + ii, ptr2[ii]);
							Mtda[gridIt*Mtdasz*3 + i*pnp+ii]=ptr2[ii] - Mtda[gridIt*Mtdasz*3 + i*pnp+ii];
							//printf("Mtda[%d] a:%f\n", gridIt*Mtdasz*3 + i*pnp+ii, Mtda[gridIt*Mtdasz*3 + i*pnp+ii]);
						}
					}
					
					for(i=0; i<3; i++)
					{
						/* compute the product (N*_i)^-1 Mtda = (N*_i)^-1 (et_i - \sum_j M_ij^T da_j).
						 * Recall that only the lower triangle of (N*_i)^-1 is stored
						 */
						ptr1 = dpt + gridIt*etranssz*3 + i*etranssz;
						
						//printf("ii:%d\n", gridIt*etranssz*3 + i*etranssz);
						ptr1[0]=0.0;
						
						ptr2_0=Ninv + gridIt*Nsz*9 + i*3*Nsz; // set ptr2 to point to (N*_i)^-1
						ptr2_1=Ninv + gridIt*Nsz*9 + i*3*Nsz + Nsz; // set ptr2 to point to (N*_i)^-1
						ptr2_2=Ninv + gridIt*Nsz*9 + i*3*Nsz + 2*Nsz; // set ptr2 to point to (N*_i)^-1
						
						ptr3_0=Mtda + gridIt*Mtdasz*3;
						ptr3_1=Mtda + gridIt*Mtdasz*3 + pnp;
						ptr3_2=Mtda + gridIt*Mtdasz*3 + pnp*2;
						
						double sum = 0.0;
						
						for(ii=0; ii<pnp; ++ii)
						{
							sum = 0.0;
							for(jj=0; jj<pnp; ++jj)
							{
								sum += ptr2_0[ii*pnp+jj]*ptr3_0[jj];
								sum += ptr2_1[ii*pnp+jj]*ptr3_1[jj];
								sum += ptr2_2[ii*pnp+jj]*ptr3_2[jj];
								//printf("ptr2_0[ii*pnp+jj]:%f, ptr3_0[jj]:%f\n",ptr2_0[ii*pnp+jj], ptr3_0[jj]);
							}
							
							//printf("ii:%d\n", gridIt*etranssz*3 + i*etranssz + ii);
							ptr1[ii]=sum;
							//printf("dpt[%d]:%f (sum:%f)\n", gridIt*etranssz*3+i*etranssz+ii, ptr1[ii], sum);
						}
					}
				}
				
				
#ifdef TIMINGS
                end = clock();
                printf("[sba_motstr_levmar_x] computing dt_i's took %0.3fs\n", 
                       (end - start) / (float) CLOCKS_PER_SEC);
#endif
				
				
                /* parameter vector updates are now in dpa, dpb, dpt */
				
                /* compute p's new estimate and ||dp||^2 */
                for(i=0, dp_L2=0.0; i<nvars; ++i)
				{
                    pdp[i]=p[i] + (tmp=dp[i]);
                    dp_L2+=tmp*tmp;
                }
                //dp_L2=sqrt(dp_L2);
				
                if(dp_L2<=eps2_sq*p_L2)
				{ 
					/* relative change in p is small, stop */
                    //if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
                    stop=2;
                    break;
                }
				
                if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ)
				{ 
					/* almost singular */
                    //if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
                    fprintf(stderr, "SBA: the matrix of the augmented normal equations is almost singular in sba_motstr_levmar_x(),\n"
                            "     minimization should be restarted from the current solution with an increased damping term\n");
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }
				
#ifdef TIMINGS
                start = clock();
#endif
				
				/*for(i=0; i<nvars; i++)
				 {
				 printf("nvar:%f\n", pdp[i]);
				 }*/
				
                (*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
                if(verbose>1) {
                    printf("mean reprojection error %g\n", sba_mean_repr_error(totalNoPoints, mnp, x, hx, &idxij, rcidxs, rcsubs));
                    fflush(stdout);
                }
				
                /* ### compute ||e(pdp)||_2 */
				if(covx==NULL)
                    pdp_eL2=nrmL2xmy(hx, x, hx, nobs); /* hx=x-hx, pdp_eL2=||hx|| */
                else
                    pdp_eL2=nrmCxmy(hx, x, hx, wght, mnp, totalNVis); /* hx=wght*(x-hx), pdp_eL2=||hx|| */
				
                if(!SBA_FINITE(pdp_eL2)){
                    if(verbose) /* identify the offending point projection */
                        sba_print_inf(hx, m, mnp, &idxij, rcidxs, rcsubs);
					
                    stop=7;
                    break;
                }
				
                /* Add in the camera constraints */
                if (use_constraints) 
				{
                    for (j = 0; j < m; j++) 
					{
                        int nnz = 0;
						
                        for(i = 0; i < totalNoPoints; i++)
                            nnz += vmask[i * m + j];
						
                        for (jj = 0; jj < cnp; jj++) 
						{
                            if (constraints[j].constrained[jj]) 
							{
                                double diff = 
								constraints[j].constraints[jj] - p[j * cnp + jj];
								
                                pdp_eL2 += 
								/*nnz **/ constraints[j].weights[jj] * diff * diff;
                            }
                        }
                    }
                }
				
                /* Add in the point constraints */
                if (use_point_constraints) 
				{
                    for (i = 0; i < totalNoPoints; i++) 
					{
                        if (point_constraints[i].constrained) {
                            for (ii=0; ii < pnp; ii++) {
                                double diff;
								if(i<n)
								{
									diff = point_constraints[i].constraints[ii] - p[m * cnp + i * pnp + ii];
								}
								else
								{
									int noCurPoints = n;
									for(j=0; j<noGrids; j++)
									{
										noCurPoints += grid_rows[j]*grid_cols[j];
										
										if(i < noCurPoints)
										{
											int noPrevPoints = 0;
											for(jj=0; jj<j; jj++)
												noPrevPoints += grid_cols[jj]*grid_rows[jj];
											
											int row = (i-noPrevPoints)/grid_cols[j];
											int col = (i-noPrevPoints)%grid_cols[j];
											
											diff = point_constraints[i].constraints[ii] - (p[m*cnp+n*pnp+j*3*3+ii]+
																						   p[m*cnp+n*pnp+j*3*3+ii+3]*col+
																						   p[m*cnp+n*pnp+j*3*3+ii+6]*row);
											break;
										}
									}
								}
								
								/* Add to the U matrix */
                                pdp_eL2 += 
								totalNVis * point_constraints[i].weight * diff * diff;
                            }
                        }
                    }
                }
				
                for(i=0, dL=0.0; i<nvars; ++i)
                    dL+=dp[i]*(mu*dp[i]+eabtrans[i]);
				
                dF=p_eL2-pdp_eL2;
				
#ifdef TIMINGS
                end = clock();
                printf("[sba_motstr_levmar_x] fn eval took %0.3fs\n", 
                       (end - start) / (float) CLOCKS_PER_SEC);
#endif
				
                if(verbose>1) {
                    printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);
                    printf("pdp_eL2: %0.3f, nvis: %d\n", pdp_eL2, nvis);
					
                    fflush(stdout);
                }
				
				
                if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
                    double max_pct_change = 0.0;
					
                    tmp=(2.0*dF/dL-1.0);
                    tmp=1.0-tmp*tmp*tmp;
                    mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
                    nu=2;
					
                    /* Check if we need to stop */
                    for (i=0; i < nobs; i++) {
                        double pct_change;
                        if (e[i] < eps5 && hx[i] < eps5) {
                            pct_change = 0.0;
                        } else {
                            double pct_change = fabs((e[i] - hx[i]) / e[i]);
                            if (pct_change > max_pct_change)
                                max_pct_change = pct_change;
                        }
                    }
					
                    printf("max_pct_change: %0.3e\n", max_pct_change);
                    fflush(stdout);
					
                    /* the test below is equivalent to the relative reduction of the RMS reprojection error: sqrt(p_eL2)-sqrt(pdp_eL2)<eps4*sqrt(p_eL2) */
                    if(pdp_eL2-2.0*sqrt(p_eL2*pdp_eL2)<(eps4_sq-1.0)*p_eL2) stop=4;
					
                    if (max_pct_change < eps5 && itno >= 4) { /* Run at least 4 iters */
                        stop=8;
                        break;
                    }
					
                    for(i=0; i<nvars; ++i) /* update p's estimate */
                        p[i]=pdp[i];
					
                    for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
                        e[i]=hx[i];
                    p_eL2=pdp_eL2;
                    break;
                }
            } /* issolved */
			
        moredamping:
            /* if this point is reached (also via an explicit goto!), either the linear system could
             * not be solved or the error did not reduce; in any case, the increment must be rejected
             */
			
            mu*=nu;
            nu2=nu<<1; // 2*nu;
            if(nu2<=nu){ /* nu has wrapped around (overflown) */
                fprintf(stderr, "SBA: too many failed attempts to increase the damping factor in sba_motstr_levmar_x()! Singular Hessian matrix?\n");
                //exit(1);
                stop=6;
                break;
            }
            nu=nu2;
			
            /* restore U, V, N diagonal entries */
            for(j=mcon; j<m; ++j)
			{
                ptr1=U + j*Usz; // set ptr1 to point to U_j
                ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
                for(i=0; i<cnp; ++i)
                    ptr1[i*cnp+i]=ptr2[i];
            }
            for(i=0; i<n; ++i)
			{
                ptr1=V + i*Vsz; // set ptr1 to point to V_i
                ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
                for(j=0; j<pnp; ++j)
                    ptr1[j*pnp+j]=ptr2[j];
            }
			for(gridIt=0; gridIt<noGrids; gridIt++)
			{
				for(i=0; i<3; ++i)
				{
					ptr1=N + gridIt*Nsz*9 + (i*4)*Nsz; // set ptr1 to point to N_i
					ptr2=diagN + gridIt*3*pnp + i*pnp; // set ptr2 to point to diagN_i
					for(j=0; j<pnp; ++j)
						ptr1[j*pnp+j]=ptr2[j];
				}
			}
        } /* inner while loop */
		
        if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
    }
	
    if(itno>=itmax) stop=3;
	
    /* restore U, V, N diagonal entries */
    for(j=mcon; j<m; ++j){
        ptr1=U + j*Usz; // set ptr1 to point to U_j
        ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
        for(i=0; i<cnp; ++i)
            ptr1[i*cnp+i]=ptr2[i];
    }
    for(i=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
        for(j=0; j<pnp; ++j)
            ptr1[j*pnp+j]=ptr2[j];
    }
	for(gridIt=0; gridIt<noGrids; gridIt++)
	{
		for(i=0; i<3; ++i)
		{
			ptr1=N + gridIt*Nsz*9 + (i*4)*Nsz; // set ptr1 to point to N_i
			ptr2=diagN + gridIt*3*pnp + i*pnp; // set ptr2 to point to diagN_i
			for(j=0; j<pnp; ++j)
				ptr1[j*pnp+j]=ptr2[j];
		}
	}
	
    if(info){
        info[0]=init_p_eL2;
        info[1]=p_eL2;
        info[2]=eabtrans_inf;
        info[3]=dp_L2;
        for(j=mcon, tmp=DBL_MIN; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            for(i=0; i<cnp; ++i)
                if(tmp<ptr1[i*cnp+i]) tmp=ptr1[i*cnp+i];
        }
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            for(j=0; j<pnp; ++j)
                if(tmp<ptr1[j*pnp+j]) tmp=ptr1[j*pnp+j];
        }
        info[4]=mu/tmp;
        info[5]=itno;
        info[6]=stop;
        info[7]=nfev;
        info[8]=njev;
        info[9]=nlss;
    }
	
    //sba_print_sol(n, m, p, cnp, pnp, x, mnp, &idxij, rcidxs, rcsubs);
    retval=(stop!=7)?  itno : SBA_ERROR;
	
freemem_and_return: /* NOTE: this point is also reached via a goto! */
	
    /* free whatever was allocated */
    free(M); free(N); free(Ninv);
    free(e);   free(eabtrans);
    free(E);   free(Yj); free(YWt); free(Zj); free(ZMt);
    free(S);   free(dp); free(Wtda); free(Mtda);
    free(rcidxs); free(rcsubs);
#ifndef SBA_DESTROY_COVS
    if(wght) free(wght);
#else
    /* nothing to do */
#endif /* SBA_DESTROY_COVS */
	
    free(hx); free(diagUVN); free(pdp);
    if(fdj_data.hxx){ // cleanup
        free(fdj_data.hxx);
        free(fdj_data.func_rcidxs);
    }
	
    sba_crsm_free(&idxij);
	
    /* free the memory allocated by the matrix inversion & linear solver routines */
    if(matinv) (*matinv)(NULL, 0);
	//sba_symat_invert_full_BK(NULL, 0);
    if(linsolver) (*linsolver)(NULL, NULL, NULL, 0, 0);
	
    return retval;
}

/* Bundle adjustment on camera and structure parameters 
 * using the sparse Levenberg-Marquardt as described in HZ p. 568
 *
 * Returns the number of iterations (>=0) if successfull, SBA_ERROR if failed
 */

int sba_motstr_levmar_x(
                        const int n,   /* number of points */
                        const int m,   /* number of images */
                        const int mcon,/* number of images (starting from the 1st) whose parameters should not be modified.
                                        * All A_ij (see below) with j<mcon are assumed to be zero
                                        */
                        char *vmask,  /* visibility mask: vmask[i, j]=1 if point i visible in image j, 0 otherwise. nxm */
                        double *p,    /* initial parameter vector p0: (a1, ..., am, b1, ..., bn).
                                       * aj are the image j parameters, bi are the i-th point parameters,
                                       * size m*cnp + n*pnp
                                       */
                        const int cnp,/* number of parameters for ONE camera; e.g. 6 for Euclidean cameras */
                        const int pnp,/* number of parameters for ONE point; e.g. 3 for Euclidean points */
                        double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
                                       * x_ij is the projection of the i-th point on the j-th image.
                                       * NOTE: some of the x_ij might be missing, if point i is not visible in image j;
                                       * see vmask[i, j], max. size n*m*mnp
                                       */
                        double *covx, /* measurements covariance matrices: (Sigma_x_11, .. Sigma_x_1m, ..., Sigma_x_n1, .. Sigma_x_nm),
                                       * where Sigma_x_ij is the mnp x mnp covariance of x_ij stored row-by-row. Set to NULL if no
                                       * covariance estimates are available (identity matrices are implicitly used in this case).
                                       * NOTE: a certain Sigma_x_ij is missing if the corresponding x_ij is also missing;
                                       * see vmask[i, j], max. size n*m*mnp*mnp
                                       */
                        const int mnp,/* number of parameters for EACH measurement; usually 2 */
                        void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
                        /* functional relation describing measurements. Given a parameter vector p,
                         * computes a prediction of the measurements \hat{x}. p is (m*cnp + n*pnp)x1,
                         * \hat{x} is (n*m*mnp)x1, maximum
                         * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                         * as working memory
                         */
                        void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
                        /* function to evaluate the sparse jacobian dX/dp.
                         * The Jacobian is returned in jac as
                         * (dx_11/da_1, ..., dx_1m/da_m, ..., dx_n1/da_1, ..., dx_nm/da_m,
                         *  dx_11/db_1, ..., dx_1m/db_1, ..., dx_n1/db_n, ..., dx_nm/db_n), or (using HZ's notation),
                         * jac=(A_11, B_11, ..., A_1m, B_1m, ..., A_n1, B_n1, ..., A_nm, B_nm)
                         * Notice that depending on idxij, some of the A_ij and B_ij might be missing.
                         * Note also that A_ij and B_ij are mnp x cnp and mnp x pnp matrices resp. and they
                         * should be stored in jac in row-major order.
                         * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                         * as working memory
                         *
                         * If NULL, the jacobian is approximated by repetitive func calls and finite
                         * differences. This is computationally inefficient and thus NOT recommended.
                         */
                        void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 

                        const int itmax,   /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
                        const int verbose, /* I: verbosity */
                        const double opts[SBA_OPTSSZ],
                        /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \epsilon4]. Respectively the scale factor for initial \mu,
                         * stopping thresholds for ||J^T e||_inf, ||dp||_2, ||e||_2 and (||e||_2-||e_new||_2)/||e||_2
                         */
                        double info[SBA_INFOSZ],
                        /* O: information regarding the minimization. Set to NULL if don't care
                         * info[0]=||e||_2 at initial p.
                         * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                         * info[5]= # iterations,
                         * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                         *                                 2 - stopped by small dp
                         *                                 3 - stopped by itmax
                         *                                 4 - stopped by small relative reduction in ||e||_2
                         *                                 5 - stopped by small ||e||_2
                         *                                 6 - too many attempts to increase damping. Restart with increased mu
                         *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
                         * info[7]= # function evaluations
                         * info[8]= # jacobian evaluations
                         * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
                         */
                        int use_constraints, camera_constraints_t *constraints,  /* Constraints on camera parameters */
                        int use_point_constraints, point_constraints_t *point_constraints,
                        double *Vout, double *Sout, double *Uout, double *Wout
                        )
{
    register int i, j, ii, jj, k, l;
    int nvis, nnz, retval;

    /* The following are work arrays that are dynamically allocated by sba_motstr_levmar_x() */
    double *jac;  /* work array for storing the jacobian, max. size n*m*(mnp*cnp + mnp*pnp) */
    double *U;    /* work array for storing the U_j in the order U_1, ..., U_m, size m*cnp*cnp */
    double *V;    /* work array for storing the *strictly upper triangles* of V_i in the order V_1, ..., V_n, size n*pnp*pnp.
                   * V also stores the lower triangles of (V*_i)^-1 in the order (V*_1)^-1, ..., (V*_n)^-1.
                   * Note that diagonal elements of V_1 are saved in diagUV
                   */

    double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_n1, ..., e_nm,
                     max. size n*m*mnp */
    double *eab;  /* work array for storing the ea_j & eb_i in the order ea_1, .. ea_m eb_1, .. eb_n size m*cnp + n*pnp */

    double *E;   /* work array for storing the e_j in the order e_1, .. e_m, size m*cnp */

    /* Notice that the blocks W_ij, Y_ij are zero iff A_ij (equivalently B_ij) is zero. This means
     * that the matrix consisting of blocks W_ij is itself sparse, similarly to the
     * block matrix made up of the A_ij and B_ij (i.e. jac)
     */
    double *W;    /* work array for storing the W_ij in the order W_11, ..., W_1m, ..., W_n1, ..., W_nm,
                     max. size n*m*cnp*pnp */
    double *Yj;   /* work array for storing the Y_ij for a *fixed* j in the order Y_1j, Y_nj,
                     max. size n*cnp*pnp */
    double *YWt;  /* work array for storing \sum_i Y_ij W_ik^T, size cnp*cnp */
    double *S;    /* work array for storing the block array S_jk, size m*m*cnp*cnp */
    double *dp;   /* work array for storing the parameter vector updates da_1, ..., da_m, db_1, ..., db_n, size m*cnp + n*pnp */
    double *Wtda; /* work array for storing \sum_j W_ij^T da_j, size pnp */
    double *wght= /* work array for storing the weights computed from the covariance inverses, max. size n*m*mnp*mnp */
        NULL;

    /* Of the above arrays, jac, e, W, Yj, wght are sparse and
     * U, V, eab, E, S, dp are dense. Sparse arrays (except Yj) are indexed
     * through idxij (see below), that is with the same mechanism as the input 
     * measurements vector x
     */

    double *pa, *pb, *ea, *eb, *dpa, *dpb; /* pointers into p, jac, eab and dp respectively */

    /* submatrices sizes */
    int Asz, Bsz, ABsz, Usz, Vsz,
        Wsz, Ysz, esz, easz, ebsz,
        YWtsz, Wtdasz, Sblsz, covsz;

    int Sdim; /* S matrix actual dimension */

    register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
    struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also
                            * the location of A_ij, B_ij in jac, etc.
                            * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                            * index k and it is used to efficiently lookup the memory locations where the non-zero
                            * blocks of a sparse matrix/vector are stored
                            */
    int maxCvis, /* max. of projections of a single point  across cameras, <=m */
        maxPvis, /* max. of projections in a single camera across points,  <=n */
        maxCPvis, /* max. of the above */
        *rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
                     column in a sparse matrix, size max(n, m) */
        *rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
                     sparse matrix, size max(n, m) */

    /* The following variables are needed by the LM algorithm */
    register int itno;  /* iteration counter */
    int issolved;
    /* temporary work arrays that are dynamically allocated */
    double *hx,         /* \hat{x}_i, max. size m*n*mnp */
        *diagUV,     /* diagonals of U_j, V_i, size m*cnp + n*pnp */
        *pdp;        /* p + dp, size m*cnp + n*pnp */

    double *diagU, *diagV; /* pointers into diagUV */

    register double mu,  /* damping constant */
        tmp; /* mainly used in matrix & vector multiplications */
    double p_eL2, eab_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
    double p_L2, dp_L2=DBL_MAX, dF, dL;
    double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2],
        eps3_sq=opts[3]*opts[3], eps4_sq=opts[4]*opts[4], eps5=opts[5];
    double init_p_eL2;
    int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;
    int nobs, nvars;
    const int mmcon=m-mcon;
    PLS linsolver=NULL;
    int (*matinv)(double *A, int m)=NULL;

    struct fdj_data_x_ fdj_data;
    void *jac_adata;

    /* Initialization */

#ifdef TIMINGS
    clock_t start = clock();
#endif

    mu=eab_inf=0.0; /* -Wall */

    /* block sizes */
    Asz=mnp * cnp; Bsz=mnp * pnp; ABsz=Asz + Bsz;
    Usz=cnp * cnp; Vsz=pnp * pnp;
    Wsz=cnp * pnp; Ysz=cnp * pnp;
    esz=mnp;
    easz=cnp; ebsz=pnp;
    YWtsz=cnp * cnp;
    Wtdasz=pnp;
    Sblsz=cnp * cnp;
    Sdim=mmcon * cnp;
    covsz=mnp * mnp;

    /* count total number of visible image points */
    for(i=nvis=0, jj=n*m; i<jj; ++i)
        nvis+=(vmask[i]!=0);

    nobs=nvis*mnp;
    nvars=m*cnp + n*pnp;
    if(nobs<nvars){
        fprintf(stderr, "SBA: sba_motstr_levmar_x() cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
        return SBA_ERROR;
    }

    /* allocate & fill up the idxij structure */
    sba_crsm_alloc(&idxij, n, m, nvis);
    for(i=k=0; i<n; ++i){
        idxij.rowptr[i]=k;
        ii=i*m;
        for(j=0; j<m; ++j)
            if(vmask[ii+j]){
                idxij.val[k]=k;
                idxij.colidx[k++]=j;
            }
    }
    idxij.rowptr[n]=nvis;

    /* find the maximum number (for all cameras) of visible image projections coming from a single 3D point */
    for(i=maxCvis=0; i<n; ++i)
        if((k=idxij.rowptr[i+1]-idxij.rowptr[i])>maxCvis) maxCvis=k;

    /* find the maximum number (for all points) of visible image projections in any single camera */
    for(j=maxPvis=0; j<m; ++j){
        for(i=ii=0; i<n; ++i)
            if(vmask[i*m+j]) ++ii;
        if(ii>maxPvis) maxPvis=ii;
    }
    maxCPvis=(maxCvis>=maxPvis)? maxCvis : maxPvis;

#if 0
    /* determine the density of matrix S */
    for(j=mcon, ii=0; j<m; ++j){
        ++ii; /* block Sjj is surely nonzero */
        for(k=j+1; k<m; ++k)
            if(sba_crsm_common_row(&idxij, j, k)) ii+=2; /* blocks Sjk & Skj are nonzero */
    }
    printf("\nS density: %.5g\n", ((double)ii)/(mmcon*mmcon)); fflush(stdout);
#endif

    /* allocate work arrays */
    /* W is big enough to hold both jac & W. Note also the extra Wsz, see the initialization of jac below for explanation */
    W=(double *)emalloc((nvis*((Wsz>=ABsz)? Wsz : ABsz) + Wsz)*sizeof(double));
    U=(double *)emalloc(m*Usz*sizeof(double));
    V=(double *)emalloc(n*Vsz*sizeof(double));
    e=(double *)emalloc(nobs*sizeof(double));
    eab=(double *)emalloc(nvars*sizeof(double));
    E=(double *)emalloc(m*cnp*sizeof(double));
    Yj=(double *)emalloc(maxPvis*Ysz*sizeof(double));
    YWt=(double *)emalloc(YWtsz*sizeof(double));
    S=(double *)emalloc(m*m*Sblsz*sizeof(double));
    dp=(double *)emalloc(nvars*sizeof(double));
    Wtda=(double *)emalloc(pnp*sizeof(double));
    rcidxs=(int *)emalloc(maxCPvis*sizeof(int));
    rcsubs=(int *)emalloc(maxCPvis*sizeof(int));
#ifndef SBA_DESTROY_COVS
    if(covx!=NULL) wght=(double *)emalloc(nvis*covsz*sizeof(double));
#else
    if(covx!=NULL) wght=covx;
#endif /* SBA_DESTROY_COVS */


    hx=(double *)emalloc(nobs*sizeof(double));
    diagUV=(double *)emalloc(nvars*sizeof(double));
    pdp=(double *)emalloc(nvars*sizeof(double));

#ifdef TIMINGS
    clock_t end = clock();
    printf("[sba_motstr_levmar_x] initialization took %0.3fs\n", 
           (end - start) / (float) CLOCKS_PER_SEC);
#endif


    /* to save resources, W and jac share the same memory: First, the jacobian
     * is computed in some working memory that is then overwritten during the
     * computation of W. To account for the case of W being larger than jac,
     * extra memory is reserved "before" jac.
     * Care must be taken, however, to ensure that storing a certain W_ij
     * does not overwrite the A_ij, B_ij used to compute it. To achieve
     * this is, note that if p1 and p2 respectively point to the first elements
     * of a certain W_ij and A_ij, B_ij pair, we should have p2-p1>=Wsz.
     * There are two cases:
     * a) Wsz>=ABsz: Then p1=W+k*Wsz and p2=jac+k*ABsz=W+Wsz+nvis*(Wsz-ABsz)+k*ABsz
     *    for some k (0<=k<nvis), thus p2-p1=(nvis-k)*(Wsz-ABsz)+Wsz. 
     *    The right side of the last equation is obviously > Wsz for all 0<=k<nvis
     *
     * b) Wsz<ABsz: Then p1=W+k*Wsz and p2=jac+k*ABsz=W+Wsz+k*ABsz and
     *    p2-p1=Wsz+k*(ABsz-Wsz), which is again > Wsz for all 0<=k<nvis
     *
     * Concluding, if jac is initialized as below, the memory allocated to all
     * W_ij is guaranteed not to overlap with that allocated to their corresponding
     * A_ij, B_ij pairs
     */
    jac=W + Wsz + ((Wsz>ABsz)? nvis*(Wsz-ABsz) : 0);

    /* set up auxiliary pointers */
    pa=p; pb=p+m*cnp;
    ea=eab; eb=eab+m*cnp;
    dpa=dp; dpb=dp+m*cnp;

    diagU=diagUV; diagV=diagUV + m*cnp;

    /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
    if(!fjac){
        fdj_data.func=func;
        fdj_data.cnp=cnp;
        fdj_data.pnp=pnp;
        fdj_data.mnp=mnp;
        fdj_data.hx=hx;
        fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
        fdj_data.func_rcidxs=(int *)emalloc(2*maxCPvis*sizeof(int));
        fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxCPvis;
        fdj_data.adata=adata;

        fjac=sba_fdjac_x;
        jac_adata=(void *)&fdj_data;
    }
    else{
        fdj_data.hxx=NULL;
        jac_adata=adata;
    }

    if(itmax==0){ /* verify jacobian */
        sba_motstr_chkjac_x(func, fjac, p, &idxij, rcidxs, rcsubs, mcon, cnp, pnp, mnp, adata, jac_adata);
        retval=0;
        goto freemem_and_return;
    }

    /* covariances Sigma_x_ij are accommodated by computing the Cholesky decompositions of their
     * inverses and using the resulting matrices w_x_ij to weigh A_ij, B_ij, and e_ij as w_x_ij A_ij,
     * w_x_ij*B_ij and w_x_ij*e_ij. In this way, auxiliary variables as U_j=\sum_i A_ij^T A_ij
     * actually become \sum_i (w_x_ij A_ij)^T w_x_ij A_ij= \sum_i A_ij^T w_x_ij^T w_x_ij A_ij =
     * A_ij^T Sigma_x_ij^-1 A_ij
     *
     * ea_j, V_i, eb_i, W_ij are weighted in a similar manner
     */
    if(covx!=NULL){
        for(i=0; i<n; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1, ptr2 to point to cov_x_ij, w_x_ij resp. */
                ptr1=covx + (k=idxij.val[rcidxs[j]]*covsz);
                ptr2=wght + k;
                if(!sba_mat_cholinv(ptr1, ptr2, mnp)){ /* compute w_x_ij s.t. w_x_ij^T w_x_ij = cov_x_ij^-1 */
                    fprintf(stderr, "SBA: invalid covariance matrix for x_ij (i=%d, j=%d) in sba_motstr_levmar_x()\n", i, rcsubs[j]);
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }
            }
        }
        sba_mat_cholinv(NULL, NULL, 0); /* cleanup */
    }

    /* compute the error vectors e_ij in hx */
    (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
    /* ### compute e=x - f(p) [e=w*(x - f(p)] and its L2 norm */
    if(covx==NULL)
        p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
    else
        p_eL2=nrmCxmy(e, x, hx, wght, mnp, nvis); /* e=wght*(x-hx), p_eL2=||e||=||x-hx||_Sigma^-1 */

    /* Add in the camera constraints */
    if (use_constraints) {
        for (j = 0; j < m; j++) {
            /* Find number of points projecting into this camera */
            int nnz = 0;

            for(i = 0; i < n; i++)
                nnz += vmask[i * m + j];

            for (jj = 0; jj < cnp; jj++) {
                if (constraints[j].constrained[jj]) {
                    double diff = 
                        constraints[j].constraints[jj] - p[j * cnp + jj];

                    p_eL2 += /*nnz **/ constraints[j].weights[jj] * diff * diff;
                }
            }
        }
    }
  
    /* Add in the point constraints */
    if (use_point_constraints) {
        for (i = 0; i < n; i++) {
            if (point_constraints[i].constrained) {
                for (ii=0; ii < pnp; ii++) {
                    double diff = 
                        point_constraints[i].constraints[ii] - 
                        p[m * cnp + i * pnp + ii];
	      
                    /* Add to the U matrix */
                    p_eL2 += nvis * point_constraints[i].weight * diff * diff;
                }
            }
        }
    }

    if(verbose) printf("initial motstr-SBA error %g [%g]\n", p_eL2, p_eL2/nvis);
    init_p_eL2=p_eL2;
    if(!SBA_FINITE(p_eL2)) stop=7;

    for(itno=0; itno<itmax && !stop; ++itno){
        /* Note that p, e and ||e||_2 have been updated at the previous iteration */

#ifdef TIMINGS
        start = clock();
#endif

        /* compute derivative submatrices A_ij, B_ij */
        (*fjac)(p, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;

#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing derivs took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif

#ifdef TIMINGS
        start = clock();
#endif

        if(covx!=NULL){
            /* compute w_x_ij A_ij and w_x_ij B_ij.
             * Since w_x_ij is upper triangular, the products can be safely saved
             * directly in A_ij, B_ij, without the need for intermediate storage
             */
            for(i=0; i<nvis; ++i){
                /* set ptr1, ptr2, ptr3 to point to w_x_ij, A_ij, B_ij, resp. */
                ptr1=wght + i*covsz;
                ptr2=jac  + i*ABsz;
                ptr3=ptr2 + Asz; // ptr3=jac  + i*ABsz + Asz;

                /* w_x_ij is mnp x mnp, A_ij is mnp x cnp, B_ij is mnp x pnp
                 * NOTE: Jamming the outter (i.e., ii) loops did not run faster!
                 */
                /* A_ij */
                for(ii=0; ii<mnp; ++ii)
                    for(jj=0; jj<cnp; ++jj){
                        for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
                            sum+=ptr1[ii*mnp+k]*ptr2[k*cnp+jj];
                        ptr2[ii*cnp+jj]=sum;
                    }

                /* B_ij */
                for(ii=0; ii<mnp; ++ii)
                    for(jj=0; jj<pnp; ++jj){
                        for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
                            sum+=ptr1[ii*mnp+k]*ptr3[k*pnp+jj];
                        ptr3[ii*pnp+jj]=sum;
                    }
            }
        }

#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing A and B took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif

        /* compute U_j = \sum_i A_ij^T A_ij */ // \Sigma here!
        /* U_j is symmetric, therefore its computation can be sped up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that A_ij is mnp x cnp
         */
        /* Also compute ea_j = \sum_i A_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */

#ifdef TIMINGS
        start = clock();
#endif

        _dblzero(U, m*Usz); /* clear all U_j */
        _dblzero(ea, m*easz); /* clear all ea_j */
        for(j=mcon; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            ptr2=ea + j*easz; // set ptr2 to point to ea_j

            nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
            for(i=0; i<nnz; ++i){
                /* set ptr3 to point to A_ij, actual row number in rcsubs[i] */
                ptr3=jac + idxij.val[rcidxs[i]]*ABsz;

                /* compute the UPPER TRIANGULAR PART of A_ij^T A_ij and add it to U_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=ii; jj<cnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*cnp+ii]*ptr3[k*cnp+jj];
                        ptr1[ii*cnp+jj]+=sum;
                    }

                    /* copy the LOWER TRIANGULAR PART of U_j from the upper one */
                    for(jj=0; jj<ii; ++jj)
                        ptr1[ii*cnp+jj]=ptr1[jj*cnp+ii];
                }

                ptr4=e + idxij.val[rcidxs[i]]*esz; /* set ptr4 to point to e_ij */
                /* compute A_ij^T e_ij and add it to ea_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*cnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }

            /* Add in the camera constraints */
            if (use_constraints) {
                for (jj=0; jj < cnp; jj++) {
                    if (constraints[j].constrained[jj]) {
                        double diff = 
                            constraints[j].constraints[jj] - p[j * cnp + jj];
                        /* Add to the U matrix */
                        ptr1[jj*cnp+jj] += /*nnz **/ constraints[j].weights[jj];
                        ptr2[jj] += /*nnz **/ constraints[j].weights[jj] * diff;
                    }
                }
            }
			
			//printf("***********U_%d***********\n", j);
			//for(ii=0; ii<cnp; ii++)
			//{
			//	for(jj=0; jj<cnp; jj++)
			//	{
			//		printf("%f ", ptr1[ii*cnp+jj]);
			//	}
			//	printf("\n");
			//}
        }
        
#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing U took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif



        /* compute V_i = \sum_j B_ij^T B_ij */ // \Sigma here!
        /* V_i is symmetric, therefore its computation can be sped up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that B_ij is mnp x pnp
         */
        /* Also compute eb_i = \sum_j B_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */

#ifdef TIMINGS
        start = clock();
#endif

        _dblzero(V, n*Vsz); /* clear all V_i */
        _dblzero(eb, n*ebsz); /* clear all eb_i */
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=eb + i*ebsz; // set ptr2 to point to eb_i

            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr3 to point to B_ij, actual column number in rcsubs[j] */
                ptr3=jac + idxij.val[rcidxs[j]]*ABsz + Asz;
      
                /* compute the UPPER TRIANGULAR PART of B_ij^T B_ij and add it to V_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=ii; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*pnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]+=sum;
                    }
                }

                ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
                /* compute B_ij^T e_ij and add it to eb_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*pnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }

            /* Add in the constraints */
            if (use_point_constraints && point_constraints[i].constrained) {
                for (ii=0; ii < pnp; ii++) {
                    double diff = 
                        point_constraints[i].constraints[ii] - 
                        p[m * cnp + i * pnp + ii];
                    
                    // printf("  p[%d][%d] diff: %0.3f\n", i, ii, diff);
                    
                    /* Add to the U matrix */
                    ptr1[ii*pnp+ii] += nvis * point_constraints[i].weight;
                    ptr2[ii] += nvis * point_constraints[i].weight * diff;
                }
            }
			//printf("***********V_%d***********\n", i);
			//for(ii=0; ii<pnp; ii++)
			//{
			//	for(jj=ii; jj<pnp; jj++)
			//	{
			//		printf("%f ", ptr1[ii*pnp+jj]);
			//	}
			//	printf("\n");
			//}
        }
        
#ifdef TIMINGS
        end = clock();
        printf("[sba_motstr_levmar_x] computing V took %0.3fs\n", 
               (end - start) / (float) CLOCKS_PER_SEC);
#endif


        /* Copy out V to Vout */
        if (Vout != NULL) {
            memcpy(Vout, V, n*Vsz*sizeof(double));
            // not quite, this is just the upper triangle
            for(i=0; i<n; ++i){
                for(ii=0; ii<pnp; ++ii){
                    for(jj=0; jj<ii; ++jj){
                        Vout[i*Vsz + ii*pnp + jj] = 
                            Vout[i*Vsz + jj*pnp + ii];
                    }
                }
            }            
        }
        
        /* compute W_ij =  A_ij^T B_ij */ // \Sigma here!
        /* Recall that A_ij is mnp x cnp and B_ij is mnp x pnp
         */
        for(i=0; i<n; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1 to point to W_ij, actual column number in rcsubs[j] */
                ptr1=W + idxij.val[rcidxs[j]]*Wsz;

                if(rcsubs[j]<mcon){ /* A_ij is zero */
                    _dblzero(ptr1, Wsz); /* clear W_ij */
                    continue;
                }

                /* set ptr2 & ptr3 to point to A_ij & B_ij resp. */
                ptr2=jac  + idxij.val[rcidxs[j]]*ABsz;
                ptr3=ptr2 + Asz;
                /* compute A_ij^T B_ij and store it in W_ij
                 * Recall that storage for A_ij, B_ij does not overlap with that for W_ij,
                 * see the comments related to the initialization of jac above
                 */
                /* assert(ptr2-ptr1>=Wsz); */
                for(ii=0; ii<cnp; ++ii)
                    for(jj=0; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr2[k*cnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]=sum;
                    }
				
				//printf("***********W%d%d***********\n", i, rcsubs[j]);
				//for(ii=0; ii<cnp; ii++)
				//{
				//	for(jj=0; jj<pnp; jj++)
				//	{
				//		printf("%f ", ptr1[ii*pnp+jj]);
				//	}
				//	printf("\n");
				//}
            }
        }

        /* Compute ||J^T e||_inf and ||p||^2 */
        for(i=0, p_L2=eab_inf=0.0; i<nvars; ++i){
            if(eab_inf < (tmp=FABS(eab[i]))) eab_inf=tmp;
            p_L2+=p[i]*p[i];
        }
        //p_L2=sqrt(p_L2);

        /* save diagonal entries so that augmentation can be later canceled.
         * Diagonal entries are in U_j and V_i
         */
        for(j=mcon; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
            for(i=0; i<cnp; ++i)
                ptr2[i]=ptr1[i*cnp+i];
        }
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
            for(j=0; j<pnp; ++j)
                ptr2[j]=ptr1[j*pnp+j];
        }

        /*
          if(!(itno%100)){
          printf("Current estimate: ");
          for(i=0; i<nvars; ++i)
          printf("%.9g ", p[i]);
          printf("-- errors %.9g %0.9g\n", eab_inf, p_eL2);
          }
        */

        /* check for convergence */
        if((eab_inf <= eps1)){
            dp_L2=0.0; /* no increment for p in this case */
            stop=1;
            break;
        }

        /* compute initial damping factor */
        if(itno==0){
            for(i=mcon*cnp, tmp=DBL_MIN; i<nvars; ++i)
                if(diagUV[i]>tmp) tmp=diagUV[i]; /* find max diagonal element */
            mu=tau*tmp;
        }

		printf("damping factor:%f\n", mu);
        /* determine increment using adaptive damping */
        while(1){

#ifdef TIMINGS
            start = clock();
#endif

            /* augment U, V */
            for(j=mcon; j<m; ++j){
                ptr1=U + j*Usz; // set ptr1 to point to U_j
                for(i=0; i<cnp; ++i)
                    ptr1[i*cnp+i]+=mu;
            }
            for(i=0; i<n; ++i){
                ptr1=V + i*Vsz; // set ptr1 to point to V_i
                for(j=0; j<pnp; ++j)
                    ptr1[j*pnp+j]+=mu;

                /* compute V*_i^-1.
                 * Recall that only the upper triangle of the symmetric pnp x pnp matrix V*_i
                 * is stored in ptr1; its (also symmetric) inverse is saved in the lower triangle of ptr1
                 */
                /* inverting V*_i with LDLT seems to result in faster overall execution compared to when using LU or Cholesky */
                //j=sba_symat_invert_LU(ptr1, pnp); matinv=sba_symat_invert_LU;
                //j=sba_symat_invert_Chol(ptr1, pnp); matinv=sba_symat_invert_Chol;
                j=sba_symat_invert_BK(ptr1, pnp); matinv=sba_symat_invert_BK;
                if(!j){
                    fprintf(stderr, "SBA: singular matrix V*_i (i=%d) in sba_motstr_levmar_x(), increasing damping\n", i);
                    goto moredamping; // increasing damping will eventually make V*_i diagonally dominant, thus nonsingular
                    //retval=SBA_ERROR;
                    //goto freemem_and_return;
                }
            }

#ifdef TIMINGS
            end = clock();
            printf("[sba_motstr_levmar_x] computing damping factor took %0.3fs\n", 
                   (end - start) / (float) CLOCKS_PER_SEC);
#endif

            _dblzero(E, m*easz); /* clear all e_j */
            /* compute the mmcon x mmcon block matrix S and e_j */

            /* Note that S is symmetric, therefore its computation can be
             * sped up by computing only the upper part and then reusing
             * it for the lower one.
             */

#ifdef TIMINGS
            start = clock();
#endif

            for(j=mcon; j<m; ++j){
                int mmconxUsz=mmcon*Usz;

                nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero Y_ij, i=0...n-1 */

                /* compute all Y_ij = W_ij (V*_i)^-1 for a *fixed* j.
                 * To save memory, the block matrix consisting of the Y_ij
                 * is not stored. Instead, only a block column of this matrix
                 * is computed & used at each time: For each j, all nonzero
                 * Y_ij are computed in Yj and then used in the calculations
                 * involving S_jk and e_j.
                 * Recall that W_ij is cnp x pnp and (V*_i) is pnp x pnp
                 */
                for(i=0; i<nnz; ++i){
                    /* set ptr3 to point to (V*_i)^-1, actual row number in rcsubs[i] */
                    ptr3=V + rcsubs[i]*Vsz;

                    /* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
                    ptr1=Yj + i*Ysz;
                    /* set ptr2 to point to W_ij resp. */
                    ptr2=W + idxij.val[rcidxs[i]]*Wsz;
                    /* compute W_ij (V*_i)^-1 and store it in Y_ij.
                     * Recall that only the lower triangle of (V*_i)^-1 is stored
                     */
                    for(ii=0; ii<cnp; ++ii){
                        ptr4=ptr2+ii*pnp;
                        for(jj=0; jj<pnp; ++jj){
                            for(k=0, sum=0.0; k<=jj; ++k)
                                sum+=ptr4[k]*ptr3[jj*pnp+k]; //ptr2[ii*pnp+k]*ptr3[jj*pnp+k];
                            for( ; k<pnp; ++k)
                                sum+=ptr4[k]*ptr3[k*pnp+jj]; //ptr2[ii*pnp+k]*ptr3[k*pnp+jj];
                            ptr1[ii*pnp+jj]=sum;
                        }
                    }
					//printf("***********Y%d%d***********\n", i, j);
					//for(ii=0; ii<cnp; ii++)
					//{
					//	for(jj=0; jj<pnp; jj++)
					//	{
					//		printf("%f ", ptr1[ii*pnp+jj]);
					//	}
					//	printf("\n");
					//}
                }

                /* compute the UPPER TRIANGULAR PART of S */
                for(k=j; k<m; ++k){ // j>=mcon
                    /* compute \sum_i Y_ij W_ik^T in YWt. Note that
                     * for an off-diagonal block defined by j, k YWt
                     * (and thus S_jk) is nonzero only if there exists
                     * a point that is visible in both the j-th and
                     * k-th images
                     */
          
                    /* Recall that Y_ij is cnp x pnp and W_ik is 
                     * cnp x pnp */ 
                    _dblzero(YWt, YWtsz); /* clear YWt */

                    for(i=0; i<nnz; ++i){
                        register double *pYWt;

                        /* find the min and max column indices of the elements in row i (actually rcsubs[i])
                         * and make sure that k falls within them. This test handles W_ik's which are
                         * certain to be zero without bothering to call sba_crsm_elmidx()
                         */
                        ii=idxij.colidx[idxij.rowptr[rcsubs[i]]];
                        jj=idxij.colidx[idxij.rowptr[rcsubs[i]+1]-1];
                        if(k<ii || k>jj) continue; /* W_ik == 0 */

                        /* set ptr2 to point to W_ik */
                        l=sba_crsm_elmidxp(&idxij, rcsubs[i], k, j, rcidxs[i]);
                        //l=sba_crsm_elmidx(&idxij, rcsubs[i], k);
                        if(l==-1) continue; /* W_ik == 0 */

                        ptr2=W + idxij.val[l]*Wsz;
                        /* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
                        ptr1=Yj + i*Ysz;

#if 0
                        matrix_product_ipp(cnp, cnp, pnp, ptr1, ptr2, YWt);
#else
                        for(ii=0; ii<cnp; ++ii){
                            ptr3=ptr1+ii*pnp;
                            pYWt=YWt+ii*cnp;

                            ptr4=ptr2;
                            for(jj=0; jj<cnp; ++jj){
                                //ptr4=ptr2+jj*pnp;
                                for(l=0, sum=0.0; l<pnp; ++l)
                                    sum+=ptr3[l]*ptr4[l]; //ptr1[ii*pnp+l]*ptr2[jj*pnp+l];
                                pYWt[jj]+=sum; //YWt[ii*cnp+jj]+=sum;
                                ptr4+=pnp;
                            }
                        }
#endif
                    }
					//printf("***********YWt***********\n");
					//for(ii=0; ii<cnp; ii++)
					//{
					//	for(jj=0; jj<cnp; jj++)
					//	{
					//		printf("%f ", YWt[ii*cnp+jj]);
					//	}
					//	printf("\n");
					//}
		  
                    /* since the linear system involving S is solved with lapack,
                     * it is preferable to store S in column major (i.e. fortran)
                     * order, so as to avoid unecessary transposing/copying.
                     */
#if MAT_STORAGE==COLUMN_MAJOR
                    ptr2=S + (k-mcon)*mmconxUsz + (j-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#else
                    ptr2=S + (j-mcon)*mmconxUsz + (k-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#endif
		  
                    if(j!=k){ /* Kronecker */
                        for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
                            for(jj=0; jj<cnp; ++jj)
                                ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
                                    -YWt[jj*cnp+ii];
#else
                        -YWt[ii*cnp+jj];
#endif
                    }
                    else{
                        ptr1=U + j*Usz; // set ptr1 to point to U_j

                        for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
                            for(jj=0; jj<cnp; ++jj)
                                ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
                                    ptr1[jj*cnp+ii] - YWt[jj*cnp+ii];
#else
                        ptr1[ii*cnp+jj] - YWt[ii*cnp+jj];
#endif
                    }
                }

                /* copy the LOWER TRIANGULAR PART of S from the upper one */
                for(k=mcon; k<j; ++k){
#if MAT_STORAGE==COLUMN_MAJOR
                    ptr1=S + (k-mcon)*mmconxUsz + (j-mcon)*cnp; // set ptr1 to point to the beginning of block j,k in S
                    ptr2=S + (j-mcon)*mmconxUsz + (k-mcon)*cnp; // set ptr2 to point to the beginning of block k,j in S
#else
                    ptr1=S + (j-mcon)*mmconxUsz + (k-mcon)*cnp; // set ptr1 to point to the beginning of block j,k in S
                    ptr2=S + (k-mcon)*mmconxUsz + (j-mcon)*cnp; // set ptr2 to point to the beginning of block k,j in S
#endif
                    for(ii=0; ii<cnp; ++ii, ptr1+=Sdim)
                        for(jj=0, ptr3=ptr2+ii; jj<cnp; ++jj, ptr3+=Sdim)
                            ptr1[jj]=*ptr3;
                }

                /* compute e_j=ea_j - \sum_i Y_ij eb_i */
                /* Recall that Y_ij is cnp x pnp and eb_i is pnp x 1 */
                ptr1=E + j*easz; // set ptr1 to point to e_j

                for(i=0; i<nnz; ++i){
                    /* set ptr2 to point to Y_ij, actual row number in rcsubs[i] */
                    ptr2=Yj + i*Ysz;

                    /* set ptr3 to point to eb_i */
                    ptr3=eb + rcsubs[i]*ebsz;
                    for(ii=0; ii<cnp; ++ii){
                        ptr4=ptr2+ii*pnp;
                        for(jj=0, sum=0.0; jj<pnp; ++jj)
                            sum+=ptr4[jj]*ptr3[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
                        ptr1[ii]+=sum;
                    }
                }

                ptr2=ea + j*easz; // set ptr2 to point to ea_j
                for(i=0; i<easz; ++i)
                    ptr1[i]=ptr2[i] - ptr1[i];
            }


#ifdef TIMINGS
            end = clock();
            printf("[sba_motstr_levmar_x] computing S and e_j took %0.3fs\n", 
                   (end - start) / (float) CLOCKS_PER_SEC);
#endif


#if 0
            if(verbose>1){ /* count the nonzeros in S */
                for(i=ii=0; i<Sdim*Sdim; ++i)
                    if(S[i]!=0.0) ++ii;
                printf("\nS density: %.5g\n", ((double)ii)/(Sdim*Sdim)); fflush(stdout);
            }
#endif

            /* solve the linear system S dpa = E to compute the da_j.
             *
             * Note that if MAT_STORAGE==1 S is modified in the following call;
             * this is OK since S is recomputed for each iteration
             */

#ifdef TIMINGS
            start = clock();
#endif

	    //issolved=sba_Axb_LU(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_LU;
            issolved=sba_Axb_Chol(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_Chol;
            //issolved=sba_Axb_BK(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_BK;
            //issolved=sba_Axb_QRnoQ(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_QRnoQ;
            //issolved=sba_Axb_QR(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_QR;
	    //issolved=sba_Axb_SVD(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, MAT_STORAGE); linsolver=sba_Axb_SVD;
	    //issolved=sba_Axb_CG(S, E+mcon*cnp, dpa+mcon*cnp, Sdim, (3*Sdim)/2, 1E-10, SBA_CG_JACOBI, MAT_STORAGE); linsolver=(PLS)sba_Axb_CG;

            ++nlss;

	    _dblzero(dpa, mcon*cnp); /* no change for the first mcon camera params */

#ifdef TIMINGS
            end = clock();
            printf("[sba_motstr_levmar_x] solving linear system took %0.3fs\n", 
                   (end - start) / (float) CLOCKS_PER_SEC);
#endif

            if(issolved){

                /* compute the db_i */

#ifdef TIMINGS
                start = clock();
#endif

                for(i=0; i<n; ++i){
                    ptr1=dpb + i*ebsz; // set ptr1 to point to db_i

                    /* compute \sum_j W_ij^T da_j */
                    /* Recall that W_ij is cnp x pnp and da_j is cnp x 1 */
                    _dblzero(Wtda, Wtdasz); /* clear Wtda */
                    nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
                    for(j=0; j<nnz; ++j){
                        /* set ptr2 to point to W_ij, actual column number in rcsubs[j] */
                        if(rcsubs[j]<mcon) continue; /* W_ij is zero */

                        ptr2=W + idxij.val[rcidxs[j]]*Wsz;

                        /* set ptr3 to point to da_j */
                        ptr3=dpa + rcsubs[j]*cnp;

                        for(ii=0; ii<pnp; ++ii){
                            ptr4=ptr2+ii;
                            for(jj=0, sum=0.0; jj<cnp; ++jj)
                                sum+=ptr4[jj*pnp]*ptr3[jj]; //ptr2[jj*pnp+ii]*ptr3[jj];
                            Wtda[ii]+=sum;
                        }
                    }

                    /* compute eb_i - \sum_j W_ij^T da_j = eb_i - Wtda in Wtda */
                    ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
                    for(ii=0; ii<pnp; ++ii)
                        Wtda[ii]=ptr2[ii] - Wtda[ii];

                    /* compute the product (V*_i)^-1 Wtda = (V*_i)^-1 (eb_i - \sum_j W_ij^T da_j).
                     * Recall that only the lower triangle of (V*_i)^-1 is stored
                     */
                    ptr2=V + i*Vsz; // set ptr2 to point to (V*_i)^-1
                    for(ii=0; ii<pnp; ++ii){
                        for(jj=0, sum=0.0; jj<=ii; ++jj)
                            sum+=ptr2[ii*pnp+jj]*Wtda[jj];
                        for( ; jj<pnp; ++jj)
                            sum+=ptr2[jj*pnp+ii]*Wtda[jj];
                        ptr1[ii]=sum;
                    }
                }

#ifdef TIMINGS
                end = clock();
                printf("[sba_motstr_levmar_x] computing db_i's took %0.3fs\n", 
                       (end - start) / (float) CLOCKS_PER_SEC);
#endif

                /* parameter vector updates are now in dpa, dpb */

                /* compute p's new estimate and ||dp||^2 */
                for(i=0, dp_L2=0.0; i<nvars; ++i){
                    pdp[i]=p[i] + (tmp=dp[i]);
                    dp_L2+=tmp*tmp;
                }
                //dp_L2=sqrt(dp_L2);

                if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
                    //if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
                    stop=2;
                    break;
                }

                if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ){ /* almost singular */
                    //if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
                    fprintf(stderr, "SBA: the matrix of the augmented normal equations is almost singular in sba_motstr_levmar_x(),\n"
                            "     minimization should be restarted from the current solution with an increased damping term\n");
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }

#ifdef TIMINGS
                start = clock();
#endif

				/*for(i=0;i<nvars;i++)
				{
					printf("nvar:%f\n", pdp[i]);
				}*/
				
                (*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
                if(verbose>1) {
                    printf("mean reprojection error %g\n", sba_mean_repr_error(n, mnp, x, hx, &idxij, rcidxs, rcsubs));
                    fflush(stdout);
                }

                /* ### compute ||e(pdp)||_2 */
                if(covx==NULL)
                    pdp_eL2=nrmL2xmy(hx, x, hx, nobs); /* hx=x-hx, pdp_eL2=||hx|| */
                else
                    pdp_eL2=nrmCxmy(hx, x, hx, wght, mnp, nvis); /* hx=wght*(x-hx), pdp_eL2=||hx|| */
                if(!SBA_FINITE(pdp_eL2)){
                    if(verbose) /* identify the offending point projection */
                        sba_print_inf(hx, m, mnp, &idxij, rcidxs, rcsubs);

                    stop=7;
                    break;
                }

                /* Add in the camera constraints */
                if (use_constraints) {
                    for (j = 0; j < m; j++) {
                        int nnz = 0;

                        for(i = 0; i < n; i++)
                            nnz += vmask[i * m + j];

                        for (jj = 0; jj < cnp; jj++) {
                            if (constraints[j].constrained[jj]) {
                                double diff = 
                                    constraints[j].constraints[jj] - p[j * cnp + jj];
			
                                pdp_eL2 += 
                                    /*nnz **/ constraints[j].weights[jj] * diff * diff;
                            }
                        }
                    }
                }

                /* Add in the point constraints */
                if (use_point_constraints) {
                    for (i = 0; i < n; i++) {
                        if (point_constraints[i].constrained) {
                            for (ii=0; ii < pnp; ii++) {
                                double diff = 
                                    point_constraints[i].constraints[ii] - 
                                    p[m * cnp + i * pnp + ii];

                                /* Add to the U matrix */
                                pdp_eL2 += 
                                    nvis * point_constraints[i].weight * diff * diff;
                            }
                        }
                    }
                }

                for(i=0, dL=0.0; i<nvars; ++i)
                    dL+=dp[i]*(mu*dp[i]+eab[i]);

                dF=p_eL2-pdp_eL2;

#ifdef TIMINGS
                end = clock();
                printf("[sba_motstr_levmar_x] fn eval took %0.3fs\n", 
                       (end - start) / (float) CLOCKS_PER_SEC);
#endif

                if(verbose>1) {
                    printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);
                    printf("pdp_eL2: %0.3f, nvis: %d\n", pdp_eL2, nvis);

                    fflush(stdout);
                }


                if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
                    double max_pct_change = 0.0;

                    tmp=(2.0*dF/dL-1.0);
                    tmp=1.0-tmp*tmp*tmp;
                    mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
                    nu=2;

                    /* Check if we need to stop */
                    for (i=0; i < nobs; i++) {
                        double pct_change;
                        if (e[i] < eps5 && hx[i] < eps5) {
                            pct_change = 0.0;
                        } else {
                            double pct_change = fabs((e[i] - hx[i]) / e[i]);
                            if (pct_change > max_pct_change)
                                max_pct_change = pct_change;
                        }
                    }

                    printf("max_pct_change: %0.3e\n", max_pct_change);
                    fflush(stdout);

                    /* the test below is equivalent to the relative reduction of the RMS reprojection error: sqrt(p_eL2)-sqrt(pdp_eL2)<eps4*sqrt(p_eL2) */
                    if(pdp_eL2-2.0*sqrt(p_eL2*pdp_eL2)<(eps4_sq-1.0)*p_eL2) stop=4;
          
                    if (max_pct_change < eps5 && itno >= 4) { /* Run at least 4 iters */
                        stop=8;
                        break;
                    }

                    for(i=0; i<nvars; ++i) /* update p's estimate */
                        p[i]=pdp[i];

                    for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
                        e[i]=hx[i];
                    p_eL2=pdp_eL2;
                    break;
                }
            } /* issolved */

        moredamping:
            /* if this point is reached (also via an explicit goto!), either the linear system could
             * not be solved or the error did not reduce; in any case, the increment must be rejected
             */

            mu*=nu;
            nu2=nu<<1; // 2*nu;
            if(nu2<=nu){ /* nu has wrapped around (overflown) */
                fprintf(stderr, "SBA: too many failed attempts to increase the damping factor in sba_motstr_levmar_x()! Singular Hessian matrix?\n");
                //exit(1);
                stop=6;
                break;
            }
            nu=nu2;

            /* restore U, V diagonal entries */
            for(j=mcon; j<m; ++j){
                ptr1=U + j*Usz; // set ptr1 to point to U_j
                ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
                for(i=0; i<cnp; ++i)
                    ptr1[i*cnp+i]=ptr2[i];
            }
            for(i=0; i<n; ++i){
                ptr1=V + i*Vsz; // set ptr1 to point to V_i
                ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
                for(j=0; j<pnp; ++j)
                    ptr1[j*pnp+j]=ptr2[j];
            }
        } /* inner while loop */

        if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
    }

    if(itno>=itmax) stop=3;

    /* restore U, V diagonal entries */
    for(j=mcon; j<m; ++j){
        ptr1=U + j*Usz; // set ptr1 to point to U_j
        ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
        for(i=0; i<cnp; ++i)
            ptr1[i*cnp+i]=ptr2[i];
    }
    for(i=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
        for(j=0; j<pnp; ++j)
            ptr1[j*pnp+j]=ptr2[j];
    }

    if (Sout != NULL) {
        /* compute derivative submatrices A_ij, B_ij */
        (*fjac)(p, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;

        /* compute U_j = \sum_i A_ij^T A_ij */ // \Sigma here!
        /* U_j is symmetric, therefore its computation can be speeded up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that A_ij is mnp x cnp
         */
        /* Also compute ea_j = \sum_i A_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */
        _dblzero(U, m*Usz); /* clear all U_j */
        _dblzero(ea, m*easz); /* clear all ea_j */
        for(j=mcon; j<m; ++j) {
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            ptr2=ea + j*easz; // set ptr2 to point to ea_j
          
            nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
            for(i=0; i<nnz; ++i){
                /* set ptr3 to point to A_ij, actual row number in rcsubs[i] */
                ptr3=jac + idxij.val[rcidxs[i]]*ABsz;

                /* compute the UPPER TRIANGULAR PART of A_ij^T A_ij and add it to U_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=ii; jj<cnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*cnp+ii]*ptr3[k*cnp+jj];
                        ptr1[ii*cnp+jj]+=sum;
                    }

                    /* copy the LOWER TRIANGULAR PART of U_j from the upper one */
                    for(jj=0; jj<ii; ++jj)
                        ptr1[ii*cnp+jj]=ptr1[jj*cnp+ii];
                }
              
                ptr4=e + idxij.val[rcidxs[i]]*esz; /* set ptr4 to point to e_ij */
                /* compute A_ij^T e_ij and add it to ea_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*cnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }
      
            /* Add in the camera constraints */
            if (use_constraints) {
                for (jj=0; jj < cnp; jj++) {
                    if (constraints[j].constrained[jj]) {
                        double diff = 
                            constraints[j].constraints[jj] - p[j * cnp + jj];
                        /* Add to the U matrix */
                        ptr1[jj*cnp+jj] += /*nnz **/ constraints[j].weights[jj];
                        ptr2[jj] += /*nnz **/ constraints[j].weights[jj] * diff;
                    }
                }
            }
        }

        if (Uout != NULL) 
            memcpy(Uout, U, m*Usz*sizeof(double));

        /* compute V_i = \sum_j B_ij^T B_ij */ // \Sigma here!
        /* V_i is symmetric, therefore its computation can be speeded up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that B_ij is mnp x pnp
         */
        /* Also compute eb_i = \sum_j B_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */
        _dblzero(V, n*Vsz); /* clear all V_i */
        _dblzero(eb, n*ebsz); /* clear all eb_i */
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
          
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr3 to point to B_ij, actual column number in rcsubs[j] */
                ptr3=jac + idxij.val[rcidxs[j]]*ABsz + Asz;
      
                /* compute the UPPER TRIANGULAR PART of B_ij^T B_ij and add it to V_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=ii; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*pnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]+=sum;
                    }
                }
              
                ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
                /* compute B_ij^T e_ij and add it to eb_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*pnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }

            /* Add in the constraints */
            if (use_point_constraints && point_constraints[i].constrained) {
                for (ii=0; ii < pnp; ii++) {
                    double diff = 
                        point_constraints[i].constraints[ii] - 
                        p[m * cnp + i * pnp + ii];

                    // printf("  p[%d][%d] diff: %0.3f\n", i, ii, diff);

                    /* Add to the U matrix */
                    ptr1[ii*pnp+ii] += nvis * point_constraints[i].weight;
                    ptr2[ii] += nvis * point_constraints[i].weight * diff;
                }
            }
        }

        /* Copy out V to Vout */
        if (Vout != NULL) {
            memcpy(Vout, V, n*Vsz*sizeof(double));
            // not quite, this is just the upper triangle
            for(i=0; i<n; ++i){
                for(ii=0; ii<pnp; ++ii){
                    for(jj=0; jj<ii; ++jj){
                        Vout[i*Vsz + ii*pnp + jj] = 
                            Vout[i*Vsz + jj*pnp + ii];
                    }
                }
            }   
        }

        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i

            /* compute V*_i^-1.
             * Recall that only the upper triangle of the symmetric pnp x pnp matrix V*_i
             * is stored in ptr1; its (also symmetric) inverse is saved in the lower triangle of ptr1
             */
            /* inverting V*_i with LDLT seems to result in faster overall execution compared to when using LU or Cholesky */
            //j=sba_symat_invert_LU(ptr1, pnp); matinv=sba_symat_invert_LU;
            //j=sba_symat_invert_Chol(ptr1, pnp); matinv=sba_symat_invert_Chol;
            j=sba_symat_invert_BK(ptr1, pnp); matinv=sba_symat_invert_BK;
          
            if(!j){
                fprintf(stderr, "Singular matrix V*_i (i=%d) in sba_motstr_levmar_x()\n", i);
                // exit(1);
                for (ii=0; ii<pnp; ii++) {
                    for (jj=0; jj<pnp; jj++) {
                        if (ii == jj) ptr2[ii*pnp+jj] = 100.0;
                        else ptr2[ii*pnp+jj] = 0.0;
                    }
                }
            }
        }

#if 0
        /* compute Y_ij = W_ij (V*_i)^-1
         * Recall that W_ij is cnp x pnp and (V*_i) is pnp x pnp
         */
        memset(Y, 0, nvis*Ysz*sizeof(double)); /* clear all Y_ij */
        for(i=0; i<n; ++i){
            /* set ptr3 to point to (V*_i)^-1 */
            ptr3=V_1 + i*Vsz;
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1 to point to Y_ij, actual column number in rcsubs[j] */
                if(rcsubs[j]<mcon) continue; /* W_ij is zero */

                ptr1=Y + idxij.val[rcidxs[j]]*Ysz;
                /* set ptr2 to point to W_ij resp. */
                ptr2=W + idxij.val[rcidxs[j]]*Wsz;
                /* compute W_ij (V*_i)^-1 and store it in Y_ij */
                for(ii=0; ii<cnp; ++ii)
                    for(jj=0; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<pnp; ++k)
                            sum+=ptr2[ii*pnp+k]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]=sum;
                    }
            }
        }
#endif

        /* compute W_ij =  A_ij^T B_ij */ // \Sigma here!
        /* Recall that A_ij is mnp x cnp and B_ij is mnp x pnp
         */
        for(i=0; i<n; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1 to point to W_ij, actual column number in rcsubs[j] */
                ptr1=W + idxij.val[rcidxs[j]]*Wsz;

                if(rcsubs[j]<mcon){ /* A_ij is zero */
                    _dblzero(ptr1, Wsz); /* clear W_ij */
                    continue;
                }

                /* set ptr2 & ptr3 to point to A_ij & B_ij resp. */
                ptr2=jac  + idxij.val[rcidxs[j]]*ABsz;
                ptr3=ptr2 + Asz;
                /* compute A_ij^T B_ij and store it in W_ij
                 * Recall that storage for A_ij, B_ij does not overlap with that for W_ij,
                 * see the comments related to the initialization of jac above
                 */
                /* assert(ptr2-ptr1>=Wsz); */
                for(ii=0; ii<cnp; ++ii) {
                    double *ptr1a = 
                        Wout + (rcsubs[j] * cnp + ii) * pnp * n + i * pnp;

                    for(jj=0; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr2[k*cnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]=sum;

                        if (Wout != NULL)
                            ptr1a[jj] = sum;
                    }
                }
            }
        }

#if 0
        if (Wout != NULL) {
            memset(Wout, 0, pnp * cnp * n * m * sizeof(double));
        }

        memset(W, 0, nvis*Wsz*sizeof(double)); /* clear all W_ij */
        for(i=0; i<n; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero W_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1 to point to W_ij, actual column number in rcsubs[j] */
                if(rcsubs[j]<mcon) continue; /* A_ij is zero */

                ptr1=W + idxij.val[rcidxs[j]]*Wsz;

                /* set ptr2 & ptr3 to point to A_ij & B_ij resp. */
                ptr2=jaca + idxij.val[rcidxs[j]]*Asz;
                ptr3=jacb + idxij.val[rcidxs[j]]*Bsz;
                /* compute A_ij^T B_ij and store it in W_ij */
                for(ii=0; ii<cnp; ++ii) {
                    double *ptr1a = 
                        Wout + (rcsubs[j] * cnp + ii) * pnp * n + i * pnp;

                    assert(rcsubs[j] >= 0);
                    assert(rcsubs[j] < m);

                    for(jj=0; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr2[k*cnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]=sum;

                        if (Wout != NULL)
                            ptr1a[jj] = sum;
                    }
                }
            }
        }
#endif

        /* Compute S with no damping */

        for(j=mcon; j<m; ++j){
            int mmconxUsz=mmcon*Usz;

            nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero Y_ij, i=0...n-1 */

            /* compute all Y_ij = W_ij (V*_i)^-1 for a *fixed* j.
             * To save memory, the block matrix consisting of the Y_ij
             * is not stored. Instead, only a block column of this matrix
             * is computed & used at each time: For each j, all nonzero
             * Y_ij are computed in Yj and then used in the calculations
             * involving S_jk and e_j.
             * Recall that W_ij is cnp x pnp and (V*_i) is pnp x pnp
             */
            for(i=0; i<nnz; ++i){
                /* set ptr3 to point to (V*_i)^-1, actual row number in rcsubs[i] */
                ptr3=V + rcsubs[i]*Vsz;
              
                /* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
                ptr1=Yj + i*Ysz;
                /* set ptr2 to point to W_ij resp. */
                ptr2=W + idxij.val[rcidxs[i]]*Wsz;
                /* compute W_ij (V*_i)^-1 and store it in Y_ij.
                 * Recall that only the lower triangle of (V*_i)^-1 is stored
                 */
                for(ii=0; ii<cnp; ++ii){
                    ptr4=ptr2+ii*pnp;
                    for(jj=0; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<=jj; ++k)
                            sum+=ptr4[k]*ptr3[jj*pnp+k]; //ptr2[ii*pnp+k]*ptr3[jj*pnp+k];
                        for( ; k<pnp; ++k)
                            sum+=ptr4[k]*ptr3[k*pnp+jj]; //ptr2[ii*pnp+k]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]=sum;
                    }
                }
            }
          
            /* compute the UPPER TRIANGULAR PART of S */
            for(k=j; k<m; ++k){ // j>=mcon
                /* compute \sum_i Y_ij W_ik^T in YWt. Note that for an off-diagonal block defined by j, k
                 * YWt (and thus S_jk) is nonzero only if there exists a point that is visible in both the
                 * j-th and k-th images
                 */
              
                /* Recall that Y_ij is cnp x pnp and W_ik is cnp x pnp */ 
                _dblzero(YWt, YWtsz); /* clear YWt */
              
                for(i=0; i<nnz; ++i){
                    register double *pYWt;
                  
                    /* find the min and max column indices of the elements in row i (actually rcsubs[i])
                     * and make sure that k falls within them. This test handles W_ik's which are
                     * certain to be zero without bothering to call sba_crsm_elmidx()
                     */
                    ii=idxij.colidx[idxij.rowptr[rcsubs[i]]];
                    jj=idxij.colidx[idxij.rowptr[rcsubs[i]+1]-1];
                    if(k<ii || k>jj) continue; /* W_ik == 0 */
                  
                    /* set ptr2 to point to W_ik */
                    l=sba_crsm_elmidxp(&idxij, rcsubs[i], k, j, rcidxs[i]);
                    //l=sba_crsm_elmidx(&idxij, rcsubs[i], k);
                    if(l==-1) continue; /* W_ik == 0 */
                  
                    ptr2=W + idxij.val[l]*Wsz;
                    /* set ptr1 to point to Y_ij, actual row number in rcsubs[i] */
                    ptr1=Yj + i*Ysz;
                    for(ii=0; ii<cnp; ++ii){
                        ptr3=ptr1+ii*pnp;
                        pYWt=YWt+ii*cnp;
                      
                        for(jj=0; jj<cnp; ++jj){
                            ptr4=ptr2+jj*pnp;
                            for(l=0, sum=0.0; l<pnp; ++l)
                                sum+=ptr3[l]*ptr4[l]; //ptr1[ii*pnp+l]*ptr2[jj*pnp+l];
                            pYWt[jj]+=sum; //YWt[ii*cnp+jj]+=sum;
                        }
                    }
                }
              
                /* since the linear system involving S is solved with lapack,
                 * it is preferable to store S in column major (i.e. fortran)
                 * order, so as to avoid unecessary transposing/copying.
                 */
#if MAT_STORAGE==COLUMN_MAJOR
                ptr2=S + (k-mcon)*mmconxUsz + (j-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#else
                ptr2=S + (j-mcon)*mmconxUsz + (k-mcon)*cnp; // set ptr2 to point to the beginning of block j,k in S
#endif
		  
                if(j!=k){ /* Kronecker */
                    for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
                        for(jj=0; jj<cnp; ++jj)
                            ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
                                -YWt[jj*cnp+ii];
#else
                    -YWt[ii*cnp+jj];
#endif
                } else {
                    ptr1=U + j*Usz; // set ptr1 to point to U_j
                  
                    for(ii=0; ii<cnp; ++ii, ptr2+=Sdim)
                        for(jj=0; jj<cnp; ++jj)
                            ptr2[jj]=
#if MAT_STORAGE==COLUMN_MAJOR
                                ptr1[jj*cnp+ii] - YWt[jj*cnp+ii];
#else
                    ptr1[ii*cnp+jj] - YWt[ii*cnp+jj];
#endif
                }
            }
          
            /* copy the LOWER TRIANGULAR PART of S from the upper one */
            for(k=mcon; k<j; ++k){
#if MAT_STORAGE==COLUMN_MAJOR
                ptr1=S + (k-mcon)*mmconxUsz + (j-mcon)*cnp; // set ptr1 to point to the beginning of block j,k in S
                ptr2=S + (j-mcon)*mmconxUsz + (k-mcon)*cnp; // set ptr2 to point to the beginning of block k,j in S
#else
                ptr1=S + (j-mcon)*mmconxUsz + (k-mcon)*cnp; // set ptr1 to point to the beginning of block j,k in S
                ptr2=S + (k-mcon)*mmconxUsz + (j-mcon)*cnp; // set ptr2 to point to the beginning of block k,j in S
#endif
                for(ii=0; ii<cnp; ++ii, ptr1+=Sdim)
                    for(jj=0, ptr3=ptr2+ii; jj<cnp; ++jj, ptr3+=Sdim)
                        ptr1[jj]=*ptr3;
            }      
        }
    
#if MAT_STORAGE==COLUMN_MAJOR
        for (i = 0; i < m * cnp; i++) {
            for (j = 0; j < m * cnp; j++) {
                Sout[i * m * cnp + j] = S[j * m * cnp + i];
            }
        }
#else
        memcpy(Sout, S, sizeof(double) * Ssz * m * m);
#endif
    }
    
    if(info){
        info[0]=init_p_eL2;
        info[1]=p_eL2;
        info[2]=eab_inf;
        info[3]=dp_L2;
        for(j=mcon, tmp=DBL_MIN; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            for(i=0; i<cnp; ++i)
                if(tmp<ptr1[i*cnp+i]) tmp=ptr1[i*cnp+i];
        }
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            for(j=0; j<pnp; ++j)
                if(tmp<ptr1[j*pnp+j]) tmp=ptr1[j*pnp+j];
        }
        info[4]=mu/tmp;
        info[5]=itno;
        info[6]=stop;
        info[7]=nfev;
        info[8]=njev;
        info[9]=nlss;
    }
                                                               
    //sba_print_sol(n, m, p, cnp, pnp, x, mnp, &idxij, rcidxs, rcsubs);
    retval=(stop!=7)?  itno : SBA_ERROR;

 freemem_and_return: /* NOTE: this point is also reached via a goto! */

    /* free whatever was allocated */
    free(W);   free(U);  free(V);
    free(e);   free(eab);
    free(E);   free(Yj); free(YWt);
    free(S);   free(dp); free(Wtda);
    free(rcidxs); free(rcsubs);
#ifndef SBA_DESTROY_COVS
    if(wght) free(wght);
#else
    /* nothing to do */
#endif /* SBA_DESTROY_COVS */

    free(hx); free(diagUV); free(pdp);
    if(fdj_data.hxx){ // cleanup
        free(fdj_data.hxx);
        free(fdj_data.func_rcidxs);
    }

    sba_crsm_free(&idxij);

    /* free the memory allocated by the matrix inversion & linear solver routines */
    if(matinv) (*matinv)(NULL, 0);
    if(linsolver) (*linsolver)(NULL, NULL, NULL, 0, 0);

    return retval;
}


/* Bundle adjustment on camera parameters only 
 * using the sparse Levenberg-Marquardt as described in HZ p. 568
 *
 * Returns the number of iterations (>=0) if successfull, SBA_ERROR if failed
 */

int sba_mot_levmar_x(
                     const int n,   /* number of points */
                     const int m,   /* number of images */
                     const int mcon,/* number of images (starting from the 1st) whose parameters should not be modified.
                                     * All A_ij (see below) with j<mcon are assumed to be zero
                                     */
                     char *vmask,  /* visibility mask: vmask[i, j]=1 if point i visible in image j, 0 otherwise. nxm */
                     double *p,    /* initial parameter vector p0: (a1, ..., am).
                                    * aj are the image j parameters, size m*cnp */
                     const int cnp,/* number of parameters for ONE camera; e.g. 6 for Euclidean cameras */
                     double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
                                    * x_ij is the projection of the i-th point on the j-th image.
                                    * NOTE: some of the x_ij might be missing, if point i is not visible in image j;
                                    * see vmask[i, j], max. size n*m*mnp
                                    */
                     double *covx, /* measurements covariance matrices: (Sigma_x_11, .. Sigma_x_1m, ..., Sigma_x_n1, .. Sigma_x_nm),
                                    * where Sigma_x_ij is the mnp x mnp covariance of x_ij stored row-by-row. Set to NULL if no
                                    * covariance estimates are available (identity matrices are implicitly used in this case).
                                    * NOTE: a certain Sigma_x_ij is missing if the corresponding x_ij is also missing;
                                    * see vmask[i, j], max. size n*m*mnp*mnp
                                    */
                     const int mnp,/* number of parameters for EACH measurement; usually 2 */
                     void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
                     /* functional relation describing measurements. Given a parameter vector p,
                      * computes a prediction of the measurements \hat{x}. p is (m*cnp)x1,
                      * \hat{x} is (n*m*mnp)x1, maximum
                      * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                      * as working memory
                      */
                     void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
                     /* function to evaluate the sparse jacobian dX/dp.
                      * The Jacobian is returned in jac as
                      * (dx_11/da_1, ..., dx_1m/da_m, ..., dx_n1/da_1, ..., dx_nm/da_m), or (using HZ's notation),
                      * jac=(A_11, ..., A_1m, ..., A_n1, ..., A_nm)
                      * Notice that depending on idxij, some of the A_ij might be missing.
                      * Note also that the A_ij are mnp x cnp matrices and they
                      * should be stored in jac in row-major order.
                      * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                      * as working memory
                      *
                      * If NULL, the jacobian is approximated by repetitive func calls and finite
                      * differences. This is computationally inefficient and thus NOT recommended.
                      */
                     void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 

                     const int itmax,   /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
                     const int verbose, /* I: verbosity */
                     const double opts[SBA_OPTSSZ],
                     /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \epsilon4]. Respectively the scale factor for initial \mu,
                      * stopping thresholds for ||J^T e||_inf, ||dp||_2, ||e||_2 and (||e||_2-||e_new||_2)/||e||_2
                      */
                     double info[SBA_INFOSZ],
                     /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]=||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small dp
                      *                                 3 - stopped by itmax
                      *                                 4 - stopped by small relative reduction in ||e||_2
                      *                                 5 - stopped by small ||e||_2
                      *                                 6 - too many attempts to increase damping. Restart with increased mu
                      *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
                      * info[7]= # function evaluations
                      * info[8]= # jacobian evaluations
                      * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
                      */
                     int use_constraints, camera_constraints_t *constraints  /* Constraints on camera parameters */
                     )
{
    register int i, j, ii, jj, k;
    int nvis, nnz, retval;

    /* The following are work arrays that are dynamically allocated by sba_mot_levmar_x() */
    double *jac; /* work array for storing the jacobian, max. size n*m*mnp*cnp */
    double *U;    /* work array for storing the U_j in the order U_1, ..., U_m, size m*cnp*cnp */

    double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_n1, ..., e_nm,
                     max. size n*m*mnp */
    double *ea;   /* work array for storing the ea_j in the order ea_1, .. ea_m, size m*cnp */

    double *dp;   /* work array for storing the parameter vector updates da_1, ..., da_m, size m*cnp */

    double *wght= /* work array for storing the weights computed from the covariance inverses, max. size n*m*mnp*mnp */
        NULL;

    /* Of the above arrays, jac, e, wght are sparse and
     * U, ea, dp are dense. Sparse arrays are indexed through
     * idxij (see below), that is with the same mechanism as the input 
     * measurements vector x
     */

    /* submatrices sizes */
    int Asz, Usz,
        esz, easz, covsz;

    register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
    struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also the location of A_ij 
                            * in jac, e_ij in e, etc.
                            * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                            * index k and it is used to efficiently lookup the memory locations where the non-zero
                            * blocks of a sparse matrix/vector are stored
                            */
    int maxCPvis, /* max. of projections across cameras & projections across points */
        *rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
                     column in a sparse matrix, size max(n, m) */
        *rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
                     sparse matrix, size max(n, m) */

    /* The following variables are needed by the LM algorithm */
    register int itno;  /* iteration counter */
    int nsolved;
    /* temporary work arrays that are dynamically allocated */
    double *hx,         /* \hat{x}_i, max. size m*n*mnp */
        *diagU,      /* diagonals of U_j, size m*cnp */
        *pdp;        /* p + dp, size m*cnp */

    register double mu,  /* damping constant */
        tmp; /* mainly used in matrix & vector multiplications */
    double p_eL2, ea_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
    double p_L2, dp_L2=DBL_MAX, dF, dL;
    double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2],
        eps3_sq=opts[3]*opts[3], eps4_sq=opts[4]*opts[4];
    double init_p_eL2;
    int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;
    int nobs, nvars;
    PLS linsolver=NULL;

    struct fdj_data_x_ fdj_data;
    void *jac_adata;

    /* Initialization */
    mu=ea_inf=0.0; /* -Wall */

    /* block sizes */
    Asz=mnp * cnp; Usz=cnp * cnp;
    esz=mnp; easz=cnp;
    covsz=mnp * mnp;
  
    /* count total number of visible image points */
    for(i=nvis=0, jj=n*m; i<jj; ++i)
        nvis+=(vmask[i]!=0);

    nobs=nvis*mnp;
    nvars=m*cnp;
    if(nobs<nvars){
        fprintf(stderr, "SBA: sba_mot_levmar_x() cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
        return SBA_ERROR;
    }

    /* allocate & fill up the idxij structure */
    sba_crsm_alloc(&idxij, n, m, nvis);
    for(i=k=0; i<n; ++i){
        idxij.rowptr[i]=k;
        ii=i*m;
        for(j=0; j<m; ++j)
            if(vmask[ii+j]){
                idxij.val[k]=k;
                idxij.colidx[k++]=j;
            }
    }
    idxij.rowptr[n]=nvis;

    /* find the maximum number of visible image points in any single camera or coming from a single 3D point */
    /* cameras */
    for(i=maxCPvis=0; i<n; ++i)
        if((k=idxij.rowptr[i+1]-idxij.rowptr[i])>maxCPvis) maxCPvis=k;

    /* points, note that maxCPvis is not reinitialized! */
    for(j=0; j<m; ++j){
        for(i=ii=0; i<n; ++i)
            if(vmask[i*m+j]) ++ii;
        if(ii>maxCPvis) maxCPvis=ii;
    }

    /* allocate work arrays */
    jac=(double *)emalloc(nvis*Asz*sizeof(double));
    U=(double *)emalloc(m*Usz*sizeof(double));
    e=(double *)emalloc(nobs*sizeof(double));
    ea=(double *)emalloc(nvars*sizeof(double));
    dp=(double *)emalloc(nvars*sizeof(double));
    rcidxs=(int *)emalloc(maxCPvis*sizeof(int));
    rcsubs=(int *)emalloc(maxCPvis*sizeof(int));
#ifndef SBA_DESTROY_COVS
    if(covx!=NULL) wght=(double *)emalloc(nvis*covsz*sizeof(double));
#else
    if(covx!=NULL) wght=covx;
#endif /* SBA_DESTROY_COVS */


    hx=(double *)emalloc(nobs*sizeof(double));
    diagU=(double *)emalloc(nvars*sizeof(double));
    pdp=(double *)emalloc(nvars*sizeof(double));

    /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
    if(!fjac){
        fdj_data.func=func;
        fdj_data.cnp=cnp;
        fdj_data.pnp=0;
        fdj_data.mnp=mnp;
        fdj_data.hx=hx;
        fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
        fdj_data.func_rcidxs=(int *)emalloc(2*maxCPvis*sizeof(int));
        fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxCPvis;
        fdj_data.adata=adata;

        fjac=sba_fdjac_x;
        jac_adata=(void *)&fdj_data;
    }
    else{
        fdj_data.hxx=NULL;
        jac_adata=adata;
    }

    if(itmax==0){ /* verify jacobian */
        sba_mot_chkjac_x(func, fjac, p, &idxij, rcidxs, rcsubs, mcon, cnp, mnp, adata, jac_adata);
        retval=0;
        goto freemem_and_return;
    }

    /* covariances Sigma_x_ij are accommodated by computing the Cholesky decompositions of their
     * inverses and using the resulting matrices w_x_ij to weigh A_ij and e_ij as w_x_ij A_ij
     * and w_x_ij*e_ij. In this way, auxiliary variables as U_j=\sum_i A_ij^T A_ij
     * actually become \sum_i (w_x_ij A_ij)^T w_x_ij A_ij= \sum_i A_ij^T w_x_ij^T w_x_ij A_ij =
     * A_ij^T Sigma_x_ij^-1 A_ij
     *
     * ea_j are weighted in a similar manner
     */
    if(covx!=NULL){
        for(i=0; i<n; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1, ptr2 to point to cov_x_ij, w_x_ij resp. */
                ptr1=covx + (k=idxij.val[rcidxs[j]]*covsz);
                ptr2=wght + k;
                if(!sba_mat_cholinv(ptr1, ptr2, mnp)){ /* compute w_x_ij s.t. w_x_ij^T w_x_ij = cov_x_ij^-1 */
                    fprintf(stderr, "SBA: invalid covariance matrix for x_ij (i=%d, j=%d) in sba_motstr_levmar_x()\n", i, rcsubs[j]);
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }
            }
        }
        sba_mat_cholinv(NULL, NULL, 0); /* cleanup */
    }

    /* compute the error vectors e_ij in hx */
    (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
    /* ### compute e=x - f(p) [e=w*(x - f(p)] and its L2 norm */
    if(covx==NULL)
        p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
    else
        p_eL2=nrmCxmy(e, x, hx, wght, mnp, nvis); /* e=wght*(x-hx), p_eL2=||e||=||x-hx||_Sigma^-1 */

    /* Add in the camera constraints */
    if (use_constraints) {
        for (j = 0; j < m; j++) {
            /* Find number of points projecting into this camera */
            int nnz = 0;

            for(i = 0; i < n; i++)
                nnz += vmask[i * m + j];

            for (jj = 0; jj < cnp; jj++) {
                if (constraints[j].constrained[jj]) {
                    double diff = 
                        constraints[j].constraints[jj] - p[j * cnp + jj];

                    p_eL2 += /*nnz **/ constraints[j].weights[jj] * diff * diff;
                }
            }
        }
    }

    if(verbose) printf("initial mot-SBA error %g [%g]\n", p_eL2, p_eL2/nvis);
    init_p_eL2=p_eL2;
    if(!SBA_FINITE(p_eL2)) stop=7;

    for(itno=0; itno<itmax && !stop; ++itno){
        /* Note that p, e and ||e||_2 have been updated at the previous iteration */

        /* compute derivative submatrices A_ij */
        (*fjac)(p, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;

        if(covx!=NULL){
            /* compute w_x_ij A_ij
             * Since w_x_ij is upper triangular, the products can be safely saved
             * directly in A_ij, without the need for intermediate storage
             */
            for(i=0; i<nvis; ++i){
                /* set ptr1, ptr2 to point to w_x_ij, A_ij, resp. */
                ptr1=wght + i*covsz;
                ptr2=jac  + i*Asz;

                /* w_x_ij is mnp x mnp, A_ij is mnp x cnp */
                for(ii=0; ii<mnp; ++ii)
                    for(jj=0; jj<cnp; ++jj){
                        for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
                            sum+=ptr1[ii*mnp+k]*ptr2[k*cnp+jj];
                        ptr2[ii*cnp+jj]=sum;
                    }
            }
        }

        /* compute U_j = \sum_i A_ij^T A_ij */ // \Sigma here!
        /* U_j is symmetric, therefore its computation can be sped up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that A_ij is mnp x cnp
         */
        /* Also compute ea_j = \sum_i A_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */
        _dblzero(U, m*Usz); /* clear all U_j */
        _dblzero(ea, m*easz); /* clear all ea_j */
        for(j=mcon; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            ptr2=ea + j*easz; // set ptr2 to point to ea_j

            nnz=sba_crsm_col_elmidxs(&idxij, j, rcidxs, rcsubs); /* find nonzero A_ij, i=0...n-1 */
            for(i=0; i<nnz; ++i){
                /* set ptr3 to point to A_ij, actual row number in rcsubs[i] */
                ptr3=jac + idxij.val[rcidxs[i]]*Asz;

                /* compute the UPPER TRIANGULAR PART of A_ij^T A_ij and add it to U_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=ii; jj<cnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*cnp+ii]*ptr3[k*cnp+jj];
                        ptr1[ii*cnp+jj]+=sum;
                    }

                    /* copy the LOWER TRIANGULAR PART of U_j from the upper one */
                    for(jj=0; jj<ii; ++jj)
                        ptr1[ii*cnp+jj]=ptr1[jj*cnp+ii];
                }

                ptr4=e + idxij.val[rcidxs[i]]*esz; /* set ptr4 to point to e_ij */
                /* compute A_ij^T e_ij and add it to ea_j */
                for(ii=0; ii<cnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*cnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }

            /* Add in the camera constraints */
            if (use_constraints) {
                for (jj=0; jj < cnp; jj++) {
                    if (constraints[j].constrained[jj]) {
                        double diff = 
                            constraints[j].constraints[jj] - p[j * cnp + jj];
                        /* Add to the U matrix */
                        ptr1[jj*cnp+jj] += /*nnz **/ constraints[j].weights[jj];
                        ptr2[jj] += /*nnz **/ constraints[j].weights[jj] * diff;
                    }
                }
            }
        }

        /* Compute ||J^T e||_inf and ||p||^2 */
        for(i=0, p_L2=ea_inf=0.0; i<nvars; ++i){
            if(ea_inf < (tmp=FABS(ea[i]))) ea_inf=tmp;
            p_L2+=p[i]*p[i];
        }
        //p_L2=sqrt(p_L2);

        /* save diagonal entries so that augmentation can be later canceled.
         * Diagonal entries are in U_j
         */
        for(j=mcon; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
            for(i=0; i<cnp; ++i)
                ptr2[i]=ptr1[i*cnp+i];
        }

        /*
          if(!(itno%100)){
          printf("Current estimate: ");
          for(i=0; i<nvars; ++i)
          printf("%.9g ", p[i]);
          printf("-- errors %.9g %0.9g\n", ea_inf, p_eL2);
          }
        */

        /* check for convergence */
        if((ea_inf <= eps1)){
            dp_L2=0.0; /* no increment for p in this case */
            stop=1;
            break;
        }

        /* compute initial damping factor */
        if(itno==0){
            for(i=mcon*cnp, tmp=DBL_MIN; i<nvars; ++i)
                if(diagU[i]>tmp) tmp=diagU[i]; /* find max diagonal element */
            mu=tau*tmp;
        }

        /* determine increment using adaptive damping */
        while(1){
            /* augment U */
            for(j=mcon; j<m; ++j){
                ptr1=U + j*Usz; // set ptr1 to point to U_j
                for(i=0; i<cnp; ++i)
                    ptr1[i*cnp+i]+=mu;
            }
 
            /* solve the linear systems U_j da_j = ea_j to compute the da_j */
            _dblzero(dp, mcon*cnp); /* no change for the first mcon camera params */
            for(j=nsolved=mcon; j<m; ++j){
                ptr1=U + j*Usz; // set ptr1 to point to U_j
                ptr2=ea + j*easz; // set ptr2 to point to ea_j
                ptr3=dp + j*cnp; // set ptr3 to point to da_j

                //nsolved+=sba_Axb_LU(ptr1, ptr2, ptr3, cnp, 0); linsolver=sba_Axb_LU;
                nsolved+=sba_Axb_Chol(ptr1, ptr2, ptr3, cnp, 0); linsolver=sba_Axb_Chol;
                //nsolved+=sba_Axb_BK(ptr1, ptr2, ptr3, cnp, 0); linsolver=sba_Axb_BK;
                //nsolved+=sba_Axb_QRnoQ(ptr1, ptr2, ptr3, cnp, 0); linsolver=sba_Axb_QRnoQ;
                //nsolved+=sba_Axb_QR(ptr1, ptr2, ptr3, cnp, 0); linsolver=sba_Axb_QR;
                //nsolved+=sba_Axb_SVD(ptr1, ptr2, ptr3, cnp, 0); linsolver=sba_Axb_SVD;
                //nsolved+=(sba_Axb_CG(ptr1, ptr2, ptr3, cnp, (3*cnp)/2, 1E-10, SBA_CG_JACOBI, 0) > 0); linsolver=(PLS)sba_Axb_CG;

                ++nlss;
            }

            if(nsolved==m){

                /* parameter vector updates are now in dp */

                /* compute p's new estimate and ||dp||^2 */
                for(i=0, dp_L2=0.0; i<nvars; ++i){
                    pdp[i]=p[i] + (tmp=dp[i]);
                    dp_L2+=tmp*tmp;
                }
                //dp_L2=sqrt(dp_L2);

                if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
                    //if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
                    stop=2;
                    break;
                }

                if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ){ /* almost singular */
                    //if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
                    fprintf(stderr, "SBA: the matrix of the augmented normal equations is almost singular in sba_mot_levmar_x(),\n"
                            "     minimization should be restarted from the current solution with an increased damping term\n");
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }

                (*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
                if(verbose>1) {
                    printf("mean reprojection error %g\n", sba_mean_repr_error(n, mnp, x, hx, &idxij, rcidxs, rcsubs));
                    fflush(stdout);
                }

                /* ### compute ||e(pdp)||_2 */
                if(covx==NULL)
                    pdp_eL2=nrmL2xmy(hx, x, hx, nobs); /* hx=x-hx, pdp_eL2=||hx|| */
                else
                    pdp_eL2=nrmCxmy(hx, x, hx, wght, mnp, nvis); /* hx=wght*(x-hx), pdp_eL2=||hx|| */
                if(!SBA_FINITE(pdp_eL2)){
                    if(verbose) /* identify the offending point projection */
                        sba_print_inf(hx, m, mnp, &idxij, rcidxs, rcsubs);

                    stop=7;
                    break;
                }

                /* Add in the camera constraints */
                if (use_constraints) {
                    for (j = 0; j < m; j++) {
                        int nnz = 0;

                        for(i = 0; i < n; i++)
                            nnz += vmask[i * m + j];

                        for (jj = 0; jj < cnp; jj++) {
                            if (constraints[j].constrained[jj]) {
                                double diff = 
                                    constraints[j].constraints[jj] - p[j * cnp + jj];
			
                                pdp_eL2 += 
                                    /*nnz **/ constraints[j].weights[jj] * diff * diff;
                            }
                        }
                    }
                }

                for(i=0, dL=0.0; i<nvars; ++i)
                    dL+=dp[i]*(mu*dp[i]+ea[i]);

                dF=p_eL2-pdp_eL2;

                if(verbose>1) {
                    printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);
                    fflush(stdout);
                }

                if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
                    tmp=(2.0*dF/dL-1.0);
                    tmp=1.0-tmp*tmp*tmp;
                    mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
                    nu=2;

                    /* the test below is equivalent to the relative reduction of the RMS reprojection error: sqrt(p_eL2)-sqrt(pdp_eL2)<eps4*sqrt(p_eL2) */
                    if(pdp_eL2-2.0*sqrt(p_eL2*pdp_eL2)<(eps4_sq-1.0)*p_eL2) stop=4;
          
                    for(i=0; i<nvars; ++i) /* update p's estimate */
                        p[i]=pdp[i];

                    for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
                        e[i]=hx[i];
                    p_eL2=pdp_eL2;
                    break;
                }
            } /* nsolved==m */

            /* if this point is reached, either at least one linear system could not be solved or
             * the error did not reduce; in any case, the increment must be rejected
             */

            mu*=nu;
            nu2=nu<<1; // 2*nu;
            if(nu2<=nu){ /* nu has wrapped around (overflown) */
                fprintf(stderr, "SBA: too many failed attempts to increase the damping factor in sba_mot_levmar_x()! Singular Hessian matrix?\n");
                //exit(1);
                stop=6;
                break;
            }
            nu=nu2;

            /* restore U diagonal entries */
            for(j=mcon; j<m; ++j){
                ptr1=U + j*Usz; // set ptr1 to point to U_j
                ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
                for(i=0; i<cnp; ++i)
                    ptr1[i*cnp+i]=ptr2[i];
            }
        } /* inner while loop */

        if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
    }

    if(itno>=itmax) stop=3;

    /* restore U diagonal entries */
    for(j=mcon; j<m; ++j){
        ptr1=U + j*Usz; // set ptr1 to point to U_j
        ptr2=diagU + j*cnp; // set ptr2 to point to diagU_j
        for(i=0; i<cnp; ++i)
            ptr1[i*cnp+i]=ptr2[i];
    }

    if(info){
        info[0]=init_p_eL2;
        info[1]=p_eL2;
        info[2]=ea_inf;
        info[3]=dp_L2;
        for(j=mcon, tmp=DBL_MIN; j<m; ++j){
            ptr1=U + j*Usz; // set ptr1 to point to U_j
            for(i=0; i<cnp; ++i)
                if(tmp<ptr1[i*cnp+i]) tmp=ptr1[i*cnp+i];
        }
        info[4]=mu/tmp;
        info[5]=itno;
        info[6]=stop;
        info[7]=nfev;
        info[8]=njev;
        info[9]=nlss;
    }
    //sba_print_sol(n, m, p, cnp, 0, x, mnp, &idxij, rcidxs, rcsubs);
    retval=(stop!=7)?  itno : SBA_ERROR;
                                                               
 freemem_and_return: /* NOTE: this point is also reached via a goto! */

    /* free whatever was allocated */
    free(jac); free(U);
    free(e); free(ea);  
    free(dp);
    free(rcidxs); free(rcsubs);
#ifndef SBA_DESTROY_COVS
    if(wght) free(wght);
#else
    /* nothing to do */
#endif /* SBA_DESTROY_COVS */

    free(hx); free(diagU); free(pdp);
    if(fdj_data.hxx){ // cleanup
        free(fdj_data.hxx);
        free(fdj_data.func_rcidxs);
    }

    sba_crsm_free(&idxij);

    /* free the memory allocated by the linear solver routine */
    if(linsolver) (*linsolver)(NULL, NULL, NULL, 0, 0);

    return retval;
}


/* Bundle adjustment on structure parameters only 
 * using the sparse Levenberg-Marquardt as described in HZ p. 568
 *
 * Returns the number of iterations (>=0) if successfull, SBA_ERROR if failed
 */

int sba_str_levmar_x(
                     const int n,   /* number of points */
                     const int m,   /* number of images */
                     char *vmask,  /* visibility mask: vmask[i, j]=1 if point i visible in image j, 0 otherwise. nxm */
                     double *p,    /* initial parameter vector p0: (b1, ..., bn).
                                    * bi are the i-th point parameters, * size n*pnp */
                     const int pnp,/* number of parameters for ONE point; e.g. 3 for Euclidean points */
                     double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
                                    * x_ij is the projection of the i-th point on the j-th image.
                                    * NOTE: some of the x_ij might be missing, if point i is not visible in image j;
                                    * see vmask[i, j], max. size n*m*mnp
                                    */
                     double *covx, /* measurements covariance matrices: (Sigma_x_11, .. Sigma_x_1m, ..., Sigma_x_n1, .. Sigma_x_nm),
                                    * where Sigma_x_ij is the mnp x mnp covariance of x_ij stored row-by-row. Set to NULL if no
                                    * covariance estimates are available (identity matrices are implicitly used in this case).
                                    * NOTE: a certain Sigma_x_ij is missing if the corresponding x_ij is also missing;
                                    * see vmask[i, j], max. size n*m*mnp*mnp
                                    */
                     const int mnp,/* number of parameters for EACH measurement; usually 2 */
                     void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
                     /* functional relation describing measurements. Given a parameter vector p,
                      * computes a prediction of the measurements \hat{x}. p is (n*pnp)x1,
                      * \hat{x} is (n*m*mnp)x1, maximum
                      * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                      * as working memory
                      */
                     void (*fjac)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
                     /* function to evaluate the sparse jacobian dX/dp.
                      * The Jacobian is returned in jac as
                      * (dx_11/db_1, ..., dx_1m/db_1, ..., dx_n1/db_n, ..., dx_nm/db_n), or (using HZ's notation),
                      * jac=(B_11, ..., B_1m, ..., B_n1, ..., B_nm)
                      * Notice that depending on idxij, some of the B_ij might be missing.
                      * Note also that B_ij are mnp x pnp matrices and they
                      * should be stored in jac in row-major order.
                      * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
                      * as working memory
                      *
                      * If NULL, the jacobian is approximated by repetitive func calls and finite
                      * differences. This is computationally inefficient and thus NOT recommended.
                      */
                     void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 

                     const int itmax,   /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
                     const int verbose, /* I: verbosity */
                     const double opts[SBA_OPTSSZ],
                     /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \epsilon4]. Respectively the scale factor for initial \mu,
                      * stopping thresholds for ||J^T e||_inf, ||dp||_2, ||e||_2 and (||e||_2-||e_new||_2)/||e||_2
                      */
                     double info[SBA_INFOSZ]
                     /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]=||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small dp
                      *                                 3 - stopped by itmax
                      *                                 4 - stopped by small relative reduction in ||e||_2
                      *                                 5 - stopped by small ||e||_2
                      *                                 6 - too many attempts to increase damping. Restart with increased mu
                      *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
                      * info[7]= # function evaluations
                      * info[8]= # jacobian evaluations
                      * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
                      */
                     )
{
    register int i, j, ii, jj, k;
    int nvis, nnz, retval;

    /* The following are work arrays that are dynamically allocated by sba_str_levmar_x() */
    double *jac;  /* work array for storing the jacobian, max. size n*m*mnp*pnp */
    double *V;    /* work array for storing the V_i in the order V_1, ..., V_n, size n*pnp*pnp */

    double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_n1, ..., e_nm,
                     max. size n*m*mnp */
    double *eb;   /* work array for storing the eb_i in the order eb_1, .. eb_n size n*pnp */

    double *dp;   /* work array for storing the parameter vector updates db_1, ..., db_n, size n*pnp */

    double *wght= /* work array for storing the weights computed from the covariance inverses, max. size n*m*mnp*mnp */
        NULL;

    /* Of the above arrays, jac, e, wght are sparse and
     * V, eb, dp are dense. Sparse arrays are indexed through
     * idxij (see below), that is with the same mechanism as the input 
     * measurements vector x
     */

    /* submatrices sizes */
    int Bsz, Vsz,
        esz, ebsz, covsz;

    register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
    struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also the location
                            * of B_ij in jac, etc.
                            * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                            * index k and it is used to efficiently lookup the memory locations where the non-zero
                            * blocks of a sparse matrix/vector are stored
                            */
    int maxCPvis, /* max. of projections across cameras & projections across points */
        *rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
                     column in a sparse matrix, size max(n, m) */
        *rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
                     sparse matrix, size max(n, m) */

    /* The following variables are needed by the LM algorithm */
    register int itno;  /* iteration counter */
    int nsolved;
    /* temporary work arrays that are dynamically allocated */
    double *hx,         /* \hat{x}_i, max. size m*n*mnp */
        *diagV,      /* diagonals of V_i, size n*pnp */
        *pdp;        /* p + dp, size n*pnp */

    register double mu,  /* damping constant */
        tmp; /* mainly used in matrix & vector multiplications */
    double p_eL2, eb_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
    double p_L2, dp_L2=DBL_MAX, dF, dL;
    double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2],
        eps3_sq=opts[3]*opts[3], eps4_sq=opts[4]*opts[4];
    double init_p_eL2;
    int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;
    int nobs, nvars;
    PLS linsolver=NULL;

    struct fdj_data_x_ fdj_data;
    void *jac_adata;

    /* Initialization */
    mu=eb_inf=tmp=0.0; /* -Wall */

    /* block sizes */
    Bsz=mnp * pnp; Vsz=pnp * pnp;
    esz=mnp; ebsz=pnp;
    covsz=mnp * mnp;

    /* count total number of visible image points */
    for(i=nvis=0, jj=n*m; i<jj; ++i)
        nvis+=(vmask[i]!=0);

    nobs=nvis*mnp;
    nvars=n*pnp;
    if(nobs<nvars){
        fprintf(stderr, "SBA: sba_str_levmar_x() cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
        return SBA_ERROR;
    }

    /* allocate & fill up the idxij structure */
    sba_crsm_alloc(&idxij, n, m, nvis);
    for(i=k=0; i<n; ++i){
        idxij.rowptr[i]=k;
        ii=i*m;
        for(j=0; j<m; ++j)
            if(vmask[ii+j]){
                idxij.val[k]=k;
                idxij.colidx[k++]=j;
            }
    }
    idxij.rowptr[n]=nvis;

    /* find the maximum number of visible image points in any single camera or coming from a single 3D point */
    /* cameras */
    for(i=maxCPvis=0; i<n; ++i)
        if((k=idxij.rowptr[i+1]-idxij.rowptr[i])>maxCPvis) maxCPvis=k;

    /* points, note that maxCPvis is not reinitialized! */
    for(j=0; j<m; ++j){
        for(i=ii=0; i<n; ++i)
            if(vmask[i*m+j]) ++ii;
        if(ii>maxCPvis) maxCPvis=ii;
    }

    /* allocate work arrays */
    jac=(double *)emalloc(nvis*Bsz*sizeof(double));
    V=(double *)emalloc(n*Vsz*sizeof(double));
    e=(double *)emalloc(nobs*sizeof(double));
    eb=(double *)emalloc(nvars*sizeof(double));
    dp=(double *)emalloc(nvars*sizeof(double));
    rcidxs=(int *)emalloc(maxCPvis*sizeof(int));
    rcsubs=(int *)emalloc(maxCPvis*sizeof(int));
#ifndef SBA_DESTROY_COVS
    if(covx!=NULL) wght=(double *)emalloc(nvis*covsz*sizeof(double));
#else
    if(covx!=NULL) wght=covx;
#endif /* SBA_DESTROY_COVS */


    hx=(double *)emalloc(nobs*sizeof(double));
    diagV=(double *)emalloc(nvars*sizeof(double));
    pdp=(double *)emalloc(nvars*sizeof(double));

    /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
    if(!fjac){
        fdj_data.func=func;
        fdj_data.cnp=0;
        fdj_data.pnp=pnp;
        fdj_data.mnp=mnp;
        fdj_data.hx=hx;
        fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
        fdj_data.func_rcidxs=(int *)emalloc(2*maxCPvis*sizeof(int));
        fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxCPvis;
        fdj_data.adata=adata;

        fjac=sba_fdjac_x;
        jac_adata=(void *)&fdj_data;
    }
    else{
        fdj_data.hxx=NULL;
        jac_adata=adata;
    }

    if(itmax==0){ /* verify jacobian */
        sba_str_chkjac_x(func, fjac, p, &idxij, rcidxs, rcsubs, pnp, mnp, adata, jac_adata);
        retval=0;
        goto freemem_and_return;
    }

    /* covariances Sigma_x_ij are accommodated by computing the Cholesky decompositions of their
     * inverses and using the resulting matrices w_x_ij to weigh B_ij and e_ij as
     * w_x_ij*B_ij and w_x_ij*e_ij. In this way, auxiliary variables as V_i=\sum_j B_ij^T B_ij
     * actually become \sum_j (w_x_ij B_ij)^T w_x_ij B_ij= \sum_j B_ij^T w_x_ij^T w_x_ij B_ij =
     * B_ij^T Sigma_x_ij^-1 B_ij
     *
     * eb_i are weighted in a similar manner
     */
    if(covx!=NULL){
        for(i=0; i<n; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1, ptr2 to point to cov_x_ij, w_x_ij resp. */
                ptr1=covx + (k=idxij.val[rcidxs[j]]*covsz);
                ptr2=wght + k;
                if(!sba_mat_cholinv(ptr1, ptr2, mnp)){ /* compute w_x_ij s.t. w_x_ij^T w_x_ij = cov_x_ij^-1 */
                    fprintf(stderr, "SBA: invalid covariance matrix for x_ij (i=%d, j=%d) in sba_motstr_levmar_x()\n", i, rcsubs[j]);
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }
            }
        }
        sba_mat_cholinv(NULL, NULL, 0); /* cleanup */
    }

    /* compute the error vectors e_ij in hx */
    (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
    /* ### compute e=x - f(p) [e=w*(x - f(p)] and its L2 norm */
    if(covx==NULL)
        p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
    else
        p_eL2=nrmCxmy(e, x, hx, wght, mnp, nvis); /* e=wght*(x-hx), p_eL2=||e||=||x-hx||_Sigma^-1 */
    if(verbose) printf("initial str-SBA error %g [%g]\n", p_eL2, p_eL2/nvis);
    init_p_eL2=p_eL2;
    if(!SBA_FINITE(p_eL2)) stop=7;

    for(itno=0; itno<itmax && !stop; ++itno){
        /* Note that p, e and ||e||_2 have been updated at the previous iteration */

        /* compute derivative submatrices B_ij */
        (*fjac)(p, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;

        if(covx!=NULL){
            /* compute w_x_ij B_ij.
             * Since w_x_ij is upper triangular, the products can be safely saved
             * directly in B_ij, without the need for intermediate storage
             */
            for(i=0; i<nvis; ++i){
                /* set ptr1, ptr2 to point to w_x_ij, B_ij, resp. */
                ptr1=wght + i*covsz;
                ptr2=jac  + i*Bsz;

                /* w_x_ij is mnp x mnp, B_ij is mnp x pnp */
                for(ii=0; ii<mnp; ++ii)
                    for(jj=0; jj<pnp; ++jj){
                        for(k=ii, sum=0.0; k<mnp; ++k) // k>=ii since w_x_ij is upper triangular
                            sum+=ptr1[ii*mnp+k]*ptr2[k*pnp+jj];
                        ptr2[ii*pnp+jj]=sum;
                    }
            }
        }

        /* compute V_i = \sum_j B_ij^T B_ij */ // \Sigma here!
        /* V_i is symmetric, therefore its computation can be sped up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that B_ij is mnp x pnp
         */
        /* Also compute eb_i = \sum_j B_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */
        _dblzero(V, n*Vsz); /* clear all V_i */
        _dblzero(eb, n*ebsz); /* clear all eb_i */
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=eb + i*ebsz; // set ptr2 to point to eb_i

            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr3 to point to B_ij, actual column number in rcsubs[j] */
                ptr3=jac + idxij.val[rcidxs[j]]*Bsz;
      
                /* compute the UPPER TRIANGULAR PART of B_ij^T B_ij and add it to V_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=ii; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*pnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]+=sum;
                    }

                    /* copy the LOWER TRIANGULAR PART of V_i from the upper one */
                    for(jj=0; jj<ii; ++jj)
                        ptr1[ii*pnp+jj]=ptr1[jj*pnp+ii];
                }

                ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
                /* compute B_ij^T e_ij and add it to eb_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*pnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }
        }

        /* Compute ||J^T e||_inf and ||p||^2 */
        for(i=0, p_L2=eb_inf=0.0; i<nvars; ++i){
            if(eb_inf < (tmp=FABS(eb[i]))) eb_inf=tmp;
            p_L2+=p[i]*p[i];
        }
        //p_L2=sqrt(p_L2);

        /* save diagonal entries so that augmentation can be later canceled.
         * Diagonal entries are in V_i
         */
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
            for(j=0; j<pnp; ++j)
                ptr2[j]=ptr1[j*pnp+j];
        }

        /*
          if(!(itno%100)){
          printf("Current estimate: ");
          for(i=0; i<nvars; ++i)
          printf("%.9g ", p[i]);
          printf("-- errors %.9g %0.9g\n", eb_inf, p_eL2);
          }
        */

        /* check for convergence */
        if((eb_inf <= eps1)){
            dp_L2=0.0; /* no increment for p in this case */
            stop=1;
            break;
        }

        /* compute initial damping factor */
        if(itno==0){
            for(i=0, tmp=DBL_MIN; i<nvars; ++i)
                if(diagV[i]>tmp) tmp=diagV[i]; /* find max diagonal element */
            mu=tau*tmp;
        }

        /* determine increment using adaptive damping */
        while(1){
            /* augment V */
            for(i=0; i<n; ++i){
                ptr1=V + i*Vsz; // set ptr1 to point to V_i
                for(j=0; j<pnp; ++j)
                    ptr1[j*pnp+j]+=mu;
            }

            /* solve the linear systems V*_i db_i = eb_i to compute the db_i */
            for(i=nsolved=0; i<n; ++i){
                ptr1=V + i*Vsz; // set ptr1 to point to V_i
                ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
                ptr3=dp + i*pnp; // set ptr3 to point to db_i

                //nsolved+=sba_Axb_LU(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_LU;
                nsolved+=sba_Axb_Chol(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_Chol;
                //nsolved+=sba_Axb_BK(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_BK;
                //nsolved+=sba_Axb_QRnoQ(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_QRnoQ;
                //nsolved+=sba_Axb_QR(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_QR;
                //nsolved+=sba_Axb_SVD(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_SVD;
                //nsolved+=(sba_Axb_CG(ptr1, ptr2, ptr3, pnp, (3*pnp)/2, 1E-10, SBA_CG_JACOBI, 0) > 0); linsolver=(PLS)sba_Axb_CG;

                ++nlss;
            }

            if(nsolved==n){

                /* parameter vector updates are now in dp */

                /* compute p's new estimate and ||dp||^2 */
                for(i=0, dp_L2=0.0; i<nvars; ++i){
                    pdp[i]=p[i] + (tmp=dp[i]);
                    dp_L2+=tmp*tmp;
                }
                //dp_L2=sqrt(dp_L2);

                if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
                    //if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
                    stop=2;
                    break;
                }

                if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ){ /* almost singular */
                    //if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
                    fprintf(stderr, "SBA: the matrix of the augmented normal equations is almost singular in sba_motstr_levmar_x(),\n"
                            "     minimization should be restarted from the current solution with an increased damping term\n");
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }

                (*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
                if(verbose>1)
                    printf("mean reprojection error %g\n", sba_mean_repr_error(n, mnp, x, hx, &idxij, rcidxs, rcsubs));
                /* ### compute ||e(pdp)||_2 */
                if(covx==NULL)
                    pdp_eL2=nrmL2xmy(hx, x, hx, nobs); /* hx=x-hx, pdp_eL2=||hx|| */
                else
                    pdp_eL2=nrmCxmy(hx, x, hx, wght, mnp, nvis); /* hx=wght*(x-hx), pdp_eL2=||hx|| */
                if(!SBA_FINITE(pdp_eL2)){
                    if(verbose) /* identify the offending point projection */
                        sba_print_inf(hx, m, mnp, &idxij, rcidxs, rcsubs);

                    stop=7;
                    break;
                }

                for(i=0, dL=0.0; i<nvars; ++i)
                    dL+=dp[i]*(mu*dp[i]+eb[i]);

                dF=p_eL2-pdp_eL2;

                if(verbose>1)
                    printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);

                if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
                    tmp=(2.0*dF/dL-1.0);
                    tmp=1.0-tmp*tmp*tmp;
                    mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
                    nu=2;

                    /* the test below is equivalent to the relative reduction of the RMS reprojection error: sqrt(p_eL2)-sqrt(pdp_eL2)<eps4*sqrt(p_eL2) */
                    if(pdp_eL2-2.0*sqrt(p_eL2*pdp_eL2)<(eps4_sq-1.0)*p_eL2) stop=4;
          
                    for(i=0; i<nvars; ++i) /* update p's estimate */
                        p[i]=pdp[i];

                    for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
                        e[i]=hx[i];
                    p_eL2=pdp_eL2;
                    break;
                }
            } /* nsolved==n */

            /* if this point is reached, either at least one linear system could not be solved or
             * the error did not reduce; in any case, the increment must be rejected
             */

            mu*=nu;
            nu2=nu<<1; // 2*nu;
            if(nu2<=nu){ /* nu has wrapped around (overflown) */
                fprintf(stderr, "SBA: too many failed attempts to increase the damping factor in sba_str_levmar_x()! Singular Hessian matrix?\n");
                //exit(1);
                stop=6;
                break;
            }
            nu=nu2;

            /* restore V diagonal entries */
            for(i=0; i<n; ++i){
                ptr1=V + i*Vsz; // set ptr1 to point to V_i
                ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
                for(j=0; j<pnp; ++j)
                    ptr1[j*pnp+j]=ptr2[j];
            }
        } /* inner while loop */

        if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
    }

    if(itno>=itmax) stop=3;

    /* restore V diagonal entries */
    for(i=0; i<n; ++i){
        ptr1=V + i*Vsz; // set ptr1 to point to V_i
        ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
        for(j=0; j<pnp; ++j)
            ptr1[j*pnp+j]=ptr2[j];
    }

    if(info){
        info[0]=init_p_eL2;
        info[1]=p_eL2;
        info[2]=eb_inf;
        info[3]=dp_L2;
        for(i=0; i<n; ++i){
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            for(j=0; j<pnp; ++j)
                if(tmp<ptr1[j*pnp+j]) tmp=ptr1[j*pnp+j];
        }
        info[4]=mu/tmp;
        info[5]=itno;
        info[6]=stop;
        info[7]=nfev;
        info[8]=njev;
        info[9]=nlss;
    }
    //sba_print_sol(n, m, p, 0, pnp, x, mnp, &idxij, rcidxs, rcsubs);
    retval=(stop!=7)?  itno : SBA_ERROR;
                                                               
 freemem_and_return: /* NOTE: this point is also reached via a goto! */

    /* free whatever was allocated */
    free(jac); free(V);
    free(e); free(eb);  
    free(dp);               
    free(rcidxs); free(rcsubs);
#ifndef SBA_DESTROY_COVS
    if(wght) free(wght);
#else
    /* nothing to do */
#endif /* SBA_DESTROY_COVS */

    free(hx); free(diagV); free(pdp);
    if(fdj_data.hxx){ // cleanup
        free(fdj_data.hxx);
        free(fdj_data.func_rcidxs);
    }

    sba_crsm_free(&idxij);

    /* free the memory allocated by the linear solver routine */
    if(linsolver) (*linsolver)(NULL, NULL, NULL, 0, 0);

    return retval;
}

int sba_str_levmar_grid_x(
						  const int n,   /* number of points */
						  const int noGrids,		/* number of grids*/
						  const int* grid_rows, /*number of rows in each grid*/
						  const int* grid_cols, /*number of cols in each grid*/
						  const int m,   /* number of images */
						  char *vmask,  /* visibility mask: vmask[i, j]=1 if point i visible in image j, 0 otherwise. nxm */
						  double *p,    /* initial parameter vector p0: (b1, ..., bn).
										 * bi are the i-th point parameters, * size n*pnp */
						  const int pnp,/* number of parameters for ONE point; e.g. 3 for Euclidean points */
						  double *x,    /* measurements vector: (x_11^T, .. x_1m^T, ..., x_n1^T, .. x_nm^T)^T where
										 * x_ij is the projection of the i-th point on the j-th image.
										 * NOTE: some of the x_ij might be missing, if point i is not visible in image j;
										 * see vmask[i, j], max. size n*m*mnp
										 */
						  double *covx, /* measurements covariance matrices: (Sigma_x_11, .. Sigma_x_1m, ..., Sigma_x_n1, .. Sigma_x_nm),
										 * where Sigma_x_ij is the mnp x mnp covariance of x_ij stored row-by-row. Set to NULL if no
										 * covariance estimates are available (identity matrices are implicitly used in this case).
										 * NOTE: a certain Sigma_x_ij is missing if the corresponding x_ij is also missing;
										 * see vmask[i, j], max. size n*m*mnp*mnp
										 */
						  const int mnp,/* number of parameters for EACH measurement; usually 2 */
						  void (*func)(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata),
						  /* functional relation describing measurements. Given a parameter vector p,
						   * computes a prediction of the measurements \hat{x}. p is (n*pnp)x1,
						   * \hat{x} is (n*m*mnp)x1, maximum
						   * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
						   * as working memory
						   */
						  void (*fjac)(double *p, int n, int nvis, int n_grid, int *grid_rows, int *grid_cols, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata),
						  /* function to evaluate the sparse jacobian dX/dp.
						   * The Jacobian is returned in jac as
						   * (dx_11/db_1, ..., dx_1m/db_1, ..., dx_n1/db_n, ..., dx_nm/db_n), or (using HZ's notation),
						   * jac=(B_11, ..., B_1m, ..., B_n1, ..., B_nm)
						   * Notice that depending on idxij, some of the B_ij might be missing.
						   * Note also that B_ij are mnp x pnp matrices and they
						   * should be stored in jac in row-major order.
						   * rcidxs, rcsubs are max(m, n) x 1, allocated by the caller and can be used
						   * as working memory
						   *
						   * If NULL, the jacobian is approximated by repetitive func calls and finite
						   * differences. This is computationally inefficient and thus NOT recommended.
						   */
						  void *adata,       /* pointer to possibly additional data, passed uninterpreted to func, fjac */ 
						  
						  const int itmax,   /* I: maximum number of iterations. itmax==0 signals jacobian verification followed by immediate return */
						  const int verbose, /* I: verbosity */
						  const double opts[SBA_OPTSSZ],
						  /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \epsilon4]. Respectively the scale factor for initial \mu,
						   * stopping thresholds for ||J^T e||_inf, ||dp||_2, ||e||_2 and (||e||_2-||e_new||_2)/||e||_2
						   */
						  double info[SBA_INFOSZ]
						  /* O: information regarding the minimization. Set to NULL if don't care
						   * info[0]=||e||_2 at initial p.
						   * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
						   * info[5]= # iterations,
						   * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
						   *                                 2 - stopped by small dp
						   *                                 3 - stopped by itmax
						   *                                 4 - stopped by small relative reduction in ||e||_2
						   *                                 5 - stopped by small ||e||_2
						   *                                 6 - too many attempts to increase damping. Restart with increased mu
						   *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
						   * info[7]= # function evaluations
						   * info[8]= # jacobian evaluations
						   * info[9]= # number of linear systems solved, i.e. number of attempts	for reducing error
						   */
						  )
{
    register int i, j, ii, jj, k;
    int nvis, nnz, retval;
	
    /* The following are work arrays that are dynamically allocated by sba_str_levmar_x() */
    double *jac;  /* work array for storing the jacobian, max. size n*m*mnp*pnp */
    double *V;    /* work array for storing the V_i in the order V_1, ..., V_n, size n*pnp*pnp */
	double *N;	  /*size 3*pnp*3*pnp*noGrids*/
	double *e;    /* work array for storing the e_ij in the order e_11, ..., e_1m, ..., e_n1, ..., e_nm,
				   max. size n*m*mnp */
	double *ebtrans;   /* work array for storing the eb_i in the order eb_1, .. eb_n size n*pnp + 3*pnp*noGrids*/
	
    double *dp;   /* work array for storing the parameter vector updates db_1, ..., db_n, size n*pnp */
	
    double *wght= /* work array for storing the weights computed from the covariance inverses, max. size n*m*mnp*mnp */
	NULL;
	
    /* Of the above arrays, jac, e, wght are sparse and
     * V, eb, dp are dense. Sparse arrays are indexed through
     * idxij (see below), that is with the same mechanism as the input 
     * measurements vector x
     */
	
    /* submatrices sizes */
    int Bsz, Vsz,
	Ptsz, THorsz, TVersz, Nsz,
	esz, ebsz, etranssz, covsz;
	
    register double *ptr1, *ptr2, *ptr3, *ptr4, sum;
    struct sba_crsm idxij; /* sparse matrix containing the location of x_ij in x. This is also the location
                            * of B_ij in jac, etc.
                            * This matrix can be thought as a map from a sparse set of pairs (i, j) to a continuous
                            * index k and it is used to efficiently lookup the memory locations where the non-zero
                            * blocks of a sparse matrix/vector are stored
                            */
    int maxCvis, maxPvis, maxCPvis, /* max. of projections across cameras & projections across points */
	*rcidxs,  /* work array for the indexes corresponding to the nonzero elements of a single row or
			   column in a sparse matrix, size max(n, m) */
	*rcsubs;  /* work array for the subscripts of nonzero elements in a single row or column of a
			   sparse matrix, size max(n, m) */
	
    /* The following variables are needed by the LM algorithm */
    register int itno;  /* iteration counter */
    int nsolved;
    /* temporary work arrays that are dynamically allocated */
    double *hx,         /* \hat{x}_i, max. size m*n*mnp */
	*diagVN,      /* diagonals of V_i, size n*pnp */
	*pdp;        /* p + dp, size n*pnp */
	
	double *diagV, *diagN; /* pointers into diagUVN */
	
    register double mu,  /* damping constant */
	tmp; /* mainly used in matrix & vector multiplications */
    double p_eL2, ebtrans_inf, pdp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+dp)||_2 */
    double p_L2, dp_L2=DBL_MAX, dF, dL;
    double tau=FABS(opts[0]), eps1=FABS(opts[1]), eps2=FABS(opts[2]), eps2_sq=opts[2]*opts[2],
	eps3_sq=opts[3]*opts[3], eps4_sq=opts[4]*opts[4];
    double init_p_eL2;
    int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;
    int nobs, nvars;
    PLS linsolver=NULL;
	double *ptr1_psum, *ptr1_THorP, *ptr1_TVerP, *ptr1_PTHor, *ptr1_hsum, *ptr1_TVerTHor, *ptr1_PTVer, *ptr1_THorTVer;
	double *ptr1_vsum, *ptr2_0, *ptr2_1, *ptr2_2, *ptr3_pt, *ptr3_thor, *ptr3_tver;
	double *eb, *et;
	double *dpt;
	
    struct fdj_data_x_ fdj_data;
    void *jac_adata;
	
    /* Initialization */
    mu=ebtrans_inf=tmp=0.0; /* -Wall */
	
    /* block sizes */
    Bsz=mnp * pnp; Vsz=pnp * pnp;
	Ptsz=mnp * pnp; THorsz=mnp * pnp; TVersz=mnp*pnp; Nsz = pnp*pnp;
    esz=mnp; ebsz=pnp; etranssz = pnp;
    covsz=mnp * mnp;
	
    /* count total number of visible image points */
    int totalNoPoints = n;
	for(i=0; i<noGrids; i++)
	{
		totalNoPoints += grid_cols[i]*grid_rows[i];
	}
	
	printf("totalPts:%d\n", totalNoPoints);
	
    /* count total number of visible image points */
    int totalNVis = 0;
	nvis = 0;
	int* nVisGrids = (int*) malloc(sizeof(int)*noGrids);
	for(i=0; i<noGrids; i++)
	{
		nVisGrids[i] = 0;
	}
	
	for(i=nvis=0, jj=totalNoPoints*m; i<jj; ++i)
	{
		if(vmask[i]!=0)
			printf("proj:%f %f\n", x[totalNVis*mnp], x[totalNVis*mnp+1]);
		
		totalNVis += (vmask[i]!=0);
		
        if(i<n*m)
		{
			nvis += (vmask[i]!=0);
			continue;
		}
		int prevPoints = n*m;
		for(j=0; j<noGrids; j++)
		{
			prevPoints += grid_cols[j]*grid_rows[j]*m;
			if(i<prevPoints)
			{
				nVisGrids[j] += (vmask[i]!=0);
				break;
			}
		}
	}
	printf("totalNVis:%d\n", totalNVis);
	
    nobs=totalNVis*mnp;
    nvars=n*pnp + noGrids*pnp*3;
    if(nobs<nvars){
        fprintf(stderr, "SBA: sba_str_levmar_x() cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n", nobs, nvars);
        return SBA_ERROR;
    }
	
	printf("allocare idxij\n");
    /* allocate & fill up the idxij structure */
    sba_crsm_alloc(&idxij, totalNoPoints, m, totalNVis);
    for(i=k=0; i<totalNoPoints; ++i){
        idxij.rowptr[i]=k;
        ii=i*m;
        for(j=0; j<m; ++j)
            if(vmask[ii+j]){
				idxij.val[k]=k;
                idxij.colidx[k++]=j;
            }
    }
    idxij.rowptr[totalNoPoints]=totalNVis;
	
    /* find the maximum number of visible image points in any single camera or coming from a single 3D point */
    /* cameras */
	printf("find maxCPvis\n");
	/* find the maximum number (for all cameras) of visible image projections coming from a single 3D point */
    for(i=maxCvis=0; i<totalNoPoints; ++i)
        if((k=idxij.rowptr[i+1]-idxij.rowptr[i])>maxCvis) maxCvis=k;
	
    /* find the maximum number (for all points) of visible image projections in any single camera */
    for(j=maxPvis=0; j<m; ++j){
        for(i=ii=0; i<totalNoPoints; ++i)
            if(vmask[i*m+j]) ++ii;
        if(ii>maxPvis) maxPvis=ii;
    }
    maxCPvis=(maxCvis>=maxPvis)? maxCvis : maxPvis;

    /* allocate work arrays */
	printf("allocare arrays\n");
    jac=(double *)emalloc((nvis*Bsz+(totalNVis-nvis)*(Ptsz + THorsz + TVersz))*sizeof(double));
    V=(double *)emalloc(n*Vsz*sizeof(double));
	N = (double *)emalloc(noGrids*3*3*Nsz*sizeof(double));
    e=(double *)emalloc(nobs*sizeof(double));
    ebtrans=(double *)emalloc(nvars*sizeof(double));
    dp=(double *)emalloc(nvars*sizeof(double));
    rcidxs=(int *)emalloc(maxCPvis*sizeof(int));
    rcsubs=(int *)emalloc(maxCPvis*sizeof(int));
	eb = ebtrans; et = ebtrans + n*pnp;
	dpt = dp + n*pnp;
	
	int gridIt;
	
#ifndef SBA_DESTROY_COVS
    if(covx!=NULL) wght=(double *)emalloc(nvis*covsz*sizeof(double));
#else
    if(covx!=NULL) wght=covx;
#endif /* SBA_DESTROY_COVS */
	
	
    hx=(double *)emalloc(nobs*sizeof(double));
    diagVN=(double *)emalloc(nvars*sizeof(double));
    pdp=(double *)emalloc(nvars*sizeof(double));
	
    /* if no jacobian function is supplied, prepare to compute jacobian with finite difference */
    if(!fjac){
        fdj_data.func=func;
        fdj_data.cnp=0;
        fdj_data.pnp=pnp;
        fdj_data.mnp=mnp;
        fdj_data.hx=hx;
        fdj_data.hxx=(double *)emalloc(nobs*sizeof(double));
        fdj_data.func_rcidxs=(int *)emalloc(2*maxCPvis*sizeof(int));
        fdj_data.func_rcsubs=fdj_data.func_rcidxs+maxCPvis;
        fdj_data.adata=adata;
		
        fjac=sba_fdjac_str_grid_x;
        jac_adata=(void *)&fdj_data;
    }
    else{
        fdj_data.hxx=NULL;
        jac_adata=adata;
    }
	
	printf("verify jacobian\n");
    if(itmax==0){ /* verify jacobian */
        sba_str_chkjac_x_grid(func, fjac, p, &idxij, rcidxs, rcsubs, n, nvis, noGrids, grid_rows, grid_cols, pnp, mnp, adata, jac_adata);
        retval=0;
        goto freemem_and_return;
    }
	
    /* covariances Sigma_x_ij are accommodated by computing the Cholesky decompositions of their
     * inverses and using the resulting matrices w_x_ij to weigh B_ij and e_ij as
     * w_x_ij*B_ij and w_x_ij*e_ij. In this way, auxiliary variables as V_i=\sum_j B_ij^T B_ij
     * actually become \sum_j (w_x_ij B_ij)^T w_x_ij B_ij= \sum_j B_ij^T w_x_ij^T w_x_ij B_ij =
     * B_ij^T Sigma_x_ij^-1 B_ij
     *
     * eb_i are weighted in a similar manner
     */
    if(covx!=NULL){
        for(i=0; i<totalNoPoints; ++i){
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero x_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr1, ptr2 to point to cov_x_ij, w_x_ij resp. */
                ptr1=covx + (k=idxij.val[rcidxs[j]]*covsz);
                ptr2=wght + k;
                if(!sba_mat_cholinv(ptr1, ptr2, mnp)){ /* compute w_x_ij s.t. w_x_ij^T w_x_ij = cov_x_ij^-1 */
                    fprintf(stderr, "SBA: invalid covariance matrix for x_ij (i=%d, j=%d) in sba_motstr_levmar_x()\n", i, rcsubs[j]);
                    retval=SBA_ERROR;
                    goto freemem_and_return;
                }
            }
        }
        sba_mat_cholinv(NULL, NULL, 0); /* cleanup */
    }
	
    /* compute the error vectors e_ij in hx */
    (*func)(p, &idxij, rcidxs, rcsubs, hx, adata); nfev=1;
	for(i=0; i<nobs; i++)
		printf("e:%f %f\n", e[i*2], e[i*2+1]);
	
    /* ### compute e=x - f(p) [e=w*(x - f(p)] and its L2 norm */
    if(covx==NULL)
        p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
    else
        p_eL2=nrmCxmy(e, x, hx, wght, mnp, totalNVis); /* e=wght*(x-hx), p_eL2=||e||=||x-hx||_Sigma^-1 */
    
	for(i=0; i<nobs; i++)
		printf("e:%f %f\n", e[i*2], e[i*2+1]);
	
	if(verbose) printf("initial str-SBA error %g [%g]\n", p_eL2, p_eL2/totalNVis);
    init_p_eL2=p_eL2;
    if(!SBA_FINITE(p_eL2)) stop=7;
	
    for(itno=0; itno<itmax && !stop; ++itno){
        /* Note that p, e and ||e||_2 have been updated at the previous iteration */
		
        /* compute derivative submatrices */
        (*fjac)(p, n, nvis, noGrids, grid_rows, grid_cols, &idxij, rcidxs, rcsubs, jac, jac_adata); ++njev;
		printf("derivative computed\n");
		
        /* compute V_i = \sum_j B_ij^T B_ij */ // \Sigma here!
        /* V_i is symmetric, therefore its computation can be sped up by
         * computing only the upper part and then reusing it for the lower one.
         * Recall that B_ij is mnp x pnp
         */
        /* Also compute eb_i = \sum_j B_ij^T e_ij */ // \Sigma here!
        /* Recall that e_ij is mnp x 1
         */
        _dblzero(V, n*Vsz); /* clear all V_i */
        _dblzero(eb, n*ebsz); /* clear all eb_i */
        for(i=0; i<n; ++i){
			printf("V for %d\n", i);
			
            ptr1=V + i*Vsz; // set ptr1 to point to V_i
            ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
			
            nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero B_ij, j=0...m-1 */
            for(j=0; j<nnz; ++j){
                /* set ptr3 to point to B_ij, actual column number in rcsubs[j] */
                ptr3=jac + idxij.val[rcidxs[j]]*Bsz;
				
                /* compute the UPPER TRIANGULAR PART of B_ij^T B_ij and add it to V_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=ii; jj<pnp; ++jj){
                        for(k=0, sum=0.0; k<mnp; ++k)
                            sum+=ptr3[k*pnp+ii]*ptr3[k*pnp+jj];
                        ptr1[ii*pnp+jj]+=sum;
                    }
					
                    /* copy the LOWER TRIANGULAR PART of V_i from the upper one */
                    for(jj=0; jj<ii; ++jj)
                        ptr1[ii*pnp+jj]=ptr1[jj*pnp+ii];
                }
				
                ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
                /* compute B_ij^T e_ij and add it to eb_i */
                for(ii=0; ii<pnp; ++ii){
                    for(jj=0, sum=0.0; jj<mnp; ++jj)
                        sum+=ptr3[jj*pnp+ii]*ptr4[jj];
                    ptr2[ii]+=sum;
                }
            }
        }
		
		/* compute N = [\sum_i \sum_j Pt_ij^T Pt_ij  \sum_i \sum_j Pt_ij^T THor_ij \sum_i \sum_j Pt_ij^T TVer_ij 
		 \sum_i \sum_j THor_ij^T Pt_ij  \sum_i \sum_j THor_ij^T THor_ij \sum_i \sum_j THor_ij^T TVer_ij
		 \sum_i \sum_j TVer_ij^T Pt_ij  \sum_i \sum_j TVer_ij^T THor_ij \sum_i \sum_j TVer_ij^T TVer_ij]
		 * Recall that Pt_ij, THor_ij, TVer_ij is mnp x pnp
         */
        /* Also compute et_i = Pt_ij^T e_ij + THor_ij^T e_ij + TVer_ij^T e_ij
         * Recall that e_ij is mnp x 1
         */
		
#ifdef TIMINGS
        start = clock();
#endif
		
        _dblzero(N, noGrids*3*3*Nsz); /* clear all N_i */
        _dblzero(et, noGrids*3*etranssz); /* clear all et_i */
		
		for(gridIt=0; gridIt<noGrids; gridIt++)
		{
			printf("N for grid %d\n");
			
			ptr1_psum=N+gridIt*Nsz*9; // set ptr1_psum to point to Pt^T*Pt
			ptr1_THorP=ptr1_psum+Nsz; // set ptr1_THorP to point to THor^T*Pt
			ptr1_TVerP=ptr1_THorP+Nsz; // set ptr1_TVerP to point to TVer^T*Pt
			ptr1_PTHor=ptr1_TVerP+Nsz; // set ptr1_PTHor to point to Pt^T*THor
			ptr1_hsum=ptr1_PTHor+Nsz; // set ptr1_hsum to point to THor^T*THor
			ptr1_TVerTHor=ptr1_hsum+Nsz; // set ptr1_TVerTHor to point to TVer^T*THor
			ptr1_PTVer=ptr1_TVerTHor+Nsz; //set ptr1_PTVer to point to Pt^T*TVer
			ptr1_THorTVer=ptr1_PTVer+Nsz; //set ptr1_THorTVer to point to THor^T*TVer
			ptr1_vsum=ptr1_THorTVer+Nsz; //set ptr1_vsum to point to TVer^T*TVer
			
			ptr2_0 = et+gridIt*etranssz*3; // set ptr2_0 to point to et_0
			ptr2_1 = ptr2_0 + etranssz; // set ptr2_1 to point to et_1
			ptr2_2 = ptr2_1 + etranssz; // set ptr2_2 to point to et_2
			
			int curPoints;
			int prevPoints = n;
			int noPrevVisPoints = nvis;
			int noCurVisPoints;
			
			for(i=0; i<gridIt; i++)
			{
				prevPoints += grid_cols[i]*grid_rows[i];
				noPrevVisPoints += nVisGrids[i];
			}
			curPoints = prevPoints;
			curPoints += grid_cols[gridIt]*grid_rows[gridIt];
			noCurVisPoints = noPrevVisPoints;
			noCurVisPoints += nVisGrids[gridIt];
			
			for(i=prevPoints; i<curPoints; ++i)
			{
				nnz=sba_crsm_row_elmidxs(&idxij, i, rcidxs, rcsubs); /* find nonzero Pt_ij,THor_ij, TVer_ij j=0...m-1 */
				
				for(j=0; j<nnz; ++j)
				{
					/* set ptr3 to point to Pt_ij, actual column number in rcsubs[j] */
					int r = idxij.val[rcidxs[j]];
					//printf("Point:%d, r:%d, prevVis:%d\n", i, r, noPrevVisPoints);
					if(r<noPrevVisPoints || r>=noCurVisPoints)
					{
						printf("Error in computation of matrix N!");
					}
					
					//printf("Computing N: vis:%d\n", r);
					
					ptr3_pt=jac + nvis*Bsz + (r-nvis)*(Ptsz+THorsz+TVersz);
					ptr3_thor=jac + nvis*Bsz + (r-nvis)*(Ptsz+THorsz+TVersz) + Ptsz;
					ptr3_tver=jac + nvis*Bsz + (r-nvis)*(Ptsz+THorsz+TVersz) + Ptsz + THorsz;
					
					/*printf("********pt%d %d********\n", i, rcsubs[j]);
					for(ii=0; ii<mnp; ii++)
					{
						for(jj=0; jj<pnp; jj++)
						{
							printf("%f ", ptr3_pt[ii*pnp+jj]);
						}
						printf("\n");
					}
					printf("********thor%d %d********\n", i, rcsubs[j]);
					for(ii=0; ii<mnp; ii++)
					{
						for(jj=0; jj<pnp; jj++)
						{
							printf("%f ", ptr3_thor[ii*pnp+jj]);
						}
						printf("\n");
					}
					printf("********tver%d %d********\n", i, rcsubs[j]);
					for(ii=0; ii<mnp; ii++)
					{
						for(jj=0; jj<pnp; jj++)
						{
							printf("%f ", ptr3_tver[ii*pnp+jj]);
						}
						printf("\n");
					}*/
					
					double sum_p, sum_ph, sum_pv, sum_h, sum_hv, sum_v;
					
					for(ii=0; ii<pnp; ++ii)
					{
						for(jj=0; jj<pnp; ++jj)
						{
							sum_p = sum_ph = sum_pv = sum_h = sum_hv = sum_v = 0.0;
							for(k=0; k<mnp; ++k)
							{
								sum_p += ptr3_pt[k*pnp+ii] * ptr3_pt[k*pnp+jj];
								sum_ph += ptr3_pt[k*pnp+ii] * ptr3_thor[k*pnp+jj];
								sum_pv += ptr3_pt[k*pnp+ii] * ptr3_tver[k*pnp+jj];
								sum_h += ptr3_thor[k*pnp+ii] * ptr3_thor[k*pnp+jj];
								sum_hv += ptr3_thor[k*pnp+ii] * ptr3_tver[k*pnp+jj];
								sum_v += ptr3_tver[k*pnp+ii] * ptr3_tver[k*pnp+jj];
							}
							ptr1_psum[ii*pnp+jj]+=sum_p;
							ptr1_THorP[ii*pnp+jj]+=sum_ph;
							ptr1_TVerP[ii*pnp+jj]+=sum_pv;
							ptr1_PTHor[ii*pnp+jj]+=sum_ph;
							ptr1_hsum[ii*pnp+jj]+=sum_h;
							ptr1_TVerTHor[ii*pnp+jj]+=sum_hv;
							ptr1_PTVer[ii*pnp+jj]+=sum_pv;
							ptr1_THorTVer[ii*pnp+jj]+=sum_hv;
							ptr1_vsum[ii*pnp+jj]+=sum_v;
						}
					}
					
					
					double sum_0, sum_1, sum_2;
					
					ptr4=e + idxij.val[rcidxs[j]]*esz; /* set ptr4 to point to e_ij */
					/* compute Pt_ij^T e_ij and add it to et_0 */
					for(ii=0; ii<pnp; ++ii)
					{
						sum_0 = sum_1 = sum_2 = 0.0;
						for(jj=0, sum=0.0; jj<mnp; ++jj)
						{
							sum_0+=ptr3_pt[jj*pnp+ii]*ptr4[jj];
							sum_1+=ptr3_thor[jj*pnp+ii]*ptr4[jj];
							sum_2+=ptr3_tver[jj*pnp+ii]*ptr4[jj];
						}
						ptr2_0[ii]+=sum_0;
						ptr2_1[ii]+=sum_1;
						ptr2_2[ii]+=sum_2;
					}
					
				}
			}
		}
			/* Compute ||J^T e||_inf and ||p||^2 */
			for(i=0, p_L2=ebtrans_inf=0.0; i<nvars; ++i){
				if(ebtrans_inf < (tmp=FABS(ebtrans[i]))) ebtrans_inf=tmp;
				p_L2+=p[i]*p[i];
			}
			//p_L2=sqrt(p_L2);
			
			/* save diagonal entries so that augmentation can be later canceled.
			 * Diagonal entries are in V_i
			 */
		printf("save diagonal\n");
		diagV = diagVN;
		diagN = diagVN + n*pnp;
		
			for(i=0; i<n; ++i){
				ptr1=V + i*Vsz; // set ptr1 to point to V_i
				ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
				for(j=0; j<pnp; ++j)
					ptr2[j]=ptr1[j*pnp+j];
			}
			
			for(gridIt=0; gridIt<noGrids; gridIt++)
			{
				for(i=0; i<3; ++i)
				{
					ptr1=N + gridIt*9*Nsz + (i*4)*Nsz; // set ptr1 to point to N_i
					ptr2=diagN + gridIt*3*pnp + i*pnp; // set ptr2 to point to diagN_i
					for(j=0; j<pnp; ++j)
					{
						ptr2[j]=ptr1[j*pnp+j];
						printf("diagN:%f\n", ptr2[j]);
					}
				}
			}
			
			/*
			 if(!(itno%100)){
			 printf("Current estimate: ");
			 for(i=0; i<nvars; ++i)
			 printf("%.9g ", p[i]);
			 printf("-- errors %.9g %0.9g\n", eb_inf, p_eL2);
			 }
			 */
			
			/* check for convergence */
			if((ebtrans_inf <= eps1)){
				dp_L2=0.0; /* no increment for p in this case */
				stop=1;
				break;
			}
			
			/* compute initial damping factor */
			if(itno==0){
				for(i=0, tmp=DBL_MIN; i<nvars; ++i)
					if(diagVN[i]>tmp) tmp=diagVN[i]; /* find max diagonal element */
				mu=tau*tmp;
			}
			
			/* determine increment using adaptive damping */
			while(1){
				/* augment V,N */
				printf("augment VN\n");
				for(i=0; i<n; ++i){
					ptr1=V + i*Vsz; // set ptr1 to point to V_i
					for(j=0; j<pnp; ++j)
						ptr1[j*pnp+j]+=mu;
				}
				
				for(gridIt=0; gridIt<noGrids; gridIt++)
				{
					for(i=0; i<3; ++i)
					{
						ptr1=N + gridIt*Nsz*9 + (i*4)*Nsz; // set ptr1 to point to N_i
						for(j=0; j<pnp; ++j)
							ptr1[j*pnp+j]+=mu;
					}
				}
				
				
				/* solve the linear systems V*_i db_i = eb_i to compute the db_i */
				for(i=nsolved=0; i<n; ++i){
					printf("solveV:%d\b", i);
					ptr1=V + i*Vsz; // set ptr1 to point to V_i
					ptr2=eb + i*ebsz; // set ptr2 to point to eb_i
					ptr3=dp + i*pnp; // set ptr3 to point to db_i
					
					//nsolved+=sba_Axb_LU(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_LU;
					nsolved+=sba_Axb_Chol(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_Chol;
					//nsolved+=sba_Axb_BK(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_BK;
					//nsolved+=sba_Axb_QRnoQ(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_QRnoQ;
					//nsolved+=sba_Axb_QR(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_QR;
					//nsolved+=sba_Axb_SVD(ptr1, ptr2, ptr3, pnp, 0); linsolver=sba_Axb_SVD;
					//nsolved+=(sba_Axb_CG(ptr1, ptr2, ptr3, pnp, (3*pnp)/2, 1E-10, SBA_CG_JACOBI, 0) > 0); linsolver=(PLS)sba_Axb_CG;
					
					++nlss;
				}
				
				for(gridIt=0; gridIt<noGrids; gridIt++)
				{
					printf("solveN:%d\n", gridIt);
					int count = 0;
					double *Ntmp = (double *) emalloc(pnp*pnp*3*3*sizeof(double));
					ptr1 = dpt + gridIt*pnp*3;
					
					ptr3=et + gridIt*etranssz*3;
					
					for(i=0; i<3; i++)
					{
						ptr2_0=N + gridIt*Nsz*9 + i*3*Nsz; // set ptr2 to point to (N*_i)^-1
						ptr2_1=N + gridIt*Nsz*9 + i*3*Nsz + Nsz; // set ptr2 to point to (N*_i)^-1
						ptr2_2=N + gridIt*Nsz*9 + i*3*Nsz + 2*Nsz; // set ptr2 to point to (N*_i)^-1
						
						for(j=0; j<pnp; j++)
						{
							for(k=0; k<pnp; k++)
							{
								Ntmp[count]=ptr2_0[k*pnp+j];
								count+=1;
							}
							for(k=0; k<pnp; k++)
							{
								Ntmp[count]=ptr2_1[k*pnp+j];
								count+=1;
							}
							for(k=0; k<pnp; k++)
							{
								Ntmp[count]=ptr2_2[k*pnp+j];
								count+=1;
							}
						}
					}
					nsolved+=sba_Axb_Chol(Ntmp, ptr3, ptr1, pnp*pnp, 0); linsolver=sba_Axb_Chol;
					
					printf("*******A*******\n");
					for(i=0; i<pnp*pnp; i++)
					{
						for(j=0; j<pnp*pnp; j++)
						{
							printf("%f ", Ntmp[i*pnp*pnp + j]);
						}
						printf("\n");
					}
					printf("********b*******");
					for(i=0; i<pnp*pnp; i++)
					{
						printf("%f\n", ptr3[i]);
					}
					
					for(i=0; i<pnp*pnp; i++)
						printf("ptr1[%d]:%f\n", i, ptr1[i]);
					
					free(Ntmp);
				}
				
				printf("nsolved:%d\n", nsolved);
				
				if(nsolved==n+noGrids){
					
					/* parameter vector updates are now in dp */
					
					/* compute p's new estimate and ||dp||^2 */
					for(i=0, dp_L2=0.0; i<nvars; ++i){
						pdp[i]=p[i] + (tmp=dp[i]);
						dp_L2+=tmp*tmp;
					}
					
					for(i=0; i<pnp*pnp; i++)
						printf("ptr1[%d]:%f\n", i, pdp[i]);
					
					printf("a\n");
					//dp_L2=sqrt(dp_L2);
					
					printf("c\n");
					(*func)(pdp, &idxij, rcidxs, rcsubs, hx, adata); ++nfev; /* evaluate function at p + dp */
					printf("d\n");
					if(verbose>1)
						printf("mean reprojection error %g\n", sba_mean_repr_error(totalNoPoints, mnp, x, hx, &idxij, rcidxs, rcsubs));
					
					if(dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
						//if(dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
						//stop=2;
						//break;
					}
					printf("b\n");
					if(dp_L2>=(p_L2+eps2)/SBA_EPSILON_SQ){ /* almost singular */
						//if(dp_L2>=(p_L2+eps2)/SBA_EPSILON){ /* almost singular */
						fprintf(stderr, "SBA: the matrix of the augmented normal equations is almost singular in sba_motstr_levmar_x(),\n"
								"     minimization should be restarted from the current solution with an increased damping term\n");
						retval=SBA_ERROR;
						goto freemem_and_return;
					}
					
					/* ### compute ||e(pdp)||_2 */
					if(covx==NULL)
						pdp_eL2=nrmL2xmy(hx, x, hx, nobs); /* hx=x-hx, pdp_eL2=||hx|| */
					else
						pdp_eL2=nrmCxmy(hx, x, hx, wght, mnp, totalNVis); /* hx=wght*(x-hx), pdp_eL2=||hx|| */
					if(!SBA_FINITE(pdp_eL2)){
						if(verbose) /* identify the offending point projection */
							sba_print_inf(hx, m, mnp, &idxij, rcidxs, rcsubs);
						
						stop=7;
						break;
					}
					
					for(i=0, dL=0.0; i<nvars; ++i)
						dL+=dp[i]*(mu*dp[i]+ebtrans[i]);
					
					dF=p_eL2-pdp_eL2;
					
					if(verbose>1)
						printf("\ndamping term %8g, gain ratio %8g, errors %8g / %8g = %g\n", mu, dL!=0.0? dF/dL : dF/DBL_EPSILON, p_eL2/nvis, pdp_eL2/nvis, p_eL2/pdp_eL2);
					
					if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
						tmp=(2.0*dF/dL-1.0);
						tmp=1.0-tmp*tmp*tmp;
						mu=mu*( (tmp>=SBA_ONE_THIRD)? tmp : SBA_ONE_THIRD );
						nu=2;
						
						/* the test below is equivalent to the relative reduction of the RMS reprojection error: sqrt(p_eL2)-sqrt(pdp_eL2)<eps4*sqrt(p_eL2) */
						if(pdp_eL2-2.0*sqrt(p_eL2*pdp_eL2)<(eps4_sq-1.0)*p_eL2) stop=4;
						
						for(i=0; i<nvars; ++i) /* update p's estimate */
							p[i]=pdp[i];
						
						for(i=0; i<nobs; ++i) /* update e and ||e||_2 */
							e[i]=hx[i];
						p_eL2=pdp_eL2;
						break;
					}
				} /* nsolved==n */
				
				/* if this point is reached, either at least one linear system could not be solved or
				 * the error did not reduce; in any case, the increment must be rejected
				 */
				
				mu*=nu;
				nu2=nu<<1; // 2*nu;
				if(nu2<=nu){ /* nu has wrapped around (overflown) */
					fprintf(stderr, "SBA: too many failed attempts to increase the damping factor in sba_str_levmar_x()! Singular Hessian matrix?\n");
					//exit(1);
					stop=6;
					break;
				}
				nu=nu2;
				
				/* restore V diagonal entries */
				for(i=0; i<n; ++i){
					ptr1=V + i*Vsz; // set ptr1 to point to V_i
					ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
					for(j=0; j<pnp; ++j)
						ptr1[j*pnp+j]=ptr2[j];
				}
				
				for(gridIt=0; gridIt<noGrids; gridIt++)
				{
					for(i=0; i<3; ++i)
					{
						ptr1=N + gridIt*Nsz*9 + (i*4)*Nsz; // set ptr1 to point to N_i
						ptr2=diagN + gridIt*3*pnp + i*pnp; // set ptr2 to point to diagN_i
						for(j=0; j<pnp; ++j)
							ptr1[j*pnp+j]=ptr2[j];
					}
				}
				
			} /* inner while loop */
			
			if(p_eL2<=eps3_sq) stop=5; // error is small, force termination of outer loop
		}
		
		if(itno>=itmax) stop=3;
		
		/* restore V diagonal entries */
		for(i=0; i<n; ++i){
			ptr1=V + i*Vsz; // set ptr1 to point to V_i
			ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
			for(j=0; j<pnp; ++j)
				ptr1[j*pnp+j]=ptr2[j];
		}
		
		for(gridIt=0; gridIt<noGrids; gridIt++)
		{
			for(i=0; i<3; ++i)
			{
				ptr1=N + gridIt*Nsz*9 + (i*4)*Nsz; // set ptr1 to point to N_i
				ptr2=diagN + gridIt*3*pnp + i*pnp; // set ptr2 to point to diagN_i
				for(j=0; j<pnp; ++j)
					ptr1[j*pnp+j]=ptr2[j];
			}
		}
		
		if(info){
			info[0]=init_p_eL2;
			info[1]=p_eL2;
			info[2]=ebtrans_inf;
			info[3]=dp_L2;
			for(i=0; i<n; ++i){
				ptr1=V + i*Vsz; // set ptr1 to point to V_i
				for(j=0; j<pnp; ++j)
					if(tmp<ptr1[j*pnp+j]) tmp=ptr1[j*pnp+j];
			}
			info[4]=mu/tmp;
			info[5]=itno;
			info[6]=stop;
			info[7]=nfev;
			info[8]=njev;
			info[9]=nlss;
		}
		//sba_print_sol(n, m, p, 0, pnp, x, mnp, &idxij, rcidxs, rcsubs);
		retval=(stop!=7)?  itno : SBA_ERROR;
		
	freemem_and_return: /* NOTE: this point is also reached via a goto! */
		
		/* free whatever was allocated */
		free(jac); free(V); free(N);
		free(e);
		free(dp);               
		free(rcidxs); free(rcsubs);
#ifndef SBA_DESTROY_COVS
		if(wght) free(wght);
#else
		/* nothing to do */
#endif /* SBA_DESTROY_COVS */
		
		free(hx); free(diagVN); free(pdp);
		if(fdj_data.hxx){ // cleanup
			free(fdj_data.hxx);
			free(fdj_data.func_rcidxs);
		}
		
		sba_crsm_free(&idxij);
		
		/* free the memory allocated by the linear solver routine */
		if(linsolver) (*linsolver)(NULL, NULL, NULL, 0, 0);
		
		return retval;
	}
