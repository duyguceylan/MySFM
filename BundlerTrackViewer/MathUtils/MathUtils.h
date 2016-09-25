/*
 *  MathUtils.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/27/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MATH_UTILS_H
#define _MATH_UTILS_H

#include "MyMatrix.h"

#include <Accelerate/Accelerate.h>
#include <numeric>
#include "minpack.h"
#include "Color.h"

typedef void (*fcn)(const int *m, const int *n, const double *x, double *fvec, int *iflag );

class MathUtils
{
private:
	
public:
	static Vec3uc generateColorFromValue(float value, float min, float max)
	{
		float v0, v1, v2, v3, v4;
		v0 = min + 0.0/4.0 * (max - min);
		v1 = min + 1.0/4.0 * (max - min);
		v2 = min + 2.0/4.0 * (max - min);
		v3 = min + 3.0/4.0 * (max - min);
		v4 = min + 4.0/4.0 * (max - min);
		
		Vec3uc col(255,255,255);
		
		unsigned char u;
		
		if (value < v0) col = Vec3uc(0, 0, 255);
		else if (value > v4) col = Vec3uc(255, 0, 0);
		
		else if (value <= v2) 
		{
			if (value <= v1) // [v0, v1]
			{
				u = (unsigned char) (255.0 * (value - v0) / (v1 - v0));
				col = Vec3uc(0, u, 255);
			}      
			else // ]v1, v2]
			{
				u = (unsigned char) (255.0 * (value - v1) / (v2 - v1));
				col = Vec3uc(0, 255, 255-u);
			}
		}
		else 
		{
			if (value <= v3) // ]v2, v3]
			{
				u = (unsigned char) (255.0 * (value - v2) / (v3 - v2));
				col = Vec3uc(u, 255, 0);
			}
			else // ]v3, v4]
			{
				u = (unsigned char) (255.0 * (value - v3) / (v4 - v3));
				col = Vec3uc(255, 255-u, 0);
			}
		}
		
		return col;
	};
	
	static void generateRandomColors(vector<Color> &colors, int noOfColors)
	{
		/* initialize random seed: */
		srand ( time(NULL) );
		
		for(int i=0; i<noOfColors; i++)
		{
			Color c(1.0, 1.0, 1.0);
			while(c[0]==1.0 && c[1]==1.0 && c[2]==1.0)
			{
				c[0] = (float)rand()/RAND_MAX;
				c[1] = (float)rand()/RAND_MAX;
				c[2] = (float)rand()/RAND_MAX;
			}
			colors.push_back(c);
		}
	};
	
	static float getDeterminentOfMatrix4f(Matrix4f m)
	{
		Vec4f firstRow(m[0][0], m[1][0], m[2][0], m[3][0]);
		Vec4f secondRow(m[0][1], m[1][1], m[2][1], m[3][1]);
		Vec4f thirdRow(m[0][2], m[1][2], m[2][2], m[3][2]);
		Vec4f fourthRow(m[0][3], m[1][3], m[2][3], m[3][3]);
		
		return dot(firstRow, crossThreeVec(secondRow, thirdRow, fourthRow));
	};

	
	static Vec4f crossThreeVec(Vec4f a, Vec4f b, Vec4f c)
	{
		float d1 = (b[2] * c[3]) - (b[3] * c[2]);
		float d2 = (b[1] * c[3]) - (b[3] * c[1]);
		float d3 = (b[1] * c[2]) - (b[2] * c[1]);
		float d4 = (b[0] * c[3]) - (b[3] * c[0]);
		float d5 = (b[0] * c[2]) - (b[2] * c[0]);
		float d6 = (b[0] * c[1]) - (b[1] * c[0]);
		
		return Vec4f(- a[1] * d1 + a[2] * d2 - a[3] * d3,
					 a[0] * d1 - a[2] * d4 + a[3] * d5,
					 - a[0] * d2 + a[1] * d4 - a[3] * d6,
					 a[0] * d3 - a[1] * d5 + a[2] * d6);
	};
	
	static void fitQuadraticPolynamial(vector<double> x, vector<double> px, float &a, float &b, float &c)
	{
		//fit a quadratic polynomial in the form of p(x) = ax^2+bx+c 
		int noX = x.size();
		double *A = new double[3*noX];
		//fill A
		for(int i=0;i<noX;i++)
		{
			A[i] = x[i]*x[i];
			A[noX+i] = x[i];
			A[noX*2+i] = 1.0;
		}
		
		//solve in least squares sense
		int noRowA = noX;
		int noColA = 3;
		int noColPx = 1;
		int lda = noRowA;
		int ldb = noRowA;
		int info;
		
		/* Query and allocate the optimal workspace */
		int lwork = -1;
		double wkopt;
        dgels_( "No transpose", &noRowA, &noColA, &noColPx, A, &lda, &(px[0]), &ldb, &wkopt, &lwork,
			   &info );
        lwork = (int)wkopt;
        double *work = (double*)malloc( lwork*sizeof(double) );
        /* Solve the equations A*X = B */
        dgels_( "No transpose", &noRowA, &noColA, &noColPx, A, &lda, &(px[0]), &ldb, work, &lwork,
			   &info );
		
		delete work;
		work = NULL;
		
		if(info==0)
		{
			a = px[0];
			b = px[1];
			c = px[2];
		}
	}
	
	static int solveLeastSquares(double *A, double *b, int *noRowA, int *noColA, int *noColB)
	{
		int lda = *noRowA; //leading dimension of A, we take the whole matrix
		int ldb = *noRowA; //leading dimension of B, equal to no of rows in A
		int info;
		
		/* Query and allocate the optimal workspace */
		int lwork = -1;
		int nRA = *noRowA;
		int nCA = *noColA;
		int nCB = *noColB;
		
		double wkopt;
        dgels_( "No transpose", &nRA, &nCA, &nCB, A, &lda, b, &ldb, &wkopt, &lwork,
			   &info );
        lwork = (int)wkopt;
        double *work = (double*)malloc( lwork*sizeof(double) );
        /* Solve the equations A*X = B */
        dgels_( "No transpose", &nRA, &nCA, &nCB, A, &lda, b, &ldb, work, &lwork,
			   &info );
		
		delete work;
		work = NULL;
		
		return info;
	}
	
	//A is a square matrix
	static int solveLinearSystem(float *A, float *b, int *noRowA, int *noColB)
	{
		int lda = *noRowA;
		int ldb = *noRowA;
		int *ipiv = new int[*noRowA];
		int info;
		int nRA = *noRowA;
		int nCB = *noColB;
		
		sgesv_( &nRA, &nCB, A, &lda, ipiv, b, &ldb, &info );
		return info;
	}
	
	//the eigenvectors are normalized
	static int eigenValueAnalysis(float *A, int *noRows, float *eigenVal, float *eigenVec)
	{
		int lda = *noRows;
		int ldvl = *noRows;
		int ldvr = *noRows;
		int info;
		int lwork;
        float wkopt;
        float* work;
        
		float *wr = new float[*noRows];
		float *wi = new float[*noRows];
		float *vl = new float[ldvl*(*noRows)];
		float *vr = new float[ldvr*(*noRows)];
		
        lwork = -1;
		sgeev_( "N", "V", noRows, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
			   &wkopt, &lwork, &info );
        lwork = (int)wkopt;
		
        work = (float*)malloc( lwork*sizeof(float) );
		
		/* Solve eigenproblem */
        sgeev_( "N", "V", noRows, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
			   work, &lwork, &info );
		
		if(info == 0)
		{
			int i,j,k;
			float p;
			
			for( i=0; i<*noRows; i++)
			{
				eigenVal[i] = wr[i];
				for(j=0;j<*noRows;j++)
				{
					eigenVec[i*(*noRows) + j] = vr[i*(*noRows) + j];
				}
			}
			
			//sort the eigenValues
			for ( i = 0; i < (*noRows)-1; i++ )
			{
				p = eigenVal[k=i];
				for( j = i+1; j < *noRows; j++ )
					if( eigenVal[j] >= p ) p = eigenVal[k=j];
				if( k != i ){
					eigenVal[k] = eigenVal[i];
					eigenVal[i] = p;
					for( j = 0; j < *noRows; j++ ){
						p = eigenVec[i*(*noRows) + j];
						eigenVec[i*(*noRows) + j] = eigenVec[k*(*noRows) + j];
						eigenVec[k*(*noRows) + j] = p;
					}
				}
			}
		}
		
		delete [] wr;
		delete [] wi;
		delete [] vl;
		delete [] vr;
		delete work;
		
		wr = NULL; wi = NULL; vl = NULL; vr = NULL;
		work = NULL;
		
		return info;
	}	
	
	static int computeSVD(double *A, int *noRows, int *noCols, double *u, double *s, double *vt)
	{
		int lda = *noRows;
		int ldu = *noRows;
		int ldvt = *noCols;
		
		int info;
		int lwork = -1;
        int *iwork = new int[8*(*noCols)];
		double wkopt;
        double* work;
		
		dgesdd_("S", noRows, noCols, A, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info);
		
		lwork = (int)wkopt;
		work = (double*)malloc( lwork*sizeof(double) );
		
		dgesdd_("S", noRows, noCols, A, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
		
		delete work;
		work = NULL;
		
		return info;
	}
	
	static int lastNullVectorT(double *AT, double * x,  int m, int n)
	{
		
		if(m<n-1) 
			return 0;
		
		int result; 
		char jobu = 'N';
		char jobvt = 'A';

		int lda=m, ldu=m, ldvt=n, lwork=8*max(m,n), info;
		int  i;
		double* work = new double[lwork+n*n+min(m,n)];
		double* vt	 = work+lwork;
		double* s	 = vt+n*n;
		dgesvd_(&jobu, &jobvt, &m, &n, AT, &lda, s, NULL, &ldu, vt, &ldvt, work, &lwork, &info);
		
		for(i = 0 ; i < n; i ++)
		{
			x[i] = vt[i*n+n-1];
		}
		
		result = (info!=0 || s[n-1] > 0.2 ||(s[n-2]<0.000001 &&  s[n-1]>0.01*s[n-2]))? 0 : 1;
		
		delete [] work;
		
		return result;
	}
	
	static void makeGaussianKernel(float sigma, float **kernel, int *windowsize)
	{
		int i, center;
		float x, fx, sum=0.0;
		
		*windowsize = 1 + 2 * ceil(2.5 * sigma);
		center = (*windowsize) / 2;
		
		(*kernel) = new float[*windowsize];
		
		for(i=0;i<(*windowsize);i++)
		{
			x = (float)(i - center);
			fx = exp(-0.5*x*x/(sigma*sigma)) / (sigma * sqrt(2.0*M_PI));
			(*kernel)[i] = fx;
			sum += fx;
		}
		
		for(i=0;i<(*windowsize);i++) 
			(*kernel)[i] /= sum;
	}
	
	static void closestRank2Matrix(Matrix3f FMat, Matrix3f &Fout)
	{
		double *u = new double[9];
		double *s = new double[9];
		double *vt = new double[9];
		memset(u, 0, sizeof(double)*9);
		memset(s, 0, sizeof(double)*9);
		memset(vt, 0, sizeof(double)*9);
		
		int noRows = 3;
		int noCols = 3;
		
		Fout = FMat;
		
		double *F = new double[9];
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				F[i*3+j] = FMat[i][j];
			}
		}
		
		int success = computeSVD(F, &noRows, &noCols, u, s, vt);
		
		if(success != 0)
		{
			printf("[closestRank2Matrix] could not compute SVD.\n");
			return ;
		}
		
		//s has diagonals a, b, c with a>=b>=c
		//Fout = u * s * vt
		Matrix3f uMat, sMat, vtMat;
		sMat.setZero();
		sMat[0][0] = s[0];
		sMat[1][1] = s[1];
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				uMat[i][j] = u[i*3+j];
				vtMat[i][j] = vt[i*3+j];
			}
		}
		
		Fout = uMat * sMat * vtMat;
		
		delete [] u;
		delete [] s;
		delete [] vt;
	}
	
	static void closestRank2MatrixSSV(Matrix3f FMat, Matrix3f &FOut)
	{
		double *u = new double[9];
		double *s = new double[9];
		double *vt = new double[9];
		memset(u, 0, sizeof(double)*9);
		memset(s, 0, sizeof(double)*9);
		memset(vt, 0, sizeof(double)*9);
		
		int noRows = 3;
		int noCols = 3;
		
		FOut = FMat;
		
		double *F = new double[9];
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				F[i*3+j] = FMat[i][j];
			}
		}
		printf("F:%f %f %f %f %f %f %f %f %f\n", F[0], F[3], F[6], F[1], F[4], F[7], F[2], F[5], F[8]);
		
		int success = computeSVD(F, &noRows, &noCols, u, s, vt);
		
		printf("u:%f %f %f %f %f %f %f %f %f\n", u[0], u[3], u[6], u[1], u[4], u[7], u[2], u[5], u[8]);
		printf("vt:%f %f %f %f %f %f %f %f %f\n", vt[0], vt[3], vt[6], vt[1], vt[4], vt[7], vt[2], vt[5], vt[8]);
		printf("s:%f %f %f %f %f %f %f %f %f\n", s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8]);
		
		if(success != 0)
		{
			printf("[closestRank2Matrix] could not compute SVD.\n");
			return ;
		}
		
		//s has diagonals a, b, c with a>=b>=c
		//Fout = u * s * vt
		Matrix3f uMat, sMat, vtMat;
		sMat.setZero();
		sMat[0][0] = (s[0]+s[1])/2.0;
		sMat[1][1] = (s[0]+s[1])/2.0;
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				uMat[i][j] = u[i*3+j];
				vtMat[i][j] = vt[i*3+j];
			}
		}
		
		FOut = uMat * sMat * vtMat;
		
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				F[i*3+j] = FOut[i][j];
			}
		}
		printf("F:%f %f %f %f %f %f %f %f %f\n", F[0], F[3], F[6], F[1], F[4], F[7], F[2], F[5], F[8]);
		
		delete [] u;
		delete [] s;
		delete [] vt;
	}
	
	static void lmdifDriver2(fcn func, int m, int n, double *xvec, double tol) 
	{
		int info;
		double *fvec;
		double gtol = 0, epsfcn = 0;
		int maxfev = 200 * (n + 1);
		double *diag;
		int mode = 1;
		double factor = 100;
		int nprint = 1;
		int nfev;
		double *fjac;
		int ldfjac = m;
		int *ipvt;
		double *qtf;
		double *wa1, *wa2, *wa3, *wa4;
		
		if (n > m) 
		{
			printf("Error: lmdif called with n > m\n");
			return;
		}
		
		fvec = (double *)malloc(sizeof(double) * m);
		diag = (double *)malloc(sizeof(double) * n);
		fjac = (double *)malloc(sizeof(double) * m * n);
		ipvt = (int *)malloc(sizeof(int) * n);
		qtf = (double *)malloc(sizeof(double) * n);
		wa1 = (double *)malloc(sizeof(double) * n);
		wa2 = (double *)malloc(sizeof(double) * n);
		wa3 = (double *)malloc(sizeof(double) * n);
		wa4 = (double *)malloc(sizeof(double) * m);
		
		lmdif_(func, &m, &n, xvec, fvec, &tol, &tol, &gtol, &maxfev, 
			   &epsfcn, diag, &mode, &factor, &nprint, &info, &nfev, 
			   fjac, &ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);
		
		switch (info) 
		{
			case 0:
				printf("Improper input parameters\n");
				break;
			case 1:
				printf("Sum of squares tolerance reached\n");
				break;
			case 2:
				printf("x is within tolerance\n");
				break;
			case 3:
				printf("Sum of squares and x are within tolerance\n");
				break;
			case 4:
				printf("fvec orthogonal\n");
				break;
			case 5:
				printf("max function calls made\n");
				break;
			case 6:
				printf("tolerance is too small (squares)\n");
				break;
			case 7:
				printf("tolerance is too small (x)\n");
				break;
			default:
				printf("???\n");
		}
		
		// Clean up
		free(fvec);
		free(diag);
		free(fjac);
		free(ipvt);
		free(qtf);
		free(wa1);
		free(wa2);
		free(wa3);
		free(wa4);
	}
	
	
};


#endif