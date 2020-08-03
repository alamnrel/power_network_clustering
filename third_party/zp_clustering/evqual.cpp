/* This function is based on the C code of L. Zelnik-Manor to estimate the nr. of clusters
based on eigenvector rotations. It computes the associated eigenvector alignment cost  
function for a number of rotation angles, which should be provided in the COLUMNS of the 
matrix theta (the other initial functions have been removed from this file and the code  
has been modified). The output is a vector containing quality value for each set of the 
rotation angles. */ 

/* The original header */
/************
 *  
 *  mex interface to compute the gradient of the eigenvectors 
 *  alignment quality
 *
 *  To mexify:   
 *               mex  evrot.cpp;
 *   
 *  [clusters,Quality,Vrot] = evrot(V,method);
 *   
 *  Input:
 *    V = eigenvecors, each column is a vector
 *    method = 1   gradient descent
 *             2   approximate gradient descent
 *
 *  Output:
 *    clusts - Resulting cluster assignment
 *    Quality = The final quality
 *    Vr = The rotated eigenvectors
 *
 * 
 *  Lihi Zelnik (Caltech) March 2005
 * 
 ************/

#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEBUG 0    /* set to 1 and mex to see print outs */
#define  EPS 2.2204e-16 


/**** DEBUG utility *****************************************/
void print_array(mxArray *A)
{
    int i,j,ind;
    
    double *p_A = mxGetPr(A);
    const size_t *idims =  mxGetDimensions(A);

    mexPrintf("\n");
    ind = 0;
    for( j=0; j < idims[1]; j++ ){
		mexPrintf("Column %d:     ", j);
        for( i=0; i < idims[0]; i++ ){
			mexPrintf("%g ", p_A[ind]);
            ind++;
        }
		mexPrintf("\n");
    }
	mexPrintf("\n");
}

/************************************************************/
/*** multiply matrices    ***/
mxArray* matrix_mult(mxArray *A, mxArray *B)
{
    mxArray  *lhs[1], *rhs[2];
                                                                                                                              
    rhs[0] = A;
    rhs[1] = B;
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=A*B */
                                                                                                          
    return lhs[0];
}

/******* compute Uab ***********/
/** Givens rotation for angles a to b **/
mxArray* build_Uab(mxArray *theta, int a, int b,
               int* ik, int* jk, int dim)
{        
    mxArray *Uab = mxCreateDoubleMatrix(dim, dim, mxREAL);
    double *p_Uab = mxGetPr(Uab);
    int ind, k,i,j,ind_ik,ind_jk;
    /* set Uab to be an identity matrix */
    for( j=0; j<dim; j++ ){
        ind = dim*j + j;
        p_Uab[ind] = 1.0;
    }        
    
    if( b < a ) {
        return Uab;
    }
    
    double *p_theta = mxGetPr(theta);
    double tt,u_ik;
    for( k=a; k<=b; k++ ){
        tt = p_theta[k];
        for( i=0; i<dim; i++ ){
            ind_ik = dim*ik[k] + i;
            ind_jk = dim*jk[k] + i;
            u_ik = p_Uab[ind_ik] * cos(tt) - p_Uab[ind_jk] * sin(tt);
            p_Uab[ind_jk] = p_Uab[ind_ik] * sin(tt) + p_Uab[ind_jk] * cos(tt);
            p_Uab[ind_ik] = u_ik;
        }                       
    }
    return Uab;
}

/** Rotate vecotrs in X with Givens rotation according to angles in theta **/
mxArray* rotate_givens(mxArray *X, mxArray *theta, 
                       int* ik, int* jk, int angle_num, int dim)
{
    mxArray *G = build_Uab(theta, 0, angle_num-1,ik,jk,dim);
    mxArray *Y = matrix_mult(X,G);        
    mxDestroyArray(G);
    return Y;
}
  
/******** alignment quality ***********/
double evqual(mxArray *X,
              int *ik, int *jk,
              int dim,int ndata)
{
    double *p_X = mxGetPr(X);
    /* take the square of all entries and find max of each row */
    double *max_values = (double*)mxCalloc(ndata,sizeof(double));
    int *max_index = (int*)mxCalloc(ndata,sizeof(int));
    int i,j,ind = 0; 
    for( j=0; j<dim; j++ ){  /* loop over all columns */
        for( i=0; i<ndata; i++ ){ /* loop over all rows */ 
            if( max_values[i] <= p_X[ind]*p_X[ind]  ){                
                max_values[i] = p_X[ind]*p_X[ind];
                max_index[i] = j;
            }
            ind++;
        }
    }
    if( DEBUG )
        mexPrintf("Found max of each row\n");
   
    /* compute cost */
    double J=0;
    ind = 0;
    for( j=0; j<dim; j++ ){  /* loop over all columns */
        for( i=0; i<ndata; i++ ){ /* loop over all rows */
            J += p_X[ind]*p_X[ind]/max_values[i];
            ind++;
        }
    }
    
    J = J/ndata;
    if( DEBUG )
        mexPrintf("Computed quality = %g\n",J);
    
    mxFree(max_values);
    mxFree(max_index);        
    return J;
}


/***************************   main   */ 
  void
  mexFunction (
  int nlhs, mxArray* plhs[],
  int nrhs, const mxArray* prhs[])
{
  /* Make sure at most two output arguments are expected */
    if (nlhs < 1) {
        mexErrMsgTxt("Too few output arguments.");
        return;
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
        return;
    }
  /* Make sure input number is sufficient */
    if (nrhs < 2) {
        mexErrMsgTxt("Too few input arguments.");
        return;
    }
    if (nrhs > 2) {
        mexErrMsgTxt("Too many input arguments.");
        return;
    }    
    if( DEBUG )
        mexPrintf("Just starting\n");
        
    /* get the number and length of eigenvectors dimensions */
    mxArray *X; X = (mxArray*)prhs[0];
    const size_t *idims =  mxGetDimensions(X);
    const int ndata = idims[0];
    const int dim = idims[1];
    if( DEBUG )
        mexPrintf("Got %d vectors of length %d\n",ndata,dim);
    
    /* get the number of angles */
    int angle_num;
    angle_num = (int)(dim*(dim-1)/2);
    
    /* get the rotation angles */
    mxArray *theta; theta = (mxArray*)prhs[1];
    idims =  mxGetDimensions(theta);
    const int n_col = idims[1];
    const int n_ang = idims[0];   
    if (n_ang != angle_num) {
        mexErrMsgTxt("Incorrect number of angles.");
        return;
    }    
    
    /* build index mapping */
    int i,j,k,ind;
    int* ik = (int*)mxCalloc(angle_num,sizeof(int));
    int* jk = (int*)mxCalloc(angle_num,sizeof(int));
    k=0;
    for( i=0; i<dim-1; i++ ){
        for( j=i+1; j<=dim-1; j++ ){
            ik[k] = i;
            jk[k] = j;
            k++;
        }
    }
    if( DEBUG )
        mexPrintf("Built index mapping for %d angles\n",k);
    
    /* definitions */
    double* Q = (double*)mxCalloc(n_col,sizeof(double));
    mxArray *Xrot, *theta0;
    double *p_theta = mxGetPr(theta);
    theta0 = mxCreateDoubleMatrix(n_ang, 1, mxREAL);
            
    for (i=0; i<n_col; i++){
        memcpy(mxGetPr(theta0), p_theta+(i)*n_ang, n_ang*sizeof(double));
        Xrot = rotate_givens(X,theta0,ik,jk,angle_num,dim);
        Q[i] = evqual(Xrot,ik,jk,dim,ndata); 
        mxDestroyArray(Xrot);
    }
	if (DEBUG){
		print_array(Xrot);		
	}
	
    /** prepare output **/
    plhs[0] = mxCreateDoubleMatrix(n_col, 1, mxREAL);
	double *p_Q = mxGetPr(plhs[0]);
	memcpy(p_Q, Q, n_col*sizeof(double));

    /* free allocated memory */
    mxFree(ik);
    mxFree(jk);
    mxDestroyArray(theta0);
    if( DEBUG )
        mexPrintf("Done evrot\n");
        
    return;
}


