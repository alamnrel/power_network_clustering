/* This function is based on the C code of L. Zelnik-Manor to estimate the nr. of clusters
based on eigenvector rotations. Its output is the rotation matrix obtained from the input
givens angles theta. */ 

/* The original header */
/**********
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
 *  Lihi Zelnik (Caltech) March.2005
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
        mexPrintf("Column %d:     ",j);
        for( i=0; i < idims[0]; i++ ){
            mexPrintf("%g ",p_A[ind]);
            ind++;
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");
}

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

 
/***************************   //////////////// main   */
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

    /* get theta */
    mxArray *theta_new; theta_new = (mxArray*)prhs[0];     
    
    /* get data dimension */
    int dim;
    dim = (int) mxGetScalar(prhs[1]);
    
    /* get the number of angles */
    int angle_num;
    angle_num = (int)(dim*(dim-1)/2);     
    if( DEBUG )
        mexPrintf("Angle number is %d\n",angle_num);
        
    /* build index mapping */
    int i,j,k;
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

  
    mxArray *Uab = build_Uab(theta_new, 0, angle_num-1, ik, jk, dim);  
    plhs[0] = mxCreateDoubleMatrix(dim, dim, mxREAL);
    memcpy(mxGetPr(plhs[0]), mxGetPr(Uab), dim*dim*sizeof(double));
    
    /* free allocated memory */
    mxFree(ik);
    mxFree(jk);

    if( DEBUG )
        mexPrintf("Done evrot\n");
        
    return;
}


