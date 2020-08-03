/* This function is based on the C code of L. Zelnik-Manor to estimate the nr. of clusters
based on eigenvector rotations. It computes the gradient of the associated eigenvector 
alignment cost function (the other initial functions have been removed from this file and 
the code has been modified). */

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
        mexPrintf("Column %d:     ",j);
        for( i=0; i < idims[0]; i++ ){
            mexPrintf("%g ",p_A[ind]);
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
        
/*** compute the matrix A given X,U1,V,U2    ***/
mxArray* buildA(mxArray *X, mxArray *U1, mxArray *Vk, mxArray *U2)
{
    mxArray  *lhs[1], *rhs[2];
                                                                                                                              
    rhs[0] = Vk;
    rhs[1] = U2;
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=Vk*U2 */
                                                                                                                              
    rhs[0] = U1;
    rhs[1] = lhs[0];
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=U1*Vk*U2 */
    mxDestroyArray(rhs[1]);
                
    rhs[0] = X;
    rhs[1] = lhs[0];
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=X*U1*Vk*U2 */
    mxDestroyArray(rhs[1]);

    return lhs[0];
}
        
/******* compute V ***********/
/** Gradient of a single Givens rotation **/
mxArray* gradU(mxArray *theta, int k,
               int* ik, int* jk, int dim)
{
    
    double *p_theta = mxGetPr(theta);
    mxArray *V = mxCreateDoubleMatrix(dim, dim, mxREAL);
    double *p_V = mxGetPr(V);

    p_V[ik[k]+dim*ik[k]] = -sin(p_theta[k]);
    p_V[ik[k]+dim*jk[k]] = cos(p_theta[k]);
    p_V[jk[k]+dim*ik[k]] = -cos(p_theta[k]);
    p_V[jk[k]+dim*jk[k]] = -sin(p_theta[k]);

    return V;
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

/****** quality gradient *******************/
double evqualitygrad(mxArray *X, mxArray* theta,
                     int *ik, int *jk,
                     int angle_num, int angle_index,
                     int dim, int ndata)
{                   
    /* build V,U,A */
    mxArray *V = gradU(theta,angle_index,ik,jk,dim);
    if( DEBUG )
        mexPrintf("Computed gradU\n");
    
    mxArray *U1 = build_Uab(theta,0,angle_index-1,ik,jk,dim);
    mxArray *U2 = build_Uab(theta,angle_index+1,angle_num-1,ik,jk,dim);
    if( DEBUG )
        mexPrintf("Computed Uab\n");
    
    mxArray *A = buildA(X, U1, V, U2);
    double *p_A = mxGetPr(A);
    if( DEBUG )
        mexPrintf("Built A\n");
  
    /* get rid of no longer needed arrays */
    mxDestroyArray(V);
    mxDestroyArray(U1);
    mxDestroyArray(U2);
    
    /* rotate vecs according to current angles */   
    mxArray *Y = rotate_givens(X,theta,ik,jk,angle_num,dim);
    double *p_Y = mxGetPr(Y);
    if( DEBUG )
        mexPrintf("Rotated according to Givens successfully\n");
       
    /* find max of each row */
    double *max_values = (double*)mxCalloc(ndata,sizeof(double));
    int *max_index = (int*)mxCalloc(ndata,sizeof(int));
    int i,j, ind = 0;
    for( j=0; j<dim; j++ ){  /* loop over all columns */
        for( i=0; i<ndata; i++ ){ /* loop over all rows */
            if( max_values[i]*max_values[i] <= p_Y[ind]*p_Y[ind]  ){
                max_values[i] = p_Y[ind];
                max_index[i] = j;
            }
            ind++;
        }
    }
    if( DEBUG )
        mexPrintf("Found max of each row\n");           
    
    /* compute gradient */
    double dJ=0, tmp1, tmp2;
    ind = 0;
    for( j=0; j<dim; j++ ){  /* loop over all columns */
        for( i=0; i<ndata; i++ ){ /* loop over all rows */
            tmp1 = p_A[ind] * p_Y[ind] / (max_values[i]*max_values[i]);
            tmp2 = p_A[ndata*max_index[i]+i] * (p_Y[ind]*p_Y[ind]) / (max_values[i]*max_values[i]*max_values[i]);
            dJ += tmp1-tmp2;
            ind++;
        }
    }
    dJ = 2*dJ/ndata;
    if( DEBUG )
        mexPrintf("Computed gradient = %g\n",dJ);
    
    mxFree(max_values);
    mxFree(max_index);
    mxDestroyArray(Y);
    mxDestroyArray(A);
    
    return dJ;
}


/***************************   //////////////// main   */
  
  void
  mexFunction (
  int nlhs, mxArray* plhs[],
  int nrhs, const mxArray* prhs[])
{
  /* Make sure at most one output argument is expected */
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
        mexPrintf("Got %d vectors of length %d\n",dim,ndata);
    
    /* get the number of angles */
    int angle_num;
    angle_num = (int)(dim*(dim-1)/2);
    
    /* get the rotation angles at which the gradient should be evaluated*/
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
    double dQ;
    int d;  // Givens angle (problem dimension)
    double* gradQ = (double*)mxCalloc(angle_num,sizeof(double));

    for( d = 0; d < angle_num; d++ ){
        dQ = evqualitygrad(X,theta,ik,jk,angle_num,d,dim,ndata);
        gradQ[d] = dQ;
    }        
    
    /** prepare output **/
    plhs[0] = mxCreateDoubleMatrix(angle_num, 1, mxREAL);
    memcpy(mxGetPr(plhs[0]), gradQ, angle_num*sizeof(double)); 
    
    /* free allocated memory */
    mxFree(ik);
    mxFree(jk);

    if( DEBUG )
        mexPrintf("Done evgradJ\n");
        
    return;
}


