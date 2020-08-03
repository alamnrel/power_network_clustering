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
 * 
 ************/

#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

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



/******* add a single Givens rotation to a previous one (calc U1*U2) ***********/
mxArray* U_add_single(mxArray *U1, 
                      mxArray *theta, 
                      int k,
                      int* ik, int* jk, int dim)
{            
    mxArray *Uab = mxDuplicateArray(U1);
    double *p_Uab = mxGetPr(Uab);
    int i,ind_ik,ind_jk;
    
    double *p_theta = mxGetPr(theta);
    double tt,u_ik;
    tt = p_theta[k];
    for( i=0; i<dim; i++ ){
        ind_ik = dim*ik[k] + i;
        ind_jk = dim*jk[k] + i;
        u_ik = p_Uab[ind_ik] * cos(tt) - p_Uab[ind_jk] * sin(tt);
        p_Uab[ind_jk] = p_Uab[ind_ik] * sin(tt) + p_Uab[ind_jk] * cos(tt);
        p_Uab[ind_ik] = u_ik;                               
    }
    return Uab;
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
                     int dim, int ndata, int siz_batch, 
                     int ibatch, bool direction)
{   
    /* Slice the ibatch(th) batch of the size siz_batch out of X */
    mxArray* Xbatch;
	Xbatch = mxCreateDoubleMatrix(siz_batch, dim, mxREAL);
    double *p_X = mxGetPr(X);
    int start, i, j;
    if(direction)
        start = siz_batch*ibatch;
    else
		start = ndata - siz_batch*(1+ibatch);      
	double *p_Xbatch = mxGetPr(Xbatch);
    for(j=0; j<dim; j++){
		memcpy(p_Xbatch + j*siz_batch, p_X + start + j*ndata, siz_batch*sizeof(double));
    }
	// const size_t *idims = mxGetDimensions(Xbatch);
	// const int ndata_batch = idims[0];
	// const int dim_batch = idims[1];
    // if( DEBUG ){
        // for (i = 0; i < dim_batch; i++){
            // for (j = 0; j < ndata_batch; j++)
                // printf("%1.4f ", *(p_Xbatch + ndata_batch*i + j));	//	
            // printf("\n ");
        // }
    // }

    /* build V,U,A */
    mxArray *V = gradU(theta,angle_index,ik,jk,dim);
    if( DEBUG )
        mexPrintf("Computed gradU\n");
    
    mxArray *U1 = build_Uab(theta,0,angle_index-1,ik,jk,dim);
    mxArray *U2 = build_Uab(theta,angle_index+1,angle_num-1,ik,jk,dim);
    if( DEBUG )
        mexPrintf("Computed Uab\n");
    
    mxArray *A = buildA(Xbatch, U1, V, U2);
    double *p_A = mxGetPr(A);
    if( DEBUG )
        mexPrintf("Built A\n");
  
    /* get rid of no longer needed arrays */
    mxDestroyArray(V);
    mxDestroyArray(U1);
    mxDestroyArray(U2);
    
    /* rotate vecs according to current angles */   
    mxArray *Y = rotate_givens(Xbatch,theta,ik,jk,angle_num,dim);
    double *p_Y = mxGetPr(Y);
    if( DEBUG )
        mexPrintf("Rotated according to Givens successfully\n");
       
    /* find max of each row */
    double *max_values = (double*)mxCalloc(siz_batch,sizeof(double));
    int *max_index = (int*)mxCalloc(siz_batch,sizeof(int));
    int ind = 0;
    for( j=0; j<dim; j++ ){  /* loop over all columns */
        for( i=0; i<siz_batch; i++ ){ /* loop over all rows */
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
        for( i=0; i<siz_batch; i++ ){ /* loop over all rows */
            tmp1 = p_A[ind] * p_Y[ind] / (max_values[i]*max_values[i]);
            tmp2 = p_A[siz_batch*max_index[i]+i] * (p_Y[ind]*p_Y[ind]) / (max_values[i]*max_values[i]*max_values[i]);
            dJ += tmp1-tmp2;
            ind++;
        }
    }
    dJ = 2*dJ/siz_batch;
    if( DEBUG )
        mexPrintf("Computed gradient = %g\n",dJ);
    
    mxFree(max_values);
    mxFree(max_index);
    mxDestroyArray(Y);
    mxDestroyArray(A);
	mxDestroyArray(Xbatch);
    return dJ;
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

  
/***************************   //////////////// main   */
  
  void
  mexFunction (
  int nlhs, mxArray* plhs[],
  int nrhs, const mxArray* prhs[])
{
  /* Make sure at most two output arguments are expected */
    if (nlhs < 3) {
        mexErrMsgTxt("Too few output arguments.");
        return;
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
        return;
    }
  /* Make sure input number is sufficient */
    if (nrhs < 1) {
        mexErrMsgTxt("Too few input arguments.");
        return;
    }
    if (nrhs > 3) {
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
    mxArray *theta = mxCreateDoubleMatrix(1,angle_num,mxREAL);
    if( DEBUG )
        mexPrintf("Angle number is %d\n",angle_num);
    
	/* get iteration scaling */
	int iter_scal;
	if (nrhs >= 2)
		iter_scal = (int)mxGetScalar(prhs[1]);
	else
		iter_scal = 1;
	if (DEBUG)
		mexPrintf("Got iteration scaling %d\n", iter_scal);

	/* get siz_batch */
	int siz_batch;
	if (nrhs >= 3)
		siz_batch = (int)mxGetScalar(prhs[2]);
	else
		siz_batch = ndata;
	if (DEBUG)
		mexPrintf("Got siz_batch %d\n", siz_batch);
    
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
    int max_iter = 1000*iter_scal; 
    int nbatch = ndata/siz_batch;
    bool direction;
    double dQ,Q,Q_new,Q_old1;
    double alpha;
    const double dmp1 = 0.9;
    const double dmp2 = 0.992;    
    const double eps = 1e-8;   //sqrt(EPS)
    int iter,d;
    mxArray *Xrot, *Xtmp1, *Xtmp2, *theta_new;
    Xtmp1 = mxCreateDoubleMatrix(dim, ndata, mxREAL);   // to contain transposed Xtmp2
    Xtmp2 = mxDuplicateArray(X);   // to contain permuted X   
    theta_new = mxDuplicateArray(theta);    
    
    double *p_theta = mxGetPr(theta);
    double *p_theta_new = mxGetPr(theta_new);
    double* m_acc = (double*)mxCalloc(angle_num,sizeof(double));
    double* v_acc = (double*)mxCalloc(angle_num,sizeof(double));
    double* m_acc_new = (double*)mxCalloc(angle_num,sizeof(double));
    double* v_acc_new = (double*)mxCalloc(angle_num,sizeof(double));    
    double* m_curr = (double*)mxCalloc(angle_num,sizeof(double));
    double* v_curr = (double*)mxCalloc(angle_num,sizeof(double));    
    double* m_corr = (double*)mxCalloc(angle_num,sizeof(double));
    double* v_corr = (double*)mxCalloc(angle_num,sizeof(double)); 
    double dang;  // for debug
    for(d = 0; d < angle_num; d++){
        m_acc[d] = 0;
        v_acc[d] = 0;
    }
    double* tmp_row = (double*)mxCalloc(dim,sizeof(double));
	srand(time(NULL));

    Q = evqual(X,ik,jk,dim,ndata); /* initial quality */
    if( DEBUG )
        mexPrintf("Q = %g\n",Q);
    Q_old1 = Q;
    alpha = 0.1;     
    iter = 0;
    direction = true;
    while( iter < max_iter ){ /* iterate to refine quality */                
        /* Simple SGD through true derivative */  
        for(i=0; i<nbatch; i++){
            iter++;
            for( d = 0; d < angle_num; d++ ){
                dQ = evqualitygrad(Xtmp2,theta,ik,jk,angle_num,d,dim,ndata,siz_batch,i,direction);  
                m_curr[d] = dmp1*m_acc[d] + (1-dmp1)*dQ;
                v_curr[d] = dmp2*v_acc[d] + (1-dmp2)*dQ*dQ;
                m_corr[d] = m_curr[d]/(1-pow(dmp1,iter));
                v_corr[d] = v_curr[d]/(1-pow(dmp2,iter));                            
            }
            for( d = 0; d < angle_num; d++ ){
                dang = alpha*m_corr[d]/(sqrt(v_corr[d])+eps);
                p_theta_new[d] = p_theta[d] - dang;
                m_acc_new[d] = m_curr[d];
                v_acc_new[d] = v_curr[d];
            }                  
            Xrot = rotate_givens(X,theta_new,ik,jk,angle_num,dim);
            Q_new = evqual(Xrot,ik,jk,dim,ndata); 
            mxDestroyArray(Xrot);   
            if( Q_new < Q ){
                for( d = 0; d < angle_num; d++ ){
                    p_theta[d] = p_theta_new[d];
                    m_acc[d] = m_acc_new[d];
                    v_acc[d] = v_acc_new[d];                
                }
                Q = Q_new;
            }                             
        }
        if (iter>2)
            if (fabs(Q_old1 - Q) < 1e-3)  // Q can only decrease!
                break;
        Q_old1 = Q;            
        direction = !direction; 
        
        // Shuffle rows of Xtmp2 (first transpose to make rows become columns)	
		mexCallMATLAB(1, &Xtmp1, 1, &Xtmp2, "transpose");
        double *p_Xtmp1 = mxGetPr(Xtmp1);
        if (ndata > 1) {
            for (size_t i = ndata - 1; i > 0; i--) {
				size_t j = (unsigned int)(rand() % (i + 1));
				memcpy(tmp_row, p_Xtmp1 + dim*j, dim*sizeof(double));
				memcpy(p_Xtmp1 + dim*j, p_Xtmp1 + dim*i, dim*sizeof(double));
				memcpy(p_Xtmp1 + dim*i, tmp_row, dim*sizeof(double));
            }
        }
		mexCallMATLAB(1, &Xtmp2, 1, &Xtmp1, "transpose");  // transpose back
        //mexCallMATLAB(0, NULL,   1, &Xtmp1, "disp");         
		//mexCallMATLAB(0, NULL,   1, &Xtmp2, "disp");  
    }
	mxFree(tmp_row); 
    mxDestroyArray(Xtmp1);   
    mxDestroyArray(Xtmp2);            
    
    if( DEBUG > 1 )
        mexPrintf("Done after %d iterations, Quality is %g\n",iter,Q);
    
    Xrot = rotate_givens(X,theta,ik,jk,angle_num,dim);
    
    /** prepare output **/
    plhs[0] = mxCreateDoubleMatrix(ndata, dim, mxREAL);
    double *p_Xout = mxGetPr(plhs[0]);
    double *p_Xrot = mxGetPr(Xrot);
    memcpy(p_Xout, p_Xrot, ndata*dim*sizeof(double));       
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *p_Q = mxGetPr(plhs[1]);
    *p_Q = Q;
    plhs[2] = mxCreateDoubleMatrix(1, angle_num, mxREAL);
    double *p_ang = mxGetPr(plhs[2]);    
    memcpy(p_ang, p_theta, angle_num*sizeof(double));       
        
    /* free allocated memory */
    mxFree(ik);
    mxFree(jk);   
    mxFree(m_acc);
    mxFree(m_curr);
    mxFree(v_acc);    
    mxFree(v_curr);
    mxFree(m_corr);    
    mxFree(v_corr);   
	mxFree(m_acc_new);
	mxFree(v_acc_new);
    mxDestroyArray(theta_new);		

    if( DEBUG )
        mexPrintf("Done evrot\n");
        
    return;
}
 

