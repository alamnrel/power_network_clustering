/*
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
*/

#include <stdio.h>
#include <math.h>       // ceil
#include <algorithm>    // std::max
#include "mex.h"
#include "matrix.h"     // mxAssert

using namespace std;
void threshold(double *, int, int, double *, int *);
void tvMultiClass(double *wValTriu, mwIndex* irs, mwIndex* jcs, double *F, double *, double *, int nVertices, int nClasses, double *tvMultiClass);
void balancesMultiClass(double *F, double *vertexWeights, int beta, bool CheegerFlag, int nVertices, int nClasses, double* eta, double *theta, double* thetaBar, double *balances);
bool comparePair(std::pair<double, double> a, std::pair<double, double> b);

// Primal infeasibility are deduced from the intermediate computations done for the primal variable update.
// Primal infeasibility is computed based on the infeasibility of the previous iterate and it is assumed that x (starting point) is primal feasible (to avoid an extra computation).
// Dual infeasibility printed is that of the previous iterate (y instead of yNew) because this is explicitely computed in the update step of the primal iterate

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    // check number of parameters.
// 	if (nrhs != 26 || nlhs != 6)
// 	{
// 		mexWarnMsgTxt("Usage: [xnew] = mexInnerProblemWithBounds(WTriu, sgBal, m, M, lambdas, Tau, Sigma, x, y, nVertices, nEdges, nClasses)");
// 		return;
// 	}
    
    if(!mxIsSparse(prhs[0])) { mexWarnMsgTxt("Weight matrix is not sparse");}
    
    // read input
    double* wValTriu = mxGetPr(prhs[0]);        // get edge weights
    mwIndex* irs = mxGetIr(prhs[0]);        // get rows of W
    mwIndex* jcs = mxGetJc(prhs[0]);        // get cols of W
    
    double *sgBal = mxGetPr(prhs[1]);
    double m = mxGetScalar(prhs[2]);  // lower bound
    double M = mxGetScalar(prhs[3]);  // upper bound
    double* lambdas = mxGetPr(prhs[4]); // current lambda values
    double* Tau = mxGetPr(prhs[5]); // Tau
    double* Sigma = mxGetPr(prhs[6]); // Sigma
    
    double *x0 = mxGetPr(prhs[7]);
    double *y0 = mxGetPr(prhs[8]);
    
    // get sizes
    int nVertices = mxGetScalar(prhs[9]);
    int nEdges = mxGetScalar(prhs[10]);
    int nClasses = mxGetScalar(prhs[11]);
    int maxIters = mxGetScalar(prhs[12]);
    int imbalance = mxGetScalar(prhs[13]);
    double cost = mxGetScalar(prhs[14]);
    double *cutEdgesPL = mxGetPr(prhs[15]);
    double *degreePL = mxGetPr(prhs[16]);
    double *sgFMax = mxGetPr(prhs[17]);
    double gamma = mxGetScalar(prhs[18]);
    int nPrintIters = mxGetScalar(prhs[19]);
    double* eta = mxGetPr(prhs[20]);
    double* theta = mxGetPr(prhs[21]);
    double* thetaBar = mxGetPr(prhs[22]);
    double *vertexWeights = mxGetPr(prhs[23]);
    mxLogical *CheegerFlagPtr = mxGetLogicals(prhs[24]);
    bool CheegerFlag = CheegerFlagPtr[0];
    double initialBCut = mxGetScalar(prhs[25]);
    double simplexProjection = mxGetScalar(prhs[26]);
    double debugMode = mxGetScalar(prhs[27]);
    //double* deltaBounds = mxGetPr(prhs[27]);
    
    int nPrimalVars = nVertices*nClasses+nEdges*nClasses+2*nClasses;
    int nEdgesTimesNClasses = nEdges*nClasses;
    int twoTimesNClasses = 2*nClasses;
    int nVerticesTimesNClasses = nVertices*nClasses;
    //int nSimplexCnstrs = nVertices;  // it will be zero if projFlag is true
    int nInEqCnstrs = nClasses*2 + 2*nEdgesTimesNClasses;
    if(!CheegerFlag){ nInEqCnstrs += nClasses; }
    int nDualVars = nInEqCnstrs; // + nVertices;
    if(simplexProjection==0)
        nDualVars += nVertices;
    
    // prepare output
    double *bestFDelta, *finalFDelta, *degeneratedF; /* output variables*/
    double *degenerateFlag, *convergenceFlag, *nIters;
    plhs[0]         = mxCreateDoubleMatrix(nVerticesTimesNClasses+twoTimesNClasses,1,mxREAL); /* create the output for x */
    bestFDelta          = mxGetPr(plhs[0]);
    plhs[1]         = mxCreateDoubleMatrix(nVerticesTimesNClasses+twoTimesNClasses,1,mxREAL); /* create the output for x */
    finalFDelta          = mxGetPr(plhs[1]);
    //plhs[2]         = mxCreateLogicalScalar(false); /* create the output for x */
    plhs[2]         = mxCreateDoubleMatrix(1,1,mxREAL); /* create the output for x */
    degenerateFlag          = mxGetPr(plhs[2]);
    degenerateFlag[0] = 0;
    plhs[3]         = mxCreateDoubleMatrix(1,1,mxREAL); /* create the output for x */
    convergenceFlag          = mxGetPr(plhs[3]);
    convergenceFlag[0] = 0;
    plhs[4]         = mxCreateDoubleMatrix(1,1,mxREAL); /* create the output for x */
    nIters = mxGetPr(plhs[4]);
    plhs[5]         = mxCreateDoubleMatrix(nVerticesTimesNClasses,1,mxREAL); /* create the output for x */
    degeneratedF = mxGetPr(plhs[5]);
    double *yPDHG;
    plhs[6] = mxCreateDoubleMatrix(nDualVars,1,mxREAL);
    yPDHG = mxGetPr(plhs[6]);
    plhs[7] = mxCreateDoubleMatrix(nVerticesTimesNClasses+twoTimesNClasses,1,mxREAL); /* create the output for x */
    double *bestFDeltaInTermsOfCut          = mxGetPr(plhs[7]);
    plhs[8] = mxCreateDoubleMatrix(1,1,mxREAL); /* create the output for x */
    double *bestLambdaPtr = mxGetPr(plhs[8]);
    
    double *x = new double[nPrimalVars];
    double *xNew = new double[nPrimalVars];
    double *xBar = new double[nPrimalVars];
    double *y = new double[nDualVars];
    double *yNew = new double[nDualVars];
    double *FNew; // = new double[nVertices*nClasses];
    double *FBar;
    double *projFNew = new double[nVertices*nClasses];
    double *Fl; // = new double[nVertices];
    double *balanceLB = new double[nClasses];
    double *mTimesDeltaPlus = new double[nClasses];
    double *MTimesDeltaMinus = new double[nClasses];
    
    double *dummyPointer;
    
    int counter, counter1, counter2, dummyCounter, primalVarCounter, dualVarCounter, yEdgesCounter1, yEdgesCounter2, sgBalCounter, FlCounter, xBarCounter, ySimplexCounter;
    int i, j, l, maxIx, rowIndex;
    int *clusterSizes = new int[nClasses];
    double dummy, tempMax, rowSum, primalObj, dualObj;
    double *dualInfeasibility = new double[nPrimalVars];
    double maxDualInfeasibility, maxPrimalInfeasibility;
    double *primalInfeasibility = new double[nDualVars];
    double maxDescentCnstrInfeasibility, maxLBCnstrInfeasibility, maxUBCnstrInfeasibility, maxEdgeCnstrInfeasibility, maxSimplexCnstrInfeasibility;
    double *primalInfeasibilityOld = new double[nDualVars];
    for(int i=0; i<nDualVars; i++){ primalInfeasibilityOld[i] = 0.0; }
    double *penalty = new double[nClasses];
    double *cutPL = new double[nClasses];
    
    double *wTransAlpha = new double[nClasses];
    double *ATransYEdges1 = new double[nVertices];
    double *ATransYEdges2 = new double[nVertices];
    double *clusters = new double[nVertices*nClasses];
    double *balances = new double[nClasses];
    double *tvs = new double[nClasses];
    double BCut, lambda = 0;
    double *feasibleF = new double[nVertices*nClasses];
    double feasibleLambda = 0.0;
    
    if(debugMode==1) mexPrintf("initial lambdas: ");
    for(l=0; l<nClasses; l++){
        lambda += lambdas[l];
        if(debugMode==1) mexPrintf("%f\t", lambdas[l]);
    }
    if(debugMode==1) mexPrintf("\n");
    for(l=0; l<nClasses; l++){
        balances[l] = 0;
    }
    
//     double eta = 0.0, thetaBar = 0.0;
//     if(degreePL[0] > 0){
//         eta = 1.0; thetaBar = imbalance;
//     }
    double *lambdaTimesEta = new double[nClasses];
    for(l=0; l<nClasses; l++){
        lambdaTimesEta[l] = lambdas[l]*eta[l];
    }
    
    for(i=0; i<nPrimalVars; i++){ x[i] = x0[i]; }
    for(i=0; i<nDualVars; i++){ y[i] = y0[i]; }
    int iter = 1;
    int bestFDeltaIteration = 0;
    
    FNew = &x[nEdgesTimesNClasses];
    for(i=0; i<nVerticesTimesNClasses+twoTimesNClasses; i++){
        bestFDelta[i] = FNew[i];
        bestFDeltaInTermsOfCut[i] = FNew[i];
    }
    threshold(FNew, nVertices, nClasses, clusters, clusterSizes);
    tvMultiClass(wValTriu, irs, jcs, clusters, cutEdgesPL, degreePL, nVertices, nClasses, tvs);
    balancesMultiClass(clusters, vertexWeights, imbalance, CheegerFlag, nVertices, nClasses, eta, theta, thetaBar, balances);
    BCut = 0.0;
    for(l=0; l<nClasses; l++){
        BCut += tvs[l]/balances[l];
        //if(debugMode==1) mexPrintf("tv: %f\t balance: %f\t eta:%f\t theta: %f\t thetaBar: %f\n", tvs[l], balances[l], eta[l], theta[l], thetaBar[l]);
    }
    //mxAssert( fabs(BCut - initialBCut) <= 1e-8, "Initial BCut and BCut computed in the mex file are not same!");
    double lambdaNew;
    double bestBCut = BCut; // better bcut than the initial bcut and corresponds to an F whose lambda is also better than initial lambda.
    double bestBCut2 = BCut; // better bcut found in the way irrespective of the corresponding F achieving better lambda (this is only for tracking best cut across all iterations).
    double bestLambda = lambda; // better lambda with the condition that the corresponding F also yields a better cut.
    //bool exitFlag = false;
    bool printFlag = false;
    
    primalObj = 0.0;
    for(l=0; l<nClasses; l++){
        primalObj += x[nPrimalVars-nClasses*2+l];
        primalObj -= x[nPrimalVars-nClasses+l];
    }
    
    if(debugMode==1) mexPrintf("iter %d\t lambda: %f\t bcut: %f\t initialBCut:%f\n", 0, lambda, BCut, initialBCut);
    if(debugMode==1) mexPrintf("\n");
    //mexPrintf("iter \t primalObj\t dualObj\t Descent\t Simplex\t Lower bound\t Edges\t\t Dual-Infeas\t lambda\t\t BCut\n");
    if(debugMode==1){ 
        if(!CheegerFlag){
            mexPrintf("iter \t primalObj\t dualObj\t DualityGap\t Descent\t Simplex\t Lower bound\t Upper bound\t Edges\t\t Dual-Infeas\t\t  lambda\t   BCut\n");
        }
        else{
            mexPrintf("iter \t primalObj\t dualObj\t DualityGap\t Descent\t Simplex\t Lower bound\t Edges\t\t Dual-Infeas\t\t  lambda\t   BCut\n");
        }
        mexPrintf("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    }
    
    bool initialLambdaHigherThanBCut = false;
    if(lambda > BCut){
        initialLambdaHigherThanBCut = true;
    }
    bool betterLambdaFound = false;
    bool betterBCutFound = false;
    bool sameBCutFound = false;
    bool bCutChanged = false;
    int nSuccessiveDegenerateIters = 0;
    int dualInfeasIndex=0;
    // needed for simplex projection
    //vector<double> dVec(nClasses, 0), FRow(nClasses, 0), bVec(nClasses, 1), lVec(nClasses, 0), projFRow(nClasses, 0);
    //set<int> switch_set;
    while( iter <= maxIters ){
        
        /****************************************
         *              PRIMAL UPDATE
         ***************************************/
        // xNew = max(0, x - Tau*(A'*y+c));
        // x update for the alpha part (edges)
        primalVarCounter = 0;
        yEdgesCounter1 = twoTimesNClasses;
        if(!CheegerFlag){
            yEdgesCounter1 += nClasses;
        }
        yEdgesCounter2 = yEdgesCounter1 + nEdgesTimesNClasses;
        for(l=0; l<nClasses; l++){
            for(i=0; i<nEdges; i++){
                dummy = y[l]*wValTriu[i] - y[yEdgesCounter1] - y[yEdgesCounter2];
                dualInfeasibility[primalVarCounter] = dummy;
                xNew[primalVarCounter] = min(max(0.0, x[primalVarCounter] - Tau[primalVarCounter]*dummy), 1.0);
                primalVarCounter++; yEdgesCounter1++; yEdgesCounter2++;
            }
        }
        
        // x update for the F part
        // First compute A'y for the edges part; this is given by \sum_{i < j} y_{ij} - \sum_{i > j} y_{ji} with (i, j) \in E
        dummyCounter = nClasses;
        sgBalCounter = 0;
        counter1 = twoTimesNClasses;
        if(!CheegerFlag){
            counter1 += nClasses;
        }
        counter2 = counter1 + nEdgesTimesNClasses;
        for(l=0;l<nClasses;l++){
            // initialze
            for(i=0;i<nVertices;i++){ ATransYEdges1[i] = 0.0; }
            dummy = 0.0;
            counter = 0;
            for(j=0; j<nVertices; j++)
            {
                for(i=0; i<jcs[j+1]-jcs[j]; i++)
                {
                    dummy = y[counter1];
                    ATransYEdges1[j] -= dummy;
                    rowIndex = irs[counter];
                    ATransYEdges1[rowIndex] += dummy;
                    counter1++;
                    counter++;
                }
            }
            
            // initialze
            for(i=0;i<nVertices;i++){ ATransYEdges2[i] = 0.0; }
            dummy=0.0;
            counter = 0;
            for(j=0; j<nVertices; j++)
            {
                for(i=0; i<jcs[j+1]-jcs[j]; i++)
                {
                    dummy = y[counter2];
                    ATransYEdges2[j] -= dummy;
                    rowIndex = irs[counter];
                    ATransYEdges2[rowIndex] += dummy;
                    counter2++;
                    counter++;
                }
            }
            
            ySimplexCounter = nInEqCnstrs;
            for(i=0; i<nVertices; i++){
                dummy = -sgBal[sgBalCounter]*(lambdas[l]*y[l] + y[dummyCounter]) - cutEdgesPL[sgBalCounter]*y[l] - gamma*sgFMax[sgBalCounter]*y[l] + ATransYEdges1[i] - ATransYEdges2[i];
                if(simplexProjection==0){
                    dummy += y[ySimplexCounter];
                }
                //dummy = -sgBal[sgBalCounter]*(lambdas[l]*y[l] + y[dummyCounter]) - gamma*sgFMax[sgBalCounter]*y[l] + ATransYEdges1[i] - ATransYEdges2[i] + y[ySimplexCounter];
                if(!CheegerFlag){
                    dummy += sgBal[sgBalCounter]*y[dummyCounter+nClasses]; // upper bound constraints
                }
                dualInfeasibility[primalVarCounter] = dummy;
                if(simplexProjection==1)
                    xNew[primalVarCounter] = x[primalVarCounter] - Tau[primalVarCounter]*dummy;
                else
                    xNew[primalVarCounter] = min(max(0.0, x[primalVarCounter] - Tau[primalVarCounter]*dummy), 1.0);
                //xNew[primalVarCounter] = max(0.0, x[primalVarCounter] + Tau[primalVarCounter]*sgBal[sgBalCounter]*(lambdas[l]*y[l] + y[dummyCounter]) - Tau[primalVarCounter]*ATransYEdges1[i] + Tau[primalVarCounter]*ATransYEdges2[i] - Tau[primalVarCounter]*y[ySimplexCounter]);
                primalVarCounter++; sgBalCounter++; ySimplexCounter++;
            }
            dummyCounter++;
        }
        
        // x update for the deltaPlus variable
        for(l=0; l<nClasses; l++){
            dummy = -m*y[l] + cost;
            dualInfeasibility[primalVarCounter] = dummy;
            //xNew[primalVarCounter] = min(max(0.0, x[primalVarCounter] - Tau[primalVarCounter]*dummy), deltaBounds[l]); // c = <deltaPlus, 1> - <deltaMinus, 1>
            xNew[primalVarCounter] = max(0.0, x[primalVarCounter] - Tau[primalVarCounter]*dummy); // c = <deltaPlus, 1> - <deltaMinus, 1>
            primalVarCounter++;
        }
        
        // x update for the deltaMinus variable
        for(l=0; l<nClasses; l++){
            dummy = M*y[l] - cost;
            dualInfeasibility[primalVarCounter] = dummy;
            //xNew[primalVarCounter] = min(max(0.0, x[primalVarCounter] - Tau[primalVarCounter]*dummy ), deltaBounds[l+nClasses]); // c = <deltaPlus, 1> - <deltaMinus, 1>
            xNew[primalVarCounter] = max(0.0, x[primalVarCounter] - Tau[primalVarCounter]*dummy ); // c = <deltaPlus, 1> - <deltaMinus, 1>
            primalVarCounter++;
        }
        
        
        /****************************************
         *              DUAL UPDATE
         ***************************************/
        // yNew = max(0, y + Sigma*(A* (2*xNew-x) - b)) All are inequality constraints; simplex constraint is taken care via projection
        double stepSize = 1;
        // y update for the descent constraint
        for(i=0; i<nPrimalVars; i++){
            xBar[i] = (1+stepSize)*xNew[i] - x[i];
            //xBar[i] = xNew[i];
        }
        
        counter = 0;
        for(l=0; l<nClasses; l++){
            wTransAlpha[l] = 0;
            for(j=0; j<nEdges; j++){ wTransAlpha[l] += wValTriu[j]*xBar[counter]; counter++; }
        }
        
        sgBalCounter = 0;
        for(l=0; l<nClasses; l++){
            balanceLB[l] = 0; penalty[l] = 1;
            cutPL[l] = degreePL[l] - lambdaTimesEta[l];
            for(i=0; i<nVertices; i++){
                balanceLB[l] += sgBal[sgBalCounter]*xBar[counter];
                penalty[l] -= sgFMax[sgBalCounter]*xBar[counter];
                cutPL[l] -= cutEdgesPL[sgBalCounter]*xBar[counter];
                counter++; sgBalCounter++;
            }
        }        
        
        for(l=0; l<nClasses; l++){ mTimesDeltaPlus[l] = m*xBar[counter]; counter++; }
        for(l=0; l<nClasses; l++){ MTimesDeltaMinus[l] = M*xBar[counter]; counter++; }
        
        for(l=0; l<nClasses; l++){
            primalInfeasibility[l] = wTransAlpha[l] + cutPL[l] + gamma* penalty[l] - lambdas[l]*balanceLB[l] - mTimesDeltaPlus[l] + MTimesDeltaMinus[l];
            yNew[l] = min(cost/m, max(cost/M, y[l] + Sigma[l]*primalInfeasibility[l] )); // this is the only dual constraint!
            //yNew[l] = max(0.0, y[l] + Sigma[l]*primalInfeasibility[l] );
        }
        dualVarCounter = nClasses;
        
        // y update for lower bound constraint
        for(l=0; l<nClasses; l++){
            primalInfeasibility[dualVarCounter] = -balanceLB[l] + m - eta[l];
            yNew[dualVarCounter] = max(0.0, y[dualVarCounter] + Sigma[dualVarCounter]*primalInfeasibility[dualVarCounter]);
            dualVarCounter++;
        }
        
        if(!CheegerFlag){
            // y update for upper bound constraint
            for(l=0; l<nClasses; l++){
                primalInfeasibility[dualVarCounter] = balanceLB[l] + eta[l] - M;
                yNew[dualVarCounter] = max(0.0, y[dualVarCounter] + Sigma[dualVarCounter]*primalInfeasibility[dualVarCounter]);
                dualVarCounter++;
            }
        }
        
        // y update for the edges constraint
        xBarCounter = 0;
        FlCounter = nEdgesTimesNClasses;
        for(l=0; l<nClasses; l++){
            Fl = &xBar[FlCounter];
            counter = 0;
            for(j=0; j<nVertices; j++){
                for(i=0; i<jcs[j+1]-jcs[j]; i++){
                    rowIndex = irs[counter];
                    primalInfeasibility[dualVarCounter] = Fl[rowIndex] - Fl[j] - xBar[xBarCounter];
                    yNew[dualVarCounter] = max(0.0, y[dualVarCounter] + Sigma[dualVarCounter]*primalInfeasibility[dualVarCounter]);
                    dualVarCounter++;
                    counter++;
                    xBarCounter++;
                }
            }
            FlCounter += nVertices;
        }
        
        xBarCounter = 0;
        FlCounter = nEdgesTimesNClasses;
        for(l=0; l<nClasses; l++){
            //Fl = &xBar[nEdgesTimesNClasses+nVertices*l];
            Fl = &xBar[FlCounter];
            counter = 0;
            for(j=0; j<nVertices; j++){
                for(i=0; i<jcs[j+1]-jcs[j]; i++){
                    rowIndex = irs[counter];
                    primalInfeasibility[dualVarCounter] = -Fl[rowIndex] + Fl[j] - xBar[xBarCounter];
                    yNew[dualVarCounter] = max(0.0, y[dualVarCounter] + Sigma[dualVarCounter]*primalInfeasibility[dualVarCounter]);
                    dualVarCounter++;
                    counter++;
                    xBarCounter++;
                }
            }
            FlCounter += nVertices;
        }
                
        // y update for the simplex constraint
        if(simplexProjection==0){
            FBar = &xBar[nEdgesTimesNClasses];
            for(i=0; i<nVertices; i++){
                rowSum = 0.0;
                counter = i;
                for(l=0; l<nClasses; l++){
                    rowSum += FBar[counter];
                    counter += nVertices;
                }                
                primalInfeasibility[dualVarCounter] = rowSum - 1;
                yNew[dualVarCounter] = y[dualVarCounter] + Sigma[dualVarCounter]*primalInfeasibility[dualVarCounter]; // this is equality constraint and hence y_i \in \R
                dualVarCounter++;
            }
        }        
        
        primalObj = 0.0;
        for(l=0; l<nClasses; l++){
            primalObj += xNew[nPrimalVars-nClasses*2+l];
            primalObj -= xNew[nPrimalVars-nClasses+l];
        }
        primalObj *= cost;
        
        // update primal infeasibilities
        for(i=0; i<nDualVars; i++){
            primalInfeasibility[i] = 0.5*(primalInfeasibility[i] + primalInfeasibilityOld[i]);
            primalInfeasibilityOld[i] = primalInfeasibility[i];
        }
        
        FNew = &xNew[nEdgesTimesNClasses];
        tvMultiClass(wValTriu, irs, jcs, FNew, cutEdgesPL, degreePL, nVertices, nClasses, tvs);
        balancesMultiClass(FNew, vertexWeights, imbalance, CheegerFlag, nVertices, nClasses, eta, theta, thetaBar, balances);
        lambdaNew = 0.0;        
        for(l=0; l<nClasses; l++){
            lambdaNew += tvs[l]/balances[l];            
        }        
        
        // Find clustering via rounding and compute the balanced cut        
        threshold(FNew, nVertices, nClasses, clusters, clusterSizes);
        tvMultiClass(wValTriu, irs, jcs, clusters, cutEdgesPL, degreePL, nVertices, nClasses, tvs);
        balancesMultiClass(clusters, vertexWeights, imbalance, CheegerFlag, nVertices, nClasses, eta, theta, thetaBar, balances);
        BCut = 0.0;
        for(l=0; l<nClasses; l++){
            BCut += tvs[l]/balances[l];
        }
        
        if(!bCutChanged && BCut!=bestBCut){
            bCutChanged = true;
        }
        
        if(BCut < bestBCut2){
            for(i=0; i<nVerticesTimesNClasses; i++){
                bestFDeltaInTermsOfCut[i] = clusters[i];
            }
            for(i=nVerticesTimesNClasses; i<nVerticesTimesNClasses+twoTimesNClasses; i++){
                bestFDeltaInTermsOfCut[i] = xNew[nEdgesTimesNClasses+i];
            }
            bestBCut2 = BCut;
        }
        printFlag = false;
        if(lambdaNew < lambda){ // We have descent in the continuous problem
            betterLambdaFound = true;
            if(isnan(BCut)){ // The continuous solution has degenerated
                degenerateFlag[0] = 1;  //true;
                if(debugMode==1) mexPrintf("iter: %d degenerated\n", iter);
                for(i=0; i<nVerticesTimesNClasses; i++){
                    degeneratedF[i] = FNew[i];
                }
            }
            else if(BCut <= bestBCut){ // The continuous solution also preserves monotonicity in the cut                
                if(BCut < bestBCut){ // The continuous solution actually yields a better cut (or the initial lambda > initial Cut), so display it and exit in the next round!
                    betterBCutFound = true;
                    //exitFlag = true;
                    if(iter%nPrintIters!=0){
                        printFlag = true;
                    }
                    // to exit early if we get a better cut.
                    // helps findBestGraphBuildingParams
                    //if (BCut<initialBCut) {
                    //    printFlag = false;
                    //}
                    bestLambda = lambdaNew;
                    bestFDeltaIteration = iter;
                    bestBCut = BCut;
                    for(i=0; i<nVerticesTimesNClasses+twoTimesNClasses; i++){
                        bestFDelta[i] = FNew[i];
                    }
                }
                else{
                    sameBCutFound = true;
                    if(lambdaNew < bestLambda){
                        bestFDeltaIteration = iter;
                        bestLambda = lambdaNew;
                        for(i=0; i<nVerticesTimesNClasses+twoTimesNClasses; i++){
                            bestFDelta[i] = FNew[i];
                        }
                    }
                }
                
            }
            else{ // The continuous solution yielded a higher cut.
                if(eta[0]!=0){ // Check if the pseudo vertices are cut-out!
                    counter = 0;
                    for(l=0; l<nClasses; l++){
                        if(clusterSizes[l] == 0){
                            degenerateFlag[0] = 1;  //true;
                            if(debugMode==1) mexPrintf("iter: %d Pseudo vertex is cut-out!\n", iter);
                            for(i=0; i<nVerticesTimesNClasses; i++){
                                degeneratedF[i] = FNew[i];
                            }
                            break;
                        }
                    }
                }
            }
        }
                
        if(degenerateFlag[0]==0){
            if(isnan(BCut)){
                nSuccessiveDegenerateIters += 1;
            }
            else{
                nSuccessiveDegenerateIters = 0;
            }
            
            if(nSuccessiveDegenerateIters>25){
                degenerateFlag[0] = 1;  //true;
                if(debugMode==1) mexPrintf("iter: %d degenerated for the last 25 iterations.. so exiting although we (might) have not seen descent in lambda...\n", iter);
                for(i=0; i<nVerticesTimesNClasses; i++){
                    degeneratedF[i] = FNew[i];
                }
            }
        }
        
        if(iter > maxIters/2 && eta[0]!=0 && degenerateFlag[0]==0 && BCut > bestBCut){ // Check if the pseudo vertices are cut-out!
            counter = 0;
            for(l=0; l<nClasses; l++){
                if(clusterSizes[l] == 0){
                    degenerateFlag[0] = 1;  //true;
                    if(debugMode==1) mexPrintf("iter: %d Pseudo vertex is cut-out! exiting although no descent in lambda\n", iter);
                    for(i=0; i<nVerticesTimesNClasses; i++){
                        degeneratedF[i] = FNew[i];
                    }
                    break;
                }
            }
        }
        
        
        //if(iter > maxIters/2 && betterLambdaFound && !(betterBCutFound || sameBCutFound)){ // The starting F is not sparse (convex splitting occured) and yielding highly imbalanced clusters (one stage before degeneration)
            //if(debugMode==1) mexPrintf("Could not improve or maintain the initial BCut (while getting strict descent in lambda) for too long... exiting\n");
            //break;
        //}
        
        //if(iter > maxIters/2 && !bCutChanged){
        //    if(debugMode==1) mexPrintf("Intial BCut stays same in spite of descent in lambda for too long... exiting so that the next inner problem has better approximation\n");
        //break;
        //}
                
        if(printFlag || degenerateFlag[0]==1 || iter%nPrintIters==0 || iter==maxIters){        
            // Compute dual objective: -yNew'*b
            dualObj = 0.0;
            if(gamma > 0){
                for(l=0; l<nClasses; l++){
                    dualObj += yNew[l]*gamma; // b is -gamma
                }
            }
            
            // b of descent constraint
            for(counter=0; counter<nClasses; counter++){
                dualObj += yNew[counter]*(degreePL[counter] - lambdaTimesEta[counter] ); // b is -degreePL+lambdaTimesEta
            }            
            
            // b of lower bound constraint
            //counter = nClasses;
            for(l=0; l<nClasses; l++){
                dualObj += yNew[counter++]*(m-eta[l]); // b is -m + eta[l]
            }
            if(!CheegerFlag){
                // b of upper bound constraint
                //counter = nClasses;
                for(l=0; l<nClasses; l++){
                    dualObj += yNew[counter++]*(eta[l]-M); // b is M - eta[l]
                }
            }
            // b of simplex constraints
            if(simplexProjection==0){
                counter = nInEqCnstrs;
                for(i=0; i<nVertices; i++){
                    dualObj += -yNew[counter++]; // b is 1
                }
            }
//             maxDualInfeasibility = fabs(min(0.0, dualInfeasibility[0]));
//             for(i=1; i<nEdgesTimesNClasses; i++){
//                 dummy = fabs(min(0.0, dualInfeasibility[i]));
//                 if(dummy > maxDualInfeasibility){
//                     maxDualInfeasibility = dummy;
//                     dualInfeasIndex =i;
//                 }
//             }
//	if simplexProjection
//            dualObj = dualObj+sum(min(0, dualInfeasibility(1:nEdgesTimesNClasses))); % constraint on alpha: [0,1]
//	    dummy = reshape(dualInfeasibility(nEdgesTimesNClasses+1:end-2*nClasses), nVertices, nClasses);
//            dualObj = dualObj+sum(min(dummy,[],2)); % constraint on F: row-wise simplex
//        else
//            dualObj = dualObj+sum(min(0,dualInfeasibility(1:end-2*nClasses)));
//        end
            
            // bug fix: compute the dual objective under the constraints: F and alpha in [0,1].
            for(i=0; i<nEdgesTimesNClasses; i++){
                dualObj += min(0.0, dualInfeasibility[i]);
            }
            
            // dual objective is different for simplexProjection.
            if(simplexProjection==1){
                for(i=nEdgesTimesNClasses; i<nEdgesTimesNClasses+nVertices; i++){
                    // find the min in each row!
                    counter = nVertices;
                    dummy = dualInfeasibility[i];
                    for(j=1;j<nClasses;j++){
                        if(dualInfeasibility[i+counter] < dummy)
                            dummy = dualInfeasibility[i+counter];
                        counter += nVertices;
                    }
                    dualObj += dummy;
                }
            }
            else{
                for(i=nEdgesTimesNClasses; i<nEdgesTimesNClasses+nVerticesTimesNClasses; i++){
                    dualObj += min(0.0, dualInfeasibility[i]);
                }
            }
            //if(deltaBounds[0] < 10000){
            //    for(i=nEdgesTimesNClasses+nVerticesTimesNClasses;i<nPrimalVars; i++){
            //        dualObj += deltaBounds[i-nEdgesTimesNClasses+nVerticesTimesNClasses]*min(0.0, dualInfeasibility[i]);
            //    }
            //}
            
            // no objective term for the delta variables: there is a corresponding dual feasibility constraint.
            // here we are checking for the maximum infeasiblity among the delta variables.
            i = nEdgesTimesNClasses+nVerticesTimesNClasses;
            maxDualInfeasibility = fabs(min(0.0, dualInfeasibility[i]));
            dualInfeasIndex = i+1;
            for(; i<nPrimalVars; i++){
                dummy = fabs(min(0.0, dualInfeasibility[i]));
                if(dummy > maxDualInfeasibility){
                    maxDualInfeasibility = dummy;
                    dualInfeasIndex =i+1; // we display the matlab index of this variable.
                }
            }
            
            counter = 0;
            maxDescentCnstrInfeasibility = max(0.0, primalInfeasibility[0]);
            while(counter < nClasses){
                dummy = max(0.0, primalInfeasibility[++counter]);
                if(dummy > maxDescentCnstrInfeasibility){
                    maxDescentCnstrInfeasibility = dummy;
                }
            }
            maxPrimalInfeasibility = maxDescentCnstrInfeasibility;
            
            maxLBCnstrInfeasibility = max(0.0, primalInfeasibility[counter]);
            while(counter < twoTimesNClasses){
                dummy = max(0.0, primalInfeasibility[++counter]);
                if(dummy > maxLBCnstrInfeasibility){
                    maxLBCnstrInfeasibility = dummy;
                }
            }
            if(maxPrimalInfeasibility < maxLBCnstrInfeasibility ){
                maxPrimalInfeasibility = maxLBCnstrInfeasibility;
            }
            
            if(!CheegerFlag){
                maxUBCnstrInfeasibility = max(0.0, primalInfeasibility[counter]);
                while(counter < twoTimesNClasses + nClasses){
                    dummy = max(0.0, primalInfeasibility[++counter]);
                    if(dummy > maxUBCnstrInfeasibility){
                        maxUBCnstrInfeasibility = dummy;
                    }
                }
                if(maxPrimalInfeasibility < maxUBCnstrInfeasibility ){
                    maxPrimalInfeasibility = maxUBCnstrInfeasibility;
                }
            }
            
            maxEdgeCnstrInfeasibility = max(0.0, primalInfeasibility[counter]);
            while(counter < nInEqCnstrs){
                dummy = max(0.0, primalInfeasibility[++counter]);
                if(dummy > maxEdgeCnstrInfeasibility){
                    maxEdgeCnstrInfeasibility = dummy;
                }
            }
            if(maxPrimalInfeasibility < maxEdgeCnstrInfeasibility ){
                maxPrimalInfeasibility = maxEdgeCnstrInfeasibility;
            }
            
            if(simplexProjection==1){
                maxSimplexCnstrInfeasibility = 0;
            }
            else{
                maxSimplexCnstrInfeasibility = fabs(primalInfeasibility[counter]);
                while(counter < nDualVars){
                    dummy = fabs(primalInfeasibility[++counter]);
                    if(dummy > maxSimplexCnstrInfeasibility){
                        maxSimplexCnstrInfeasibility = dummy;
                    }
                }
            }
            if(maxPrimalInfeasibility < maxSimplexCnstrInfeasibility ){
                maxPrimalInfeasibility = maxSimplexCnstrInfeasibility;
            }
            
            if(debugMode==1){
                if(!CheegerFlag){
                    mexPrintf("%d\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.8f\t %1.8f\n", iter, primalObj, dualObj, primalObj-dualObj, maxDescentCnstrInfeasibility, maxSimplexCnstrInfeasibility, maxLBCnstrInfeasibility, maxUBCnstrInfeasibility, maxEdgeCnstrInfeasibility, maxDualInfeasibility, lambdaNew, BCut);
                }
                else{
                    mexPrintf("%d\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t (%d) %1.8f\t %1.8f\n", iter, primalObj/cost, dualObj/cost, (primalObj-dualObj)/cost, maxDescentCnstrInfeasibility, maxSimplexCnstrInfeasibility, maxLBCnstrInfeasibility, maxEdgeCnstrInfeasibility, maxDualInfeasibility, dualInfeasIndex,lambdaNew, BCut);
                }
            }
            
            if(degenerateFlag[0]){ // Continuous solution degenerated
                break;
            }
            
            if(betterBCutFound && !printFlag){ // We got better cut
                break;
//                 if(bestLambda > bestBCut){
//                     nVerticesMoved = 0;
//                     for(l=0; l<nClasses; l++){
//                         for(i=0; i<nVertices; i++){
//                             if(initialClusters[counter] != best
//                         }
//                     }
//                 }
//                 else{
//                     break;
//                 }
            }
            
            if(!betterBCutFound && !printFlag){
                nPrintIters *= 2;
            }
                        
            if(maxDualInfeasibility <= 1e-3 && maxPrimalInfeasibility <= 1e-3 && primalObj ==0 && dualObj > -0.5) {
                convergenceFlag[0] = 1;
                if(debugMode==1) mexPrintf("Solution converged upto numerical accuracy of dualObj\t");
                break;
            }       
        }
        if (iter > 800 && bestLambda >= lambda){
        //if (iter%nPrintIters==0 && bestLambda >= lambda){
            if(debugMode==1) mexPrintf("No descent in Lambda.. exiting\n"); // This is needed because PDHG is not a monotonic descent method and so we don't want to exhaust all iterations just to find out that optimal is zero.
            break;
        }
        
        dummyPointer = xNew; xNew = x; x = dummyPointer;
        dummyPointer = yNew; yNew = y; y = dummyPointer;
        
        iter++;
        
    }
    
    if(debugMode==1){
        mexPrintf("Returning with lambda: %1.8f\t BCut: %1.8f (bestFDeltaIteration = %d)\n", bestLambda, bestBCut, bestFDeltaIteration);
        mexPrintf("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        if(!CheegerFlag){
            mexPrintf("iter \t primalObj\t dualObj\t DualityGap\t Descent\t Simplex\t Lower bound\t Upper bound\t Edges\t\t Dual-Infeas\t\t  lambda\t   BCut\n");
        }
        else{
            mexPrintf("iter \t primalObj\t dualObj\t DualityGap\t Descent\t Simplex\t Lower bound\t Edges\t\t Dual-Infeas\t\t  lambda\t   BCut\n");
        }
    }
    
    bestLambdaPtr[0] = bestLambda;
    nIters[0] = iter-1;
    FNew = &xNew[nEdgesTimesNClasses];
    for(i=0; i<nVerticesTimesNClasses+twoTimesNClasses; i++){ finalFDelta[i] = FNew[i]; }
    for(i=0; i<nDualVars; i++){ yPDHG[i] = yNew[i]; }
    delete x; delete y;
    delete xNew; delete yNew; delete xBar;
    delete projFNew; delete balanceLB; delete dualInfeasibility; delete primalInfeasibility; delete primalInfeasibilityOld;
    delete penalty; delete cutPL;
    delete wTransAlpha; delete ATransYEdges1; delete ATransYEdges2; delete clusters; delete clusterSizes;
    delete mTimesDeltaPlus; delete MTimesDeltaMinus;
    delete balances; delete tvs; delete lambdaTimesEta; delete feasibleF;
    
}

void threshold(double *F, int nVertices, int nClasses, double *clusters, int *clusterSizes){
    
    int i, l, maxIx, counter;
    double tempMax, dummy;
    
    for(l=0; l<nClasses; l++){
        clusterSizes[l] = 0;
    }
    
    for(i=0; i<nVertices; i++){
        maxIx = i; // first column is the maximum by default for any vertex i.
        counter = i; // counter goes over the ith row
        tempMax = F[counter];
        for(l=1;l<nClasses;l++){
            counter += nVertices; // start from the second element in the ith row
            dummy = F[counter];
            if(dummy > tempMax)
            {
                clusters[maxIx] = 0;
                maxIx = counter; // the current column has a higher value than all the previous columns
                tempMax = dummy;
            }
            else{
                clusters[counter] = 0; // the current column does not achieve maximum for the ith row.
            }
        }
        clusters[maxIx] = 1; // by now maxIx is the linear index for the column achieving maximum in ith row.
        clusterSizes[maxIx/nVertices] += 1; // the column index is the quotient
    }
    
    
}

void tvMultiClass(double *wValTriu, mwIndex* irs, mwIndex* jcs, double *F, double *cutEdgesPL, double *degreePL, int nVertices, int nClasses, double *tvMultiClass){
    
    int i, j, l, counter, len, rowIndex, columnCounter = 0;
    double *Fl;
    
    for(l=0; l<nClasses; l++){
        tvMultiClass[l] = 0.0;
        counter = 0;
        
        Fl = &F[columnCounter];
        for(j=0; j<nVertices; j++){
            len = jcs[j+1]-jcs[j];
            for(i=0; i<len; i++){
                rowIndex = irs[counter];
                tvMultiClass[l] += wValTriu[counter]*fabs(Fl[rowIndex]-Fl[j]);
                counter++;
            }
        }
        columnCounter += nVertices;
    }
    
    double *cutPL = new double[nClasses];
    counter = 0;
    for(l=0; l<nClasses; l++){
        cutPL[l] = degreePL[l];
        for(i=0; i<nVertices; i++){
            cutPL[l] -= cutEdgesPL[counter]*F[counter];
            counter++;
        }
    }
    
    for(l=0; l<nClasses; l++){
        tvMultiClass[l] += cutPL[l];
    }
    delete cutPL;
}

bool comparePair(std::pair<double, double> a, std::pair<double, double> b){
    return a.first < b. first;
}

void balancesMultiClass(double *F, double *vertexWeights, int beta, bool CheegerFlag, int nVertices, int nClasses, double* eta, double *theta, double* thetaBar, double *balances){
    
    int i, l, counter = 0;
    
    if(!CheegerFlag){
        for(l=0; l<nClasses; l++){
            balances[l] = eta[l];
            for(i=0; i<nVertices; i++){ balances[l] += vertexWeights[i]*F[counter]; counter++; }
        }
        return;
    }
    
    double volV = 0.0;
    for(i=0; i<nVertices; i++){ volV += vertexWeights[i]; }
    
    double cumSortedVolume = 0.0;
    double medianThreshold; // (beta*volV - thetaBar)/(beta+1.0);
    //double dummy = beta*volV /(beta+1.0);
    
    //double *sortedFl = new double[nVertices]; // Fl: l^th column of F
    std::pair<double, double>* sortedFlVertexWeightsPair = new std::pair<double, double>[nVertices];
    
    for(l=0;l<nClasses;l++){
        
        for(i=0; i<nVertices; i++){ sortedFlVertexWeightsPair[i].first = F[counter]; sortedFlVertexWeightsPair[i].second = vertexWeights[i]; counter++; }
        std::sort(sortedFlVertexWeightsPair, sortedFlVertexWeightsPair+nVertices, comparePair); // sorts the array in the range: [start, last);
        
        //dummy -= thetaBar[l]/(beta+1.0);
        medianThreshold = (beta*volV + theta[l] - thetaBar[l])/(beta+1.0);
        balances[l] = eta[l];
        i = 0;
        cumSortedVolume = sortedFlVertexWeightsPair[i].second;
        // ixMinus
        while(i<nVertices && cumSortedVolume <= medianThreshold){
            balances[l] -= sortedFlVertexWeightsPair[i].second * sortedFlVertexWeightsPair[i].first;
            i++;
            if(i<nVertices)
                cumSortedVolume += sortedFlVertexWeightsPair[i].second;
        }
        // ixZero
        while(i< nVertices && cumSortedVolume < medianThreshold + sortedFlVertexWeightsPair[i].second){
            balances[l] += ((beta+1)*cumSortedVolume - beta*volV - sortedFlVertexWeightsPair[i].second + thetaBar[l] - theta[l]) * sortedFlVertexWeightsPair[i].first;
            i++;
            if(i<nVertices)
                cumSortedVolume += sortedFlVertexWeightsPair[i].second;
        }
        // ixPlus
        while(i<nVertices){
            balances[l] += beta*sortedFlVertexWeightsPair[i].second*sortedFlVertexWeightsPair[i].first;
            i++;
        }
    }
    delete sortedFlVertexWeightsPair;
}