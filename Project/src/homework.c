#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)



double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
            
        
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            double r = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
                r += x[i]*phi[i];
                }
                double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }      

            // Partie axisymétrique :  
            
            if (theProblem -> planarStrainStress  == AXISYM){

                for (i = 0; i < theSpace->n; i++) { 
                    for(j = 0; j < theSpace->n; j++) {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * r * dphidx[j] +
                                                dphidy[i] * c * r * dphidy[j] +
                                                phi[i] * ((b * dphidx[j]) + (a * phi[j] / r)) +
                                                dphidx[i] * b * phi[j]) * jac * weight;
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * r * dphidy[j] +
                                                dphidy[i] * c * r * dphidx[j] +
                                                phi[i] * b * dphidy[j]) * jac * weight;
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * r * dphidx[j] +
                                                dphidx[i] * c * r * dphidy[j] +
                                                dphidy[i] * b * phi[j]) * jac * weight;
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * r * dphidy[j] +
                                                dphidx[i] * c * r * dphidx[j]) * jac * weight;
                    }
                }
                
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * g * r * rho * jac * weight; } 

            }


            // Partie plan :
            else {
                for (i = 0; i < theSpace->n; i++) { 
                    for(j = 0; j < theSpace->n; j++) {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                                dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                                dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                                dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                                dphidx[i] * c * dphidx[j]) * jac * weight; }}
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
    
    }

    // Je dois créer un nouveau the Rule d'intégration pour le cas axisymmétrique des conditions de Neumann
    // Sur les segments, pour intégrer on doit multiplier r et les fonctions de forme. On va intégrer pour les 
    // points -sqrt(3)/3 et sqrt(3)/3. On va donc créer un nouveau theRule avec ces points et les poids qui vont avec. 
    // Ces poids valent 1. Le jacobien est déja géré par le sqrt(pow(x2-x1,2) + pow(y2-y1,2))/2.
    //


    // Neumann
    for (j = 0; j < theProblem->nBoundaryConditions; j++) {
        femBoundaryCondition *theCondition = theProblem->conditions[j];
        femDomain *theDomain = theCondition->domain;
        if (theCondition -> type == NEUMANN_X || theCondition -> type == NEUMANN_Y){
            for (iEdge = 0; iEdge < theDomain -> nElem; iEdge++) { 
                
                int iSegment = theDomain->elem[iEdge]; 
                int iNode1 = theEdges->elem[2*iSegment];
                int iNode2 = theEdges->elem[2*iSegment+1];

                int iNodes[2] = {iNode1, iNode2};
                int iNode;

                double x1 = theNodes->X[iNode1];
                double x2 = theNodes->X[iNode2];

                double y1 = theNodes->Y[iNode1];
                double y2 = theNodes->Y[iNode2];
                
                double integral;
                double phi1, phi2;
                double phi[2];

                double r;
                double weight;

                double xsi[2];
                double xsik;

                // Neumann, cas axisymétrique :
                if (theProblem -> planarStrainStress  == AXISYM){
                    // on va devoir réitérer pour chaque point d'intégration afin de calculer la valeur de r
                    // Pour un triangles, les fonctions de forme associée au segment frontière sont (1-xsi)/2 et (1+xsi)/2
                        xsi[0] = -sqrt(3)/3;
                        xsi[1] = sqrt(3)/3;
                        weight = 1.0;

                        for (int i = 0; i < 2; i++){
                            iNode = iNodes[i];

                            for (int k = 0; k < 2; k++){
                                xsik = xsi[k];
                                phi[0] = (1-xsik)/2;
                                phi[1] = (1+xsik)/2;
                                r = ((1-xsik)/2)*x1 + ((1+xsik)/2)*x2;

                                integral = sqrt(pow(x2-x1,2) + pow(y2-y1,2))/2 * r * weight * phi[i];
                                
                                if (theCondition -> type == NEUMANN_X){
                                    B[2*iNode] += integral*theCondition->value;
                                }
                                else if (theCondition -> type == NEUMANN_Y) {
                                    B[2*iNode+1] += integral*theCondition->value; // weight = 1 donc pas esoin et jac est deja h/2 car dx/dxsi = h/2 et c'est en 1D 
                                }
                                else {
                                    Error("The boundary condition is not a Neumann X, Y condition");
                                    }
                            }
                        }
                }
                // Neumann, cas plan :
                else{
                    integral = sqrt(pow(x2-x1,2) + pow(y2-y1,2))/2;
                    }  

                if (theCondition -> type == NEUMANN_X){
                    B[2*iNode1] += integral*theCondition->value;
                    B[2*iNode2] += integral*theCondition->value; 
                }
                else if (theCondition -> type == NEUMANN_Y){
                    B[2*iNode1+1] += integral*theCondition->value;
                    B[2*iNode2+1] += integral*theCondition->value; // weight = 1 donc pas esoin et jac est deja h/2 car dx/dxsi = h/2 et c'est en 1D
                    
                }
                else {
                    Error("The boundary condition is not a Neumann X, Y condition");
                }
            }
        }
    }   
    
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value);
            // addboundarycondition crée un domaine avec tous les points où il y a une condition dessus
            // grace à cette liste on sait prendre les edges qui nous intéressent et ensuite on prendre
            // les noeuds de ces arrêtes qui nous intéressent et on applique l'intégrale et on le remplace
            // dans le second membre, attention u = x et v = y

             }
    }
                            
    return femFullSystemEliminate(theSystem);
}

