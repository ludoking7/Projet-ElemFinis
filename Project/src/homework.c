#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)



double *femElasticitySolve(femProblem *theProblem, femSolverType typeSolver)
{
    femFullSystem *theFullSystem;
    femBandSystem *theBandSystem;

    double **A;
    double *B;
    int size;

    if (theProblem->theSolver->type == FEM_BAND){
        theBandSystem = (femBandSystem*) theProblem->system;
        A = theBandSystem->A;
        B = theBandSystem->B;
        size = theBandSystem->size;
    }
    else if(theProblem->theSolver->type == FEM_BAND){
        theFullSystem = (femBandSystem*) theProblem->system;
        A = theFullSystem->A;
        B = theFullSystem->B;
        size = theFullSystem->size;
    }
    else {
        Error("Unexpected solver type");
    }

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
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            // renumber
            int renumberIndex = theMesh->elem[iElem*nLocal+j];
            map[j]  = theMesh->number[renumberIndex];

            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[renumberIndex];
            y[j]    = theNodes->Y[renumberIndex];}

        
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
                fprintf(stdout, "Axisymétrique\n");
                for (i = 0; i < theSpace->n; i++) { 
                    for(j = 0; j < theSpace->n; j++) {
                        if(typeSolver == FEM_BAND){
                            if(mapX[i] <= mapX[j]){
                                A[mapX[i]][mapX[j]] += (dphidx[i] * a * r * dphidx[j] +
                                                        dphidy[i] * c * r * dphidy[j] +
                                                        phi[i] * ((b * dphidx[j]) + (a * phi[j] / r)) +
                                                        dphidx[i] * b * phi[j]) * jac * weight;
                            }
                            else{
                                A[mapX[i]][mapX[j]] = 0.0;
                            }

                            if(mapX[i] <= mapY[j]){
                                A[mapX[i]][mapY[j]] += (dphidx[i] * b * r * dphidy[j] +
                                                        dphidy[i] * c * r * dphidx[j] +
                                                        phi[i] * b * dphidy[j]) * jac * weight;
                            }
                            else{
                                A[mapX[i]][mapY[j]] = 0.0;
                            }

                            if(mapY[i] <= mapX[j]){
                                A[mapY[i]][mapX[j]] += (dphidy[i] * b * r * dphidx[j] +
                                                        dphidx[i] * c * r * dphidy[j] +
                                                        dphidy[i] * b * phi[j]) * jac * weight;
                            }
                            else{
                                A[mapY[i]][mapX[j]] = 0.0;
                            }

                            if (mapY[i] <= mapY[j]){(
                                    A[mapY[i]][mapY[j]] +=  dphidy[i] * a * r * dphidy[j] +
                                                            dphidx[i] * c * r * dphidx[j]) * jac * weight;
                            }
                            else{
                                A[mapY[i]][mapY[j]] = 0.0;
                            }

                        }
                        else {
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
                }
                
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * g * r * rho * jac * weight; } 

            }


            // Partie plan :
            else if (theProblem -> planarStrainStress  != AXISYM){
                //fprintf(stdout, "Plan\n");
                for (i = 0; i < theSpace->n; i++) { 
                    for(j = 0; j < theSpace->n; j++) {
                        if(typeSolver == FEM_BAND){
                            if(mapX[i] <= mapX[j]){
                                A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                                                        dphidy[i] * c * dphidy[j]) * jac * weight;
                            }
                            else{
                                A[mapX[i]][mapX[j]] = 0.0;
                            }

                            if(mapX[i] <= mapY[j]){
                                A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                                                        dphidy[i] * c * dphidx[j]) * jac * weight;
                            }
                            else{
                                A[mapX[i]][mapY[j]] = 0.0;
                            }

                            if(mapY[i] <= mapX[j]){
                                A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                                                        dphidx[i] * c * dphidy[j]) * jac * weight;
                            }
                            else{
                                A[mapY[i]][mapX[j]] = 0.0;
                            }

                            if (mapY[i] <= mapY[j]){
                                A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                                                        dphidx[i] * c * dphidx[j]) * jac * weight;
                            }
                            else{
                                A[mapY[i]][mapY[j]] = 0.0;
                            }
                        }
                    }
                }
                
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight; }
            }

            else{
                fprintf(stdout, "Erreur : Wrong planarStrainStress value\n");
                exit(1);
            }
        }
    
    }


    // Neumann
    for (j = 0; j < theProblem->nBoundaryConditions; j++) {
        femBoundaryCondition *theCondition = theProblem->conditions[j];
        femDomain *theDomain = theCondition->domain;
        if (theCondition -> type == NEUMANN_X || theCondition -> type == NEUMANN_Y){
            for (iEdge = 0; iEdge < theDomain -> nElem; iEdge++) {
                int iSegment = theDomain->elem[iEdge]; 
                int iNode1 = theEdges->elem[2*iSegment];
                int iNode2 = theEdges->elem[2*iSegment+1];

                int iNodeRenum1 = theNodes->number[iNode1];
                int iNodeRenum2 = theNodes->number[iNode2];

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


                        for (int i = 0; i < 2; i++){
                            iNode = iNodes[i];
                            weight = theRule->weight[i];

                            for (int k = 0; k < 2; k++){
                                xsik = xsi[k];
                                phi[0] = (1-xsik)/2;
                                phi[1] = (1+xsik)/2;
                                r = ((1-xsik)/2)*x1 + ((1+xsik)/2)*x2;

                                integral = sqrt(pow(x2-x1,2) + pow(y2-y1,2))/2 * r * weight * phi[i];
                                
                                if (theCondition -> type == NEUMANN_X){
                                    B[2*iNodeRenum1] += integral*theCondition->value;
                                }
                                else if (theCondition -> type == NEUMANN_Y) {
                                    B[2*iNodeRenum2+1] += integral*theCondition->value; // weight = 1 donc pas esoin et jac est deja h/2 car dx/dxsi = h/2 et c'est en 1D
                                }
                                else {
                                    Error("The boundary condition is not a Neumann X, Y condition");
                                    }
                            }
                        }
                }
                // Neumann, cas plan :
                else if(theProblem -> planarStrainStress  != AXISYM){
                    integral = sqrt(pow(x2-x1,2) + pow(y2-y1,2))/2;

                    if (theCondition -> type == NEUMANN_X){
                        B[2*iNodeRenum1] += integral*theCondition->value;
                        B[2*iNodeRenum2] += integral*theCondition->value;
                        }
                    else if (theCondition -> type == NEUMANN_Y){
                        B[2*iNodeRenum1+1] += integral*theCondition->value;
                        B[2*iNodeRenum2+1] += integral*theCondition->value; // weight = 1 donc pas esoin et jac est deja h/2 car dx/dxsi = h/2 et c'est en 1D
                        
                    }
                    else {
                        Error("The boundary condition is not a Neumann X, Y condition");
                    }
                }
                else{
                    fprintf(stdout, "Erreur : Wrong planarStrainStress value\n");
                    exit(1);
                }

            }
        }
    }   

    if(theProblem->theSolver->type == FEM_BAND){
        int *theConstrainedNodes = theProblem->constrainedNodes;
        for (int i=0; i < theBandSystem->size; i++) {
            if (theConstrainedNodes[i] != -1) {
                double value = theProblem->conditions[theConstrainedNodes[i]]->value;
                femFullSystemConstrain(theBandSystem,i,value);
    }
    

            // addboundarycondition crée un domaine avec tous les points où il y a une condition dessus
            // grace à cette liste on sait prendre les edges qui nous intéressent et ensuite on prendre
            // les noeuds de ces arrêtes qui nous intéressent et on applique l'intégrale et on le remplace
            // dans le second membre, attention u = x et v = y

             }

    double* solution = femBandSystemEliminate(theBandSystem);
    femFullSystemPrint(theBandSystem);

    femDenumber(theNodes, theBandSystem->size, solution);

    return solution;
    }

    else{
        int *theConstrainedNodes = theProblem->constrainedNodes;
        for (int i=0; i < theFullSystem->size; i++) {
            if (theConstrainedNodes[i] != -1) {
                double value = theProblem->conditions[theConstrainedNodes[i]]->value;
                femFullSystemConstrain(theFullSystem, i, value);
            }
        }

        double* solution = femFullSystemEliminate(theFullSystem);

        return solution;

    }




    
}

