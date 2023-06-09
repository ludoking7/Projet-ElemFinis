/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"

int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n\n\n");


//
//  -1- Lecture des donnees
//
    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3; 
    double g   = 9.81;
    double q = 100;
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../../connexion/mesh.txt");
    
    femProblem* theProblem = femElasticityRead(theGeometry,"../../connexion/problem.txt");
    double *theSoluce = theProblem->system->B;
    int n = theGeometry->theNodes->nNodes;

    femFieldRead(&n,2,&theSoluce[0],"../../connexion/U.txt");
    femFieldRead(&n,2,&theSoluce[1],"../../connexion/V.txt");
    femElasticityPrint(theProblem);

    double *X = theProblem->geometry->theNodes->X;
    double *Y = theProblem->geometry->theNodes->Y;

    printf(" ==== Soluce : \n");
    for (int i=0; i<theProblem->geometry->theNodes->nNodes*2; i++) {
        if (i % 2 == 0) printf("U%d = %14.7e", i/2, theSoluce[i]);
        else printf("V%d = %14.7e", i/2, theSoluce[i]); 
        
        if (theProblem->constrainedNodes[i] == -1) {
            if (i%2==0) printf("\t expected : %14.7e\n", q*nu*(X[i/2])/E);
            else        printf("\t expected : %14.7e\n", -q*Y[i/2]/E);
        }
        else printf("\t exepcted :  0.0\n");
    }
    
//
//  -2- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i=0; i<n; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]); }
  
    double hMin = femMin(normDisplacement,n);  
    double hMax = femMax(normDisplacement,n);  
    printf(" ==== Minimum displacement          : %14.7e \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e \n",hMax);
        
//
//  -3- Visualisation 
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Project 2022-23 ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
             glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);  }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    femElasticityFree(theProblem) ; 
    geoFree();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
