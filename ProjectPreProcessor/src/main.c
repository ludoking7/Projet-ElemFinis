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

int main(int argc, char* argv[])
{

//
//  -1- Construction de la geometrie
//

    double Lx = 1.0;
    double Ly = 1.0;
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();

    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly;
    theGeometry->h           =  Lx * 0.05;
    theGeometry->elementType = FEM_TRIANGLE;

    //geoBasicElasticityProblem();

    designHouse(20.0, 20.0, 4.0, 6.0, 4.0, 5.0, 7.2, 0.5);

    //geoMeshGenerate();      // Utilisation de OpenCascade

//    geoMeshGenerateGeo();   // Utilisation de outils de GMSH
    // Attention : les entit�s sont diff�rentes !
    // On a aussi invers� la g�omtrie pour rire !

//    geoMeshGenerateGeoFile("../data/mesh.geo");   // Lecture fichier geo
    /*
      geoMeshImport();
      //geoMeshPrint();
      //geoSetDomainName(1, "Top");
      //geoSetDomainName(3, "Bottom");
      geoSetDomainName(0,"Symetry");
      geoSetDomainName(7,"Bottom");
      geoMeshWrite("../../data/mesh.txt");

  //
  //  -2- Definition du probleme
  //

      double E   = 211.e9;
      double nu  = 0.3;
      double rho = 7.85e3;
      double g   = 9.81;
      femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN, FEM_BAND);
      femElasticityAddBoundaryCondition(theProblem,"Symetry",DIRICHLET_X,0.0);
      femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
      femElasticityPrint(theProblem);
      femElasticityWrite(theProblem,"../../data/problem.txt");*/

//
//  -2.1- Definition du probleme basic elasticity
//
/*
    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3;
    double g   = 0;
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS, FEM_BAND);
    theProblem->constrainedNodes[0] = 0;  //A enlever !
    double q = 100;
    femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_Y, -q);
    femElasticityPrint(theProblem);
    femElasticityWrite(theProblem, "../../data/problem.txt");
*/
// -2.2- Definition du probleme maison

    geoMeshImport();
    //geoMeshPrint();
    geoSetDomainName(2,"BottomL");
    //geoSetDomainName(4,"LeftDoor");
    //geoSetDomainName(6,"BottomDoor");
    //geoSetDomainName(8,"RightDoor");
    geoSetDomainName(10,"BottomR");
    geoMeshWrite("../../data/mesh.txt");
    double E   = 14e11;
    double nu  = 0.22;
    double rho = 1.8e3;
    double g   = 9.81;
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN,FEM_BAND, FEM_YNUM);
    femElasticityAddBoundaryCondition(theProblem,"BottomL",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"BottomR",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"BottomL",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"BottomR",DIRICHLET_X,0.0);
    //femElasticityAddBoundaryCondition(theProblem,"BottomDoor",DIRICHLET_Y,0.0);
    //femElasticityAddBoundaryCondition(theProblem,"RightDoor",DIRICHLET_X,0.0);
    //femElasticityAddBoundaryCondition(theProblem,"LeftDoor",DIRICHLET_X,0.0);
    femElasticityPrint(theProblem);
    femElasticityWrite(theProblem,"../../data/problem.txt");


//
//  -3- Champ de la taille de r�f�rence du maillage
//

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);
    double hMax = femMax(meshSizeField,theNodes->nNodes);
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);

//
//  -4- Visualisation
//

    if(argc >= 2 && strcmp(argv[1], "view") == 0){
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
                glfemPlotField(theGeometry->theElements,meshSizeField);
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
    }


    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);
    femElasticityFree(theProblem) ;
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return 0;
}


