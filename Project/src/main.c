/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.h"

int main(void)
{  
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../../connexion/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../../connexion/problem.txt");
    femElasticityPrint(theProblem);
    theProblem->constrainedNodes[0] = 0;
    double *theSoluce = femElasticitySolve(theProblem); 
    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../../connexion/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../../connexion/V.txt");
    femElasticityFree(theProblem); 
    geoFree();
    return 0;  
}

 
