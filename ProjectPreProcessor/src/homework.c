#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie 🙂
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade 
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH


double geoSize(double x, double y) {

    femGeo* theGeometry = geoGetGeometry();
    return theGeometry->h * (1.0 - 0.5 * x);
}


void geoBasicElasticityProblem() {

    //    
    //              q_y force
    //        +-----------------+ ^ 
    //        |                 | |
    //        |                 | | y = LyPlate
    //        |                 | |
    // ------ *-----------------* _ ------
    //        |-----------------> 
    //            x = LxPlate
    //

    femGeo *theGeometry = geoGetGeometry();

    int ierr;
    gmshClear(&ierr);

    int L = 3;
    int H = 1;

    gmshModelOccAddPoint(0,0,0,4,1,&ierr);
    gmshModelOccAddPoint(0,H,0,4,2,&ierr);
    gmshModelOccAddPoint(L,H,0,4,3,&ierr);
    gmshModelOccAddPoint(L,0,0,4,4,&ierr);
    gmshModelOccAddLine(1,2,1,&ierr);
    gmshModelOccAddLine(2,3,2,&ierr);
    gmshModelOccAddLine(3,4,3,&ierr);
    gmshModelOccAddLine(4,1,4,&ierr);
    int curves[4] = {1,2,3,4};
    gmshModelOccAddWire(curves, 4, 1, 1, &ierr);
    int wiretags[1] = {1};
    gmshModelOccAddPlaneSurface(wiretags, 1, 1, &ierr);

    gmshModelOccSynchronize(&ierr);
    gmshModelMeshGenerate(2, &ierr);

    return;
}



void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}


void geoMeshGenerateGeo() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(-w/2, -h/2, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( w/2, -h/2, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( w/2,  h/2, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(-w/2,  h/2, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(-w/2,    r, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0.,      r, 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(0.,     -r, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(-w/2,   -r, 0., lc, 8, &ierr);
    int p9 = gmshModelGeoAddPoint(0.,     0., 0., lc, 9, &ierr); // center of circle


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

 //   gmshFltkRun(&ierr);
}
void geoMeshGenerateGeobon() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(0., -0.32, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( 0.6, 0., 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( 1,  0., 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(1, 0.1, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(0, 0.1, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0,0., 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(0.4, 0., 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(0,   -0.2, 0., lc, 8, &ierr);

    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    //int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    int l6 = gmshModelGeoAddLine(p6, p7, 6, &ierr);
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    int lTags[] = {l1, l2, l3, l4, l5, l6, l7, l8}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

 //   gmshFltkRun(&ierr);
}

void geoMeshGenerateGeo2bon() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(0., -0.55, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( 0.55, 0., 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( 1,  0., 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(1, 0.1, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(0, 0.1, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0,0., 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(0.4, 0., 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(0,   -0.4, 0., lc, 8, &ierr);
    int p9 = gmshModelGeoAddPoint(0,   0.0, 0., lc, 9, &ierr); //centre cercle 1
    int p10 = gmshModelGeoAddPoint(1.5,   0.05, 0., lc, 10, &ierr); //centre cercle 2

    //int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l1= gmshModelGeoAddCircleArc(p2, p9, p1, 1, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddCircleArc(p4, p10, p3, 3, 0., 0., 0., &ierr);
    //int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddLine(p6, p7, 6, &ierr);
    int l7 = gmshModelGeoAddCircleArc(p8, p9, p7, 7, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    //int l8 = gmshModelGeoAddLine(p7, p8, 8, &ierr);
    int l9 = gmshModelGeoAddLine(p8, p1, 9, &ierr);

    int lTags[] = {-l1, l2, -l3, l4, l5, l6, -l7, l9}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

 //   gmshFltkRun(&ierr);
}

void geoMeshGenerateGeo2() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(0.7, h/-3., 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint(1., h/-3., 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( 1.,0., 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint( 1.5, 0., 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(1.5, h/2., 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0.5,    h/2., 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(0.5,  0.  , 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(w/2.,  0., 0., lc, 8, &ierr);
    
    //int p9 = gmshModelGeoAddPoint(0.,     0., 0., lc, 9, &ierr); // center of circle


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    //int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    int l6 = gmshModelGeoAddLine(p6, p7, 6, &ierr);
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    int lTags[] = {l1, l2, l3, l4, l5, l6, l7, l8}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

 //   gmshFltkRun(&ierr);

}
void geoMeshGenerateGeo3() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(0., 0., 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint(1.2, 0., 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( 1.2,1.2, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint( 0., 1.2, 0., lc, 4, &ierr);
    
    
    //int p9 = gmshModelGeoAddPoint(0.,     0., 0., lc, 9, &ierr); // center of circle


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p1, 4, &ierr);
    

    int lTags[] = {l1, l2, l3, l4}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 4, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",4,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

 //   gmshFltkRun(&ierr);
}

void geoMeshGenerateGeoFile(const char *filename){
    femGeo* theGeometry = geoGetGeometry();
    int ierr;
    gmshOpen(filename, &ierr); ErrorGmsh(ierr);
    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr); 
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr); 
        gmshModelMeshGenerate(2,&ierr);  }
 
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
    return;
}