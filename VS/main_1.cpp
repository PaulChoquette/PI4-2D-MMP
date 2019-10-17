#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional> // std::divides 
#include "main.h" 
 

int main()
{
////////////////////////////////////////////////////////////////////
    connect maillage;
    matrix mat;
    maillage.npoin = 9;
    maillage.nelem = 8;
    maillage.nnode = 3;
    maillage.poinperFace = 2;
    maillage.nfael = 3;
    maillage.NbNdPerFace = 2;
    maillage.NbFacePerELEM = 3;
    maillage.ddl = 2;
////////////////////////////////////////////////////////////////////
    int inpoel_[8][3] = {{1, 2, 4}, {2, 5, 4}, {2, 3, 5}, {3, 6, 5}, {4, 5, 7}, {5, 8, 7}, {5, 6, 8}, {6, 9, 8}};
    int lpofa_[2][3] = {{1, 2, 3},{2, 3, 1}};
    int lnofa_[3]= {2, 2, 2};
    double coord_[9][2] = {{0, 0}, {1, 0},{2, 0},{0, 1}, {1, 1},{2, 1},{0, 2},{1, 2},{2, 2}};
////////////////////////////////////////////////////////////////////
    int **inpoel = mat.generateMatrix(maillage.nelem, maillage.nnode); 
    int **lnofa = mat.generateMatrix(maillage.nfael,1); 
    int **lpofa = mat.generateMatrix(maillage.poinperFace, maillage.nfael); 
    double **coord = mat.generateMatrix_double(maillage.npoin, 2); 
////////////////////////////////////////////////////////////////////
    for(int i = 0; i < maillage.nelem; i++) 
    {
        for(int j = 0; j < maillage.nnode; j++) 
        {
            inpoel[i][j] = inpoel_[i][j];
        }
    } 
    for(int i=0; i<3; i++)
    {
        *(int*)lnofa[i] = lnofa_[i];
    }
    for(int i=0; i<maillage.poinperFace; i++)
    {
        for(int j=0; j<maillage.nfael; j++)
        {
            lpofa[i][j] = lpofa_[i][j];
        }
    }
    for(int i = 0; i < maillage.npoin; i++) 
    {
        for(int j = 0; j < 2; j++) 
        {
            coord[i][j] = coord_[i][j];
        }
    } 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int **esuel = maillage.init(inpoel,lpofa,lnofa,"esuel");
    cout << "esuel (ou elem2elem)";cout << "\n";
    mat.printMatrix(esuel,maillage.nelem,maillage.nfael);
    
    int **elem2face = maillage.init(inpoel,lpofa,lnofa,"elem2face");
    cout << "elem2face";cout << "\n";
    mat.printMatrix(elem2face,maillage.nelem,maillage.nfael);

    int **face2elem = maillage.init(inpoel,lpofa,lnofa,"face2elem");
    cout << "face2elem";cout << "\n";
    mat.printMatrix(face2elem,maillage.FaceNumber,2);

    int **face2node = maillage.init(inpoel,lpofa,lnofa,"face2node");
    cout << "face2node";cout << "\n";
    mat.printMatrix(face2node,maillage.FaceNumber,2);
    
    cout << "coord";cout << "\n";
    mat.printMatrix_double(coord,maillage.npoin,2);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int **Face2Vec = maillage.Face2Vec(elem2face, face2node, coord);
    cout << "Face2Vec";cout << "\n";
    mat.printMatrix(Face2Vec,maillage.FaceNumber, maillage.NbNdPerFace);

    int **Elem2Node = maillage.Elem2Node(elem2face,face2node);
    cout << "Elem2Node";cout << "\n";
    mat.printMatrix(Elem2Node,maillage.nelem, maillage.nfael);

    double **Elem2Vec_x = maillage.Elem2Vec_x(Elem2Node,coord);
    cout << "Elem2Vec_x";cout << "\n";
    mat.printMatrix_double(Elem2Vec_x,maillage.nelem, maillage.NbNdPerElem);

    double **Elem2Vec_y = maillage.Elem2Vec_y(Elem2Node,coord);
    cout << "Elem2Vec_y";cout << "\n";
    mat.printMatrix_double(Elem2Vec_y,maillage.nelem, maillage.NbNdPerElem);  

    cout << " --- Elem2Vec --- ";cout << "\n";
    mat.printXY_double(Elem2Vec_x,Elem2Vec_y,maillage.nelem, maillage.NbNdPerElem);    

    
    double **Elem2Area = maillage.Elem2Area(Elem2Node,coord);
    cout << "Elem2Area";cout << "\n";
    mat.printMatrix_double(Elem2Area,maillage.nelem, 1); 


    double **Elem2Normal_x = maillage.Elem2Normal(Elem2Node,coord,Elem2Vec_x,Elem2Vec_y, "x");
    cout << "Elem2Normal_x";cout << "\n";
    mat.printMatrix_double(Elem2Normal_x,maillage.nelem, maillage.NbNdPerElem);

    double **Elem2Normal_y = maillage.Elem2Normal(Elem2Node,coord,Elem2Vec_x,Elem2Vec_y, "y");
    cout << "Elem2Normal_y";cout << "\n";
    mat.printMatrix_double(Elem2Normal_y,maillage.nelem, maillage.NbNdPerElem);  

    cout << " --- Elem2Normal --- ";cout << "\n";
    mat.printXY_double(Elem2Normal_x,Elem2Normal_y,maillage.nelem, maillage.NbNdPerElem);

}


