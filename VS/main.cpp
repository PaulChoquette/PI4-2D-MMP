#include <iostream>
#include <sstream> 
#include <fstream>
#include <string>
#include <cctype>
#include <algorithm>
#include "main.h"
 
using namespace std;


int main()
{
	Reader FileContents;
	matrix matrix;
	FileContents.nnodemax = 4;

	//ifstream meshfile("naca0012_129x129_1B_JAMESON.su2");
	string File_Name = "square.su2";

	if (FileContents.OpenFile(File_Name))
	{
		FileContents.read_1(File_Name);
		FileContents.file.close();
	}


	cout << "--- INPOEL ---";cout << "\n";
	matrix.printMatrix(FileContents.inpoel1,FileContents.nelem,FileContents.nnodemax);

	cout << "--- Coord ---";cout << "\n";
	matrix.printMatrix_double(FileContents.coord,FileContents.npoin,FileContents.ndime);

    cout << "--- NELEM ";cout << "  :  ";
    cout << FileContents.nelem;cout << "			";
	cout << "--- Position ";cout << "  :  ";
    cout << FileContents.lineNelem;cout << "\n";

    cout << "--- NPOIN ";cout << "  :  ";
    cout << FileContents.npoin;cout << "			";
	cout << "--- Position ";cout << "  :  ";
    cout << FileContents.poinlinen;cout << "\n";

	cout << "--- NNODE ";cout << "   :  ";
	cout << FileContents.nnode;cout << "\n";

	cout << "--- ndime ";cout << "   :  ";
	cout << FileContents.ndime;cout << "\n";

	
//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
    connect maillage;
    maillage.npoin = FileContents.npoin;
    maillage.nelem = FileContents.nelem;
    //maillage.nnode = FileContents.nnode;
	maillage.vtk = FileContents.vtk;
    maillage.poinperFace = 2;
    maillage.nfaelmax = 4;
	maillage.nnodemax = FileContents.nnodemax;
	maillage.NbNdPerFace = 2; //OK parce qu'on est en 2D pour l'instant...
    maillage.ddl = 2;
////////////////////////////////////////////////////////////////////
   // int inpoel_[8][3] = {{1, 2, 4}, {2, 5, 4}, {2, 3, 5}, {3, 6, 5}, {4, 5, 7}, {5, 8, 7}, {5, 6, 8}, {6, 9, 8}};
//    int lpofa_[2][3] = {{1, 2, 3},{2, 3, 1}};
 //   int lnofa_[3]= {2, 2, 2};
    double coord_[9][2] = {{0, 0}, {1, 0},{2, 0},{0, 1}, {1, 1},{2, 1},{0, 2},{1, 2},{2, 2}};
////////////////////////////////////////////////////////////////////
    int **inpoel = FileContents.inpoel1; // matrix.generateMatrix(maillage.nelem, maillage.nnode); 
 //   int **lnofa = matrix.generateMatrix(maillage.nfael,1); 
//    int **lpofa = matrix.generateMatrix(maillage.poinperFace, maillage.nfael); 
    double **coord = FileContents.coord; // matrix.generateMatrix_double(maillage.npoin, 2); 
////////////////////////////////////////////////////////////////////
/*
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
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << ".....";cout << "\n";
    int **esuel = maillage.init(inpoel,"esuel");
    cout << ".....";cout << "\n";
    cout << "esuel (ou elem2elem)";cout << "\n";
    cout << ".....";cout << "\n";
    matrix.printMatrix(esuel,maillage.nelem,maillage.nfaelmax);
    cout << ".....";cout << "\n";
    
    int **elem2face = maillage.init(inpoel,"elem2face");
    cout << "elem2face";cout << "\n";
    matrix.printMatrix(elem2face,maillage.nelem,maillage.nfaelmax);

    int **face2elem = maillage.init(inpoel,"face2elem");
    cout << "face2elem";cout << "\n";
    matrix.printMatrix(face2elem,maillage.FaceNumber,2);

    int **face2node = maillage.init(inpoel,"face2node");
    cout << "face2node";cout << "\n";
    matrix.printMatrix(face2node,maillage.FaceNumber,2);
    
    cout << "coord";cout << "\n";
    matrix.printMatrix_double(coord,maillage.npoin,2);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int **Face2Vec = maillage.Face2Vec(elem2face, face2node, coord);
    cout << "Face2Vec";cout << "\n";
    matrix.printMatrix(Face2Vec,maillage.FaceNumber, maillage.NbNdPerFace);

    int **Elem2Node = maillage.Elem2Node(elem2face,face2node);
    cout << "Elem2Node";cout << "\n";
    matrix.printMatrix(Elem2Node,maillage.nelem, maillage.nfaelmax);

    double **Elem2Vec_x = maillage.Elem2Vec_x(Elem2Node,coord);
    cout << "Elem2Vec_x";cout << "\n";
    matrix.printMatrix_double(Elem2Vec_x,maillage.nelem, maillage.nnodemax);

    double **Elem2Vec_y = maillage.Elem2Vec_y(Elem2Node,coord);
    cout << "Elem2Vec_y";cout << "\n";
    matrix.printMatrix_double(Elem2Vec_y,maillage.nelem, maillage.nnodemax);  

    cout << " --- Elem2Vec --- ";cout << "\n";
    matrix.printXY_double(Elem2Vec_x,Elem2Vec_y,maillage.nelem, maillage.nnodemax);    

    
    double **Elem2Area = maillage.Elem2Area(Elem2Node,coord);
    cout << "Elem2Area";cout << "\n";
    matrix.printMatrix_double(Elem2Area,maillage.nelem, 1); 


    double **Elem2Normal_x = maillage.Elem2Normal(Elem2Node,coord,Elem2Vec_x,Elem2Vec_y, "x");
    cout << "Elem2Normal_x";cout << "\n";
    matrix.printMatrix_double(Elem2Normal_x,maillage.nelem, maillage.nnodemax);

    double **Elem2Normal_y = maillage.Elem2Normal(Elem2Node,coord,Elem2Vec_x,Elem2Vec_y, "y");
    cout << "Elem2Normal_y";cout << "\n";
    matrix.printMatrix_double(Elem2Normal_y,maillage.nelem, maillage.nnodemax);  

    cout << "elem2face";cout << "\n";
    matrix.printMatrix(elem2face,maillage.nelem,maillage.nfaelmax);

    cout << " --- Elem2Normal --- ";cout << "\n";
    matrix.printXY_double(Elem2Normal_x,Elem2Normal_y,maillage.nelem, maillage.nnodemax);

//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
    primitive prim;
    prim.NbrPrimitive = 4;
    prim.rho_0 = 0.01;
    prim.P_0 = 1000;
    prim.u_0 = 1;
    prim.v_0 = 1;
    double **PrimVec = prim.init_0(prim.rho_0,prim.P_0,prim.u_0,prim.v_0,prim.NbrPrimitive,maillage.nelem);
    cout << " --- PrimVec --- ";cout << "\n";
    matrix.printMatrix_double(PrimVec,maillage.nelem, prim.NbrPrimitive);  


//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
 
    return 0;

}