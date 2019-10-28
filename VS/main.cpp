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

	string File_Name = "square_tube.su2";
	//string File_Name = "naca0012_129x129_1B_JAMESON.su2";

	if (FileContents.OpenFile(File_Name))
	{
		FileContents.read_file(File_Name);
		FileContents.file.close();
	}

//Debug Printout, removing accelerates code immensely
	cout << "--- INPOEL ---";cout << "\n";
	matrix.printMatrix_unsigned(FileContents.inpoel1_wf,FileContents.nelem+FileContents.totalmarken,FileContents.nnodemax);

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

	cout << "--- ndime ";cout << "   :  ";
	cout << FileContents.ndime;cout << "\n";

	
//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
    connect maillage;
    maillage.npoin = FileContents.npoin;
    maillage.nelem = FileContents.nelem+FileContents.totalmarken;
	maillage.vtk = FileContents.vtk_wf;
    maillage.poinperFace = 2;
    maillage.nfaelmax = 4;
	maillage.nnodemax = FileContents.nnodemax;
	maillage.NbNdPerFace = 2; //OK parce qu'on est en 2D pour l'instant...
    maillage.ddl = 2;
////////////////////////////////////////////////////////////////////
    unsigned **inpoel = FileContents.inpoel1_wf; // matrix.generateMatrix(maillage.nelem, maillage.nnode); 											//   int **lnofa = matrix.generateMatrix(maillage.nfael,1); 
//    int **lpofa = matrix.generateMatrix(maillage.poinperFace, maillage.nfael); 
    double **coord = FileContents.coord; // matrix.generateMatrix_double(maillage.npoin, 2); 
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

    //int **inpoel = maillage.inpoel(elem2face,face2node);
    cout << "inpoel";cout << "\n";
    //matrix.printMatrix(inpoel,maillage.nelem, maillage.nfaelmax);

    double **Elem2Vec_x = maillage.Elem2Vec_x(inpoel,coord);
    cout << "Elem2Vec_x";cout << "\n";
    matrix.printMatrix_double(Elem2Vec_x,maillage.nelem, maillage.nnodemax);

    double **Elem2Vec_y = maillage.Elem2Vec_y(inpoel,coord);
    cout << "Elem2Vec_y";cout << "\n";
    matrix.printMatrix_double(Elem2Vec_y,maillage.nelem, maillage.nnodemax);  

    cout << " --- Elem2Vec --- ";cout << "\n";
    matrix.printXY_double(Elem2Vec_x,Elem2Vec_y,maillage.nelem, maillage.nnodemax);    

    
    double **Elem2Area = maillage.Elem2Area(inpoel,coord);
    cout << "Elem2Area";cout << "\n";
    matrix.printMatrix_double(Elem2Area,maillage.nelem, 1); 


    double **Elem2Normal_x = maillage.Elem2Normal(inpoel,coord,Elem2Vec_x,Elem2Vec_y, "x");
    cout << "Elem2Normal_x";cout << "\n";
    matrix.printMatrix_double(Elem2Normal_x,maillage.nelem, maillage.nnodemax);

    double **Elem2Normal_y = maillage.Elem2Normal(inpoel,coord,Elem2Vec_x,Elem2Vec_y, "y");
    cout << "Elem2Normal_y";cout << "\n";
    matrix.printMatrix_double(Elem2Normal_y,maillage.nelem, maillage.nnodemax);  

    cout << "elem2face";cout << "\n";
    matrix.printMatrix(elem2face,maillage.nelem,maillage.nfaelmax);

    cout << " --- Elem2Normal --- ";cout << "\n";
    matrix.printXY_double(Elem2Normal_x,Elem2Normal_y,maillage.nelem, maillage.nnodemax);

    double **Face2Normal_u = maillage.Face2Normal(Elem2Normal_x,Elem2Normal_y,elem2face,"x");
    double **Face2Normal_v = maillage.Face2Normal(Elem2Normal_x,Elem2Normal_y,elem2face,"y");
    cout << " --- Face2Normal --- ";cout << "\n";
    matrix.printXY_double(Face2Normal_u,Face2Normal_v,maillage.FaceNumber, 1);

    //cout << " --- maillage.vtk --- ";cout << "\n";
    //matrix.printMatrix(maillage.vtk,maillage.nelem, 1);

//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
    mesh mesh_;
    mesh_.elem2elem = esuel;
    mesh_.nbFace = maillage.FaceNumber;
    mesh_.Face2Normal_u = Face2Normal_u;
    mesh_.Face2Normal_v = Face2Normal_v;
    mesh_.face2elem = face2elem;
    mesh_.NbrPrimitive = 4;
    mesh_.nelem = maillage.nelem;
    mesh_.nbNode = FileContents.npoin;
	mesh_.celltype = FileContents.celltype;
    mesh_.rho_0 = 0.01;
    mesh_.P_0 = 1000;
    mesh_.u_0 = 1;
    mesh_.v_0 = 1;
//
    mesh_.primitive_init(mesh_.rho_0,mesh_.P_0,mesh_.u_0,mesh_.v_0,mesh_.NbrPrimitive);
//
    cout << "mesh_.nelem : ";cout << mesh_.nelem;cout << "\n";
    //cout << "mesh_.NbrPrimitive : ";cout << mesh_.NbrPrimitive;cout << "\n";
    //cout << "primitive_0";cout << "\n";
    //cout << "- Rho - P - U - V -";cout << "\n";
    //matrix.printMatrix_double(mesh_.primitive_0,mesh_.nelem,mesh_.NbrPrimitive);
//  
    mesh_.saveConservative();
    mesh_.roe_compute();

    mesh_.savePrimitive();
//
    cout << " --- massFlux_ --- ";cout << "\n";
    matrix.printMatrix_double(mesh_.massFlux_,mesh_.nbFace, 1); 

    //cout << " --- momentumFlux_x --- ";cout << "\n";
    //matrix.printMatrix_double(mesh_.momentumFlux_x,mesh_.nbFace, 1); 

    //cout << " --- momentumFlux_y --- ";cout << "\n";
    //matrix.printMatrix_double(mesh_.momentumFlux_y,mesh_.nbFace, 1); 

    //cout << " --- energyFlux_ --- ";cout << "\n";
    //matrix.printMatrix_double(mesh_.energyFlux_,mesh_.nbFace, 1);  

    //cout << " --- saveConservative --- ";cout << "\n";
    //matrix.printMatrix_double(mesh_.conservative1_,mesh_.nelem, 1);

    //cout << " --- savePrimitive --- ";cout << "\n";
    //matrix.printMatrix_double(mesh_.density_,mesh_.nelem, 1);

    cout << " --- Elem_flux --- ";cout << "\n";
    matrix.printMatrix_double(mesh_.Elem_flux,mesh_.nelem, 4);
    
//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
	double** SumVector = matrix.generateMatrix_double(FileContents.nelem, 1);
	for (int i = 0; i < FileContents.nelem; i++) {
		int znfael = maillage.Get_nnode(FileContents.vtk, i);
		for (int j = 0; j <znfael ; j++) {
			if (j == 0) {
				*(double*)SumVector[i] = 0.0;
			}
			*(double*)SumVector[i] = *(double*)SumVector[i] + Elem2Vec_x[i][j] + Elem2Vec_y[i][j];
		}
	}

	Writer Output;
	Output.npoin = FileContents.npoin;
	Output.nelem = FileContents.nelem;
	Output.nnodemax = FileContents.nnodemax;
	string varname[3] = {"X","Y","SUmVec"};
	Output.Write_Output("test.dat", 3,varname, FileContents.inpoel1, coord, SumVector);
//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
 
    return 0;

}