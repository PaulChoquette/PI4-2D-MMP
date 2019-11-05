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
    double Gamma = 1.4;

	Reader FileContents;
	matrix matrix;
	FileContents.nnodemax = 4;

    //string File_Name = "naca0012_129x129_1B_JAMESON.su2";
	//string File_Name = "square_tube.su2";
	string File_Name = "33_33.su2";

	if (FileContents.OpenFile(File_Name))
	{
		FileContents.read_file(File_Name);
		FileContents.file.close();
	}
	//matrix.printMatrix_unsigned(FileContents.inpoel1_wf,FileContents.nelem+FileContents.totalmarken,FileContents.nnodemax);
	//matrix.printMatrix_double(FileContents.coord,FileContents.npoin,FileContents.ndime);

	
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
    unsigned **inpoel = FileContents.inpoel1_wf; 
    double **coord = FileContents.coord; 
    int **esuel = maillage.init(inpoel,"esuel");
    int **elem2face = maillage.init(inpoel,"elem2face");
    int **face2elem = maillage.init(inpoel,"face2elem");
    int **face2node = maillage.init(inpoel,"face2node");
    double **Face2Vec = maillage.Face2Vec(elem2face, face2node, coord);
    double **Elem2Vec_x = maillage.Elem2Vec_x(inpoel,coord);
    double **Elem2Vec_y = maillage.Elem2Vec_y(inpoel,coord);
    double **Elem2Area = maillage.Elem2Area(inpoel,coord);
    double **Elem2Normal_x = maillage.Elem2Normal(inpoel,coord,Elem2Vec_x,Elem2Vec_y, "x");
    double **Elem2Normal_y = maillage.Elem2Normal(inpoel,coord,Elem2Vec_x,Elem2Vec_y, "y");
    double **Face2Normal_u = maillage.Face2Normal(Elem2Normal_x,Elem2Normal_y,elem2face,"x");
    double **Face2Normal_v = maillage.Face2Normal(Elem2Normal_x,Elem2Normal_y,elem2face,"y");

    maillage.Face2Area_compute(Face2Vec);

//////////////////////-------------------------------------------------------------------------//////////////////////
//////////////////////-------------------------------------------------------------------------//////////////////////
    mesh mesh_;
    mesh_.elem2elem = esuel;
    mesh_.nbFace = maillage.FaceNumber;
    mesh_.Face2Area = maillage.Face2Area;
    mesh_.Elem2Normal_u = Elem2Normal_x;
    mesh_.Elem2Normal_v = Elem2Normal_y;
	mesh_.Face2Normal_u = Face2Normal_u;
	mesh_.Face2Normal_v = Face2Normal_v;
    mesh_.face2elem = face2elem;
    mesh_.NbrPrimitive = 4;
    mesh_.nelem = maillage.nelem;
    mesh_.nelem_REAL = FileContents.nelem;
    mesh_.nbNode = FileContents.npoin;
	mesh_.celltype = FileContents.celltype;
    mesh_.rho_0 = 1.0;
    mesh_.P_0 = 1.0;
    mesh_.u_0 = 0.5*sqrt(1.4);
    mesh_.v_0 = 0.0;
    cout << " --- Elem2Area --- ";cout << "\n";
    matrix.printMatrix_double(Elem2Area,mesh_.nelem,1);

cout << "-------------------------------------------------------------------------------------------------------------------\n";
cout << "-------------------------------------------------------------------------------------------------------------------\n";
cout << "----------------------------------        DEBUT RESOLUTION        -------------------------------------------------\n";
cout << "-------------------------------------------------------------------------------------------------------------------\n";
cout << "-------------------------------------------------------------------------------------------------------------------\n";
    double cfl = 0.9;
    double somme_sX, somme_sY, lam_x, lam_y;
    int nfael;
    long double dt;
    

    mesh_.primitive_init(mesh_.rho_0,mesh_.P_0,mesh_.u_0,mesh_.v_0,mesh_.NbrPrimitive);
    double Critere = 1;
    int it = 1;
    //while (Critere > 0.00000000001)
	double** soundSpeed = matrix.generateMatrix_double(mesh_.nelem, 1);
	double** dts = matrix.generateMatrix_double(mesh_.nelem, 1);
    while (it < 1000)
    {
        cout << "-------------------------------------------------------------------------------------------------------------------\n";
        cout << "-------------------------------------------------------------------------------------------------------------------\n";
        cout << " IT # : ";cout << it;cout<<"\n";
        it +=1;
        cout << "-------------------------------------------------------------------------------------------------------------------\n";
        cout << "-------------------------------------------------------------------------------------------------------------------\n";

        mesh_.saveConservative();
        mesh_.roe_compute();
    //////////////////////-------------------------------------------------------------------------//////////////////////
        for (int k = 0; k < mesh_.nelem_REAL; ++k) // 0 to 8
        {
            nfael = maillage.Get_nfael(maillage.vtk,k);
            somme_sX = 0;
            somme_sY = 0;
            for (int i=0; i<nfael; i++)
            {
                int face_ = elem2face[k][i] - 1;
                somme_sX += abs( Face2Normal_u[face_][0]*mesh_.Face2Area[face_][0] );
                somme_sY += abs( Face2Normal_v[face_][0]*mesh_.Face2Area[face_][0] );
            }
            double delta_Sx = 0.5*somme_sX;
            double delta_Sy = 0.5*somme_sY;

            soundSpeed[k][0] = sqrt((mesh_.pressure_[k][0]/ mesh_.density_[k][0])*Gamma);

            lam_x = (abs(mesh_.velocity_x[k][0]) + soundSpeed[k][0]) * delta_Sx;
            lam_y = (abs(mesh_.velocity_y[k][0]) + soundSpeed[k][0]) * delta_Sy;


            dts[k][0] = cfl * Elem2Area[k][0] / (lam_x + lam_y);
            //cout << dts[k][0];cout << "\n";
        }
    //////////////////////-------------------------------------------------------------------------//////////////////////
        mesh_.ExpliciteTime_euler(dts, Elem2Area);
        mesh_.savePrimitive();
		/*
        cout << " --- Elem_flux --- ";cout << "\n";
        matrix.printMatrix_double(mesh_.Elem_flux,mesh_.nelem, 4);
        cout << " --- density_ --- ";cout << "\n";
        matrix.printMatrix_double(mesh_.density_,mesh_.nelem, 1);
        cout << " --- velocity_x --- ";cout << "\n";
        matrix.printMatrix_double(mesh_.velocity_x,mesh_.nelem, 1);
        cout << " --- velocity_y --- ";cout << "\n";
        matrix.printMatrix_double(mesh_.velocity_y,mesh_.nelem, 1);
        cout << " --- pressure_ --- ";cout << "\n";
        matrix.printMatrix_double(mesh_.pressure_,mesh_.nelem, 1);
		*/
    //////////////////////-------------------------------------------------------------------------//////////////////////
        Critere = mesh_.Elem_flux[0][0]/mesh_.density_[0][0];
        for (int i=1; i<mesh_.nelem_REAL; i++)
        {
            if (mesh_.Elem_flux[i][0]/mesh_.density_[i][0] > Critere)
            {
                Critere = mesh_.Elem_flux[i][0]/mesh_.density_[i][0];
            }
        }
        cout << " ------------------------------------------------------------------";cout << "\n";
        cout << "Critere : ";cout << Critere;cout << "\n";
        cout << " ------------------------------------------------------------------";cout << "\n";
		matrix.deleteMatrix(mesh_.Elem_flux, mesh_.nelem, 4);
		matrix.deleteMatrix(mesh_.massFlux_, mesh_.nbFace, 1);
		matrix.deleteMatrix(mesh_.momentumFlux_x, mesh_.nbFace, 1);
		matrix.deleteMatrix(mesh_.momentumFlux_y, mesh_.nbFace, 1);
		matrix.deleteMatrix(mesh_.energyFlux_, mesh_.nbFace, 1);
    }
    
	double** vararr = matrix.generateMatrix_double(mesh_.nelem_REAL,4);

	for (int i = 0; i < mesh_.nelem_REAL; i++) {
		vararr[i][0] = mesh_.density_[i][0];
		vararr[i][1] = mesh_.velocity_x[i][0];
		vararr[i][2] = mesh_.velocity_y[i][0];
		vararr[i][3] = mesh_.pressure_[i][0];
	}


    //////////////////////-------------------------------------------------------------------------//////////////////////
    //////////////////////-------------------------------------------------------------------------//////////////////////
        double** SumVector = matrix.generateMatrix_double(FileContents.nelem, 1);
        for (int i = 0; i < FileContents.nelem; i++) 
        {
            int znfael = maillage.Get_nnode(FileContents.vtk, i);
            for (int j = 0; j <znfael ; j++)
            {
                if (j == 0)
                {
                    *(double*)SumVector[i] = 0.0;
                }
                *(double*)SumVector[i] = *(double*)SumVector[i] + Elem2Vec_x[i][j] + Elem2Vec_y[i][j];
            }
        }

        Writer Output;
        Output.npoin = FileContents.npoin;
        Output.nelem = FileContents.nelem;
        Output.nnodemax = FileContents.nnodemax;
        string varname[6] = {"X","Y","Density","Velocity_x","Velocity_y","Pressure" };
        Output.Write_Output("test.dat", 6,varname, FileContents.inpoel1, coord, vararr);
    //////////////////////-------------------------------------------------------------------------//////////////////////
    //////////////////////-------------------------------------------------------------------------//////////////////////

    return 0;

}