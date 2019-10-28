#include <iostream>
#include <sstream> 
#include <fstream>
#include <string>
#include <cctype>
#include <algorithm>
#include "main.h"
 
using namespace std;

void mesh::primitive_init(double rho_0,double P_0,double u_0,double v_0, int NbrPrimitive)
{
    matrix matrice;
    //pointer to pointer, to the element at row 1, column 1
    double ** mat = new double*[nelem];
    //for each row, create an array with size equal to number of elements in a column
    for(int i = 0; i < nelem; i++) 
    {
        mat[i] = new double[NbrPrimitive];
    }
    //return the pointer to pointer
    for(int i = 0; i < nelem; i++)
    {
        mat[i][0] = rho_0;
        mat[i][1] = P_0;
        mat[i][2] = u_0;
        mat[i][3] = v_0;
        /*
        
        string cell_type = Get_cell_type(i);
        if (cell_type == "center")
        {
            mat[i][0] = rho_0;
            mat[i][1] = P_0;
            mat[i][2] = u_0;
            mat[i][3] = v_0;
        }
        if (cell_type == "FarField")
        {
            mat[i][0] = rho_0;
            mat[i][1] = P_0;
            mat[i][2] = u_0;
            mat[i][3] = v_0;
        }
        if (cell_type == "NoSlipWall")
        {
            double n_u = abs(Face2Normal_u[i][0]);
            double n_v = abs(Face2Normal_v[i][0]);
            double n_u = ;
            double n_v = ;
            int Elem_voisin = elem2elem[i][0];
            mat[i][0] = mat[Elem_voisin][0];
            mat[i][1] = mat[Elem_voisin][1];
            mat[i][2] = mat[Elem_voisin][2] - 2*(mat[Elem_voisin][2]*n_u + mat[Elem_voisin][3]*n_v)*n_u;
            mat[i][3] = mat[Elem_voisin][3] - 2*(mat[Elem_voisin][2]*n_u + mat[Elem_voisin][3]*n_v)*n_v;
        }
        */
    }
    primitive_0 = mat;
    density_ = matrice.generateMatrix_double(nelem, 1); 
    pressure_ = matrice.generateMatrix_double(nelem, 1);
    velocity_x = matrice.generateMatrix_double(nelem, 1); 
    velocity_y = matrice.generateMatrix_double(nelem, 1);  
    for (int i=0;i<nelem;i++)
    {
        density_[i][0] = primitive_0[i][0];
        pressure_[i][0] = primitive_0[i][1];
        velocity_x[i][0] = primitive_0[i][2];
        velocity_y[i][0] = primitive_0[i][3];
    }
    conservative1_ = matrice.generateMatrix_double(nelem, 1); 
    conservative2_ = matrice.generateMatrix_double(nelem, 1);
    conservative3_ = matrice.generateMatrix_double(nelem, 1); 
    conservative4_ = matrice.generateMatrix_double(nelem, 1);
}



void mesh::roe_compute()
{
    double const Gamma = 1.4;
    matrix mat;
// IN & OUT //
//
    Elem_flux = mat.generateMatrix_double(nelem, 4);
    massFlux_ = mat.generateMatrix_double(nbFace, 1); 
    momentumFlux_x = mat.generateMatrix_double(nbFace, 1); 
    momentumFlux_y = mat.generateMatrix_double(nbFace, 1); 
    energyFlux_ = mat.generateMatrix_double(nbFace, 1); 
//
// IN & OUT //
/*
    cout << " --- massFlux_ --- ";cout << "\n";
    mat.printMatrix_double(massFlux_,nbFace, 1); 
    cout << " --- momentumFlux_x --- ";cout << "\n";
    mat.printMatrix_double(momentumFlux_x,nbFace, 1); 
    cout << "density_";cout << "\n";
    mat.printMatrix_double(density_,nelem, 1);
    cout << "pressure_";cout << "\n";
    mat.printMatrix_double(pressure_,nelem, 1);
    cout << "velocity_x";cout << "\n";
    mat.printMatrix_double(velocity_x,nelem, 1);
    cout << "velocity_y";cout << "\n";
    mat.printMatrix_double(velocity_y,nelem, 1);  
*/  
///////////////////////
// Boucle sur faces //
///////////////////////
    for (int i=0; i<nbFace; i++)
    {
        int Elem_L = face2elem[i][0] - 1;
        int Elem_R = face2elem[i][1] - 1;

        double n_u = (Face2Normal_u[i][0]);
        double n_v = (Face2Normal_v[i][0]);
        
        double roL = density_[Elem_L][0];
        double roR = density_[Elem_R][0];
        double uL = velocity_x[Elem_L][0];
        double uR = velocity_x[Elem_R][0];
        double vL = velocity_y[Elem_L][0];
        double vR = velocity_y[Elem_R][0];
        double pL = pressure_[Elem_L][0];
        double pR = pressure_[Elem_R][0];
        double eL = 0.5 * (uL * uL + vL * vL);
        double eR = 0.5 * (uR * uR + vR * vR);
        double hL = eL + pL / roL;
        double hR = eR + pR / roR;
        sqrtRoL = sqrt(roL);
        sqrtRoR = sqrt(roR);
// ROE_val
        roRoe = sqrt(roL*roR);


        uRoe = (uL*sqrtRoL + uR*sqrtRoR)/(sqrtRoL + sqrtRoR);
        vRoe = (vL*sqrtRoL + vR*sqrtRoR)/(sqrtRoL + sqrtRoR);

        hRoe = (hL*sqrtRoL + hR*sqrtRoR)/(sqrtRoL + sqrtRoR);

        qRoe2 = uRoe*uRoe + vRoe*vRoe;
        cRoe2 = (Gamma-1.0)*(hRoe-qRoe2/2.0);
        cRoe = sqrt(cRoe2);
        harten = 0.1*sqrt(Gamma*(pL+pR)/(roL+roR));

        deltaRo = roR-roL;
        deltaU  = uR-uL;
        deltaV  = vR-vL;
        
        V_Roe = uRoe*n_u + vRoe*n_v;
        double V_R = uR*n_u + vR*n_v;
        double V_L = uL*n_u + vL*n_v;
        double roAvg = 0.5*(roR + roL);
        double pAVG = 0.5*(pR + pL);
        double E_AVG = 0.5*(eR + eL);
        double uAVG = 0.5*(uL + uR); 
        double vAVG = 0.5*(vL + vR); 
        //double V_avg = 0.5*(V_R + V_L);
        double V_avg = uAVG*n_u + vAVG*n_v;
        deltaV_  = V_R-V_L;
        deltaP = pR-pL;

        lambda1 = abs(V_Roe-cRoe);
        if (lambda1<=harten)
        {
            lambda1 = (lambda1*lambda1 + harten*harten)/(2.0*harten);
        }
        lambda2 = abs(V_Roe);
        if (lambda2<=harten)
        {
            lambda2 = (lambda2*lambda2 + harten*harten)/(2.0*harten);
        }
        lambda3 = abs(V_Roe+cRoe);
        if (lambda3<harten)
        {
            lambda3 = (lambda3*lambda3 + harten*harten)/(2.0*harten);
        }

        coef1 = lambda1*(deltaP-roRoe*cRoe*deltaV_)/(2.0*cRoe2);
        f11 = coef1;
        f12 = coef1 * (uRoe-cRoe*n_u);
        f13 = coef1 * (vRoe-cRoe*n_v);
        f14 = coef1 * (hRoe-cRoe*V_Roe);

        coef2 = deltaRo-deltaP/cRoe2;
        f21 = lambda2*(coef2);
        f22 = lambda2*(coef2*uRoe + roRoe*(deltaU-deltaV_*n_u));
        f23 = lambda2*(coef2*vRoe + roRoe*(deltaV-deltaV_*n_v));
        f24 = lambda2*(coef2*qRoe2/2.0 + roRoe*(uRoe*deltaU + vRoe*deltaV - V_Roe*deltaV_));

        coef3 = lambda3 * (deltaP+roRoe*cRoe*deltaV_)/(2.0*cRoe2);
        f31 = coef3;
        f32 = coef3*(uRoe + cRoe*n_u);
        f33 = coef3*(vRoe + cRoe*n_v); 
        f34 = coef3*(hRoe + cRoe*V_Roe);

// Là c'est vraiment le shcéma de Roe
        ////////////////////////////
/*
        massFlux = uAvg*roAvg
        momentumFlux = uAvg**2*roAvg+pAvg
        energyFlux = uAvg*(0.5 * uAvg**2 * roAvg + pAvg/(gamma -1.0) + pAvg)
*/
/////////////////////////////////////////////////////
////// HHHMMMMMMMMMMMMMMM.....  /////////////////////
/////////////////////////////////////////////////////
/*
        massFlux_dob = 0.5* ( roL*V_L + roR*V_R );
        momentumFlux_X_dob = 0.5* ( (roL*uL*V_L + n_u*pL) + (roR*uR*V_R + n_u*pR) );     //V_avg*uAvg*roAvg + n_u*pAvg;
        momentumFlux_Y_dob = 0.5* ( (roL*vL*V_L + n_v*pL) + (roR*vR*V_R + n_v*pR) );     //V_avg*vAvg*roAvg + n_v*pAvg;
        energyFlux_dob = 0.5* ( roL*hL*V_L + roR*hR*V_R );  //V_avg*H_avg*roAvg; //(0.5 * uAvg*uAvg * roAvg + pAvg/(Gamma -1.0) + pAvg);
*/
        massFlux_dob = roAvg*V_avg;
        momentumFlux_X_dob = roAvg*uAVG*V_avg + n_u*pAVG;     //V_avg*uAvg*roAvg + n_u*pAvg;
        momentumFlux_Y_dob = roAvg*vAVG*V_avg + n_v*pAVG;     //V_avg*vAvg*roAvg + n_v*pAvg;
        energyFlux_dob = roAvg*E_AVG*V_avg + V_avg*pAVG;  //V_avg*H_avg*roAvg; //(0.5 * uAvg*uAvg * roAvg + pAvg/(Gamma -1.0) + pAvg);
/*
        cout << "massFlux_dob : ";cout << massFlux_dob;cout<< "\n";
        cout << "momentumFlux_X_dob : ";cout << momentumFlux_X_dob;cout<< "\n";
        cout << "momentumFlux_Y_dob : ";cout << momentumFlux_Y_dob;cout<< "\n";
        cout << "energyFlux_dob : ";cout << energyFlux_dob;cout<< "\n";
*/ 
        massFlux_dob -= 0.5*(f11+f21+f31);
        momentumFlux_X_dob -= 0.5*(f12+f22+f32);
        momentumFlux_Y_dob -= 0.5*(f13+f23+f33);
        energyFlux_dob -= 0.5*(f14+f24+f34);
/*
        cout << "massFlux_dob : ";cout << massFlux_dob;cout<< "\n";
        cout << "momentumFlux_X_dob : ";cout << momentumFlux_X_dob;cout<< "\n";
        cout << "momentumFlux_Y_dob : ";cout << momentumFlux_Y_dob;cout<< "\n";
        cout << "energyFlux_dob : ";cout << energyFlux_dob;cout<< "\n";
        ////////////////////////////
*/        
        massFlux_[i][0] = massFlux_dob;
        momentumFlux_x[i][0] = momentumFlux_X_dob;
        momentumFlux_y[i][0] = momentumFlux_Y_dob;
        energyFlux_[i][0] = energyFlux_dob;
    }
    for (int j=0; j<nbFace; j++)
    {
        int Elem_L = face2elem[j][0] - 1;
        int Elem_R = face2elem[j][1] - 1;
        //cout << "Elem_L : ";cout << Elem_L;cout << "\n";
        //cout << "Elem_R : ";cout << Elem_R;cout << "\n";
        double n_u = abs(Face2Normal_u[j][0]);
        double n_v = abs(Face2Normal_v[j][0]); 
        int Face_num = j;
        double face_vec_x = (Face2Normal_u[j][0]);
        double face_vec_y = (Face2Normal_v[j][0]);
        //cout << "face_vec_x : ";cout << face_vec_x;cout << "\n";
        //cout << "face_vec_y : ";cout << face_vec_y;cout << "\n";
        Elem_flux[Elem_L][0] += (n_u + n_v)*massFlux_[Face_num][0];
        Elem_flux[Elem_L][1] += (n_u + n_v)*momentumFlux_x[Face_num][0];
        Elem_flux[Elem_L][2] += (n_u + n_v)*momentumFlux_y[Face_num][0];
        Elem_flux[Elem_L][3] += (n_u + n_v)*energyFlux_[Face_num][0];

        Elem_flux[Elem_R][0] -= (n_u + n_v)*massFlux_[Face_num][0];
        Elem_flux[Elem_R][1] -= (n_u + n_v)*momentumFlux_x[Face_num][0];
        Elem_flux[Elem_R][2] -= (n_u + n_v)*momentumFlux_y[Face_num][0];
        Elem_flux[Elem_R][3] -= (n_u + n_v)*energyFlux_[Face_num][0];
        //mat.printMatrix_double(Elem_flux,8, 4);
    }

}

double mesh::Get_LocalMach(int Face_number,string choix)
{
    matrix mat;

    double const Gamma = 1.4;
    int Elem_L = face2elem[Face_number][0] - 1;
    int Elem_R = face2elem[Face_number][1] - 1;
    double **density_ = mat.generateMatrix_double(nelem, 1); 
    double **pressure_ = mat.generateMatrix_double(nelem, 1);
    double **velocity_x = mat.generateMatrix_double(nelem, 1); 
    double **velocity_y = mat.generateMatrix_double(nelem, 1);  
    for (int i=0;i<nelem;i++)
    {
        density_[i][0] = primitive_0[i][0];
        pressure_[i][0] = primitive_0[i][1];
        velocity_x[i][0] = primitive_0[i][2];
        velocity_y[i][0] = primitive_0[i][3];
    }
    double roL = density_[Elem_L][0];
    double roR = density_[Elem_R][0];
    double pL = pressure_[Elem_L][0];
    double pR = pressure_[Elem_R][0];
    double a_speed = sqrt(Gamma*(pL+pR)/(roL+roR));
    
    return 0;
}

void mesh::saveConservative()
{
    for (int i=0; i < nelem; i++) // 0 to 9
    {
        conservative1_[i][0] = density_[i][0];
        conservative2_[i][0] = density_[i][0] * velocity_x[i][0];
        conservative3_[i][0] = density_[i][0] * velocity_y[i][0];
        conservative4_[i][0] = density_[i][0]* (0.5*(velocity_x[i][0]*velocity_x[i][0] + velocity_y[i][0]*velocity_y[i][0])  +   pressure_[i][0] / density_[i][0]);
    }
}

void mesh::savePrimitive()
{
    for (int i = 0; i < nelem; i++)
    {
        density_[i][0] = conservative1_[i][0];
        velocity_x[i][0] = conservative2_[i][0] / conservative1_[i][0];
        velocity_y[i][0] = conservative3_[i][0] / conservative1_[i][0];
        pressure_[i][0] = (conservative4_[i][0]/conservative1_[i][0] - 0.5*((velocity_x[i][0]*velocity_x[i][0]) + velocity_y[i][0]*velocity_y[i][0]))/density_[i][0];
    }
}

