 #include <iostream>
#include <sstream> 
#include <fstream>
#include <string>
#include <cctype>
#include <algorithm>
#include "main.h"
 
using namespace std;

double **primitive::init_0(double rho_0,double P_0,double u_0,double v_0, int NbrPrimitive, int nelem)
{
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
    }
    return mat;
}




