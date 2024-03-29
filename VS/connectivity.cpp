#include <iostream>
#include <vector>  
#include <algorithm>
#include <cmath> 
#include <functional> // std::divides 
#include "main.h"

int** connect::Get_lpofa(int** vtk, int ielem)
{
	int vtk_ind = *(int*)vtk[ielem];
	matrix mat;
	int** lpofa;
	if (*(int*)vtk[ielem] == 5) {
		lpofa = mat.generateMatrix(2, 3);
		int lpofa_[2][3] = { {1, 2, 3},{2, 3, 1} };
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				lpofa[i][j] = lpofa_[i][j];
			}
		}
		return lpofa;
	}
	else if (*(int*)vtk[ielem] == 9) {
		lpofa = mat.generateMatrix(2, 4);
		int lpofa_[2][4] = { {1, 2, 3, 4},{2, 3, 4, 1} };
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				lpofa[i][j] = lpofa_[i][j];
			}
		}
		return lpofa;
	}
	else if (*(int*)vtk[ielem] == 3) {
		lpofa = mat.generateMatrix(2, 1);
		int lpofa_[2][1] = { {1},{2} };
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 1; j++)
			{
				lpofa[i][j] = lpofa_[i][j];
			}
		}
		return lpofa;
	}
	else {
		return 0;
	}
}

int** connect::Get_lnofa(int** vtk, int ielem)
{
	int vtk_ind = *(int*)vtk[ielem];
	matrix mat;
	int** lnofa;
	if (*(int*)vtk[ielem] == 5) {
		lnofa = mat.generateMatrix(3, 1);
		int lnofa_[3] = { 2, 2, 2 };
		for (int i = 0; i < 3; i++)
		{
			*(int*)lnofa[i] = lnofa_[i];
		}
		return lnofa;
	}
	else if (*(int*)vtk[ielem] == 9) {
		lnofa = mat.generateMatrix(4, 1);
		int lnofa_[4] = { 2, 2, 2, 2 };
		for (int i = 0; i < 4; i++)
		{
			*(int*)lnofa[i] = lnofa_[i];
		}
		return lnofa;

	}
	else if (*(int*)vtk[ielem] == 3) {
		lnofa = mat.generateMatrix(1, 1);
		int lnofa_[1] = { 2 };
		for (int i = 0; i < 1; i++)
		{
			*(int*)lnofa[i] = lnofa_[i];
		}
		return lnofa;

	}
	else {
		return 0;
	}
}

int connect::Get_nnode(int ** vtk, int ielem)
{
	int vtk_ind = *(int*)vtk[ielem];
	matrix mat;
	int innode;
	int** lnofa;
	if (*(int*)vtk[ielem] == 5) {
		innode = 3;
		return innode;
	}
	else if (*(int*)vtk[ielem] == 9) {
		innode = 4;
		return innode;
	}
	else if (*(int*)vtk[ielem] == 3) {
		innode = 2;
		return innode;
	}
	else {
		return 0;
	}
}

int connect::Get_nfael(int** vtk, int ielem)
{
	int vtk_ind = *(int*)vtk[ielem];
	matrix mat;
	int infael;
	int** lnofa;
	if (*(int*)vtk[ielem] == 5) 
    {
		infael = 3;
		return infael;
	}
	else if (*(int*)vtk[ielem] == 9) 
    {
		infael = 4;
		return infael;
	}
	else if (*(int*)vtk[ielem] == 3) 
    {
		infael = 1;
		return infael;
	}
	else {
		return 0;
	}
}

int ** connect::init(unsigned **inpoel, string choix)
{
    matrix mat;
    int ipoil;
    int **esup2 = mat.generateMatrix(npoin+1, 1);
    int **esup2_test = mat.generateMatrix(npoin+1, 1); 
	for (int ipoin = 0; ipoin < npoin+1; ipoin++)
	{
		*(int*)esup2[ipoin] = 0;
		*(int*)esup2_test[ipoin] = 0;
	}
    for (int ielem=0;  ielem < nelem; ielem++)
    {
		nnode = Get_nnode(vtk, ielem);
        for (int inode = 0; inode < nnode; inode++ )
        {
            ipoil = inpoel[ielem][inode] ;
            *(int*)esup2[ipoil] += 1;
            *(int*)esup2_test[ipoil] += 1;
        }
    }
////////////////////////////////////////////////////////////////////
    for (int ipoin = 1;  ipoin<npoin+1; ipoin++)
    {
        *(int*)esup2[ipoin] += *(int*)esup2[ipoin-1];
        *(int*)esup2_test[ipoin] += *(int*)esup2_test[ipoin-1];
    }
////////////////////////////////////////////////////////////////////
    int size = 0;
    int istor = 0;
    int ipoin = 0;
    // sizing de esup1
    for (int ielem = 0; ielem < nelem; ielem++)
    {
		nnode = Get_nnode(vtk, ielem);
        for (int inode = 0; inode < nnode; inode++)
        {
            ipoin = inpoel[ielem][inode] - 1;
            istor = *(int*)esup2_test[ipoin] + 1;
            *(int*)esup2_test[ipoin]=istor;
            //cout << (ielem + 1);cout << "   ";
            //cout << (istor - 1);cout << "\n";
            if (size < *(int*)esup2_test[ipoin])
            {
                size = (*(int*)esup2_test[ipoin]);
            }
        }
    }
    //cout << "size : ";cout << size;cout << "\n";
    int **esup1 = mat.generateMatrix(size, 1); 
    for (int ielem = 0; ielem < nelem; ielem++)
    {
		nnode = Get_nnode(vtk, ielem);
        for (int inode = 0; inode < nnode; inode++)
        {
            ipoin = inpoel[ielem][inode] - 1;
            istor = *(int*)esup2[ipoin] + 1;
            *(int*)esup2[ipoin]=istor;
            *(int*)esup1[istor-1]=ielem + 1;
            
        }
    }
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
    for (int ipoin = npoin; ipoin >= 1; ipoin--)
    {
        *(int*)esup2[ipoin] = *(int*)esup2[ipoin-1];
    } 
    *(int*)esup2[0] = 0;
////////////////////////////////////////////////////////////////////
    int **lpoin = mat.generateMatrix(npoin, 1); 
    int **psup2 = mat.generateMatrix(npoin+1, 1); 
    istor = 0;
    *(int*)psup2[0]=0;
    size = 0;
    int ielem, jpoin;
    for (int ipoin = 1 ; ipoin <= npoin; ipoin++)
    {
        for (int iesup = *(int*)esup2[ipoin - 1] + 1; iesup <= *(int*)esup2[ipoin + 1 - 1]; iesup++ )
        {
            ielem = *(int*)esup1[iesup-1];
			nnode = Get_nnode(vtk, ielem-1);
            for (int inode = 1; inode <= nnode; inode++ )
            {
                jpoin = inpoel[ielem - 1][inode -1 ];
                if ( (jpoin != ipoin) && (*(int*)lpoin[jpoin-1] != ipoin) )
                {
                    istor += 1;
                    if (size < istor)
                    {
                        size = istor;
                    }
                }
            }
        }
    }
    istor = 0;
    int **psup1 = mat.generateMatrix(size, 1); 
    //cout << "psup2 : ";cout <<"\n";
    for (int ipoin = 1 ; ipoin <= npoin; ipoin++)
    {
        for (int iesup = *(int*)esup2[ipoin-1]+1; iesup <= *(int*)esup2[ipoin+1 - 1]; iesup++ )
        {
            ielem = *(int*)esup1[iesup-1];
			nnode = Get_nnode(vtk, ielem-1);
            for (int inode = 1; inode <= nnode; inode++ )
            {
                jpoin = inpoel[ielem - 1][inode -1 ];
                if ( (jpoin != ipoin) && (*(int*)lpoin[jpoin-1] != ipoin) )
                {
                    istor += 1;
                    *(int*)psup1[istor - 1] = jpoin;
                    *(int*)lpoin[jpoin-1] = ipoin;
                }
            }
        }
        *(int*)psup2[ipoin] = istor;
        //cout << *(int*)psup2[ipoin];cout <<"\n";
    }
////////////////////////////////////////////////////////////////////
    //------------------------------------------------------------------------------
    // ÉLÉMENTS - ÉLÉMENTS

    int **rpoin = mat.generateMatrix(npoin, 1); 

    int **esuel = mat.generateMatrix(nelem, nfaelmax); 
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		for (int ifael = 0; ifael < nfaelmax; ifael++)
		{
			esuel[ielem][ifael] = 0;
		}
	}

    int nnofa, ifael, jelem, nnofj, jfael, jnofa, icoun, inofa, ypoin;
    int **lhelp = mat.generateMatrix(poinperFace,1); 

    for (ielem = 0; ielem < nelem ; ielem++)
    {
		int** ilnofa = Get_lnofa(vtk, ielem);
		int infael = Get_nfael(vtk, ielem);
        //cout << " ..... ";cout << "\n";
        for (ifael = 0 ; ifael < infael ; ifael++)
        {
            //cout << " ..... ";cout << "\n";
			int** ilpofa = Get_lpofa(vtk, ielem);
            nnofa =  *(int*)ilnofa[ifael];
            //cout << " nnofa : ";cout << nnofa;cout << "\n";
            // lhelp(1:nnofa)=inpoel(lpofa(1:nnofa,ifael),ielem)
            for (inofa = 0 ; inofa < nnofa; inofa++)
            {
                //cout << " ..... ";cout << "\n";
                *(int*)lhelp[inofa] = inpoel[ielem][ilpofa[inofa][ifael]-1];
                *(int*)rpoin[*(int*)lhelp[inofa] - 1] = 1;
                ipoin = *(int*)lhelp[0];
            }
            for (istor =  *(int*)esup2[ipoin - 1]; istor <  *(int*)esup2[ipoin]  ; istor++)
            {
                //cout << " ..... ";cout << "\n";
                jelem =  *(int*)esup1[istor] - 1;
                if (jelem != ielem)
                {
					int** jlnofa = Get_lnofa(vtk, jelem);
					int jnfael = Get_nfael(vtk, jelem);
                    for (jfael = 0 ; jfael < jnfael; jfael++)
                    {
						int** jlpofa = Get_lpofa(vtk, jelem);
                        nnofj =  *(int*)jlnofa[jfael];
                        if (nnofj == nnofa)
                        {
                            icoun = 0;
                            for (jnofa = 0 ; jnofa < nnofa; jnofa++)
                            {
                                jpoin = inpoel[jelem][jlpofa[jnofa][jfael]-1];
                                icoun +=  *(int*)rpoin[jpoin-1];
                            }
                            if (icoun == nnofa)
                            {
                                esuel[ielem][ifael] = jelem + 1;
                            }
                        }
                    }
                }
            }
            for (ypoin = 0 ; ypoin < nnofa ; ypoin++ )
            {
                *(int*)rpoin[*(int*)lhelp[ypoin] - 1] = 0;
            }
        }
    }
    
////////////////////////////////////////////////////////////////////
    //------------------------------------------------------------------------------
    // FACE2ELEM, ELEM2FACE, FACE2NODE

    int faceIndex = 0;
    FaceNumber = 0;
    int yelem,yfael;
    int numface = 0;

    int **elem2face = mat.generateMatrix(nelem,nfaelmax); 
    
    for (int ielem = 0 ; ielem <= nelem - 1; ielem++)
    {
		nfael = Get_nfael(vtk, ielem);
        for (int iface = 0; iface < nfael; iface++)
        {
            jelem = esuel[ielem][iface];
            if (ielem < jelem)
            {
                FaceNumber += 1;
            }
            if (jelem == 0)
            {
                FaceNumber += 1;
            }
        }
    }
    int **face2elem = mat.generateMatrix(FaceNumber,2); 

    int **face2node = mat.generateMatrix(FaceNumber,2); 

    int **mmIndex = mat.generateMatrix(nelem,nelem); 

    for (int ielem = 0 ; ielem < nelem; ielem++)
    {
		int** lnofa = Get_lnofa(vtk, ielem);
		int** lpofa = Get_lpofa(vtk, ielem);
		nfael = Get_nfael(vtk, ielem);
        for (int iface = 0; iface < nfael; iface++)
        {
            jelem = esuel[ielem][iface] ;
            if (ielem < jelem-1)
            {
                face2elem[faceIndex][0] = ielem+1;
                face2elem[faceIndex][1] = jelem;
                elem2face[ielem][iface] = faceIndex+1;

                face2node[numface][0]=inpoel[ielem][lpofa[0][iface]-1];
                face2node[numface][1]=inpoel[ielem][lpofa[1][iface]-1];
                numface += 1;

                for (yelem = 0; yelem <= nelem -1 ; yelem++)
                {
					int ynfael = Get_nfael(vtk, yelem);
                    for (yfael = 0; yfael <= ynfael - 1; yfael++)
                    {
                        if ((esuel[yelem][yfael]-1 == ielem) && (yelem==jelem-1))
                        {
                            mmIndex[esuel[yelem][yfael]-1][yelem] = faceIndex+1;
                            mmIndex[yelem][esuel[yelem][yfael]-1] = faceIndex+1;
                        }
                    }
                }
                faceIndex += 1;
            }
            if (jelem == 0)
            {
                face2elem[faceIndex][0] = ielem + 1;
                face2elem[faceIndex][1] = -1;
                elem2face[ielem][iface] = faceIndex+1;
                faceIndex += 1;

                face2node[numface][0]=inpoel[ielem][lpofa[0][iface]-1];
                face2node[numface][1]=inpoel[ielem][lpofa[1][iface]-1];
                numface += 1;
            }
            if (jelem != 0 && ielem >= jelem - 1)
            {
                elem2face[ielem][iface] = mmIndex[jelem-1][ielem];
            }
        }
    }
    if (choix == "esuel")
    {
        return esuel;
    }
    else if (choix == "elem2face")
    {
        return elem2face;
    }
    else if (choix == "face2elem")
    {
        return face2elem;
    }
    else if (choix == "face2node")
    {
        return face2node;
    }
    return 0;

}



double ** connect::Face2Vec(int **elem2face, int **face2node, double **coord) 
{
    matrix mat;
    double **Face2Vec = mat.generateMatrix_double(FaceNumber, nnodemax); 
    for(int i = 0; i < nelem; i++) 
    {
		nfael = Get_nfael(vtk, i);
        for(int k = 0; k < nfael; k++)
        {
            for(int j = 0; j < ddl; j++)
            {
                Face2Vec[elem2face[i][k] - 1][j] = coord[face2node[elem2face[i][k] - 1][1] -1][j] - coord[face2node[elem2face[i][k] - 1][0] - 1][j];
            }
        }
    }
    return Face2Vec;
}


double ** connect::Elem2Vec_x(unsigned **Elem2Node,double **coord)
{
    nfael = 3; // aussi = Nbr face par elem
    matrix mat;
    double **Elem2Vec_x = mat.generateMatrix_double(nelem, nfaelmax); 
    for(int i = 0; i < nelem; i++) 
    {
		nfael = Get_nfael(vtk, i);
		if (nfael == 1) {
			Elem2Vec_x[i][0] = coord[Elem2Node[i][1] - 1][0] - coord[Elem2Node[i][0] - 1][0];
		}
		else {
			for (int k = 0; k < nfael; k++)
			{
				if (k == nfael - 1)
				{
					Elem2Vec_x[i][k] = coord[Elem2Node[i][0] - 1][0] - coord[Elem2Node[i][k] - 1][0];
				}
				else
				{
					Elem2Vec_x[i][k] = coord[Elem2Node[i][k + 1] - 1][0] - coord[Elem2Node[i][k] - 1][0];
				}
			}
		}
    }
    return Elem2Vec_x;
}

void connect::Face2Area_compute(double **Face2Vec)
{
    matrix mat;
    Face2Area = mat.generateMatrix_double(FaceNumber, 1); 
    for(int i = 0; i < FaceNumber; i++) 
    {
        
        Face2Area[i][0] = sqrt(Face2Vec[i][0]*Face2Vec[i][0] + Face2Vec[i][1]*Face2Vec[i][1]);
    }
}


double ** connect::Elem2Vec_y(unsigned **Elem2Node,double **coord)
{
    matrix mat;
    double **Elem2Vec_y = mat.generateMatrix_double(nelem, nfaelmax); 
    for(int i = 0; i < nelem; i++) 
    {
		nfael = Get_nfael(vtk, i);
		if (nfael == 1) {
			Elem2Vec_y[i][0] = coord[Elem2Node[i][1] - 1][1] - coord[Elem2Node[i][0] - 1][1];
		}
		else {
			for (int k = 0; k < nfael; k++)
			{
				if (k == nfael - 1)
				{
					Elem2Vec_y[i][k] = coord[Elem2Node[i][0] - 1][1] - coord[Elem2Node[i][k] - 1][1];
				}
				else
				{
					Elem2Vec_y[i][k] = coord[Elem2Node[i][k + 1] - 1][1] - coord[Elem2Node[i][k] - 1][1];
				}
			}
		}
    }
    return Elem2Vec_y;
}


double ** connect::Elem2Area(unsigned **Elem2Node, double **coord)
{
    matrix mat;
    double **Elem2Area = mat.generateMatrix_double(nelem, 1); 
    for(int i = 0; i < nelem; i++) 
    {
		double sumnode=0;
		nnode = Get_nnode(vtk, i);
		// Area definition according to https://www.mathopenref.com/coordpolygonarea.html
		for (int j = 0; j < nnode; j++) {
			if (j != nnode - 1) {
				sumnode += coord[Elem2Node[i][j] - 1][0] * coord[Elem2Node[i][j + 1] - 1][1] - coord[Elem2Node[i][j] - 1][1] * coord[Elem2Node[i][j + 1] - 1][0];
			}
			else {
				sumnode += coord[Elem2Node[i][j] - 1][0] * coord[Elem2Node[i][0] - 1][1] - coord[Elem2Node[i][j] - 1][1] * coord[Elem2Node[i][0] - 1][0];
			}
		}
            double dArea = fabs(sumnode/2.0);
            Elem2Area[i][0] = dArea;
    }
    return Elem2Area;
}


double ** connect::Elem2Normal(unsigned **Elem2Node,double **coord,double **Elem2Vec_x,double **Elem2Vec_y, string choix)
{
    nfael = 3; // aussi = Nbr face par elem
    matrix mat;
    double **Elem2Normal_x = mat.generateMatrix_double(nelem, nfaelmax);
    double **Elem2Normal_y = mat.generateMatrix_double(nelem, nfaelmax); 

    for(int i = 0; i < nelem; i++)
    {   
		nfael = Get_nfael(vtk, i);
        for(int k = 0; k < nfael; k++)
        {
            double Y = 1.0;
            double X = 1.0;
            /*
            if (Elem2Vec_y[i][k] == 0)
            {
                Y = 1.0;
                X = 0.0;
            }
            else
            {
                Y = -(Elem2Vec_x[i][k]*X)/Elem2Vec_y[i][k];
            }
            double norme = sqrt(X*X + Y*Y);
            Elem2Normal_x[i][k] = X/norme;
            Elem2Normal_y[i][k] = Y/norme;
            */

            
            long double phi = -	1.57079632679489661923;
            X = Elem2Vec_x[i][k]*cos(phi) - Elem2Vec_y[i][k]*sin(phi);
            Y = Elem2Vec_x[i][k]*sin(phi) + Elem2Vec_y[i][k]*cos(phi);

            double norme = sqrt(X*X + Y*Y);
            Elem2Normal_x[i][k] = X/norme;
            Elem2Normal_y[i][k] = Y/norme;
            
        }
    }
    
    if (choix == "x")
    {
        return Elem2Normal_x;
    }
    else if (choix == "y")
    {
        return Elem2Normal_y;
    }
    return 0;
}

double ** connect::Face2Normal(double **Elem2Normal_x, double **Elem2Normal_y, int **elem2face, string choix)
{
    matrix mat;
    //cout << "\nFaceNumber";cout << FaceNumber;cout << "\n";  OKK
    double ** Face2Normal_u = mat.generateMatrix_double(FaceNumber, 1);
    double ** Face2Normal_v = mat.generateMatrix_double(FaceNumber, 1);
    //mat.printMatrix_double(Face2Normal_u,FaceNumber,1);
    for(int i=0; i<nelem; i++)
    {   
        NbNdPerElem = Get_nfael(vtk, i);
        for (int j=0; j<NbNdPerElem; j++)
        {
            int Elem_Number = i;
            int Face_Number = elem2face[i][j] - 1;
            if (Face2Normal_u[Face_Number][0] == 0)
            {
                Face2Normal_u[Face_Number][0] = Elem2Normal_x[Elem_Number][j];
                Face2Normal_v[Face_Number][0] = Elem2Normal_y[Elem_Number][j];
            }
        }
    }
    if (choix == "x")
    {
        return Face2Normal_u;
    }
    else if (choix == "y")
    {
        return Face2Normal_v;
    }
    return 0;
}



