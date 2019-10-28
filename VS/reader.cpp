//#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <string> 
#include "main.h" 
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Reader::read_file(string File_Name) 
{
	elem = 0, poin = 0;
	
	OpenFile(File_Name);

	ifstream f (File_Name);
	ndime = 0;
	nelem = 0;
	npoin = 0;
	nmark = 0;


	linen = 0;
	lineNelem =0;
	poinlinen =0;

	if (file.is_open())
	{
		while (getline(f, line))
		{
			const char* cline = line.c_str();
			linen++;
	///////////////////
	//// ndime ////////
	///////////////////
			//If number of dimensions not specified yet, check
			if (ndime == 0)
			{
				ndime = Readndime(line);
			}
	///////////////////
	//// nelem ////////
	///////////////////
			//If number of elements not specified yet, check
			if (nelem == 0)
			{
				nelem = Readnelem(line);
				//If nelem is defined, store line number and allocate memory for elem 2 poin array
				if (nelem != 0) 
				{
	///////////////////
	//// inpoel1 //////
	///////////////////
					lineNelem = linen;
					inpoel1 = matrix.generateMatrix_unsigned(nelem, nnodemax); 
					vtk = matrix.generateMatrix(nelem,1);
					//matrix.printMatrix(inpoel1, nelem,nnode);
					continue;
				}
			}
			//If line describes the nodes of each element, store in elemsv
			if (lineNelem != 0) 
			{
				if (linen > lineNelem && linen <= lineNelem + nelem)
				{
					FillE2P_VTK(cline);
					//Restart line reading loop because line has no other purpose
					continue;
				}
			}
			if (npoin == 0) 
			{
	///////////////////
	//// npoin ////////
	///////////////////
				npoin = Readnpoin(line);
				//If npoin is defined, store line number and allocate memory for elem 2 poin array
				if (npoin) 
				{
					coord = matrix.generateMatrix_double(npoin,ndime);
					poinlinen = linen;
					
					/*
					for (unsigned i = 0; i < nelem; i++) 
					{
						for (unsigned j = 0; j < ndime; j++)
						{
							coord[i][j] = -1;
						}
					}
					*/
					continue;
				}
			}
			if (poinlinen != 0)
			{
				if (linen > poinlinen && linen <= poinlinen + npoin)
				{
					FillCoord(cline);
					//Restart line reading loop because line has no other purpose
					continue;
				}
			}
	///////////////////
	//// nmark ////////
	///////////////////
			if (nmark == 0) {
				nmark = Readnmark(line);
				if (nmark != 0) {
					marklinen = linen;
					marknl = 2 * nmark; //2 lines per marker for marker name and elemn
					markername = new string[nmark];
					markerdata = new unsigned**[nmark];
					melemnv = matrix.generateMatrix(nmark, 1);
					markn = 0;
					step = 0;
					continue;
				}
			}
			if (marklinen != 0) {
				if (linen > marklinen&& linen <= marklinen + marknl) {
					//3 steps : step 0 is reading marker tag, step 1 is reading marker elemn and final step is filling markerdata
					if (step == 0) {
						markername[markn] = FillMarkTag(line);
						step++;
					}
					else if (step == 1) {
						markelemn = Readmarkelemn(line);
						*(int*)melemnv[markn] = markelemn;
						markerdata[markn] = matrix.generateMatrix_unsigned(markelemn,2);  // 2 since we are in 2D and boundaries will always be lines
						marknl += markelemn; //add nelem since each elem has 1 line
						imen = 0; //counter for marker element number used in FillMarker
						step++;
					}
					else if (step == 2) {
						FillMarker(cline, markn);
						if (imen == markelemn) {
							markn++;
							step = 0;
						}
					}
					if (markn == nmark) {
						totalmarken = marknl - 2 * nmark;
						celltype = new string[nelem + totalmarken];
						inpoel1_wf = matrix.generateMatrix_unsigned(nelem + totalmarken, nnodemax);
						vtk_wf = matrix.generateMatrix(nelem + totalmarken, 1);
						for (int ielem = 0; ielem < nelem; ielem++) {
							*(int*)vtk_wf[ielem] = *(int*)vtk[ielem];
							celltype[ielem] = "center";
							for (int inode = 0; inode < nnodemax; inode++) {
								inpoel1_wf[ielem][inode] = inpoel1[ielem][inode];
							}
						}
						int imelem = 0;
						for (int k = 0; k < nmark; k++) {
							for (int ifc = 0; ifc < *(int*)melemnv[k]; ifc++) {
								*(int*)vtk_wf[nelem + imelem] = 3;
								string farfield = "farfield";
								string slipwall = "slipwall";
								if (markername[k] == "farfield") {
									celltype[nelem + imelem] = farfield;
								}
								else if (markername[k] == "airfoil") {
									celltype[nelem + imelem] = slipwall;
								}
								
								for (int j = 0; j < 2; j++) {
									inpoel1_wf[nelem + imelem][j] = markerdata[k][ifc][j];
								}
								imelem++;
							}
						}
					}
				}
			}
			

		}
	}
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////



unsigned Reader::Readndime(const string& line) 
{
	size_t nelempos = line.find("NDIME= ");
	unsigned dim = 0;
	if (nelempos != string::npos) {
		char ndimech;
		const char* pdimech = &ndimech;
		ndimech = line.back();
		dim = atoi(pdimech);
	}
	return dim;
}

unsigned Reader::Readnelem(const string& line) {

	unsigned nel = 0;

	size_t nelempos = line.find("NELEM= ");
	if (nelempos != string::npos) {
		nelempos = line.find_last_of(" ") + 1;

		string nelemch;

		nelemch = line.substr(nelempos, line.length() - nelempos);
		nel = stoul(nelemch);
	}
	return nel;
}

unsigned Reader::Readnmark(const string& line) {

	unsigned short nma = 0;

	size_t nelempos = line.find("NMARK= ");
	if (nelempos != string::npos) {
		nelempos = line.find_last_of(" ") + 1;

		string nmarkch;

		nmarkch = line.substr(nelempos, line.length() - nelempos);
		nma = stoul(nmarkch);
	}
	return nma;
}

unsigned Reader::Readmarkelemn(const string& line) {

	unsigned short men = 0;

	size_t nelempos = line.find("MARKER_ELEMS= ");
	if (nelempos != string::npos) {
		nelempos = line.find_last_of(" ") + 1;

		string nmarkch;

		nmarkch = line.substr(nelempos, line.length() - nelempos);
		men = stoul(nmarkch);
	}
	return men;
}
//	inpoel1 = matrix.generateMatrix(nelem, nnode);
void Reader::FillE2P_VTK(const char* cline)
{
	char* end;
	unsigned j = 0;
	for (unsigned c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;
		//cout << c;cout << "\n";
		if (j == 0) {
			*(int*)vtk[elem] = c;
		}
		else 
		{
			//cout << "j";cout << j;cout << "\n";
			inpoel1[elem][j - 1] = c + 1;
		}
		if (*(int*)vtk[elem] == 5) {
			if (j == 4) {
				inpoel1[elem][j - 1] = inpoel1[elem][0];
			}
		}
		j++;
		
		
	}
	elem++;
}

string Reader::FillMarkTag(const string& line)
{
	string marker_tag;
	size_t mtpos = line.find("MARKER_TAG= ");
	if (mtpos != string::npos) {
		mtpos = line.find_last_of(" ") + 1;
		marker_tag = line.substr(mtpos, line.length() - mtpos);
	}

	return marker_tag;
}


void Reader::FillMarker(const char* cline, int markn)
{
	char* end;
	unsigned j = 0;
	for (unsigned c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;
		if (j != 0) {
			markerdata[markn][imen][j - 1] = c + 1;
		}
		j++;
	}
	imen++;
}

unsigned Reader::Readnpoin(const string& line)
{
	unsigned npo = 0;

	size_t npoinpos = line.find("NPOIN= ");
	if (npoinpos != string::npos) {
		npoinpos = line.find_last_of(" ") + 1;

		string npoinch;

		npoinch = line.substr(npoinpos, line.length() - npoinpos);
		npo = stoul(npoinch);
	}

	return npo;
}

double ** Reader::FillCoord(const char* cline)
{
	char* end;
	unsigned j = 0;
	for (long double c = strtod(cline, &end); cline != end; c = strtod(cline, &end))
	{
		cline = end;
		coord[poin][j] = c;
		j++;
	}
	poin++;

	return coord;
}


bool Reader::IsDefined(int numb)
{
	return !numb;
}

bool Reader::OpenFile(string filename)
{
	file.open(filename);

	if (file.is_open()) {
		return true;
	}
	else {
		return false;
	}
}
