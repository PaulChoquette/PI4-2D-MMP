//#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "main.h" 
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Reader::read_1(string File_Name) 
{
	elem = 0, poin = 0;
	
	OpenFile(File_Name);
	unsigned linen = 0;
	string line;
	ifstream f (File_Name);
	ndime = 0;
	nelem = 0;
	npoin = 0;

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
					inpoel1 = matrix.generateMatrix(nelem, nnode); 
					//matrix.printMatrix(inpoel1, nelem,nnode);
					continue;
				}
			}
			//If line describes the nodes of each element, store in elemsv
			if (lineNelem != 0) 
			{
				if (linen > lineNelem && linen <= lineNelem + nelem)
				{
					FillE2P(cline);
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
//	inpoel1 = matrix.generateMatrix(nelem, nnode);
void Reader::FillE2P(const char* cline)
{
	char* end;
	unsigned j = 0;
	for (unsigned c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;
		//cout << c;cout << "\n";
		if (j != 0) 
		{
			//cout << "j";cout << j;cout << "\n";
			inpoel1[elem][j - 1] = c + 1;
		}
		j++;
	}
	elem++;
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
		cout << c; cout << "\n";
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
