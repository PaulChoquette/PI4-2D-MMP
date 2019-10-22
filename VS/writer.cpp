//#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "main.h" 
using namespace std;

void Writer::Write_Output(string filename, int nvar,string* varname, unsigned** inpoel, double** coord, double ** vararr) {
	outfile.open(filename);
	if (outfile.is_open()) {
		outfile << "VARIABLES=";
		for (int i = 0; i < nvar; i++) {
			outfile << "\"" << varname[i] << "\"";
			if (i != nvar - 1) {
				outfile << ",";
			}
			else {
				outfile << "\n";
			}
		}
		outfile << "ZONE T=\"Yeethaw\"\n";
		outfile << "Nodes=" << npoin << ", " << "Elements=" << nelem << ", " << "ZONETYPE=FEQuadrilateral\n";
		outfile << "DATAPACKING=BLOCK\n";
		outfile << "VARLOCATION=([3]=CELLCENTERED) \n";
		//write out variable data
		
		for (int ivar = 0; ivar < nvar; ivar++) {
			if (ivar < 2) {
				for (int ipoin = 0; ipoin < npoin; ipoin++) {
					//First two variables are X and Y, thus writing the coord matrix

					outfile << coord[ipoin][ivar] << "\n";

				}
			}
			else {
				for (int ielem = 0; ielem < nelem;ielem++)
				outfile << vararr[ielem][ivar-2] << "\n";
			}
		}
		//write out inpoel connectivity matrix
		for (int ielem = 0; ielem < nelem; ielem++) {
			for (int inode = 0; inode < nnodemax; inode++) {
				outfile << inpoel[ielem][inode] << " ";
			}
			outfile << "\n";
		}
	}

}