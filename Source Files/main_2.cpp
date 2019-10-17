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
	FileContents.nnode = 3;

	//ifstream meshfile("naca0012_129x129_1B_JAMESON.su2");
	string File_Name = "square.su2";

	if (FileContents.OpenFile(File_Name))
	{
		FileContents.read_1(File_Name);
		FileContents.file.close();
	}


	cout << "--- INPOEL ---";cout << "\n";
	matrix.printMatrix(FileContents.inpoel1,FileContents.nelem,FileContents.nnode);

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

	return 0;

}