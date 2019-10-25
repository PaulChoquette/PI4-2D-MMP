
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm> 
#include <cmath> 
#include <functional> // std::divides 
using namespace std;

class matrix 
{
    public: 
        void printMatrix(int ** mat, int rows, int cols);
        void printMatrix_double(double ** mat, int rows, int cols) ;
        void printMatrix_unsigned(unsigned ** mat, int rows, int cols) ;
        void printXY_double(double ** matX, double ** matY, int rows, int cols);
        void deleteMatrix(int ** mat, int rows, int cols);
        int ** generateMatrix(int rows, int cols);
        double ** generateMatrix_double(int rows, int cols);
        unsigned ** generateMatrix_unsigned(unsigned rows, unsigned cols);
};

template <class T>
class Array {
    public:
        Array();
        ~Array();

        void push(T value);
        T pop();
        void set(size_t index, T value);
        T get(size_t index);

        size_t size();
        size_t capacity();
        void print();
        int is_empty();

    private:
        size_t m_size;
        size_t m_capacity;
        T* m_data;
        void resize();
};

class Reader {
	public:
		ifstream file;
		matrix matrix;

		unsigned linen;
		string line;

		void read_file(string File_Name);
		unsigned ndime, nelem, npoin;
		unsigned short nmark;
		unsigned nnode; // Nbre noeud par element
		unsigned nnodemax;
		unsigned **inpoel1;
		unsigned** inpoel1_wf; //inpoel1 with fantom cells
		double **coord;
		int** vtk;
		int** vtk_wf;
		unsigned elem, poin;
		unsigned lineNelem;
		unsigned poinlinen;
		unsigned marklinen;
		unsigned marknl, melem, markn, markelemn;
		int** melemnv;
		string* markername;
		unsigned step;
		unsigned*** markerdata;
		unsigned imen, totalmarken;

		bool OpenFile(string filename);
		bool IsDefined(int numb);

		unsigned Readndime(const string& line);
		unsigned Readnelem(const string& line);
		unsigned Readnmark(const string& line);
		unsigned Readmarkelemn(const string& line);
		void FillE2P_VTK(const char* cline);
		void FillMarker(const char* cline, int markn);
		string FillMarkTag(const string& line);

		unsigned Readnpoin(const string& line);
		double ** FillCoord(const char* cline);


};

class connect 
{
	public:
        // Attributs
        int npoin, nelem, nnode;
        int nfael;
		int nfaelmax,nnodemax;
        int poinperFace;
        int FaceNumber;
        int NbNdPerFace;
        int NbNdPerElem;
        int ddl;
		int** vtk;
        // Pointeur
        

        // Methods 
        int ** init(unsigned **inpoel, string choix);
		int ** Get_lpofa(int** vtk, int ielem);
		int ** Get_lnofa(int** vtk, int ielem);
		int Get_nnode(int** vtk, int ielem);
		int Get_nfael(int** vtk, int ielem);
        int ** Face2Vec(int **elem2face, int **face2node, double **coord);
        int ** Elem2Node(int **elem2face,int **face2node);
        double ** Elem2Vec_x(unsigned **Elem2Node,double **coord);
        double ** Elem2Vec_y(unsigned **Elem2Node,double **coord);
        double ** Elem2Area(unsigned **Elem2Node,double **coord);
        double ** Elem2Normal(unsigned **Elem2Node,double **coord,double **Elem2Vec_x,double **Elem2Vec_Y, string choix);
        

};

class primitive
{
    public:
        // Attributs
        double rho_0, P_0, u_0, v_0;
        int NbrPrimitive;
        // Method
        double ** init_0(double rho_0,double P_0,double u_0,double v_0, int NbrPrimitive, int nelem);



};


class roe
{
    public:
        // Attributs

        // Method
        void compute();



};

class Writer
{
	public:
		//Attributs
		ofstream outfile;
		int npoin, nelem, nnodemax;
		//Methods
		void Write_Output(string filename, int nvar,string* varname, unsigned** inpoel,double ** coord, double ** vararr);
};

