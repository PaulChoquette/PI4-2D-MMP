
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
        unsigned ** generateMatrix_unsigned(int rows, int cols);
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

		void read_1(string File_Name);
		unsigned ndime, nelem, npoin;
		int nnode; // Nbre noeud par element
		int **inpoel1;
		double **coord;
		unsigned elem, poin;
		int lineNelem;
		int poinlinen;

		bool OpenFile(string filename);
		bool IsDefined(int numb);

		unsigned Readndime(const string& line);

		unsigned Readnelem(const string& line);
		unsigned** AlloE2P(unsigned**& elem2poin, unsigned ndime, unsigned nelem);
		void FillE2P(const char* cline);

		unsigned Readnpoin(const string& line);
		long double** AlloCoord(long double**& oldcoord, unsigned ndime, unsigned npoin);
		double ** FillCoord(const char* cline);


};

class connect 
{
	public:
        // Attributs
        int npoin, nelem, nnode;
        int nfael;
        int poinperFace;
        int FaceNumber;
        int NbNdPerFace;
        int NbNdPerElem;
        int NbFacePerELEM;
        int ddl;
        // Pointeur
        

        // Methods 
        int ** init(int **inpoel,int **lpofa,int **lnofa, string choix);
        int ** Face2Vec(int **elem2face, int **face2node, double **coord);
        int ** Elem2Node(int **elem2face,int **face2node);
        double ** Elem2Vec_x(int **Elem2Node,double **coord);
        double ** Elem2Vec_y(int **Elem2Node,double **coord);
        double ** Elem2Area(int **Elem2Node,double **coord);
        double ** Elem2Normal(int **Elem2Node,double **coord,double **Elem2Vec_x,double **Elem2Vec_Y, string choix);
        

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


