#include "main.h"
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void matrix::printMatrix(int ** mat, int rows, int cols) 
{
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            std::cout << mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
void matrix::printXY_double(double ** matX, double ** matY, int rows, int cols) 
{
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            std::cout << "[ " << matX[i][j] << " " << matY[i][j] << " ]";
        }
        std::cout << std::endl;
    }
}

void matrix::printMatrix_double(double ** mat, int rows, int cols) 
{
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            std::cout << mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void matrix::printMatrix_unsigned(unsigned** mat, int rows, int cols)
{
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			std::cout << mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void matrix::deleteMatrix(int ** mat, int rows, int cols) {
    //delete rows
    for (int i = 0; i < rows; ++i) {
        delete[] mat[i];
    }
    //delete the pointer 
    delete[] mat;    
}




int ** matrix::generateMatrix(int rows, int cols) 
{
    //pointer to pointer, to the element at row 1, column 1
    int ** mat = new int*[rows];
    //for each row, create an array with size equal to number of elements in a column
    for(int i = 0; i < rows; i++) 
    {
        mat[i] = new int[cols];
    }

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			mat[i][j] = 0;
		}
	}
    //return the pointer to pointer
    return mat;
}

double ** matrix::generateMatrix_double(int rows, int cols) 
{
    //pointer to pointer, to the element at row 1, column 1
    double ** mat = new double*[rows];
    //for each row, create an array with size equal to number of elements in a column
    for(int i = 0; i < rows; i++) 
    {
        mat[i] = new double[cols];
    }
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			mat[i][j] = 0.0;
		}
	}
    //return the pointer to pointer
    return mat;
}

unsigned** matrix::generateMatrix_unsigned(unsigned rows, unsigned cols)
{
	//pointer to pointer, to the element at row 1, column 1
	unsigned** mat = new unsigned* [rows];
	//for each row, create an array with size equal to number of elements in a column
	for (int i = 0; i < rows; i++)
	{
		mat[i] = new unsigned[cols];
	}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			mat[i][j] = 0;
		}
	}
	//return the pointer to pointer
	return mat;
}

#define MIN_CAPACITY 10
#define GROWTH_FACTOR 2
 

template <class T>
Array<T>::Array() 
{
    m_capacity = MIN_CAPACITY;
    m_size = 0;
    m_data = (T*)malloc(m_capacity * sizeof(*m_data));
    if (!m_data)
        throw std::bad_alloc(); 
}

template <class T>
Array<T>::~Array() 
{
    free(m_data);
}

template <class T>
void Array<T>::push(T value) 
{
    if (m_size >= m_capacity)
        resize();
    *(m_data + m_size++) = value;
}

template <class T>
T Array<T>::pop() 
{
    return *(m_data + --m_size);
}

template <class T>
void Array<T>::set(size_t index, T value) 
{
    while (index >= m_size)
        this->push(0);
    *(m_data + index) = value;
}


template <class T>
T Array<T>::get(size_t index) {
    return *(m_data + index);
}

template <class T>
size_t Array<T>::size() 
{ 
    return m_size;
}

template <class T>
size_t Array<T>::capacity() 
{
    return m_capacity;
}

template <class T>
void Array<T>::resize() 
{
    size_t capacity = m_capacity*GROWTH_FACTOR;
    T *tmp = (T*)realloc(m_data, capacity * sizeof(*m_data));
    if (!tmp)
        throw std::bad_alloc();
    m_data = tmp;
    m_capacity = capacity;
}

template <class T>
void Array<T>::print() 
{
    size_t i = 0;
    while (i < m_size) {
        std::cout << "Array[" << i << "] = " << this->get(i) << "\n";
        i++;
    }
}

template <class T>
int Array<T>::is_empty() 
{
    return m_size == 0;
}

template class Array<int>;
template class Array<float>;
template class Array<std::string>;


