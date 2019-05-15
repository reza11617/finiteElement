#include "Sparse.h"

Sparse::Sparse(unsigned int* index_i, unsigned int* index_j, double* value_x, unsigned int x_size, unsigned int rows, unsigned int cols )
  : i(index_i), j(index_j), value(value_x), valueSize(x_size), numberOfRows(rows), numberOfColumns(cols)
{
  if (numberOfRows == numberOfColumns) //considers symmetry if its a qubic matrix 
    symmetry = true;
  
};

Sparse::Sparse()
{
};

void Sparse::iSet(unsigned int* index_i)
{
  i = index_i;
};
void Sparse::jSet(unsigned int* index_j)
{
  j = index_j;
};
void Sparse::xSet(double* x)
{
  value = x;
};

void Sparse::sizeSet(unsigned int x_size, unsigned int rows, unsigned int cols)
{
  valueSize = x_size;
  numberOfRows = rows;
  numberOfColumns = cols;
  if (numberOfRows == numberOfColumns) //considers symmetry if its a qubic matrix
    symmetry = true;
}

Sparse::~Sparse()
{
  delete[] i;
  delete[] j;
  delete[] value;
}

void Sparse::Print()
{
  std::cout << "\033[1;34m[SparseMatrix]: \033[0m" << numberOfRows  << " x " << numberOfColumns << " ";
  if (symmetry)
    std::cout << "\033[32msymmetry\033[0m " ;
  unsigned int size = valueSize <20 ? valueSize : 20; 
  std::cout << "print size: " << size << std::endl;
  for (unsigned int counter = 0; counter < size; counter++)
    std::cout << i[counter] << "\t"<< j[counter] <<"\t: " << value[counter] <<std::endl;
};


unsigned int Sparse::Printer(Sparse & s)
{
  s.Print();
  return 0;
};
