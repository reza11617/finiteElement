#include "Sparse.h"

Sparse::Sparse(unsigned int* index_i, unsigned int* index_j, float* value_x, unsigned int x_size, unsigned int rows, unsigned int cols )
  : i(index_i), j(index_j), value(value_x), valueSize(x_size), numberOfRows(rows), numberOfColumns(cols)
{
  Log::Logger().Info("Sparse Created");
  if (numberOfRows == numberOfColumns) //considers symmetry if its a qubic matrix 
    symmetry = true;  
};

Sparse::Sparse()
{
  Log::Logger().Info("Sparse Created");
};

Sparse::Sparse(unsigned int x_size, unsigned int sizeOfMatrix, int deviceID)
  : valueSize(x_size), numberOfRows(sizeOfMatrix), numberOfColumns(sizeOfMatrix), device(deviceID)
{
  symmetry = true;
  if (device == 0) //CPU
    {
      Log::Logger().Info("Sparse Created in cpu");
      temp_i = new unsigned int[valueSize];
      temp_j = new unsigned int[valueSize];
      temp_value = new float[valueSize];
    } else if (device == 1) //GPU
    {
      Log::Logger().Info("Sparse Created in GPU");
      cudaMallocManaged(&temp_i, valueSize*sizeof(float));
      cudaMallocManaged(&temp_j, valueSize*sizeof(float));
      cudaMallocManaged(&temp_value, valueSize*sizeof(float));      
    } else
    {
    Log::Logger().Error("Sparse Matrix error: 0 for cpu and 1 for gpu");
    }
  i = temp_i;
  j = temp_j;
  value = temp_value;
};

Sparse::~Sparse()
{
  if (device == 0) //CPU
    {
      Log::Logger().Info("Sparse deleted by CPU");
      
      delete[] i;
      delete[] j;
      delete[] value;
    } else if (device == 1) //GPU
    {
      Log::Logger().Info("Sparse deleted by GPU");
      cudaFree(i);
      cudaFree(j);
      cudaFree(value);
    } 
}

void Sparse::Assemble()
{
  
}

void Sparse::Assemble(unsigned int* dofFree, unsigned int dofFreeSize)
// very poor algorithem that does the job but needs to changedofFreeSize;
{
  // build DOFs
  unsigned int sizeOfAssembeledMatrix = dofFreeSize*(dofFreeSize+1)/2;
  if (device == 0) { // cpu
    a_temp_i = new unsigned int[sizeOfAssembeledMatrix];
    a_temp_j = new unsigned int[sizeOfAssembeledMatrix];
    a_temp_value = new float[sizeOfAssembeledMatrix];
  } else {
    cudaMallocManaged(&a_temp_i, valueSize*sizeof(float));
    cudaMallocManaged(&a_temp_j, valueSize*sizeof(float));
    cudaMallocManaged(&a_temp_value, valueSize*sizeof(float));      
  }
  unsigned int c = 0;
  for (unsigned int ii = 0; ii < dofFreeSize; ii++) {
    for (unsigned int jj = 0; jj <= ii; jj++) {
      a_temp_i[c] = dofFree[ii]; 
      a_temp_j[c++] = dofFree[jj];
    }
  }
  // use DOF free to assemble
  for (unsigned int ii = 0; ii < sizeOfAssembeledMatrix; ii++) {
    a_temp_value[ii] = 0;
    for (unsigned int jj = 0; jj < valueSize; jj++) {
      if (i[jj] == a_temp_i[ii] && j[jj] == a_temp_j[ii])
	a_temp_value[ii] = a_temp_value[ii] + value[jj];
    }
  }
  // delete previous sparse
  if (device == 0) { //CPU 
    delete[] temp_i;
    delete[] temp_j;
    delete[] temp_value;
  } else  { //gpu
    cudaFree(temp_i);
    cudaFree(temp_j);
    cudaFree(temp_value);
  } 
  // build the new sparse
  numberOfRows = dofFreeSize;
  numberOfColumns = dofFreeSize;
  valueSize = sizeOfAssembeledMatrix;
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;
}

void Sparse::iSet(unsigned int* index_i)
{
  i = index_i;
};
void Sparse::jSet(unsigned int* index_j)
{
  j = index_j;
};
void Sparse::xSet(float* x)
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


unsigned int Sparse::Printer(Sparse& s)
{
  s.Print();
  return 0;
};

float Sparse::Compare(Sparse& A, Sparse& B)
// return the max Error between the two matrices
{
  float MaxError = 0.00;
  for (int i = 0; i < A.valueSize; i++)
    MaxError = fmax(MaxError, fabs(A.value[i] - B.value[i]));
  return MaxError;
}
