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

Sparse::Sparse(unsigned int x_size, unsigned int rowSize, unsigned int columnSize, int deviceID)
  : valueSize(x_size), numberOfRows(rowSize), numberOfColumns(columnSize), device(deviceID)
{
  if (device == 0) //CPU
    {
      Log::Logger().Info("Sparse Created in cpu");
      temp_i = new unsigned int[valueSize];
      temp_j = new unsigned int[valueSize];
      temp_value = new float[valueSize];
    } else if (device == 1) //GPU
    {
      Log::Logger().Info("Sparse Created in GPU");
      cudaMallocManaged(&temp_i, valueSize*sizeof(unsigned int));
      cudaMallocManaged(&temp_j, valueSize*sizeof(unsigned int));
      cudaMallocManaged(&temp_value, valueSize*sizeof(float));      
    } else
    {
    Log::Logger().Error("Sparse Matrix error: 0 for cpu and 1 for gpu");
    }
  i = temp_i;
  j = temp_j;
  value = temp_value;  
}
Sparse::Sparse(unsigned int x_size, unsigned int sizeOfMatrix, int deviceID)
  : Sparse(x_size, sizeOfMatrix, sizeOfMatrix, deviceID)
{
  symmetry = true;
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
// very poor algorithem that does the job but needs to change
{
  Timer timer("Time spend in CPU for assembling: ");
  // build the new sparse
  numberOfRows = dofFreeSize;
  if (numberOfColumns != 1) // for vectors
    numberOfColumns = dofFreeSize;

  // build DOFs
  unsigned int sizeOfAssembeledMatrix = numberOfRows*numberOfColumns;
  if (device == 0) { // cpu
    a_temp_i = new unsigned int[sizeOfAssembeledMatrix];
    a_temp_j = new unsigned int[sizeOfAssembeledMatrix];
    a_temp_value = new float[sizeOfAssembeledMatrix];
  } else {
    cudaMallocManaged(&a_temp_i, sizeOfAssembeledMatrix*sizeof(unsigned int));
    cudaMallocManaged(&a_temp_j, sizeOfAssembeledMatrix*sizeof(unsigned int));
    cudaMallocManaged(&a_temp_value, sizeOfAssembeledMatrix*sizeof(float));      
  }
  unsigned int dofFreeSize_j = dofFreeSize;
  if (numberOfColumns == 1) // for load vectors
    dofFreeSize_j = 1;
      
  unsigned int c = 0;
  for (unsigned int ii = 0; ii < dofFreeSize; ii++) {
    for (unsigned int jj = 0; jj < dofFreeSize_j; jj++) {
      a_temp_i[c] = dofFree[ii]; 
      if (numberOfColumns != 1) // for vectors
	a_temp_j[c++] = dofFree[jj];
      else
	a_temp_j[c++] = 1; // column should be one due to being vector
    }
  }
  // use DOF free to assemble
  for (unsigned int ii = 0; ii < sizeOfAssembeledMatrix; ii++) {
    a_temp_value[ii] = 0;
    for (unsigned int jj = 0; jj < valueSize; jj++) {
      if (i[jj] == a_temp_i[ii] && j[jj] == a_temp_j[ii])
	a_temp_value[ii] = a_temp_value[ii] + value[jj];
      else if (i[jj] == a_temp_j[ii] && j[jj] == a_temp_i[ii])
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
  valueSize = sizeOfAssembeledMatrix;
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;
}

float Sparse::Compare(Sparse& A, Sparse& B)
// return the max Error between the two matrices
{
  float MaxError = 0.00;
  for (int i = 0; i < A.valueSize; i++)
    MaxError = fmax(MaxError, fabs(A.value[i] - B.value[i]));
  return MaxError;
}

void Sparse::solver(Sparse& matrix, Sparse& vector) {
  Timer timer("Time spend in GPU for inverting: ");
  cudaError cudaStatus;
  cusolverStatus_t  cusolverStatus;
  cusolverDnHandle_t  handle;
  float   *Work;                                         //   workspace
  int   *info , Lwork;                      //   info , workspace  size
  cudaStatus = cudaGetDevice (0);
  cusolverStatus = cusolverDnCreate (& handle ); //  create  handle
  cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;
  cudaMallocManaged (&info ,sizeof(int ));// unified  mem. for  info
  //  compute  workspace  size  and  prepare  workspace
  cusolverStatus = cusolverDnSpotrf_bufferSize(handle , uplo ,matrix.numberOfRows,matrix.value,matrix.numberOfRows,& Lwork );
  cudaMallocManaged (&Work ,Lwork*sizeof(float )); //mem.for  Work
  //  Cholesky  decomposition   d_A=L*L^T, lower  triangle  of d_A is
  //  replaced  by the  factor L
  cusolverStatus = cusolverDnSpotrf(handle,uplo,matrix.numberOfRows,matrix.value,matrix.numberOfRows,Work,Lwork,info);
  cudaStatus = cudaDeviceSynchronize ();
  //  solve A*X=B,   where A is  factorized  by  potrf  function
  // B is  overwritten  by the  solution
  cusolverStatus = cusolverDnSpotrs(handle,uplo,matrix.numberOfRows,1,matrix.value,matrix.numberOfRows,vector.value,matrix.numberOfRows,info);
  cudaStatus = cudaDeviceSynchronize ();
  cudaStatus = cudaFree(info);
  cudaStatus = cudaFree(Work);
  cusolverStatus = cusolverDnDestroy(handle);
  //cudaStatus = cudaDeviceReset();
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

void Sparse::set_numberOfRows(unsigned int x) {
  numberOfRows = x;
}

void Sparse::set_numberOfColumns(unsigned int x) {
  numberOfColumns = x;
  Log::Logger().Info(valueSize);
  if (x == 1 ) { // if we are dealing with vecotrs 
    j = new unsigned int[valueSize];
    for (int ci = 0; ci < valueSize; ci++) {
      j[ci] = 1;
    }
  }
}
  

void Sparse::set_valueSize(unsigned int x) {
  valueSize = x;
}

void Sparse::set_device(int x) {
  device = x;
}

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
  std::cout << "print size: " << size << ", nnz = " << valueSize << std::endl;
  for (unsigned int counter = 0; counter < size; counter++)
    std::cout << i[counter] << "\t"<< j[counter] <<"\t: " << value[counter] <<std::endl;
};


unsigned int Sparse::Printer(Sparse& s)
{
  s.Print();
  return 0;
};
