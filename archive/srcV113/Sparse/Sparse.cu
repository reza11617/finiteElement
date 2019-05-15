#include "Sparse.h"

Sparse::Sparse() {};

Sparse::Sparse(unsigned int x_size, unsigned int rowSize, unsigned int columnSize)
  : valueSize(x_size), numberOfRows(rowSize), numberOfColumns(columnSize)
{
  Log::Logger().Info("Sparse created");
  cudaMallocManaged(&i, valueSize*sizeof(unsigned int));
  cudaMallocManaged(&j, valueSize*sizeof(unsigned int));  
  cudaMallocManaged(&value, valueSize*sizeof(float));
  cudaMallocManaged(&nnz_inRow, numberOfRows*sizeof(unsigned int));
  cudaMemset(j,0,valueSize*sizeof(unsigned int));
  cudaMemset(i,0,valueSize*sizeof(unsigned int));
  cudaMemset(value, 0 , valueSize*sizeof(float));
}

Sparse::Sparse(unsigned int x_size, unsigned int sizeOfMatrix)
  : Sparse(x_size, sizeOfMatrix, sizeOfMatrix)
{
  symmetry = true;
};

Sparse::~Sparse() {
  Log::Logger().Info("Sparse deleted");
  cudaFree(i);
  cudaFree(j);
  cudaFree(value);
  cudaFree(nnz_inRow);
}

// This function return a CSC matrix
void Sparse::Assemble()
// coo -> a triplet sparse format
{
  unsigned int * a_temp_i;  unsigned int * a_temp_j; float * a_temp_value;
  // build the indices vector
  unsigned int* indices = new unsigned int[valueSize];
  for ( unsigned int c = 0; c < valueSize; c++) { indices[c] = c;};
  // timer
  Timer timer("Time spend in assembley: ");
  // sorting i
  std::sort(indices,indices+valueSize, sort_indices(i));
  // sorting j
  unsigned int a = 0, b =0;
  for (unsigned int c = 0; c < numberOfRows; c++) {
    a = a + nnz_inRow[c]; b = a + nnz_inRow[c+1];
    std::sort(indices + a, indices+b, sort_indices_j(j,value));
  }
  

  // copy to new array
  cudaMallocManaged(&a_temp_i, valueSize*sizeof(unsigned int));
  cudaMallocManaged(&a_temp_j, valueSize*sizeof(unsigned int)); 
  cudaMallocManaged(&a_temp_value,valueSize*sizeof(float));
  unsigned int vs = 0; // new value size
  for (unsigned int c = 0; c < valueSize; c++) {
    //if (abs(value[indices[c]]) > 0.001) {
    a_temp_i[vs] = i[indices[c]];
    a_temp_j[vs] = j[indices[c]];
    a_temp_value[vs] = value[indices[c]];
    vs++;
    //}
  }
  //valueSize = vs;
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;
  
  //
  delete[] indices;
}
  /*  
  Log::Logger().Info(rowPtr,numberOfRows+1); 
  unsigned int counter = 0;
  for (unsigned int r = 1; r < numberOfRows; r++) {
    a = 0;
    for (unsigned int c = rowPtr[r]; c < rowPtr[r+1]; c++) {
      if (abs(value[indices[counter]]) > 0.001) {
	a++;
      }
      counter++;
    }
    rowPtr[r+1] = a;
  }
  // just change the variable for test
  unsigned int tempValueSize = valueSize; 
  //valueSize = rowPtr[numberOfRows];
  cudaDeviceSynchronize;
  */

float Sparse::Compare(Sparse& A, Sparse& B)
// return the max Error between the two matrices
{
  float MaxError = 0.00;
  for (unsigned int counter = 0; counter < A.valueSize; counter++)
    MaxError = fmax(MaxError, fabs(A.value[counter] - B.value[counter]));
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

// Setters
void Sparse::set_i(unsigned int* index_i) { i = index_i;}
void Sparse::set_j(unsigned int* index_j) { j = index_j;}
void Sparse::set_x(float* x) { value = x;}
void Sparse::set_numberOfRows(unsigned int x) { numberOfRows = x;}
void Sparse::set_valueSize(unsigned int x) { valueSize = x;}

void Sparse::set_numberOfColumns(unsigned int x) {
  numberOfColumns = x;
  if (x == 1 ) { // if we are dealing with vecotrs 
    j = new unsigned int[valueSize];
    for (int ci = 0; ci < valueSize; ci++) {
      j[ci] = 1;
    }
  }
}
  
// Geters
unsigned int Sparse::get_valueSize() const {return valueSize;};
unsigned int Sparse::get_numberOfRows() const { return numberOfRows;}
unsigned int Sparse::get_numberOfColumns() const {return numberOfColumns;}
float * Sparse::get_value() const {return value;}
unsigned int * Sparse::get_i() const {return i;}
unsigned int * Sparse::get_j() const {return j;}

// -- override the cout << oprator 
std::ostream& operator<< (std::ostream &out, Sparse const& sp) {
  const float* x = sp.get_value();
  const unsigned int* i = sp.get_i();
  const unsigned int* j = sp.get_j();
  for (int c = 0; c < sp.get_valueSize(); c++) {
    out << '\t' << i[c] << '\t' << j[c] << ':' << '\t' << x[c] << '\n';
  }
  return out;
}

// printers
void Sparse::Print() {
  std::cout << "\033[1;34m[SparseMatrix]: \033[0m" << numberOfRows  << " x " << numberOfColumns << " ";
  if (symmetry)
    std::cout << "\033[32msymmetry\033[0m " ;
  unsigned int size = valueSize <72 ? valueSize : 72; 
  std::cout << "print size: " << size << ", nnz = " << valueSize << std::endl;
  for (unsigned int counter = 0; counter < size; counter++)
    std::cout << i[counter] << "\t"<< j[counter] <<"\t: " << value[counter] <<std::endl;
};

unsigned int Sparse::Printer(Sparse& s)
{
  s.Print();
  return 0;
};

// ---------------------- sort_indices struct -------------------------
sort_indices::sort_indices(unsigned int* var)
  : dofSorted(var) {};

bool sort_indices::operator()(unsigned int i, unsigned int j) const { return dofSorted[i] < dofSorted[j];};

// ---------------------- sort_indices_j struct ------------------------- 

sort_indices_j::sort_indices_j(unsigned int* var, float* value)
  : dofSorted(var), x(value) {};

bool sort_indices_j::operator()(unsigned int i, unsigned int j) {
  if (dofSorted[i] == dofSorted[j]) {
    x[i] = x[i] + x[j]; x[j] = 0;
    return true;
  } else {return dofSorted[i] < dofSorted[j];}
};
















