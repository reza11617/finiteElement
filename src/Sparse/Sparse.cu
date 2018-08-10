#include "Sparse.h"

Sparse::Sparse() {};

Sparse::Sparse(unsigned int x_size, unsigned int rowSize, unsigned int columnSize)
  : valueSize(x_size), numberOfRows(rowSize), numberOfColumns(columnSize)
{
  Log::Logger().Info("Sparse created");
  cudaMallocManaged(&i, valueSize*sizeof(unsigned int));
  cudaMallocManaged(&j, valueSize*sizeof(unsigned int));  
  cudaMallocManaged(&value, valueSize*sizeof(float));
  cudaMallocManaged(&nnz_inRow, (numberOfRows+1)*sizeof(unsigned int));
  cudaMemset(j,0,valueSize*sizeof(unsigned int));
  cudaMemset(i,0,valueSize*sizeof(unsigned int));
  cudaMemset(value, 0 , valueSize*sizeof(float));
  type = 0; // means the sparse matrix is COO format
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
void Sparse::Assemble(float tol)
// coo -> a triplet sparse format
{
  type = 1; // type is CSC
  unsigned int * a_temp_i;  unsigned int * a_temp_j; float * a_temp_value;
  // build the indices vector
  unsigned int* indices = new unsigned int[valueSize];
  for ( unsigned int c = 0; c < valueSize; c++) { indices[c] = c;};
  // timer
  Timer timer("Time spend in assembley: ");
  //Log::Logger().Info(nnz_inRow,numberOfRows+1);
  // sorting i
  std::sort(indices,indices+valueSize, sort_indices(i));
  // sorting j 
  cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(unsigned int));
  unsigned int a = nnz_inRow[0], b =a;
  a_temp_i[0] = 0; // rowPtr
  for (unsigned int c = 1; c <= numberOfRows; c++) {
    b = b + nnz_inRow[c];
    std::sort(indices + a, indices+b, sort_indices_j(j,value));
    // -- remove the zero entries (counting the number of nnz in each row)
    a_temp_i[c] = nnz_inRow[c];
    for (unsigned int count = a; count < b; count++) 
      if (abs(value[indices[count]]) < tol) {a_temp_i[c]--;}
    a_temp_i[c] += a_temp_i[c-1];
    //
    a = b;
  }
  
  // copy to new array
  cudaMallocManaged(&a_temp_j, (a_temp_i[numberOfRows])*sizeof(unsigned int)); 
  cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(float));
  unsigned int vs = 0; // new value size
  for (unsigned int c = 0; c < valueSize; c++) {
    if (abs(value[indices[c]]) > tol) {
      a_temp_j[vs] = j[indices[c]] - 1;
      a_temp_value[vs] = value[indices[c]];
      vs++;
    }
  }
  valueSize = a_temp_i[numberOfRows];
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;
  
  delete[] indices;
}

float Sparse::Compare(Sparse& A, Sparse& B)
// return the max Error between the two matrices
{
  float MaxError = 0.00;
  for (unsigned int counter = 0; counter < A.valueSize; counter++)
    MaxError = fmax(MaxError, fabs(A.value[counter] - B.value[counter]));
  return MaxError;
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
unsigned int Sparse::get_type() const {return type;}
// -- override the cout << oprator 
std::ostream& operator<< (std::ostream &out, Sparse const& sp) {
  const float* x = sp.get_value();
  const unsigned int* j = sp.get_j();
  const unsigned int* i = sp.get_i();
  unsigned int iCounter = 0;
  for (int c = 0; c < sp.get_valueSize(); c++) {
    out << '\t' ;
    if (!sp.get_type()) { // if the format is COO
      out << i[c];
	} else if (c == i[iCounter+1]-1) { // if the format is CSC
      out << i[iCounter+1];
      iCounter++;
    }
    out << '\t' << j[c] << ':' << '\t' << x[c] << '\n';
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
















