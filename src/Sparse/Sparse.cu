#include "Sparse.h"

Sparse::Sparse() {};

Sparse::Sparse(unsigned int x_size, unsigned int rowSize, unsigned int columnSize)
  : valueSize(x_size), numberOfRows(rowSize), numberOfColumns(columnSize)
{
  Log::Logger().Info("Sparse created");
  cudaMallocManaged(&i, valueSize*sizeof(unsigned int));
  cudaMallocManaged(&j, valueSize*sizeof(unsigned int));  
  cudaMallocManaged(&value, valueSize*sizeof(double));
  cudaMemset(j,0,valueSize*sizeof(unsigned int));
  cudaMemset(i,0,valueSize*sizeof(unsigned int));
  cudaMemset(value, 0 , valueSize*sizeof(double));
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
}

// This function return a CSC matrix
void Sparse::STLAssemble(double tol)
// coo -> a triplet sparse format
{
  type = 1; // type is CSC
  unsigned int * a_temp_i;  unsigned int * a_temp_j; double * a_temp_value;
  // build the indices vector
  unsigned int* indices = new unsigned int[valueSize];
  for ( unsigned int c = 0; c < valueSize; c++) { indices[c] = c;};
  // timer
    Timer *timer = new Timer("Time spend in sorting for STL Assembley: ");
  // sorting i
  std::sort(indices,indices+valueSize, sort_indices(i));
  // number on non-zeros in each row
  unsigned int* nnz_inRow = new unsigned int[numberOfRows+1](); 
  unsigned int rowCounter = 0;
  for (unsigned int c = 0; c < valueSize; c++) {
    if (i[indices[c]] == rowCounter) {
      nnz_inRow[rowCounter] = c+1;
    } else {
      rowCounter++;
    }
  }
  // sorting j 
  cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(unsigned int));
  a_temp_i[0] = 0;// rowPtr
  for (unsigned int c = 1; c <= numberOfRows; c++) {
    std::sort(indices + nnz_inRow[c-1], indices + nnz_inRow[c], sort_indices_j(j,value));
    // -- remove the zero entries (counting the number of nnz in each row)
    a_temp_i[c] = 0;
    for (unsigned int count = nnz_inRow[c-1]; count < nnz_inRow[c]; count++) 
      if (abs(value[indices[count]]) > tol) {a_temp_i[c]++;}
    a_temp_i[c] += a_temp_i[c-1];
  }
  delete timer; Timer t("Time for rest of STL Assembley: ");
  // copy to new array
  cudaMallocManaged(&a_temp_j, (a_temp_i[numberOfRows])*sizeof(unsigned int)); 
  cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  unsigned int vs = 0; // new value size
  for (unsigned int c = 0; c < valueSize; c++) {
    if (abs(value[indices[c]]) > tol) {
      //a_temp_i[vs] = i[indices[c]];
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
  delete[] nnz_inRow;
}

// This function return a CSC matrix using thrust algorithm
void Sparse::ThrustAssemble(double tol)
// coo -> a triplet sparse format
{
  cudaDeviceSynchronize();
  type = 1; // type is CSC
  unsigned int * a_temp_i;  unsigned int * a_temp_j; double * a_temp_value;
  // build the indices vector
  unsigned int* indices;  cudaMallocManaged(&indices, valueSize*sizeof(unsigned int));
  for ( unsigned int c = 0; c < valueSize; c++) { indices[c] = c;};
  // timer
    Timer *timer = new Timer("Time spend in sorting for Thrust Assembley: ");
  // sorting i
  thrust::sort(thrust::device,indices,indices+valueSize, sort_indices(i));
  // number on non-zeros in each row
  unsigned int* nnz_inRow;
  cudaMallocManaged(&nnz_inRow, (numberOfRows+1)*sizeof(unsigned int));
  cudaMemset(nnz_inRow, 0, (numberOfRows+1)*sizeof(unsigned int));
  unsigned int rowCounter = 0;
  for (unsigned int c = 0; c < valueSize; c++) {
    if (i[indices[c]] == rowCounter) {
      nnz_inRow[rowCounter] = c+1;
    } else {
      rowCounter++;
    }
  }
  // sorting j 
  cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(unsigned int));
  a_temp_i[0] = 0;// rowPtr
  for (unsigned int c = 1; c <= numberOfRows; c++) {
    thrust::sort(indices + nnz_inRow[c-1], indices + nnz_inRow[c], sort_indices_j(j,value));
    // -- remove the zero entries (counting the number of nnz in each row)
    a_temp_i[c] = 0;
    for (unsigned int count = nnz_inRow[c-1]; count < nnz_inRow[c]; count++) 
      if (abs(value[indices[count]]) > tol) {a_temp_i[c]++;}
    a_temp_i[c] += a_temp_i[c-1];
  }
  delete timer; Timer t("Time for rest of  Assembley: ");
  // copy to new array
  cudaMallocManaged(&a_temp_j, (a_temp_i[numberOfRows])*sizeof(unsigned int)); 
  cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  unsigned int vs = 0; // new value size
  for (unsigned int c = 0; c < valueSize; c++) {
    if (abs(value[indices[c]]) > tol) {
      //a_temp_i[vs] = i[indices[c]];
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
  cudaFree(indices);
  cudaFree(nnz_inRow);
}

double Sparse::Compare(Sparse& A, Sparse& B)
// return the max Error between the two matrices
{
  double MaxError = 0.00;
  for (unsigned int counter = 0; counter < A.valueSize; counter++)
    MaxError = fmax(MaxError, fabs(A.value[counter] - B.value[counter]));
  return MaxError;
}
// -- new assemblers
void Sparse::STLAssemble2(double tol) {
  
  type = 1; // type is CSC
  // build the indices vector
  unsigned int* indices;  cudaMallocManaged(&indices, valueSize*sizeof(unsigned int));
  for ( unsigned int c = 0; c < valueSize; c++) { indices[c] = c;};
  // timer
  Timer *timer = new Timer("Time spend in sorting for STL Assembley2: ");
  // sorting i
  std::sort(indices,indices+valueSize, newSort(i,j,numberOfRows));
  //thrust::sort(thrust::device,indices,indices+valueSize, newSort(i,j,numberOfRows));
  delete timer; timer = new Timer("Time for rest of STL Assembley2: ");
  // -- define rowPtr
  unsigned int* a_temp_i; cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(unsigned int));
  cudaMemset(a_temp_i,0,(numberOfRows+1)*sizeof(unsigned int));
  // -- find how many zeros are present
  size_t zeroCounter = 0;
  while (!i[indices[zeroCounter]]) {zeroCounter++;}
  // -- add duplicates and creat the rowPtr
  size_t count;
  for (unsigned int c = 0; c < valueSize-1; c++) {
    count = c+1;
    while (i[indices[c]] == i[indices[count]] && j[indices[c]] == j[indices[count]]) {
      value[indices[c]] += value[indices[count]];
      i[indices[count]] = 0; j[indices[count]] = 0; value[indices[count]] = 0;
	count++;
    }
    if (abs(value[indices[c]]) > tol) {
      a_temp_i[i[indices[c]]]++;
    } 
    c = count-1;
  }
  if (abs(value[indices[valueSize-1]]) > tol) {a_temp_i[i[indices[valueSize-1]]]++;}
  // -- print to screen
  //for (unsigned int c = 0; c < valueSize; c++)
  //  std::cout<< c << '\t' << i[indices[c]] << '\t' << j[indices[c]] << '\t' << value[indices[c]] << '\n';
  // -- sum rowPtr
  for (unsigned int c = 1; c<=numberOfRows; c++)
    a_temp_i[c] += a_temp_i[c-1];
  delete timer; Timer t("Time for copy of STL Assembley2: ");
  // copy to new variables
  unsigned int* a_temp_j; cudaMallocManaged(&a_temp_j,     (a_temp_i[numberOfRows])*sizeof(unsigned int)); 
  double*   a_temp_value; cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  unsigned int valueCounter = 0;
  for (unsigned int c = 0; c < valueSize; c++) {
    if(abs(value[indices[c]]) > tol) {
      a_temp_j[valueCounter] = j[indices[c]]-1;
      a_temp_value[valueCounter] = value[indices[c]];
      valueCounter++;
    }
  }
  // change variables
  valueSize = a_temp_i[numberOfRows];
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;  
  cudaFree(indices);
}

// -- new assemblers thrust
void Sparse::ThrustAssemble2(double tol) {
  type = 1; // type is CSC
  // build the indices vector
  unsigned int* indices;  cudaMallocManaged(&indices, valueSize*sizeof(unsigned int));
  for ( unsigned int c = 0; c < valueSize; c++) { indices[c] = c;};
  // timer
  Timer *timer = new Timer("Time spend in sorting for STL Assembley2: ");
  // sorting i
  //std::sort(indices,indices+valueSize, newSort(i,j,numberOfRows));
  cudaDeviceSynchronize();
  thrust::sort(thrust::device,indices,indices+valueSize, newSort(i,j,numberOfRows));
  delete timer; Timer t("Time for rest of STL Assembley2: ");
  // -- define rowPtr
  unsigned int* a_temp_i; cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(unsigned int));
  cudaMemset(a_temp_i,0,(numberOfRows+1)*sizeof(unsigned int));
  // -- find how many zeros are present
  size_t zeroCounter = 0;
  while (!i[indices[zeroCounter]]) {zeroCounter++;}
  // -- add duplicates and creat the rowPtr
  size_t count;
  for (unsigned int c = zeroCounter; c < valueSize-1; c++) {
    count = c+1;
    while (i[indices[c]] == i[indices[count]] && j[indices[c]] == j[indices[count]]) {
      value[indices[c]] += value[indices[count]];
      i[indices[count]] = 0; j[indices[count]] = 0; value[indices[count]] = 0;
	count++;
    }
    if (abs(value[indices[c]]) > tol) {
      a_temp_i[i[indices[c]]]++;
    } else {
      i[indices[c]] = 0; 
    }
    c = count-1;
  }
  if (abs(value[indices[valueSize-1]]) > tol) {a_temp_i[i[indices[valueSize-1]]]++;}
  // -- sum rowPtr
  for (unsigned int c = 1; c<=numberOfRows; c++)
    a_temp_i[c] += a_temp_i[c-1];
  // copy to new variables
  unsigned int* a_temp_j; cudaMallocManaged(&a_temp_j,     (a_temp_i[numberOfRows])*sizeof(unsigned int)); 
  double*   a_temp_value; cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  size_t valueCounter = 0;
  for (unsigned int c = zeroCounter; c < valueSize; c++) {
    if(i[indices[c]]) {
      a_temp_j[valueCounter] = j[indices[c]]-1;
      a_temp_value[valueCounter] = value[indices[c]];
      valueCounter++;
    }
  }
  // change variables
  valueSize = a_temp_i[numberOfRows];
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;  
  cudaFree(indices);
}

// Setters
void Sparse::set_i(unsigned int* index_i) { i = index_i;}
void Sparse::set_j(unsigned int* index_j) { j = index_j;}
void Sparse::set_x(double* x) { value = x;}
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
double * Sparse::get_value() const {return value;}
unsigned int * Sparse::get_i() const {return i;}
unsigned int * Sparse::get_j() const {return j;}
unsigned int Sparse::get_type() const {return type;}

// -- printers
void Sparse::Print() {
  std::cout << "\033[1;34m[SparseMatrix]: \033[0m" << numberOfRows  << " x " << numberOfColumns << " ";
  if (symmetry)
    std::cout << "\033[32msymmetry\033[0m " ;
  unsigned int size = valueSize <72 ? valueSize : 72; 
  std::cout << "print size: " << size << ", nnz = " << valueSize << std::endl;
  for (unsigned int counter = 0; counter < size; counter++)
    std::cout << i[counter] << "\t"<< j[counter] <<"\t: " << value[counter] <<std::endl;
};

unsigned int Sparse::Printer(Sparse& s) {
  s.Print();
  return 0;
};

// -- override the cout << oprator 
std::ostream& operator<< (std::ostream &out, Sparse const& sp) {
  const double* x = sp.get_value();
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

// ---------------------- sort_indices struct -------------------------
sort_indices::sort_indices(unsigned int* var)
  : dofSorted(var) {};

bool sort_indices::operator()(unsigned int i, unsigned int j) const { return dofSorted[i] < dofSorted[j];};

// ---------------------- sort_indices_j struct ------------------------- 

sort_indices_j::sort_indices_j(unsigned int* j, double* value)
  : dofSorted(j), x(value) { };

bool sort_indices_j::operator()(unsigned int i, unsigned int j) {
  if (dofSorted[i] == dofSorted[j]) {
    x[i] = x[i] + x[j]; x[j] = 0;
    return false;
  } else {return dofSorted[i] < dofSorted[j];}
};

// ----------------------- new Sort struct ---------------------
newSort::newSort(unsigned int* i_index, unsigned int* j_index, unsigned int numberOfRows)
  : i(i_index), j(j_index), nRow(numberOfRows) {
};

bool newSort::operator()(unsigned int first, unsigned int second) {
  // compares the actual place in big stiffness matrix
  a = (unsigned long long int)(i[first]) *nRow + (unsigned long long int)(j[first]);
  b = (unsigned long long int)(i[second])*nRow + (unsigned long long int)(j[second]);
  return a < b;
}
