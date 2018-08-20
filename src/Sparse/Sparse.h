#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <algorithm> 
#include <thrust/sort.h>

#include "../Log/Log.h"
#include "../Timer/Timer.h"

//macro
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif


//#include <cuda_runtime_api.h>
//#include <cuda.h>
//#include  <cusolverDn.h>


class Sparse
{
  //variables
public:
  double* value;
  unsigned int* i;
  unsigned int* j;
  unsigned int valueSize; //how many nonzeros
  unsigned int type;
private:
  unsigned int numberOfRows; //number of rows
  unsigned int numberOfColumns; //number of columns
  bool symmetry = false; //by default matrix is not symmetry
  //cuda inverse related
  
  
  // Methods
public:
  Sparse(unsigned int, unsigned int, unsigned int);
  Sparse(unsigned int, unsigned int);
  Sparse();
  ~Sparse();
  void STLAssemble(double);
  void STLAssemble2(double);
  void ThrustAssemble(double);
  void ThrustAssemble2(double);
  
  void set_numberOfRows(unsigned int);
  void set_numberOfColumns(unsigned int);
  void set_valueSize(unsigned int);
  void set_i(unsigned int*);
  void set_j(unsigned int*);
  void set_x(double*);
  unsigned int get_valueSize() const;
  unsigned int get_numberOfRows() const;
  unsigned int get_numberOfColumns() const;
  double * get_value() const;
  unsigned int * get_i() const;
  unsigned int * get_j() const;
  unsigned int get_type() const;
  
private:
  void Print();

  // staic methods
public:
  static unsigned int Printer(Sparse & s);
  static double Compare(Sparse& , Sparse&);
};

std::ostream& operator<< (std::ostream &, Sparse const& );



struct sort_indices {
private:
  unsigned int* dofSorted; // the variable that index is going to sorted base of
public:
  CUDA_HOSTDEV sort_indices(unsigned int*);
  CUDA_HOSTDEV bool operator()(unsigned int i, unsigned int j) const;
};

struct sort_indices_j {
private:
  double* x; // holds values of stiffness matrix
  unsigned int* dofSorted; // the variable that index is going to sorted base of
public:
  CUDA_HOSTDEV sort_indices_j(unsigned int*, double*); 
  CUDA_HOSTDEV bool operator()(unsigned int i, unsigned int j);
};


// -- new Sort
struct newSort {
  unsigned int nRow;
  unsigned int* i; unsigned int* j;
  unsigned long long int a;
  unsigned long long int b;
  CUDA_HOSTDEV newSort(unsigned int*, unsigned int*, unsigned int);
  CUDA_HOSTDEV bool operator()(unsigned int, unsigned int);
};
#endif
