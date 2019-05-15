#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <algorithm> 

#include "../Log/Log.h"
#include "../Timer/Timer.h"

//#include <cuda_runtime_api.h>
//#include <cuda.h>
#include  <cusolverDn.h>

struct sort_indices {
private:
  unsigned int* dofSorted; // the variable that index is going to sorted base of
public:
  sort_indices(unsigned int*);
  bool operator()(unsigned int i, unsigned int j) const;
};

struct sort_indices_j {
private:
  unsigned int* dofSorted; // the variable that index is going to sorted base of
  float* x; // holds values of stiffness matrix
public:
  sort_indices_j(unsigned int*, float*); 
  bool operator()(unsigned int i, unsigned int j);
};

class Sparse
{
  //variables
public:
  float* value;
  unsigned int* i;
  unsigned int* j;
  unsigned int valueSize; //how many nonzeros
  unsigned int* nnz_inRow; // number of non-zero in each row;
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
  void Assemble();
  void set_numberOfRows(unsigned int);
  void set_numberOfColumns(unsigned int);
  void set_valueSize(unsigned int);
  void set_i(unsigned int*);
  void set_j(unsigned int*);
  void set_x(float*);
  unsigned int get_valueSize() const;
  unsigned int get_numberOfRows() const;
  unsigned int get_numberOfColumns() const;
  float * get_value() const;
  unsigned int * get_i() const;
  unsigned int * get_j() const;
  
  
private:
  void Print();

  // staic methods
public:
  static unsigned int Printer(Sparse & s);
  static float Compare(Sparse& , Sparse&);
  static  void solver(Sparse& , Sparse&);
};

std::ostream& operator<< (std::ostream &, Sparse const& );


#endif
