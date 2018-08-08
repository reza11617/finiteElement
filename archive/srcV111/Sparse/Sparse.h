#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include "../Log/Log.h"
#include "../Timer/Timer.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include  <cusolverDn.h>

class Sparse
{
  //variables
public:
  float* value;
  unsigned int* i;
  unsigned int* j;
  unsigned int valueSize; //how many nonzeros
private:
  unsigned int* temp_i;
  unsigned int* temp_j;
  float* temp_value;
  unsigned int* a_temp_i;
  unsigned int* a_temp_j;
  float* a_temp_value;
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
  void Assemble(unsigned int*, unsigned int);
  void set_numberOfRows(unsigned int);
  void set_numberOfColumns(unsigned int);
  void set_valueSize(unsigned int);
  void set_i(unsigned int*);
  void set_j(unsigned int*);
  void set_x(float*);
private:
  void Print();

  // staic methods
public:
  static unsigned int Printer(Sparse & s);
  static float Compare(Sparse& , Sparse&);
  static  void solver(Sparse& , Sparse&);
};
#endif
