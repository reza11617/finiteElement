#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include "../Log/Log.h"

class Sparse
{
private:
  unsigned int numberOfRows; //number of rows
  unsigned int numberOfColumns; //number of columns
  unsigned int valueSize; //how many nonzeros
  unsigned int* i;
  unsigned int* j;
  bool symmetry = false; //by default matrix is not symmetry 
  float* value;
  // Methods
public:
  Sparse(unsigned int*, unsigned int*, double*, unsigned int, unsigned int, unsigned int);
  Sparse();
  ~Sparse();
  void iSet(unsigned int*);
  void jSet(unsigned int*);
  void xSet(double*);
  void sizeSet(unsigned int, unsigned int, unsigned int);
  ~Sparse();
  static unsigned int Printer(Sparse & s);
  void Print();
};
#endif
