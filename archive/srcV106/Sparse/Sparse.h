#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include "../Log/Log.h"

class Sparse
{
  //variables
public:
  float* value;
  unsigned int* i;
  unsigned int* j;
private:
  unsigned int numberOfRows; //number of rows
  unsigned int numberOfColumns; //number of columns
  unsigned int valueSize; //how many nonzeros
  int device;
  bool symmetry = false; //by default matrix is not symmetry 
  // Methods
public:
  Sparse(unsigned int*, unsigned int*, float*, unsigned int, unsigned int, unsigned int);
  Sparse(unsigned int, unsigned int, int);
  Sparse();
  ~Sparse();
  void iSet(unsigned int*);
  void jSet(unsigned int*);
  void xSet(float*);
  void sizeSet(unsigned int, unsigned int, unsigned int);
  void Print();

  // staic methods
public:
  static unsigned int Printer(Sparse & s);
  static float Compare(Sparse& , Sparse&);
};
#endif
