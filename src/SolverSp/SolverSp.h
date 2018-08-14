#ifndef SOLVERSP_H
#define SOLVERSP_H

#include <iostream>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cusolverSp.h>

#include "../Timer/Timer.h"

// This class solves the symmetric sparse systems A(N*N) * x(N) = F(N)
class SolverSp {
private:
  int* rowPtr;
  int* colIndices;
  double* sparseMatrix;  // matrix has to be CSR
  double* rightHandSideVector;
  double* leftHandSideVector;
  int nnz; // number of non-zeros
  int N; // size of the matrix and vectors
  double tolerance; // tolerance to decide singularity
  int reorder;
public:
  SolverSp(double*, unsigned int*, unsigned int*, unsigned int, unsigned int, double*, double*);
  //SolverSp(double*, unsigned int*, unsigned int*, unsigned int, unsigned int,double*, double*, double);
  ~SolverSp();
private:
  int SolverSpChol(double*, int*, int*);
  void SolverSpQR(double*, int*, int*);
  void fixMatrices();
};
#endif
