#ifndef STIFFNESSMATRIXFIRSTORDER_H
#define STIFFNESSMATRIXFIRSTORDER_H

#include <iostream>
#include <cuda.h>

#include "../StiffnessMatrix.h"
#include "../../Log/Log.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixFirstOrder : public StiffnessMatrix
{
  //VARIABLES
public:

protected:
  
private:

  //METHODS
public:
  StiffnessMatrixFirstOrder();
  StiffnessMatrixFirstOrder(Material&, Geometry&, unsigned int);
  virtual ~StiffnessMatrixFirstOrder();
  virtual Sparse& GetStiffnessMatrix() = 0;
  int GetStiffnessMatrixSize() override;
  CUDA_HOSTDEV void stiffnessMatrixCalculation(unsigned int, unsigned int, double*, unsigned int*, double*, double*, double*, unsigned int* ,double*, unsigned int*, unsigned int*, unsigned int*, const double*);
  CUDA_HOSTDEV void constantCreator(unsigned int numberElement, double* c, double* x, double* y, unsigned int* mesh);
protected:
private:


};
#endif


