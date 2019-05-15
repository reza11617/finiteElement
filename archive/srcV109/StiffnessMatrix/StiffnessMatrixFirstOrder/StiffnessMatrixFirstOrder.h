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
  CUDA_HOSTDEV StiffnessMatrixFirstOrder();
  CUDA_HOSTDEV StiffnessMatrixFirstOrder(Material&, Geometry&, unsigned int);
  CUDA_HOSTDEV virtual ~StiffnessMatrixFirstOrder();
  virtual Sparse& GetStiffnessMatrix() = 0;
  int GetStiffnessMatrixSize() override;
  CUDA_HOSTDEV void stiffnessMatrixCalculation(unsigned int, unsigned int, unsigned int, float*, unsigned int*, float*, float*, float*, float*);
  CUDA_HOSTDEV void constantCreator(unsigned int numberElement, float* c, float* x, float* y, unsigned int* mesh);
  CUDA_HOSTDEV void DOFCreator(unsigned int,unsigned int, unsigned int, unsigned int, unsigned int *, unsigned int *, unsigned int *);
protected:
private:


};
#endif


