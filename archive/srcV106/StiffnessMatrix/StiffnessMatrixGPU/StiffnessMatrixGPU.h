#ifndef STIFFNESSMATRIXGPU_H
#define STIFFNESSMATRIXGPU_H

#include <iostream>
#include "../StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../../Log/Log.h"
#include "../../Timer/Timer.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixGPU : public StiffnessMatrixFirstOrder
{
  // variables
public:
  float* x_d; float* y_d; unsigned int* mesh_d; float* D_d; // required arrays to copy to GPU
  float* c_d;
  Sparse* stiffMat_d;
  float* integrationNode_d; unsigned int* integrationPos_d; float* integrationWeight_d;  
  //unsigned int simulationSize;
private:
  unsigned int blockSize;
  // methods
public:
  CUDA_HOSTDEV StiffnessMatrixGPU(Material&, Geometry&, unsigned int);
  CUDA_HOSTDEV ~StiffnessMatrixGPU();
  Sparse& GetStiffnessMatrix() override;
private:

};
#endif
