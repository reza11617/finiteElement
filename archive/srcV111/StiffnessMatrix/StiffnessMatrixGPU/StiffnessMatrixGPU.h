#ifndef STIFFNESSMATRIXGPU_H
#define STIFFNESSMATRIXGPU_H

#include <iostream>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "../StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../../Log/Log.h"
#include "../../Timer/Timer.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixGPU : public StiffnessMatrixFirstOrder
{
  // variables
public:
  float* D_d; // required arrays to copy to GPU
private:
  unsigned int blockSize;
  // methods
public:
  StiffnessMatrixGPU(Material&, Geometry&, unsigned int);
  ~StiffnessMatrixGPU();
  Sparse& GetStiffnessMatrix() override;
private:

};
#endif
