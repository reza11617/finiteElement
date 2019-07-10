#ifndef STIFFNESSMATRIXSINGLECPU
#define STIFFNESSMATRIXSINGLECPU

#include <iostream>
#include <thread>
#include <string>
#include "../StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../StiffnessMatrix.h"
#include "../../Log/Log.h"
#include "../../Timer/Timer.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixSingleCPU : public StiffnessMatrixFirstOrder
{
  //variables
private:
  unsigned concurentThreadsSupported;
  unsigned int* simulationPerThread;
  unsigned int threadNumber;
  //methods
public:
  StiffnessMatrixSingleCPU(Material&, Geometry&, unsigned int);
  ~StiffnessMatrixSingleCPU();
  Sparse* GetStiffnessMatrix() override;
private:  
};
#endif
