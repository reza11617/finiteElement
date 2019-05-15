#ifndef STIFFNESSMATRIXPARALLELCPU
#define STIFFNESSMATRIXPARALLELCPU

#include <iostream>
#include <thread>
#include <string>
#include "../StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../StiffnessMatrix.h"
#include "../../Log/Log.h"
#include "../../Timer/Timer.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixParallelCPU : public StiffnessMatrixFirstOrder
{
  //variables
private:
  unsigned concurentThreadsSupported;
  unsigned int* simulationPerThread;
  unsigned int threadNumber;
  //methods
public:
  StiffnessMatrixParallelCPU(Material&, Geometry&, unsigned int);
  StiffnessMatrixParallelCPU(Material&, Geometry&, unsigned int, unsigned int); // takes number of the desired cores
  ~StiffnessMatrixParallelCPU();
  Sparse& GetStiffnessMatrix() override;
private:  
  void GetStiffnessMatrixForEachThread(unsigned int);
};
#endif
