#ifndef STIFFNESSMATRIXPARALLELCPU
#define STIFFNESSMATRIXPARALLELCPU

#include <iostream>
#include <thread>
#include "/StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"

class StiffnessMatrixParallelCPU : public StiffnessMatrixFirstOrder
{
  //variables
private:
  static unsigned concurentThreadsSupported;
  unsigned int* simulationPerThread;
  unsigned int threadNumber;
  //methods
public:
  StiffnessMatrixFirstOrder(Material&, Geometry&, unsigned int);
  ~StiffnessMatrixFirstOrder();
  float* GetStiffnessMatrix();
private:  
  void stiffnessMatrixFirstOrderForEachThread(unsigned int);
  void threadInitializer();
};
#endif
