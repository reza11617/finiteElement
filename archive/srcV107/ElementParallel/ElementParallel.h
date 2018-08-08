#ifndef ELEMENTPARALLEL_H
#define ELEMENTPARALLEL_H

#include <iostream>
#include <thread>
#include <string>
#include "../Material/Material.h"
#include "../Geometry/Geometry.h"
#include "../Element/Element.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"

class ElementParallel : public Element
{
  // variables
private:
  static unsigned concurentThreadsSupported;
  unsigned int* simulationPerThread;
  unsigned int threadNumber;
public:
  ElementParallel(Material&, Geometry&, unsigned int);
  ~ElementParallel();
  void MultiCoreSimulation () override;
  virtual void GpuSimulations();
protected:  
  void stiffnessMatrixFirstOrderForEachThread(unsigned int);
  void threadInitializer();
};
#endif
