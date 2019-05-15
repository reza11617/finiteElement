#ifndef ELEMENTCUDA_H
#define ELEMENTCUDA_H

#include <iostream>
#include <thread>
#include <string>
#include "../Material/Material.h"
#include "../Geometry/Geometry.h"
#include "../ElementParallel/ElementParallel.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"

class ElementCuda : public ElementParallel
{
  // variables
public:
  static const unsigned int blockSize = 256;
  // constants for stiffness matrix
  float* x_d; float* y_d; unsigned int* mesh_d; float* D_d; // required arrays to copy to GPU
  // required variables
  unsigned int simulationSize;
protected:

private:
  //methods
public:
  ElementCuda(Material&, Geometry&, unsigned int);
  ~ElementCuda();
  void InitializeGpu();
  void constantCreator(unsigned int) override;
  void GpuSimulations() override;
  
protected:
private:
};
#endif
