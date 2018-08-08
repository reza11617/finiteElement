#ifndef ELEMENTCUDA_H
#define ELEMENTCUDA_H

#include <iostream>
#include <thread>
#include <string>
#include "../Material/Material.h"
#include "../Geometry/Geometry.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"

class ElementCuda
{
  // variables
private:
  static const unsigned int dimention = 2;
  static const unsigned int sizeStiffMatPerEle = 36;
  static const unsigned int blockSize = 256;
  static unsigned concurentThreadsSupported;
  // integration point variables
  float* integrationNode;
  unsigned int* integrationPos;
  float* integrationWeight;
  // constants for stiffness matrix
  float* c;
  float* x_d; float* y_d; unsigned int* mesh_d; float* D_d; // required arrays to copy to GPU
  // required variables
  unsigned int simulationSize;
  unsigned int numberOfIntegrationPoint;
  unsigned int numberOfElements; //number of elements
  unsigned int stiffMatSize;
  unsigned int nipSquared; // number of interation point squared
  Material* material;
  Geometry* geometry;
public:
  float* stiffMat; // stiffness matrix a 2D vector rows: 36*NOE and columns: N$
  // Methods
private:

public:
  ElementCuda(Material&, Geometry&, unsigned int);
  ~ElementCuda();
  void integrationPoint(); // creates integration points and weight
  //void stiffnessMatrixCreator(unsigned int);
  //void DOFCreator(unsigned int,unsigned int);
  void constantCreator();
  void stiffnessMatrixFirstOrder();
  //void stiffnessMatrixFirstOrderForEachThread(unsigned int);
  //void threadInitializer();
  //void MultiCoreSimulation();
  //void SingleCoreSimulation();
  // Static Methods
private:
public:


};
#endif
