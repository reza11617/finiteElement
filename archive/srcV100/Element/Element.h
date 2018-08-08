#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include <thread>
#include <string>
#include "../Material/Material.h"
#include "../Geometry/Geometry.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"


class Element
{
  // variables
private:
  static const unsigned int dimention = 2;
  static const unsigned int sizeStiffMatPerEle = 36;
  static unsigned concurentThreadsSupported;
  unsigned int simulationSize;
  unsigned int* simulationPerThread;
  unsigned int threadNumber;
  unsigned int numberOfIntegrationPoint;
  unsigned int numberOfElements; //number of elements
  unsigned int stiffMatSize;
  double* integrationNode; unsigned int* integrationPos; double* integrationWeight;
  unsigned int nipSquared; // number of interation point squared
  Material* material;
  Geometry* geometry;
  unsigned int* DOF_i; // DOF_i a vector rows: 36*NOE  
  unsigned int* DOF_j; // (DOF_i, DOF_j) this shows the dofs related to stiffness matrix for the project
  double* c; // c1x to x3y -> builds vectors of constants required in calculation of stiffness matrix each row 6 constant for each element
public:
  double* stiffMat; // stiffness matrix a 2D vector rows: 36*NOE and columns: NIP^2
  // Methods
private:

public:
  Element(Material&, Geometry&, unsigned int);
  ~Element();
  void integrationPoint(); // creates integration points and weight
  void stiffnessMatrixCreator(unsigned int);
  void DOFCreator(unsigned int,unsigned int);
  void constantCreator(unsigned int);
  void stiffnessMatrixFirstOrder(unsigned int, unsigned int);
  void stiffnessMatrixFirstOrderForEachThread(unsigned int);
  void threadInitializer();
  void MultiCoreSimulation();
  void SingleCoreSimulation();
  // Static Methods
private:
public:

  
};
#endif
