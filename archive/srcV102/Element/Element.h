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
protected:
  static const unsigned int dimention = 2;
  static const unsigned int sizeStiffMatPerEle = 36;
  unsigned int simulationSize;
  unsigned int numberOfIntegrationPoint;
  unsigned int numberOfElements; //number of elements
  unsigned int nipSquared; // number of interation point squared
  float* integrationNode; unsigned int* integrationPos; float* integrationWeight;
  float* c; // c1x to x3y -> builds vectors of constants required in calculation of stiffness matrix each row 6 constant for each element
  unsigned int stiffMatSize;
  Material* material;
  Geometry* geometry;
public:
  unsigned int* DOF_i; // DOF_i a vector rows: 36*NOE  
  unsigned int* DOF_j; // (DOF_i, DOF_j) this shows the dofs related to stiffness matrix for the project
  float* stiffMat; // stiffness matrix a 2D vector rows: 36*NOE and columns: NIP^2
  // Methods
private:

public:
  Element(Material&, Geometry&, unsigned int);
  virtual ~Element();
  //setters

  //getters
  unsigned int getStiffMatSize();
protected:
  void Initialize();
  void integrationPoint(); // creates integration points and weight
  void stiffnessMatrixCreator(unsigned int);
  void DOFCreator(unsigned int,unsigned int);
  virtual void constantCreator(unsigned int);
  void stiffnessMatrixFirstOrder(unsigned int, unsigned int);
public:
  void SingleCoreSimulation(); 
  virtual void MultiCoreSimulation() = 0;
  virtual void GpuSimulations() = 0;
};

#endif
