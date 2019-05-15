#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H

#include <iostream>
#include <string>

#include "../Material/Material.h"
#include "../Geometry/Geometry.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"

class StiffnessMatrix
{
  //VARIABLES
public:
  unsigned int* DOF_i; // DOF_i a vector rows: 36*NOE  
  unsigned int* DOF_j; // (DOF_i, DOF_j) this shows the dofs related to stiffness matrix for the project
  float* stiffMat; // stiffness matrix a 2D vector rows: 36*NOE and columns: NIP^2  
protected:
  static const unsigned int dimention = 2;
  unsigned int numberOfIntegrationPoint;
  unsigned int numberOfElements; //number of elements
  unsigned int nipSquared; // number of interation point squared
  float* integrationNode; unsigned int* integrationPos; float* integrationWeight;
  Material* material;
  Geometry* geometry;
  float* c;
  unsigned int stiffMatSize;
  unsigned int sizeStiffMatPerEle;
  unsigned int simulationSize;
  //std::stirng hardwareType;
private:
  //METHODS
public:
  StiffnessMatrix();
  StiffnessMatrix(Material&, Geometry&, unsigned int);
  virtual ~StiffnessMatrix();
  //void SetHardwareType(std::string&);
  virtual float* GetStiffnessMatrix() = 0;
  //virtual float* GetStiffnessMatrixCPU() = 0;
  //virtual float* GetStiffnessMatrixParallelCPU() = 0;
  //virtual float* GetStiffnessMatrixGPU() = 0;
protected:
  
private:
  void integrationPoint();
};
#endif
