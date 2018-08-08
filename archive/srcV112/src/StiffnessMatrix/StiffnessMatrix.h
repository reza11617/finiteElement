#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H

//macro
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#include <iostream>
#include <string>
#include <cuda_runtime_api.h>
#include <cuda.h>

#include "../Material/Material.h"
#include "../Geometry/Geometry.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"
#include "../Sparse/Sparse.h"

class StiffnessMatrix
{
  //VARIABLES
public:
  Sparse* stiffMat;
  //unsigned int* DOF_i; // DOF_i a vector rows: 36*NOE  
  //unsigned int* DOF_j; // (DOF_i, DOF_j) this shows the dofs related to stiffness matrix for the project
  //float* stiffMat; // stiffness matrix a 2D vector rows: 36*NOE and columns: NIP^2  
protected:
  static const unsigned int dimention = 2;
  unsigned int numberOfIntegrationPoint;
  unsigned int numberOfElements; //number of elements
  unsigned int nipSquared; // number of interation point squared
  float* integrationNode; unsigned int* integrationPos; float* integrationWeight;
  Material* material;
  Geometry* geometry;
  float* c;
  unsigned int sizeStiffMatPerEle;
  unsigned int simulationSize;
  unsigned int stiffMatSize;
  unsigned int globalStiffMatSize;
  //std::stirng hardwareType;
private:
  //METHODS
public:
  StiffnessMatrix(Material&, Geometry&, unsigned int);
  virtual ~StiffnessMatrix();
  virtual int GetStiffnessMatrixSize () = 0;
  virtual Sparse& GetStiffnessMatrix() = 0;
protected:
  unsigned int globalStiffMatSizeCalculator(unsigned int*, unsigned int, unsigned int);
private:
  CUDA_HOSTDEV void integrationPoint();
};
#endif
