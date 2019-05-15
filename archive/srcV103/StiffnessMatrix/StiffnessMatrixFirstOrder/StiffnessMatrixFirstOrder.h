#ifndef STIFFNESSMATRIXFIRSTORDER_H
#define STIFFNESSMATRIXFIRSTORDER_H

#include <iostream>

#include "../StiffnessMatrix.h"
#include "../../Log/Log.h"

class StiffnessMatrixFirstOrder : public StiffnessMatrix
{
  //VARIABLES
public:

protected:
  
private:

  //METHODS
public:
  StiffnessMatrixFirstOrder();
  StiffnessMatrixFirstOrder(Material&, Geometry&, unsigned int);
  virtual ~StiffnessMatrixFirstOrder();
  virtual float* GetStiffnessMatrix();
  void stiffnessMatrixCalculation(unsigned int, unsigned int, unsigned int, float*, unsigned int*, float*, float*, float*, float*);
  void constantCreator(unsigned int, float*, float*, float*, unsigned int*);
protected:
  void DOFCreator(unsigned int,unsigned int);
private:


};
#endif


