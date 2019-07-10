#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "../Log/Log.h"

class Material {
private:
  double ElasticModulus;
  double poissonRatio;
public:
  double *materialMatrix;
  Material(double, double);
  ~Material();
  void elasticMaterial();
};

#endif
