#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include "../Log/Log.h"

class Material {
private:
  float ElasticModulus;
  float poissonRatio;
public:
  float *materialMatrix;
  Material(float, float);
  ~Material();
  void elasticMaterial();
};

#endif
