#include "Material.h"

Material::Material(float E, float nu)
  : ElasticModulus(E), poissonRatio(nu)
{
  elasticMaterial();
};
Material::~Material()
{
  Log::Logger().Info("Material deleted");
  delete[] materialMatrix;
};
void Material::elasticMaterial()
{
  // this function will create a vector for material matrix
  // D11 D22 D33 D12 D13 D23
  materialMatrix = new float[6];
  float E1 = ElasticModulus*(1-poissonRatio)/((1-poissonRatio*2)*(1+poissonRatio));
  float E2 = poissonRatio*E1/(1-poissonRatio);
  float G  = ElasticModulus/(2*(1+poissonRatio));
  materialMatrix[0] = E1;
  materialMatrix[1] = E1;
  materialMatrix[2] =  G;
  materialMatrix[3] = E2;
  materialMatrix[4] = 0;
  materialMatrix[5] = 0;
}
