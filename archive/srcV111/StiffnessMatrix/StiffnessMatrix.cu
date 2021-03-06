#include "StiffnessMatrix.h"

StiffnessMatrix::StiffnessMatrix(Material& mat, Geometry& geo, unsigned int n)
  : material(&mat), geometry(&geo), numberOfIntegrationPoint(n)
{
  numberOfElements = geometry->get_numberOfElementsG();
  nipSquared = numberOfIntegrationPoint*numberOfIntegrationPoint;
  simulationSize = numberOfElements*nipSquared;
  // integration points
  cudaMallocManaged(&integrationNode, numberOfIntegrationPoint*sizeof(float));
  cudaMallocManaged(&integrationPos, numberOfIntegrationPoint*dimention*numberOfIntegrationPoint*sizeof(unsigned int));
  cudaMallocManaged(&integrationWeight, numberOfIntegrationPoint*sizeof(float));
  integrationPoint();
  Log::Logger().Info("StiffnessMatrix Created by CPU");
};

StiffnessMatrix::~StiffnessMatrix()
{
  Log::Logger().Info("StiffnessMatrix Deleted by CPU");
  cudaFree(integrationNode);
  cudaFree(integrationPos);
  cudaFree(integrationWeight);
};

void StiffnessMatrix::integrationPoint()
// Creats the integration points
// XI = integrationNode[integrationPos[i]] YI = integrationNode[integrationPos[i+1]] 
{

  unsigned int counter = 0;
  for (unsigned int i = 0; i < numberOfIntegrationPoint; i++)
    for (unsigned int j = 0; j < numberOfIntegrationPoint; j++)
      {
	integrationPos[counter++] = i;
	integrationPos[counter++] = j;
      };
  
  if (numberOfIntegrationPoint == 1) {
    integrationNode[0] = 0; integrationWeight[0] = 4;
  } else if (numberOfIntegrationPoint == 2) {
    integrationNode[0] = -0.57735; integrationWeight[0] = 1.0;
    integrationNode[1] =  0.57735; integrationWeight[1] = 1.0;
  } else if (numberOfIntegrationPoint == 3) {
    integrationNode[0] = -0.774596; integrationWeight[0] = 0.555556;
    integrationNode[1] =  0.0     ; integrationWeight[1] = 0.888889;
    integrationNode[2] =  0.774596; integrationWeight[2] = 0.555556;
  } else if (numberOfIntegrationPoint == 4) {
    integrationNode[0] = -0.861136; integrationWeight[0] = 0.347855;
    integrationNode[1] = -0.339981; integrationWeight[1] = 0.652145;
    integrationNode[2] =  0.339981; integrationWeight[2] = 0.652145;
    integrationNode[3] =  0.861136; integrationWeight[3] = 0.347855;
  } else if (numberOfIntegrationPoint == 5) {
    integrationNode[0] = -0.90618;  integrationWeight[0] = 0.236927;
    integrationNode[1] = -0.538469; integrationWeight[1] = 0.478629;
    integrationNode[2] =  0.0     ; integrationWeight[2] = 0.568889;
    integrationNode[3] =  0.538469; integrationWeight[3] = 0.478629;
    integrationNode[4] =  0.90618;  integrationWeight[4] = 0.236927;
  } else {
    printf("Integration points more than five is under construction");
  }
};
