#include "StiffnessMatrixSingleCPU.h"

StiffnessMatrixSingleCPU::StiffnessMatrixSingleCPU(Material& mat, Geometry &geo, unsigned int n)
  : StiffnessMatrixFirstOrder(mat, geo, n)
{ 
  Log::Logger().Info("StiffnessMatrixSingleCPU created");
  stiffMat = new Sparse(stiffMatSize,geometry->numberOfNodes*2,0);
};

StiffnessMatrixSingleCPU::~StiffnessMatrixSingleCPU()
{
  Log::Logger().Info("StiffnessMatrixSingleCPU Deleted");
  delete stiffMat;
}

Sparse* StiffnessMatrixSingleCPU::GetStiffnessMatrix()
{
  for (unsigned int i = 0; i<numberOfElements; i++)
    constantCreator(i, c, geometry->x, geometry->y, geometry->mesh);
  Timer timer("Time spend in CPU using signle core: ");
  for (unsigned int i = 0; i<numberOfElements; i++)
    {
      for (unsigned int j = 0; j<nipSquared; j++)
	{
	  DOFCreator(i,j);
	  stiffnessMatrixCalculation(i, j,  nipSquared,integrationNode, integrationPos, integrationWeight, c, material->materialMatrix, stiffMat->value);
	}
    }
  return stiffMat;
}
