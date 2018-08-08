#include "StiffnessMatrixSingleCPU.h"

StiffnessMatrixSingleCPU::StiffnessMatrixSingleCPU(Material& mat, Geometry &geo, unsigned int n)
  : StiffnessMatrixFirstOrder(mat, geo, n)
{ 
  Log::Logger().Info("StiffnessMatrixSingleCPU created");
};

StiffnessMatrixSingleCPU::~StiffnessMatrixSingleCPU()
{
  Log::Logger().Info("StiffnessMatrixSingleCPU Deleted");
}

Sparse& StiffnessMatrixSingleCPU::GetStiffnessMatrix()
{
  for (unsigned int i = 0; i<numberOfElements; i++)
    constantCreator(i, c, geometry->get_x(), geometry->get_y(), geometry->get_mesh());
  Timer timer("Time spend in CPU using signle core: ");
  for (unsigned int i = 0; i<numberOfElements; i++)
    stiffnessMatrixCalculation(i, nipSquared,integrationNode, integrationPos, integrationWeight, c, material->materialMatrix, geometry->get_mesh(), stiffMat->value, stiffMat->i, stiffMat->j, stiffMat->nnz_inRow, geometry->get_Dof().get_free());
  return *stiffMat;
}
