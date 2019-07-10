#include "StiffnessMatrixParallelCPU.h"

StiffnessMatrixParallelCPU::StiffnessMatrixParallelCPU(Material& mat, Geometry &geo, unsigned int n)
  : StiffnessMatrixParallelCPU(mat, geo, n,(std::thread::hardware_concurrency()-1))
{ 
};

StiffnessMatrixParallelCPU::StiffnessMatrixParallelCPU(Material& mat, Geometry &geo, unsigned int n, unsigned int numberOfCores)
  : StiffnessMatrixFirstOrder(mat, geo, n), concurentThreadsSupported(numberOfCores)
{
  Log::Logger().Info("StiffnessMatrixParallelCPU created with " + std::to_string(concurentThreadsSupported) + " threads");
  simulationPerThread = new unsigned int[concurentThreadsSupported];
  for (unsigned int i = 0; i < concurentThreadsSupported; i++)
    {
      if (i == concurentThreadsSupported-1)
        simulationPerThread[i] = simulationSize/concurentThreadsSupported + simulationSize%concurentThreadsSupported;
      else
        simulationPerThread[i] = simulationSize/concurentThreadsSupported;
    };
  
}

StiffnessMatrixParallelCPU::~StiffnessMatrixParallelCPU()
{
  Log::Logger().Info("StiffnessMatrixParallelCPU Deleted");
  delete[] simulationPerThread;
}

float* StiffnessMatrixParallelCPU::GetStiffnessMatrix()
{
  for (unsigned int i = 0; i<numberOfElements; i++)
    constantCreator(i, c, geometry->x, geometry->y, geometry->mesh);
  std::cout << "Time spend in MultiCoreSimulation function running on " << concurentThreadsSupported << " cores is: " << std::endl;
  Timer timer;
  std::thread t[concurentThreadsSupported]; // number of threads being used in this program
  for (unsigned int i = 0; i<concurentThreadsSupported; i++)
    {
      t[i] = std::thread(&StiffnessMatrixParallelCPU::GetStiffnessMatrixForEachThread,this,i);
    }
  for (unsigned int i = 0; i<concurentThreadsSupported; i++)
    {
      t[i].join();
    }
  return stiffMat;
}

void StiffnessMatrixParallelCPU::GetStiffnessMatrixForEachThread(unsigned int threadId)
{
  unsigned int counter = threadId*(simulationSize/concurentThreadsSupported);
  for (unsigned int i = 0; i<simulationPerThread[threadId]; i++)
    {
      int numberElement = (counter+i)/nipSquared;
      int noIP = (counter+i)%nipSquared;
      DOFCreator(numberElement,noIP);
      stiffnessMatrixCalculation(numberElement, noIP,  nipSquared,integrationNode, integrationPos, integrationWeight, c, material->materialMatrix, stiffMat);

    }

}
