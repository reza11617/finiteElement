#include "ElementParallel.h"
#define D material->materialMatrix

unsigned int ElementParallel::concurentThreadsSupported = 5; //(std::thread::hardware_concurrency()-1);

ElementParallel::ElementParallel(Material& mat, Geometry &geo, unsigned int n)
  : Element(mat, geo, n)
{
  // Initialize();
  // Multi-core processing
  
};

ElementParallel::~ElementParallel()
{
  Log::Logger().Info("ElementParallel Deleted");
  delete[] simulationPerThread;
}

void ElementParallel::MultiCoreSimulation()
{
  Initialize();
  threadInitializer();
  std::cout << "Time spend in MultiCoreSimulation function running on " << concurentThreadsSupported << " cores is: " << std::endl;
  Timer timer;
  std::thread t[concurentThreadsSupported]; // number of threads being used in this program
  for (unsigned int i = 0; i<concurentThreadsSupported; i++)
    {
      t[i] = std::thread(&ElementParallel::stiffnessMatrixFirstOrderForEachThread,this,i);
    }
  for (unsigned int i = 0; i<concurentThreadsSupported; i++)
    {
      t[i].join();
    }
}

void ElementParallel::threadInitializer()
{
  simulationSize = numberOfElements*nipSquared;
  simulationPerThread = new unsigned int[concurentThreadsSupported];
  for (unsigned int i = 0; i < concurentThreadsSupported; i++)
    {
      if (i == concurentThreadsSupported-1)
        simulationPerThread[i] = simulationSize/concurentThreadsSupported + simulationSize%concurentThreadsSupported;
      else
        simulationPerThread[i] = simulationSize/concurentThreadsSupported;
    }
  //Log::Logger().Info(simulationPerThread,concurentThreadsSupported);
}

void ElementParallel::stiffnessMatrixFirstOrderForEachThread(unsigned int threadId)
{
  unsigned int counter = threadId*(simulationSize/concurentThreadsSupported);
  for (unsigned int i = 0; i<simulationPerThread[threadId]; i++)
    {
      //std::cout<<(counter+i)/nipSquared <<" , " << (counter+i)%nipSquared << std::endl;
      stiffnessMatrixFirstOrder((counter+i)/nipSquared, (counter+i)%nipSquared);
    }
}

void ElementParallel::GpuSimulations() {};
