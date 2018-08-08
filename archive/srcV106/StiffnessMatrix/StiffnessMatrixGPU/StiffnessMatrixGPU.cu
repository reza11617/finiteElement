#include "StiffnessMatrixGPU.h"

StiffnessMatrixGPU::StiffnessMatrixGPU(Material& mat, Geometry &geo, unsigned int n)
  : StiffnessMatrixFirstOrder(mat, geo, n)
{
#if __CUDA_ARCH__
  printf("Stiffness Matrix GPU Created by GPU");
#elif !defined(__CUDA_ARCH__)
  int device = -1;
  cudaGetDevice(&device);
  // define array in GPU
  cudaMallocManaged(&integrationNode_d, numberOfIntegrationPoint*sizeof(float));
  cudaMallocManaged(&integrationPos_d, numberOfIntegrationPoint*dimention*numberOfIntegrationPoint*sizeof(unsigned int));
  cudaMallocManaged(&integrationWeight_d, numberOfIntegrationPoint*sizeof(float));
  cudaMemcpy(integrationNode_d, integrationNode, numberOfIntegrationPoint*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(integrationPos_d, integrationPos, numberOfIntegrationPoint*dimention*numberOfIntegrationPoint*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy(integrationWeight_d,integrationWeight, numberOfIntegrationPoint*sizeof(float), cudaMemcpyHostToDevice);
  // copying arrays in geometry class to the gpu
  cudaMallocManaged(&x_d, geometry->numberOfNodes*sizeof(float));
  cudaMallocManaged(&y_d, geometry->numberOfNodes*sizeof(float));
  cudaMallocManaged(&mesh_d, geometry->numberOfElementsG*4*sizeof(unsigned int));
  cudaMemcpy(x_d, geometry->x, geometry->numberOfNodes*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(y_d, geometry->y, geometry->numberOfNodes*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(mesh_d, geometry->mesh, geometry->numberOfElementsG*4*sizeof(unsigned int), cudaMemcpyHostToDevice);
  // constant stiffMat and DOFs
  cudaMallocManaged(&c_d, numberOfElements*6*sizeof(float));
  stiffMat_d = new Sparse(stiffMatSize,geometry->numberOfNodes*2,1);
  // copy from the material matarix
  cudaMallocManaged(&D_d, 6*sizeof(float));
  cudaMemcpy(D_d, material->materialMatrix, 6*sizeof(float), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  Log::Logger().Info("StiffnessMatrixGPU created by CPU");
#endif
};

StiffnessMatrixGPU::~StiffnessMatrixGPU()
{
#if __CUDA_ARCH__
  printf("Stiffness Matrix GPU deleted by GPU");
#elif !defined(__CUDA_ARCH__)
  Log::Logger().Info("StiffnessMatrixGPU deleted by CPU");
  cudaFree(integrationPos_d);
  cudaFree(integrationNode_d);
  cudaFree(integrationWeight_d);
  cudaFree(c_d);
  cudaFree(x_d);
  cudaFree(y_d);
  cudaFree(mesh_d);
  cudaFree(D_d);
#endif
}

__global__ void constantCreatorKernel(int n, float* c, float* x, float* y, unsigned int* mesh, StiffnessMatrixGPU *s)
{
  //printf("in the function\n blockDim.x = %d, gridDim.x = %d, blockIdx.x = %d\n", blockDim.x,gridDim.x, blockIdx.x);
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    {
      //printf("i is %d stride is %d threadID = %d\n",i,stride,threadIdx.x);
      s->constantCreator(i, c, x, y, mesh);
    }
};


__global__ void StiffnessMatrixKernel(unsigned int n, unsigned int nip, float* in, unsigned int* ip, float* iw, float* c, float* D, float* k, StiffnessMatrixGPU *obj)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    {
      //printf("i is: %d \n", i);
      unsigned int noIP = i%nip;
      unsigned int numberElement = i/nip;
      obj->stiffnessMatrixCalculation(numberElement, noIP, nip, in, ip, iw, c, D, k);
    }
}

Sparse& StiffnessMatrixGPU::GetStiffnessMatrix()
{
  blockSize = 1024;
  //numberOfElements=33;
  int numBlocks = (numberOfElements + blockSize-1)/blockSize;
  constantCreatorKernel<<<numBlocks, blockSize>>>(numberOfElements, c_d, x_d, y_d, mesh_d, this);
  cudaDeviceSynchronize();
  numBlocks = (simulationSize + blockSize-1)/blockSize;
  Timer timer("Time spend in GPU: ");
  StiffnessMatrixKernel<<<numBlocks, blockSize>>>(simulationSize, nipSquared, integrationNode_d, integrationPos_d, integrationWeight_d, c_d, D_d, stiffMat_d->value, this);
  cudaDeviceSynchronize();
  return *stiffMat_d;
}
