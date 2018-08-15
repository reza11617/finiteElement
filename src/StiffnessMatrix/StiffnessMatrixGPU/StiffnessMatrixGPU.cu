#include "StiffnessMatrixGPU.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
     fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

StiffnessMatrixGPU::StiffnessMatrixGPU(Material& mat, Geometry &geo, unsigned int n)
  : StiffnessMatrixFirstOrder(mat, geo, n)
{
  int device = -1;
  cudaGetDevice(&device);
  // copy from the material matarix
  cudaMallocManaged(&D_d, 6*sizeof(double));
  cudaMemcpy(D_d, material->materialMatrix, 6*sizeof(double), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  Log::Logger().Info("StiffnessMatrixGPU created by CPU");
};

StiffnessMatrixGPU::~StiffnessMatrixGPU()
{
  Log::Logger().Info("StiffnessMatrixGPU deleted by CPU");
  cudaFree(D_d);
}

__global__ void constantCreatorKernel(int n, double* c, double* x, double* y, unsigned int* mesh, StiffnessMatrixGPU *s)
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


__global__ void StiffnessMatrixKernel(unsigned int n, unsigned int nip, double* in, unsigned int* ip, double* iw, double* c, double* D, unsigned int* mesh, double* k, unsigned int* i_index, unsigned int *j_index, unsigned int* dofFree,StiffnessMatrixGPU *obj)
{
  int index  = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    {
      obj->stiffnessMatrixCalculation(i, nip, in, ip, iw, c, D, mesh, k, i_index, j_index, dofFree);
    }
}

Sparse& StiffnessMatrixGPU::GetStiffnessMatrix()
{
  blockSize = 32;
  //numberOfElements=33;
  int numBlocks = (numberOfElements + blockSize-1)/blockSize;
  constantCreatorKernel<<<numBlocks, blockSize>>>(numberOfElements, c, geometry->get_x(), geometry->get_y(), geometry->get_mesh(), this);
  cudaDeviceSynchronize();
  numBlocks = (simulationSize + blockSize-1)/blockSize;
  Timer timer("Time spend in GPU: ");
  StiffnessMatrixKernel<<<numBlocks, blockSize>>>(numberOfElements, nipSquared, integrationNode, integrationPos, integrationWeight, c, D_d, geometry->get_mesh(), stiffMat->value, stiffMat->i, stiffMat->j , geometry->get_Dof().get_free(),this);
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
  cudaDeviceSynchronize();
  return *stiffMat;
}
