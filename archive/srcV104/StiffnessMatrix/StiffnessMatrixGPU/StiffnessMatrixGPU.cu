#include "StiffnessMatrixGPU.h"

__device__ void stiffnessMatrixCalculation2(unsigned int numberElement, unsigned int noIP, unsigned int nip ,float* in, unsigned int* ip, float* iw, float* c, float* D, float* k)
// numberElement -> the element number needed to be calculated
// noIP -> the integration point number needed to be calculated
// nip is the number of integration point squared.
// in is the integrationNode
// ip -> integrationPos
// iw -> integrationWeight
// c -> constants
// D -> material matrix
// k -> stiffness matrix
{
  double XI = in[ip[2*noIP]]; double YI = in[ip[2*noIP+1]];
  // Jacobian
  double J11 = c[numberElement*6+0]*YI-c[numberElement*6+1]; double J12 = c[numberElement*6+3]*YI-c[numberElement*6+4];
  double J21 = c[numberElement*6+0]*XI-c[numberElement*6+2]; double J22 = c[numberElement*6+3]*XI-c[numberElement*6+5];
  double detJ = J11*J22-J12*J21;
  double WeightPerDetJ = (iw[ip[2*noIP]]*iw[ip[2*noIP+1]])/detJ;
  // derveativs of the shape function N1x N2x ... N1y N2y ...
  double Ni[8] = {J22*( YI-1)/4 -  J12*( XI-1)/4, J22*(-YI+1)/4 -  J12*(-XI-1)/4, \
		  J22*( YI+1)/4 -  J12*( XI+1)/4, J22*(-YI-1)/4 -  J12*(-XI+1)/4, \
		  J11*( XI-1)/4 -  J21*( YI-1)/4, J11*(-XI-1)/4 -  J21*(-YI+1)/4, \
		  J11*( XI+1)/4 -  J21*( YI+1)/4, J11*(-XI+1)/4 -  J21*(-YI-1)/4};
  // multiplication of shape functions N1x^2 N1x*N2x ....
  double N[36];
  unsigned int counterN = 0;
  for (unsigned int i = 0; i < 8; i++)
    {
      for (unsigned int j = i; j < 8 ; j++)
	N[counterN++] = Ni[i]*Ni[j];
    };
  // find the position to start filling the stiffness matrix
  unsigned int counter = 36*(numberElement*nip+noIP);
  // writes all 36 components of the 8 by 8 stiffness Matrix considering symmetry
  k[counter+0]  = WeightPerDetJ*(D[0]*N[0] + 2*D[4]*N[4] + D[2]*N[26]);
  k[counter+1]  = WeightPerDetJ*(D[4]*N[0] + D[5]*N[26] + D[3]*N[4] + D[2]*N[4]);
  k[counter+2]  = WeightPerDetJ*(D[2]*N[0] + 2*D[5]*N[4] + D[1]*N[26]);
  k[counter+3]  = WeightPerDetJ*(D[0]*N[1] + D[4]*N[5] + D[4]*N[11] + D[2]*N[27]);
  k[counter+4]  = WeightPerDetJ*(D[4]*N[1] + D[3]*N[11] + D[2]*N[5] + D[5]*N[27]);
  k[counter+5]  = WeightPerDetJ*(D[0]*N[8] + 2*D[4]*N[12] + D[2]*N[30]);
  k[counter+6]  = WeightPerDetJ*(D[4]*N[1] + D[3]*N[5] + D[2]*N[11] + D[5]*N[27]);
  k[counter+7]  = WeightPerDetJ*(D[2]*N[1] + D[5]*N[5] + D[5]*N[11] + D[1]*N[27]);
  k[counter+8]  = WeightPerDetJ*(D[4]*N[8] + D[5]*N[30] + D[3]*N[12] + D[2]*N[12]);
  k[counter+9]  = WeightPerDetJ*(D[2]*N[8] + 2*D[5]*N[12] + D[1]*N[30]);
  k[counter+10] = WeightPerDetJ*(D[0]*N[2] + D[4]*N[6] + D[4]*N[17] + D[2]*N[28]);
  k[counter+11] = WeightPerDetJ*(D[4]*N[2] + D[3]*N[17] + D[2]*N[6] + D[5]*N[28]);
  k[counter+12] = WeightPerDetJ*(D[0]*N[9] + D[4]*N[13] + D[4]*N[18] + D[2]*N[31]);
  k[counter+13] = WeightPerDetJ*(D[4]*N[9] + D[3]*N[18] + D[2]*N[13] + D[5]*N[31]);
  k[counter+14] = WeightPerDetJ*(D[0]*N[15] + 2*D[4]*N[19] + D[2]*N[33]);
  k[counter+15] = WeightPerDetJ*(D[4]*N[2] + D[3]*N[6] + D[2]*N[17] + D[5]*N[28]);
  k[counter+16] = WeightPerDetJ*(D[2]*N[2] + D[5]*N[6] + D[5]*N[17] + D[1]*N[28]);
  k[counter+17] = WeightPerDetJ*(D[4]*N[9] + D[3]*N[13] + D[2]*N[18] + D[5]*N[31]);
  k[counter+18] = WeightPerDetJ*(D[2]*N[9] + D[5]*N[13] + D[5]*N[18] + D[1]*N[31]);
  k[counter+19] = WeightPerDetJ*(D[4]*N[15] + D[5]*N[33] + D[3]*N[19] + D[2]*N[19]);
  k[counter+20] = WeightPerDetJ*(D[2]*N[15] + 2*D[5]*N[19] + D[1]*N[33]);
  k[counter+21] = WeightPerDetJ*(D[0]*N[3] + D[4]*N[7] + D[4]*N[22] + D[2]*N[29]);
  k[counter+22] = WeightPerDetJ*(D[4]*N[3] + D[3]*N[22] + D[2]*N[7] + D[5]*N[29]);
  k[counter+23] = WeightPerDetJ*(D[0]*N[10] + D[4]*N[14] + D[4]*N[23] + D[2]*N[32]);
  k[counter+24] = WeightPerDetJ*(D[4]*N[10] + D[3]*N[23] + D[2]*N[14] + D[5]*N[32]);
  k[counter+25] = WeightPerDetJ*(D[0]*N[16] + D[4]*N[20] + D[4]*N[24] + D[2]*N[34]);
  k[counter+26] = WeightPerDetJ*(D[4]*N[16] + D[3]*N[24] + D[2]*N[20] + D[5]*N[34]);
  k[counter+27] = WeightPerDetJ*(D[0]*N[21] + 2*D[4]*N[25] + D[2]*N[35]);
  k[counter+28] = WeightPerDetJ*(D[4]*N[3] + D[3]*N[7] + D[2]*N[22] + D[5]*N[29]);
  k[counter+29] = WeightPerDetJ*(D[2]*N[3] + D[5]*N[7] + D[5]*N[22] + D[1]*N[29]);
  k[counter+30] = WeightPerDetJ*(D[4]*N[10] + D[3]*N[14] + D[2]*N[23] + D[5]*N[32]);
  k[counter+31] = WeightPerDetJ*(D[2]*N[10] + D[5]*N[14] + D[5]*N[23] + D[1]*N[32]);
  k[counter+32] = WeightPerDetJ*(D[4]*N[16] + D[3]*N[20] + D[2]*N[24] + D[5]*N[34]);
  k[counter+33] = WeightPerDetJ*(D[2]*N[16] + D[5]*N[20] + D[5]*N[24] + D[1]*N[34]);
  k[counter+34] = WeightPerDetJ*(D[4]*N[21] + D[5]*N[35] + D[3]*N[25] + D[2]*N[25]);
  k[counter+35] = WeightPerDetJ*(D[2]*N[21] + 2*D[5]*N[25] + D[1]*N[35]);
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
  cudaMallocManaged(&stiffMat_d, stiffMatSize*sizeof(float));
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
  cudaFree(stiffMat_d);
  cudaFree(x_d);
  cudaFree(y_d);
  cudaFree(mesh_d);
  cudaFree(D_d);
#endif
}
__global__ void StiffnessMatrixKernel(unsigned int n, unsigned int nip, float* in, unsigned int* ip, float* iw, float* c, float* D, float* k, StiffnessMatrixGPU *obj)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    {
      unsigned int noIP = i%nip;
      unsigned int numberElement = i/nip;
      obj->stiffnessMatrixCalculation(numberElement, noIP, nip, in, ip, iw, c, D, k);
      //stiffnessMatrixCalculation2(numberElement, noIP, nip, in, ip, iw, c, D, k);
    }
}

float* StiffnessMatrixGPU::GetStiffnessMatrix()
{
  blockSize = 256;
  //numberOfElements=33;
  int numBlocks = (numberOfElements + blockSize-1)/blockSize;
  constantCreatorKernel<<<numBlocks, blockSize>>>(numberOfElements, c_d, x_d, y_d, mesh_d, this);
  cudaDeviceSynchronize();
  Timer timer("Time spend in GPU: ");
  StiffnessMatrixKernel<<<numBlocks, blockSize>>>(simulationSize, nipSquared, integrationNode_d, integrationPos_d, integrationWeight_d, c_d, D_d, stiffMat_d, this);
  cudaDeviceSynchronize();
  return stiffMat_d;
}
