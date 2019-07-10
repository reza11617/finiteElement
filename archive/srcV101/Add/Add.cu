#include <iostream>
#include <math.h>
#include "Add.h"

__global__ void init(int n, float *x, float *y) {
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }
}

// Kernel function to add the elements of two arrays
__global__ void add(int n, float *x, float *y)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    y[i] = x[i] + y[i];
  
  printf("number is %d\n",index);
  //printf("index = %d   blockIdx.x = %d    blockDim.x = %d   threadIdx.x = %d\n",index,blockIdx.x,blockDim.x,threadIdx.x);
}

int runCuda()
{
  std::cout << "here "<< std::endl;
  int N = 100;
  float *x, *y;

  // set your gpu
  cudaSetDevice(0);

  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged(&x, N*sizeof(float));
  cudaMallocManaged(&y, N*sizeof(float));

  /*
  // initialize x and y arrays on the host
  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }
  */
  int blockSize = 128;
  int numBlocks = (N + blockSize-1)/blockSize;


  //initialize x and y on the device
  init<<<numBlocks, blockSize>>>(N, x, y);

  // Run kernel on 1M elements on the GPU
  add<<<numBlocks, blockSize>>>(N, x, y);

  // Wait for GPU to finish before accessing on host
  cudaDeviceSynchronize();

  // Check for errors (all values should be 3.0f)
  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = fmax(maxError, fabs(y[i]-3.0f));
  std::cout << "Max error: " << maxError << std::endl;

  // Free memory
  cudaFree(x);
  cudaFree(y);
  return 0;
}
