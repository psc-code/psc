
#include "psc_cuda.h"
#include "math.h"

__constant__ static float _dt;

static void set_constants()
{
  _dt = psc.dt;
}

__device__ static void
push_part_yz_a_one(int n, struct d_part d_part)
{
  float4 xi4  = d_part.xi4[n];
  float4 pxi4 = d_part.pxi4[n];
  
  float root = 1. / sqrt(1. + sqr(pxi4.x) + sqr(pxi4.y) + sqr(pxi4.z));
  float vyi = pxi4.y * root;
  float vzi = pxi4.z * root;
  
  xi4.y += vyi * .5 * _dt;
  xi4.z += vzi * .5 * _dt;

  d_part.xi4[n] = xi4;
}

__global__ static void
push_part_yz_a(int n_part, struct d_part d_part, int stride)
{
  int n = threadIdx.x + blockDim.x * blockIdx.x;

  while (n < n_part) {
    push_part_yz_a_one(n, d_part);
    n += stride;
  }
}

void
cuda_push_part_yz_a()
{
  struct psc_cuda *cuda = psc.c_ctx;

  set_constants();

  const int threadsPerBlock = 128;
  const int gridSize = 256;
  int dimGrid[2]  = { gridSize, 1 };
  int dimBlock[2] = { threadsPerBlock, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_yz_a, (psc.n_part, cuda->d_part,
			      gridSize * threadsPerBlock));

  for (int n = 0; n < psc.n_part; n++) {
  }
}

