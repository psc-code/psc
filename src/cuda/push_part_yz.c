
#include "psc_cuda.h"
#include "math.h"
#include "profile/profile.h"

__constant__ static float _dt;

static void set_constants()
{
  float __dt = psc.dt;
  check(cudaMemcpyToSymbol(_dt, &__dt, sizeof(_dt)));
}

__device__ static void
push_part_yz_a_one(int n, struct d_part d_part)
{
  float4 xi4  = d_part.xi4[n];
  float4 pxi4 = d_part.pxi4[n];
  
  float root = 1.f / sqrt(1.f + sqr(pxi4.x) + sqr(pxi4.y) + sqr(pxi4.z));
  float vyi = pxi4.y * root;
  float vzi = pxi4.z * root;
  
  xi4.y += vyi * .5f * _dt;
  xi4.z += vzi * .5f * _dt;

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

EXTERN_C void
cuda_push_part_yz_a()
{
  static int pr;
  if (!pr) {
    pr = prof_register("cuda_part_yz_a", 1., 0, psc.n_part * 12 * sizeof(double));
  }
  prof_start(pr);

  struct psc_cuda *cuda = (struct psc_cuda *) psc.c_ctx;

  set_constants();

  const int threadsPerBlock = 128;
  const int gridSize = 256;
  int dimGrid[2]  = { gridSize, 1 };
  int dimBlock[2] = { threadsPerBlock, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_yz_a, (psc.n_part, cuda->d_part,
			      gridSize * threadsPerBlock));

  prof_stop(pr);
}

EXTERN_C void
__cuda_particles_from_fortran(struct psc_cuda *cuda)
{
  check(cudaMalloc((void **) &cuda->d_part.xi4,  psc.n_part * sizeof(float4)));
  check(cudaMalloc((void **) &cuda->d_part.pxi4, psc.n_part * sizeof(float4)));

  check(cudaMemcpy(cuda->d_part.xi4, cuda->xi4, psc.n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(cuda->d_part.pxi4, cuda->pxi4, psc.n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__cuda_particles_to_fortran(struct psc_cuda *cuda)
{
  check(cudaMemcpy(cuda->xi4, cuda->d_part.xi4, psc.n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(cuda->pxi4, cuda->d_part.pxi4, psc.n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaFree(cuda->d_part.xi4));
  check(cudaFree(cuda->d_part.pxi4));
}
