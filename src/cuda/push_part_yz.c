
#include "psc_cuda.h"
#include "math.h"
#include "profile/profile.h"

__constant__ static real _dt;

static void set_constants()
{
  real __dt = psc.dt;
  check(cudaMemcpyToSymbol(_dt, &__dt, sizeof(_dt)));
}

__device__ static void
calc_vxi(real vxi[3], struct d_particle p)
{
  real root = rsqrtr(real(1.) + sqr(p.pxi[0]) + sqr(p.pxi[1]) + sqr(p.pxi[2]));

  for (int d = 1; d < 3; d++) {
    vxi[d] = p.pxi[d] * root;
  }
}

__device__ static void
push_xi_halfdt(struct d_particle *p, const real vxi[3])
{
  for (int d = 1; d < 3; d++) {
    p->xi[d] += real(.5) * _dt * vxi[d];
  }
}

__device__ static void
push_part_yz_a_one(int n, struct d_part d_particles)
{
  struct d_particle p;
  LOAD_PARTICLE(p, d_particles, n);
  real vxi[3];

  // x^n, p^n -> x^(n+0.5), p^n
  
  calc_vxi(vxi, p);
  push_xi_halfdt(&p, vxi);

  STORE_PARTICLE_POS(p, d_particles, n);
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
    pr = prof_register("cuda_part_yz_a", 1., 0, psc.n_part * 12 * sizeof(float));
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
