
#include "cuda_iface.h"
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_mparticles_const.h"
#include "cuda_bits.h"
#include "psc_bits.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <curand_kernel.h>

#include <cstdio>

#define THREADS_PER_BLOCK 256

static cuda_heating_foil foil;
__device__ static cuda_heating_foil d_foil; // FIXME, could use const memory

// ----------------------------------------------------------------------
// cuda_heating_params

struct cuda_heating_params {
  float_3 *d_xb_by_patch;
};

static cuda_heating_params h_prm;

// ----------------------------------------------------------------------
// cuda_heating_params_set

static void
cuda_heating_params_set(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  ierr = cudaMalloc(&h_prm.d_xb_by_patch, cmprts->n_patches * sizeof(float_3));
  cudaCheck(ierr);
  ierr = cudaMemcpy(h_prm.d_xb_by_patch, cmprts->xb_by_patch,
		    cmprts->n_patches * sizeof(float_3), cudaMemcpyHostToDevice);
  cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_heating_params_free

static void
cuda_heating_params_free()
{
  cudaError_t ierr;

  ierr = cudaFree(&h_prm.d_xb_by_patch);
  cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// foil_get_H

static float
foil_get_H(float *xx)
{
  if (xx[2] <= foil.zl || xx[2] >= foil.zh) {
    return 0;
  }

  return foil.fac * exp(-(sqr(xx[0] - foil.xc) +
			  sqr(xx[1] - foil.yc)) / sqr(foil.rH));
}

// ----------------------------------------------------------------------
// d_foil_get_H

__device__ static float
d_foil_get_H(float *xx)
{
  if (xx[2] <= d_foil.zl || xx[2] >= d_foil.zh) {
    return 0;
  }

  return d_foil.fac * exp(-(sqr(xx[0] - d_foil.xc) +
			    sqr(xx[1] - d_foil.yc)) / sqr(d_foil.rH));
}

// ----------------------------------------------------------------------
// bm_normal2

static inline float2
bm_normal2(void)
{
  float u1, u2;
  do {
    u1 = random() * (1.f / RAND_MAX);
    u2 = random() * (1.f / RAND_MAX);
  } while (u1 <= 0.f);

  float2 rv;
  rv.x = sqrtf(-2.f * logf(u1)) * cosf(2.f * M_PI * u2);
  rv.y = sqrtf(-2.f * logf(u1)) * sinf(2.f * M_PI * u2);
  return rv;
}

// ----------------------------------------------------------------------
// particle_kick

static void
particle_kick(float4 *pxi4, float H)
{
  float2 r01 = bm_normal2();
  float2 r23 = bm_normal2();

  float Dp = sqrtf(H * foil.heating_dt);

  pxi4->x += Dp * r01.x;
  pxi4->y += Dp * r01.y;
  pxi4->z += Dp * r23.x;
}

// ----------------------------------------------------------------------
// d_particle_kick

__device__ static void
d_particle_kick(float4 *pxi4, float H, curandState *state)
{
  float2 r01 = curand_normal2(state);
  float r2 = curand_normal(state);

  float Dp = sqrtf(H * d_foil.heating_dt);

  pxi4->x += Dp * r01.x;
  pxi4->y += Dp * r01.y;
  pxi4->z += Dp * r2;
}

// ----------------------------------------------------------------------
// cuda_heating_run_foil_gold

void
cuda_heating_run_foil_gold(struct cuda_mparticles *cmprts)
{
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);
  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);

  for (int b = 0; b < cmprts->n_blocks; b++) {
    int p = b / cmprts->n_blocks_per_patch;
    for (int n = d_off[b]; n < d_off[b+1]; n++) {
      float4 xi4 = d_xi4[n];

      int prt_kind = cuda_float_as_int(xi4.w);
      if (prt_kind != foil.kind) {
	continue;
      }

      float *xb = &cmprts->xb_by_patch[p][0];
      float xx[3] = {
	xi4.x + xb[0],
	xi4.y + xb[1],
	xi4.z + xb[2],
      };

      float H = foil_get_H(xx);
      // float4 pxi4 = d_pxi4[n];
      // printf("%s xx = %g %g %g H = %g px = %g %g %g\n", (H > 0) ? "H" : " ",
      // 	     xx[0], xx[1], xx[2], H,
      // 	     pxi4.x, pxi4.y, pxi4.z);
      // pxi4.w = H;
      // d_pxi4[n] = pxi4;
      if (H > 0) {
	float4 pxi4 = d_pxi4[n];
	particle_kick(&pxi4, H);
	d_pxi4[n] = pxi4;
	// printf("H xx = %g %g %g H = %g px = %g %g %g\n", xx[0], xx[1], xx[2], H,
	//        pxi4.x, pxi4.y, pxi4.z);
      }
    }
  }
}

// ----------------------------------------------------------------------
// k_heating_run_foil

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
k_heating_run_foil(struct cuda_heating_params prm, float4 *d_xi4, float4 *d_pxi4,
		   unsigned int *d_off, curandState *d_curand_states)
{
  int block_pos[3], ci0[3];
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(block_pos, ci0);
  int bid = find_bid();
  int id = threadIdx.x + bid * THREADS_PER_BLOCK;
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  float_3 xb; // __shared__
  xb[0] = prm.d_xb_by_patch[p][0];
  xb[1] = prm.d_xb_by_patch[p][1];
  xb[2] = prm.d_xb_by_patch[p][2];

  /* Copy state to local memory for efficiency */
  curandState local_state = d_curand_states[id];

  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    float4 xi4 = d_xi4[n];
    
    int prt_kind = __float_as_int(xi4.w);
    if (prt_kind != d_foil.kind) {
      continue;
    }

    float xx[3] = {
      xi4.x + xb[0],
      xi4.y + xb[1],
      xi4.z + xb[2],
    };
    float H = d_foil_get_H(xx);
    //d_pxi4[n].w = H;
    if (H > 0) {
      float4 pxi4 = d_pxi4[n];
      d_particle_kick(&pxi4, H, &local_state);
      d_pxi4[n] = pxi4;
    }
  }

  d_curand_states[id] = local_state;
}

// ----------------------------------------------------------------------
// k_curand_setup

__global__ static void
k_curand_setup(curandState *d_curand_states, int b_my)
{
  int bid = blockIdx.y * b_my + blockIdx.x;
  int id = threadIdx.x + bid * THREADS_PER_BLOCK;

  curand_init(1234, id % 1024, 0, &d_curand_states[id]); // FIXME, % 1024 hack
}

// ----------------------------------------------------------------------
// heating_run_foil

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
heating_run_foil(struct cuda_mparticles *cmprts, curandState *d_curand_states)
{
  dim3 dimGrid(cmprts->b_mx[1], cmprts->b_mx[2] * cmprts->n_patches);

  k_heating_run_foil<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (h_prm, cmprts->d_xi4, cmprts->d_pxi4, cmprts->d_off,
       d_curand_states);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// cuda_heating_setup_foil

void
cuda_heating_setup_foil(struct cuda_heating_foil *_foil)
{
  foil = *_foil;

  float width = foil.zh - foil.zl;
  foil.fac = (8.f * pow(foil.T, 1.5)) / (sqrt(foil.Mi) * width);

  cudaError_t ierr;
  ierr = cudaMemcpyToSymbol(d_foil, &foil, sizeof(d_foil), 0,
			    cudaMemcpyHostToDevice);
  cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_heating_run_foil

void
cuda_heating_run_foil(struct cuda_mparticles *cmprts)
{
  printf("cuda_heating_run_foil\n");
  //return cuda_heating_run_foil_gold(cmprts);

  static bool first_time = true;
  static curandState *d_curand_states;
  if (first_time) {
    cuda_mparticles_const_set(cmprts);
    cuda_heating_params_set(cmprts);

    dim3 dimGrid(cmprts->b_mx[1], cmprts->b_mx[2] * cmprts->n_patches);
    
    cudaError_t ierr;
    int n_threads = dimGrid.x * dimGrid. y * THREADS_PER_BLOCK;
    ierr = cudaMalloc(&d_curand_states, n_threads * sizeof(*d_curand_states));
    cudaCheck(ierr);

    k_curand_setup<<<dimGrid, THREADS_PER_BLOCK>>>(d_curand_states, cmprts->b_mx[1]);
    cuda_sync_if_enabled();

    first_time = false;
  }

  if (cmprts->need_reorder) { 
    cuda_mparticles_reorder(cmprts);
    cmprts->need_reorder = false;
  }

  if (cmprts->bs[0] == 1 && cmprts->bs[1] == 2 && cmprts->bs[2] == 2) {
    heating_run_foil<1, 2, 2>(cmprts, d_curand_states);
  } else if (cmprts->bs[0] == 1 && cmprts->bs[1] == 4 && cmprts->bs[2] == 4) {
    heating_run_foil<1, 4, 4>(cmprts, d_curand_states);
  } else {
    assert(0);
  }
  
  if (0) {
    cuda_heating_params_free();
    
    cudaError_t ierr;
    ierr = cudaFree(d_curand_states);
    cudaCheck(ierr);
  }
}

