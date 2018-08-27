
#include "cuda_iface.h"
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"
#include "psc_bits.h"
#include "heating_spot_foil.hxx"
#include "heating_cuda_impl.hxx"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <curand_kernel.h>

#include <cstdio>

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// cuda_heating_params

struct cuda_heating_params {
  float_3 *d_xb_by_patch;
};

// ----------------------------------------------------------------------
// cuda_heating_params_set

template<typename BS>
static void cuda_heating_params_set(cuda_heating_params& h_prm, cuda_mparticles<BS>* cmprts)
{
  cudaError_t ierr;

  ierr = cudaMalloc(&h_prm.d_xb_by_patch, cmprts->n_patches * sizeof(float_3));
  cudaCheck(ierr);
  ierr = cudaMemcpy(h_prm.d_xb_by_patch, cmprts->xb_by_patch.data(),
		    cmprts->n_patches * sizeof(float_3), cudaMemcpyHostToDevice);
  cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_heating_params_free

static void cuda_heating_params_free(cuda_heating_params& h_prm)
{
  cudaError_t ierr;

  ierr = cudaFree(&h_prm.d_xb_by_patch);
  cudaCheck(ierr);
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
// k_curand_setup

__global__ static void
k_curand_setup(curandState *d_curand_states, int b_my)
{
  int bid = blockIdx.y * b_my + blockIdx.x;
  int id = threadIdx.x + bid * THREADS_PER_BLOCK;

  curand_init(1234, id % 1024, 0, &d_curand_states[id]); // FIXME, % 1024 hack
}


struct cuda_heating_foil;

template<typename BS>
__global__ static void
k_heating_run_foil(cuda_heating_foil d_foil, DMparticlesCuda<BS> dmprts, struct cuda_heating_params prm,
		   curandState *d_curand_states);

// ======================================================================
// cuda_heating_foil

struct cuda_heating_foil : HeatingSpotFoilParams
{
  cuda_heating_foil(const HeatingSpotFoilParams& params, int kind, double heating_dt)
    : HeatingSpotFoilParams(params), kind(kind), heating_dt(heating_dt),
      h_prm_{}
  {
    float width = zh - zl;
    fac = (8.f * pow(T, 1.5)) / (sqrt(Mi) * width);
  }

  ~cuda_heating_foil()
  {
#if 0 // FIXME
    cuda_heating_params_free(h_prm_);
    
    cudaError_t ierr;
    ierr = cudaFree(d_curand_states_);
    cudaCheck(ierr);
#endif
  }
  
  __host__  __device__ float get_H(float *xx)
  {
    if (xx[2] <= zl || xx[2] >= zh) {
      return 0;
    }
    
    return fac * exp(-(sqr(xx[0] - xc) +
		       sqr(xx[1] - yc)) / sqr(rH));
  }
  
  // ----------------------------------------------------------------------
  // particle_kick
  
  __host__ void particle_kick(float4 *pxi4, float H)
  {
    float2 r01 = bm_normal2();
    float2 r23 = bm_normal2();
    
    float Dp = sqrtf(H * heating_dt);
    
    pxi4->x += Dp * r01.x;
    pxi4->y += Dp * r01.y;
    pxi4->z += Dp * r23.x;
  }

  // ----------------------------------------------------------------------
  // d_particle_kick
  
  __device__ void d_particle_kick(float4 *pxi4, float H, curandState *state)
  {
    float2 r01 = curand_normal2(state);
    float r2 = curand_normal(state);
    
    float Dp = sqrtf(H * heating_dt);
    
    pxi4->x += Dp * r01.x;
    pxi4->y += Dp * r01.y;
    pxi4->z += Dp * r2;
  }

  // ----------------------------------------------------------------------
  // run_foil

  template<typename BS>
  void run_foil(cuda_mparticles<BS>* cmprts, curandState *d_curand_states)
  {
    dim3 dimGrid = BlockSimple<BS, dim_yz>::dimGrid(*cmprts);
    
    k_heating_run_foil<BS>
      <<<dimGrid, THREADS_PER_BLOCK>>>(*this, *cmprts, h_prm_, d_curand_states);
    cuda_sync_if_enabled();
  }
  
  // ----------------------------------------------------------------------
  // operator()
  
  template<typename BS>
  void operator()(cuda_mparticles<BS>* cmprts)
  {
    //return cuda_heating_run_foil_gold(cmprts);
    
    static bool first_time = true;
    if (first_time) {
      cuda_heating_params_set(h_prm_, cmprts);
      
      dim3 dimGrid = BlockSimple<BS, dim_yz>::dimGrid(*cmprts);
      int n_threads = dimGrid.x * dimGrid.y * THREADS_PER_BLOCK;
      
      cudaError_t ierr;
      ierr = cudaMalloc(&d_curand_states_, n_threads * sizeof(*d_curand_states_));
      cudaCheck(ierr);
      
      k_curand_setup<<<dimGrid, THREADS_PER_BLOCK>>>(d_curand_states_, cmprts->b_mx()[1]);
      cuda_sync_if_enabled();
      
      first_time = false;
    }
    
    if (cmprts->need_reorder) { 
      cmprts->reorder();
    }
    
    run_foil<BS>(cmprts, d_curand_states_);
  }

  // params
  int kind;

  // state (FIXME, shouldn't be part of the interface)
  float fac;
  float heating_dt;

  cuda_heating_params h_prm_;
  curandState* d_curand_states_;
};

// ----------------------------------------------------------------------
// cuda_heating_run_foil_gold

template<typename BS>
void cuda_heating_run_foil_gold(cuda_heating_foil& foil, cuda_mparticles<BS>* cmprts)
{
  for (int b = 0; b < cmprts->n_blocks; b++) {
    int p = b / cmprts->n_blocks_per_patch;
    for (int n = cmprts->d_off[b]; n < cmprts->d_off[b+1]; n++) {
      float4 xi4 = cmprts->d_xi4[n];

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

      float H = foil.get_H(xx);
      // float4 pxi4 = d_pxi4[n];
      // printf("%s xx = %g %g %g H = %g px = %g %g %g\n", (H > 0) ? "H" : " ",
      // 	     xx[0], xx[1], xx[2], H,
      // 	     pxi4.x, pxi4.y, pxi4.z);
      // pxi4.w = H;
      // d_pxi4[n] = pxi4;
      if (H > 0) {
	float4 pxi4 = cmprts->d_pxi4[n];
	foil.particle_kick(&pxi4, H);
	cmprts->d_pxi4[n] = pxi4;
	// printf("H xx = %g %g %g H = %g px = %g %g %g\n", xx[0], xx[1], xx[2], H,
	//        pxi4.x, pxi4.y, pxi4.z);
      }
    }
  }
}

// ----------------------------------------------------------------------
// k_heating_run_foil

template<typename BS>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
k_heating_run_foil(cuda_heating_foil d_foil, DMparticlesCuda<BS> dmprts, struct cuda_heating_params prm,
		   curandState *d_curand_states)
{
  BlockSimple<BS, dim_yz> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  float_3 xb; // __shared__
  xb[0] = prm.d_xb_by_patch[current_block.p][0];
  xb[1] = prm.d_xb_by_patch[current_block.p][1];
  xb[2] = prm.d_xb_by_patch[current_block.p][2];


  int id = threadIdx.x + current_block.bid * THREADS_PER_BLOCK;
  /* Copy state to local memory for efficiency */
  curandState local_state = d_curand_states[id];

  int block_begin = dmprts.off_[current_block.bid];
  int block_end = dmprts.off_[current_block.bid + 1];
  for (int n : in_block_loop(block_begin, block_end)) {
    if (n < block_begin) {
      continue;
    }
    float4 xi4 = dmprts.xi4_[n];
    
    int prt_kind = __float_as_int(xi4.w);
    if (prt_kind != d_foil.kind) {
      continue;
    }

    float xx[3] = {
      xi4.x + xb[0],
      xi4.y + xb[1],
      xi4.z + xb[2],
    };
    float H = d_foil.get_H(xx);
    //d_pxi4[n].w = H;
    if (H > 0) {
      float4 pxi4 = dmprts.pxi4_[n];
      d_foil.d_particle_kick(&pxi4, H, &local_state);
      dmprts.pxi4_[n] = pxi4;
    }
  }

  d_curand_states[id] = local_state;
}

// ======================================================================

template<typename BS>
template<typename FUNC>
HeatingCuda<BS>::HeatingCuda(const Grid_t& grid, int interval, int kind, FUNC get_H)
{
  foil_ = new cuda_heating_foil{get_H, kind, interval * grid.dt};
}

template<typename BS>
HeatingCuda<BS>::~HeatingCuda()
{
  delete foil_;
}

template<typename BS>
void HeatingCuda<BS>::operator()(MparticlesCuda<BS>& mprts)
{
  (*foil_)(mprts.cmprts());
}
  
// ======================================================================

template struct HeatingCuda<BS144>;
template HeatingCuda<BS144>::HeatingCuda(const Grid_t& grid, int interval, int kind, HeatingSpotFoil get_H);

template struct HeatingCuda<BS444>;
template HeatingCuda<BS444>::HeatingCuda(const Grid_t& grid, int interval, int kind, HeatingSpotFoil get_H);
