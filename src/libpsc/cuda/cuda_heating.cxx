
#include "cuda_iface.h"
#include "cuda_mparticles.hxx"
#include "cuda_bits.h"
#include "cuda_base.hxx"
#include "psc_bits.h"
#include "heating_spot_foil.hxx"
#include "heating_cuda_impl.hxx"
#include "balance.hxx"
#include "rng_state.hxx"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <cstdio>

#define THREADS_PER_BLOCK 128

extern std::size_t mem_heating;

using Float3 = Vec3<float>;

// ----------------------------------------------------------------------
// bm_normal2

static inline float2 bm_normal2(void)
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
// d_particle_kick

__device__ void d_particle_kick(DParticleCuda& prt, float H, float heating_dt,
                                RngStateCuda::Rng& rng)
{
  float2 r01 = rng.normal2();
  float r2 = rng.normal();

  float Dp = sqrtf(H * heating_dt);

  prt.u[0] += Dp * r01.x;
  prt.u[1] += Dp * r01.y;
  prt.u[2] += Dp * r2;
}

// ----------------------------------------------------------------------
// k_heating_run_foil

template <typename BS, typename HS>
__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
  k_heating_run_foil(HS foil, DMparticlesCuda<BS> dmprts, float heating_dt,
                     Float3* d_xb_by_patch, RngStateCuda::Device rng_state,
                     int n_blocks)
{
  BlockSimple2<BS, typename HS::dim> current_block;

  /* Copy state to local memory for efficiency */
  int bid = blockIdx.x;
  int id = threadIdx.x + bid * blockDim.x;

  auto local_rng = rng_state[id];

  for (; bid < n_blocks; bid += gridDim.x) {
    current_block.init(dmprts, bid);

    Float3 xb; // __shared__
    xb[0] = d_xb_by_patch[current_block.p][0];
    xb[1] = d_xb_by_patch[current_block.p][1];
    xb[2] = d_xb_by_patch[current_block.p][2];

    int block_begin = dmprts.off_[current_block.bid];
    int block_end = dmprts.off_[current_block.bid + 1];
    for (int n : in_block_loop(block_begin, block_end)) {
      if (n < block_begin) {
        continue;
      }
      auto prt = dmprts.storage[n];

      float xx[3] = {
        prt.x[0] + xb[0],
        prt.x[1] + xb[1],
        prt.x[2] + xb[2],
      };
      float H = foil(xx, prt.kind);
      if (H > 0.f) {
        d_particle_kick(prt, H, heating_dt, local_rng);
        dmprts.storage.store_momentum(prt, n);
      }
    }
  }

  rng_state[id] = local_rng;
}

// ======================================================================
// cuda_heating_foil

template <typename HS>
struct cuda_heating_foil
{
  cuda_heating_foil(const Grid_t& grid, const HS& heating_spot,
                    double heating_dt)
    : heating_dt_(heating_dt),
      heating_spot_{heating_spot},
      first_time_{true},
      rng_state_{get_rng_state()}
  {}

  // no copy constructor / assign, to catch performance issues
  cuda_heating_foil(const cuda_heating_foil&) = delete;
  cuda_heating_foil& operator=(const cuda_heating_foil&) = delete;

  ~cuda_heating_foil() { mem_heating -= allocated_bytes(d_xb_by_patch_); }

  void reset() { first_time_ = true; }

  // ----------------------------------------------------------------------
  // operator()

  template <typename BS>
  void operator()(cuda_mparticles<BS>* cmprts)
  {
    prof_barrier("heating start");

    // return cuda_heating_run_foil_gold(cmprts);
    if (cmprts->n_prts == 0) {
      return;
    }

    dim3 dimGrid = BlockSimple2<BS, typename HS::dim>::dimGrid(*cmprts);

    if (first_time_) { // FIXME
      mem_heating -= allocated_bytes(d_xb_by_patch_);
      d_xb_by_patch_ = cmprts->xb_by_patch;
      mem_heating += allocated_bytes(d_xb_by_patch_);

      rng_state_.resize(dimGrid.x * dimGrid.y * dimGrid.z * THREADS_PER_BLOCK);

      prof_barrier("heating first");
      first_time_ = false;
    }

    if (cmprts->need_reorder) {
      cmprts->reorder();
    }

    int n_blocks = cmprts->b_mx()[0] * cmprts->b_mx()[1] * cmprts->b_mx()[2] *
                   cmprts->n_patches();
    k_heating_run_foil<BS><<<dimGrid, THREADS_PER_BLOCK>>>(
      heating_spot_, *cmprts, heating_dt_, d_xb_by_patch_.data().get(),
      rng_state_, n_blocks);
    cuda_sync_if_enabled();
  }

  // state (FIXME, shouldn't be part of the interface)
  bool first_time_;
  float heating_dt_;
  HS heating_spot_;

  psc::device_vector<Float3> d_xb_by_patch_;
  RngStateCuda& rng_state_;
};

// ----------------------------------------------------------------------
// particle_kick

__host__ void particle_kick(DParticleCuda& prt, float H, float heating_dt)
{
  float2 r01 = bm_normal2();
  float2 r23 = bm_normal2();

  float Dp = sqrtf(H * heating_dt);

  prt.x[0] += Dp * r01.x;
  prt.x[1] += Dp * r01.y;
  prt.x[2] += Dp * r23.x;
}

// ----------------------------------------------------------------------
// cuda_heating_run_foil_gold

template <typename BS, typename HS>
void cuda_heating_run_foil_gold(HS& foil, float heating_dt,
                                cuda_mparticles<BS>* cmprts)
{
  for (int b = 0; b < cmprts->n_blocks; b++) {
    int p = b / cmprts->n_blocks_per_patch;
    for (int n = cmprts->d_off[b]; n < cmprts->d_off[b + 1]; n++) {
      auto prt = cmprts->storage[n];

      float* xb = &cmprts->xb_by_patch[p][0];
      float xx[3] = {
        prt.x[0] + xb[0],
        prt.x[1] + xb[1],
        prt.x[2] + xb[2],
      };

      float H = foil(xx, prt.kind);
      // float4 pxi4 = d_pxi4[n];
      // printf("%s xx = %g %g %g H = %g px = %g %g %g\n", (H > 0) ? "H" : " ",
      // 	     xx[0], xx[1], xx[2], H,
      // 	     pxi4.x, pxi4.y, pxi4.z);
      // pxi4.w = H;
      // d_pxi4[n] = pxi4;
      if (H > 0) {
        auto prt = cmprts->storage[n];
        particle_kick(prt, H, heating_dt);
        cmprts->storage.store_momentum(prt, n);
        // printf("H xx = %g %g %g H = %g px = %g %g %g\n", xx[0], xx[1], xx[2],
        // H,
        //        pxi4.x, pxi4.y, pxi4.z);
      }
    }
  }
}

// ======================================================================

template <typename HS, typename MP>
HeatingCuda<HS, MP>::HeatingCuda(const Grid_t& grid, int interval,
                                 HS heating_spot)
  : foil_{new cuda_heating_foil<HS>{grid, heating_spot, interval * grid.dt}},
    balance_generation_cnt_{-1}
{}

template <typename HS, typename MP>
HeatingCuda<HS, MP>::~HeatingCuda()
{
  delete foil_;
}

template <typename HS, typename MP>
void HeatingCuda<HS, MP>::reset(const MP& mprts)
{
  foil_->reset();
}

template <typename HS, typename MP>
void HeatingCuda<HS, MP>::operator()(MP& mprts)
{
  if (psc_balance_generation_cnt > this->balance_generation_cnt_) {
    balance_generation_cnt_ = psc_balance_generation_cnt;
    reset(mprts);
  }

  (*foil_)(mprts.cmprts());
}

// ======================================================================

template struct HeatingCuda<HeatingSpotFoil<dim_yz>, MparticlesCuda<BS144>>;
template struct HeatingCuda<HeatingSpotFoil<dim_xyz>, MparticlesCuda<BS444>>;
