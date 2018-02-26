
#pragma once

#if (DIM & DIM_X)
#define IF_DIM_X(s) s do{} while(0)
#define IF_NOT_DIM_X(s) do{} while(0)
#else
#define IF_DIM_X(s) do{} while(0)
#define IF_NOT_DIM_X(s) s do{} while(0)
#endif

#if (DIM & DIM_Y)
#define IF_DIM_Y(s) s do{} while(0)
#define IF_NOT_DIM_Y(s) do{} while(0)
#else
#define IF_DIM_Y(s) do{} while(0)
#define IF_NOT_DIM_Y(s) s do{} while(0)
#endif

#if (DIM & DIM_Z)
#define IF_DIM_Z(s) s do{} while(0)
#define IF_NOT_DIM_Z(s) do{} while(0)
#else
#define IF_DIM_Z(s) do{} while(0)
#define IF_NOT_DIM_Z(s) s do{} while(0)
#endif

#if DIM == DIM_1
using opt_dim = dim_1;
#elif DIM == DIM_X
using opt_dim = dim_x;
#elif DIM == DIM_Y
using opt_dim = dim_y;
#elif DIM == DIM_Z
using opt_dim = dim_z;
#elif DIM == DIM_XY
using opt_dim = dim_xy;
#elif DIM == DIM_XZ
using opt_dim = dim_xz;
#elif DIM == DIM_YZ
using opt_dim = dim_yz;
#elif DIM == DIM_XYZ
using opt_dim = dim_xyz;
#else
#error DIM must be defined
#endif

#if ORDER == ORDER_1ST
using opt_order = opt_order_1st;
#elif ORDER == ORDER_2ND
using opt_order = opt_order_2nd;
#else
#error unknown ORDER
#endif

#if ORDER == ORDER_1ST
#if IP_VARIANT == IP_VARIANT_EC
using opt_ip = opt_ip_1st_ec;
#else
using opt_ip = opt_ip_1st;
#endif
#elif ORDER == ORDER_2ND
using opt_ip = opt_ip_2nd;
#endif

#if CALC_J == CALC_J_1VB_SPLIT
using opt_calcj = opt_calcj_1vb_split;
#elif CALC_J == CALC_J_1VB_VAR1
using opt_calcj = opt_calcj_1vb_var1;
#elif CALC_J == CALC_J_1VB_2D
using opt_calcj = opt_calcj_1vb_2d;
#endif

#define CUDA_CONSTANT
#define CUDA_DEVICE
#define __forceinline__
#define atomicAdd(addr, val) \
  do { *(addr) += (val); } while (0)


// ----------------------------------------------------------------------
// c_prm: constant parameters

#define MAX_NR_KINDS (10)

struct const_params {
  real_t dt; // FIXME, do we need both dt and dqs? or maybe get rid of xl/yl/zl
  real_t dqs;
  real_t fnqs;
  real_t fnqxs, fnqys, fnqzs;
  real_t dxi[3];
};

struct params_1vb {
  // particle-related
  real_t dq_kind[MAX_NR_KINDS];
  int b_mx[3];

  // field-related
  int mx[3];
  int ilg[3];
};

CUDA_CONSTANT static struct const_params c_prm;
CUDA_CONSTANT static struct params_1vb prm;

// ----------------------------------------------------------------------
// c_prm_set

static void
c_prm_set(const Grid_t& grid)
{
  struct const_params prm;

  prm.dt = grid.dt;
  prm.dqs = .5f * grid.eta * prm.dt;
  prm.fnqs = grid.fnqs;

  assert(grid.n_patches() > 0);

  for (int d = 0; d < 3; d++) {
    prm.dxi[d] = 1.f / grid.dx[d];
  }

  prm.fnqxs = grid.dx[0] * grid.fnqs / grid.dt;
  prm.fnqys = grid.dx[1] * grid.fnqs / grid.dt;
  prm.fnqzs = grid.dx[2] * grid.fnqs / grid.dt;

#ifndef __CUDACC__
  c_prm = prm;
#else
  check(cudaMemcpyToSymbol(c_prm, &prm, sizeof(prm)));
#endif
}

// ----------------------------------------------------------------------
// params_1vb

static void _mrc_unused
params_1vb_set(const Grid_t& grid)
{
  auto& kinds = grid.kinds;
  struct params_1vb params;

#if CALC_J == CALC_J_1VB_2D && DIM != DIM_YZ
#error inc_params.c: CALC_J_1VB_2D only works for DIM_YZ
#endif

  assert(kinds.size() <= MAX_NR_KINDS);
  for (int k = 0; k < kinds.size(); k++) {
    params.dq_kind[k] = .5f * grid.eta * grid.dt * kinds[k].q / kinds[k].m;
  }

#ifndef __CUDACC__
  prm = params;
#else
  check(cudaMemcpyToSymbol(prm, &params, sizeof(prm)));
#endif
}

