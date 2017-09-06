
// ----------------------------------------------------------------------
// c_prm: constant parameters

struct const_params {
  particle_real_t dqs;
  particle_real_t fnqs;
#if (DIM & DIM_X)
  particle_real_t xl;
  particle_real_t dxi;
  particle_real_t fnqxs;
#endif
#if (DIM & DIM_Y)
  particle_real_t yl;
  particle_real_t dyi;
  particle_real_t fnqys;
#endif
#if (DIM & DIM_Z)
  particle_real_t zl;
  particle_real_t dzi;
  particle_real_t fnqzs;
#endif
};

struct const_params c_prm;

static void
c_prm_set(struct psc *psc)
{
  particle_real_t dt = psc->dt;

  c_prm.dqs = .5f * ppsc->coeff.eta * dt;
  c_prm.fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;

  IF_DIM_X( c_prm.xl = .5f * dt; );
  IF_DIM_Y( c_prm.yl = .5f * dt; );
  IF_DIM_Z( c_prm.zl = .5f * dt; );

  assert(psc->nr_patches > 0);

  IF_DIM_X( c_prm.dxi = 1.f / psc->patch[0].dx[0]; );
  IF_DIM_Y( c_prm.dyi = 1.f / psc->patch[0].dx[1]; );
  IF_DIM_Z( c_prm.dzi = 1.f / psc->patch[0].dx[2]; );

  IF_DIM_X( c_prm.fnqxs = ppsc->patch[0].dx[0] * c_prm.fnqs / dt; );
  IF_DIM_Y( c_prm.fnqys = ppsc->patch[0].dx[1] * c_prm.fnqs / dt; );
  IF_DIM_Z( c_prm.fnqzs = ppsc->patch[0].dx[2] * c_prm.fnqs / dt; );
}

// ----------------------------------------------------------------------
// params_1vb

#include "psc.h"
#include "cuda_wrap.h"

#define MAX_NR_KINDS (10)

struct params_1vb {
  // particle-related
  particle_real_t dt;
  particle_real_t fnqs, fnqxs, fnqys, fnqzs;
  particle_real_t dxi[3];
  particle_real_t dq_kind[MAX_NR_KINDS];
  int b_mx[3];

  // field-related
  int mx[3];
  int ilg[3];
};

CUDA_CONSTANT static struct params_1vb prm;

static void _mrc_unused
params_1vb_set(struct psc *psc,
	       struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct params_1vb params;

  params.dt = psc->dt;
  for (int d = 0; d < 3; d++) {
    params.dxi[d] = 1.f / ppsc->patch[0].dx[d];
  }

  params.fnqs   = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
#if CALC_J == CALC_J_1VB_2D

#if !(DIM == DIM_YZ)
#error inc_params.c: CALC_J_1VB_2D only works for DIM_YZ
#endif
  params.fnqxs = params.fnqs;
#else
  params.fnqxs = ppsc->patch[0].dx[0] * params.fnqs / params.dt;
#endif
  params.fnqys = ppsc->patch[0].dx[1] * params.fnqs / params.dt;
  params.fnqzs = ppsc->patch[0].dx[2] * params.fnqs / params.dt;

  assert(psc->nr_kinds <= MAX_NR_KINDS);
  for (int k = 0; k < ppsc->nr_kinds; k++) {
    params.dq_kind[k] = .5f * ppsc->coeff.eta * params.dt * ppsc->kinds[k].q / ppsc->kinds[k].m;
  }

  if (mprts && mprts->nr_patches > 0) {
#if PSC_PARTICLES_AS_CUDA2
    struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);
    for (int d = 0; d < 3; d++) {
      params.b_mx[d] = mprts_sub->b_mx[d];
    }
#else
    assert(0);
#endif
  }

  if (mflds) {
#if PSC_FIELDS_AS_CUDA2
    struct psc_mfields_cuda2 * mflds_sub = psc_mfields_cuda2(mflds);
    for (int d = 0; d < 3; d++) {
      params.mx[d] = mflds_sub->im[d];
      params.ilg[d] = mflds_sub->ib[d];
      assert(mflds_sub->ib[d] == -2 || mflds_sub->im[d] == 1); // assumes BND == 2
    }
#else
    assert(0);
#endif
  }

#ifndef __CUDACC__
  prm = params;
#else
  check(cudaMemcpyToSymbol(prm, &params, sizeof(prm)));
#endif
}

