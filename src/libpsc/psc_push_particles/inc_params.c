
#include "cuda_wrap.h"

#define MAX_NR_KINDS (10)

struct params_1vb {
  particle_real_t dt;
  particle_real_t fnqs, fnqxs, fnqys, fnqzs;
  particle_real_t dxi[3];
  particle_real_t dq_kind[MAX_NR_KINDS];
#ifdef __CUDACC__
  int mx[3];
  int ilg[3];
  int b_mx[3];
#endif
};

CUDA_CONSTANT static struct params_1vb prm;

static void
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

#ifdef __CUDACC__
  if (mprts && mprts->nr_patches > 0) {
    struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);
    for (int d = 0; d < 3; d++) {
      params.b_mx[d] = mprts_sub->b_mx[d];
    }
  }

  if (mflds) {
    struct psc_mfields_cuda2 * mflds_sub = psc_mfields_cuda2(mflds);
    for (int d = 0; d < 3; d++) {
      params.mx[d] = mflds_sub->im[d];
      params.ilg[d] = mflds_sub->ib[d];
      assert(mflds_sub->ib[d] == -2 || mflds_sub->im[d] == 1); // assumes BND == 2
    }
  }
#endif

#ifndef __CUDACC__
  prm = params;
#else
  check(cudaMemcpyToSymbol(prm, &params, sizeof(prm)));
#endif
}

