
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "psc_cuda.h"
#include "particles_cuda.h"

void
cuda_mfields_params_set(struct cuda_mfields_params *mflds_prm,
			struct cuda_mfields *cmflds)
{
  for (int d = 0; d < 3; d++) {
    mflds_prm->mx[d] = cmflds->im[d];
    mflds_prm->ilg[d] = cmflds->ib[d];
    if (d != 0) {
      assert(cmflds->ib[d] == -BND);
    } else {
      assert(cmflds->im[d] == 1);
    }
  }
}

void
cuda_mparticles_params_set(struct cuda_mparticles_params *mprts_prm,
			   struct cuda_mparticles *cmprts, struct psc *psc)
{
  mprts_prm->fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;

  for (int d = 0; d < 3; d++) {
    mprts_prm->b_mx[d] = cmprts->b_mx[d];
    mprts_prm->b_dxi[d] = cmprts->b_dxi[d];
    mprts_prm->dxi[d] = 1.f / cmprts->dx[d];
  }
}

