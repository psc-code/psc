
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
set_params(struct cuda_params *prm, struct psc *psc,
	   struct cuda_mparticles *cmprts)
{
  prm->dt = psc->dt;
  for (int d = 0; d < 3; d++) {
    prm->dxi[d] = 1.f / psc->patch[0].dx[d];
  }

  prm->dqs    = .5f * psc->coeff.eta * psc->dt;
  prm->fnqs   = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
  prm->fnqxs  = psc->patch[0].dx[0] * prm->fnqs / psc->dt;
  prm->fnqys  = psc->patch[0].dx[1] * prm->fnqs / psc->dt;
  prm->fnqzs  = psc->patch[0].dx[2] * prm->fnqs / psc->dt;
  assert(psc->nr_kinds <= MAX_KINDS);
  for (int k = 0; k < psc->nr_kinds; k++) {
    prm->dq[k] = prm->dqs * psc->kinds[k].q / psc->kinds[k].m;
  }

  if (cmprts) {
    for (int d = 0; d < 3; d++) {
      prm->b_mx[d] = cmprts->b_mx[d];
      prm->b_dxi[d] = cmprts->b_dxi[d];
    }
  }
}

void
free_params(struct cuda_params *prm)
{
}

