
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

// ----------------------------------------------------------------------
// psc_push_fields_cuda_push_mflds_E

static void
psc_push_fields_cuda_push_mflds_E(struct psc_push_fields *push, struct psc_mfields *mflds_base,
				  double dt_fac)
{
  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(JXI, HX + 3);
  
  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_E_yz(mf->cmflds, dt_fac * ppsc->dt);
  } else {
    assert(0);
  }

  mf.put_as(mflds_base, EX, EX + 3);
}

// ----------------------------------------------------------------------
// psc_push_fields_cuda_push_mflds_H

static void
psc_push_fields_cuda_push_mflds_H(struct psc_push_fields *push, struct psc_mfields *mflds_base,
				  double dt_fac)
{
  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(JXI, HX + 3);

  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_H_yz(mf->cmflds, dt_fac * ppsc->dt);
  } else {
    assert(0);
  }

  mf.put_as(mflds_base, HX, HX + 3);
}

// ======================================================================
// psc_push_fields: subclass "cuda"

struct psc_push_fields_ops_cuda : psc_push_fields_ops {
  psc_push_fields_ops_cuda() {
    name                  = "cuda";
    push_mflds_E          = psc_push_fields_cuda_push_mflds_E;
    push_mflds_H          = psc_push_fields_cuda_push_mflds_H;
  }
} psc_push_fields_cuda_ops;
