
#include "psc_push_fields_private.h"

#include "psc_cuda2.h"

// ----------------------------------------------------------------------
// psc_push_fields_cuda2_push_mflds_E

static void
psc_push_fields_cuda2_push_mflds_E(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda2", JXI, HX + 3);
  if (ppsc->domain.gdims[0] == 1) {
    //    cuda2_push_mflds_E_yz_gold(mflds_cuda);
    cuda2_push_mflds_E_yz(mflds);
  } else {
    assert(0);
  }
  psc_mfields_put_as(mflds, mflds_base, EX, EX + 3);
}

// ----------------------------------------------------------------------
// psc_push_fields_cuda2_push_mflds_H

static void
psc_push_fields_cuda2_push_mflds_H(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda2", EX, HX + 3);
  if (ppsc->domain.gdims[0] == 1) {
    cuda2_push_mflds_H_yz(mflds);
  } else {
    assert(0);
  }
  psc_mfields_put_as(mflds, mflds_base, HX, HX + 3);
}

// ======================================================================
// psc_push_fields: subclass "cuda2"

struct psc_push_fields_ops psc_push_fields_cuda2_ops = {
  .name                  = "cuda2",
  .push_mflds_E          = psc_push_fields_cuda2_push_mflds_E,
  .push_mflds_H          = psc_push_fields_cuda2_push_mflds_H,
};
