
#include "psc_push_fields_private.h"

#include "psc_cuda2.h"

#define GOLD

// ----------------------------------------------------------------------
// psc_push_fields_cuda2_push_mflds_E

static void
psc_push_fields_cuda2_push_mflds_E(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda2", JXI, HX + 3);
  int *gdims = ppsc->domain.gdims;
  if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
#ifdef GOLD
    cuda2_push_mflds_E_yz_gold(mflds);
#else
    psc_mfields_cuda2_copy_to_device(mflds);
    cuda2_push_mflds_E_yz(mflds);
    psc_mfields_cuda2_copy_to_host(mflds);
#endif
  } else if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
#ifdef GOLD
    cuda2_push_mflds_E_xyz_gold(mflds);
#else
    assert(0);
#endif
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
  int *gdims = ppsc->domain.gdims;
  if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
#ifdef GOLD
    cuda2_push_mflds_H_yz_gold(mflds);
#else
    psc_mfields_cuda2_copy_to_device(mflds);
    cuda2_push_mflds_H_yz(mflds);
    psc_mfields_cuda2_copy_to_host(mflds);
#endif
  } else if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
#ifdef GOLD
    cuda2_push_mflds_H_xyz_gold(mflds);
#else
    assert(0);
#endif
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
