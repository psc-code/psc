
#include "psc_push_fields_private.h"

#include "psc_acc.h"

// ----------------------------------------------------------------------
// psc_push_fields_acc_push_mflds_E

static void
psc_push_fields_acc_push_mflds_E(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "acc", JXI, HX + 3);
  int *gdims = ppsc->domain.gdims;
  if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
    acc_push_mflds_E_yz(mflds);
  } else if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
    acc_push_mflds_E_xyz(mflds);
  } else {
    assert(0);
  }
  psc_mfields_put_as(mflds, mflds_base, EX, EX + 3);
}

// ----------------------------------------------------------------------
// psc_push_fields_acc_push_mflds_H

static void
psc_push_fields_acc_push_mflds_H(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "acc", EX, HX + 3);
  int *gdims = ppsc->domain.gdims;
  if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
    acc_push_mflds_H_yz(mflds);
  } else if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
    acc_push_mflds_H_xyz(mflds);
  } else {
    assert(0);
  }
  psc_mfields_put_as(mflds, mflds_base, HX, HX + 3);
}

// ======================================================================
// psc_push_fields: subclass "acc"

struct psc_push_fields_ops psc_push_fields_acc_ops = {
  .name                  = "acc",
  .push_mflds_E          = psc_push_fields_acc_push_mflds_E,
  .push_mflds_H          = psc_push_fields_acc_push_mflds_H,
};
