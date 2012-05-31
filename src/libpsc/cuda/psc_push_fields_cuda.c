
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "psc.h"
#include "psc_bnd_fields.h"

EXTERN_C void cuda_push_fields_E_yz(int p, struct psc_fields *flds);
EXTERN_C void cuda_push_fields_H_yz(int p, struct psc_fields *flds);

// FIXME, the logic is duplicated from push_fields "single2",
// and somehow should be generalized

// ----------------------------------------------------------------------
// step_E

static void
step_E(struct psc_mfields *mflds)
{
  psc_stats_start(st_time_field);
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds_base = psc_mfields_get_patch(mflds, p);
    struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", JXI, HX + 3);
    if (ppsc->domain.gdims[0] == 1) {
      cuda_push_fields_E_yz(flds->p, flds);
    } else {
      assert(0);
    }
    psc_fields_put_as(flds, flds_base, EX, EX + 3);
  }
  psc_stats_stop(st_time_field);
}

// ----------------------------------------------------------------------
// step_H

static void
step_H(struct psc_mfields *mflds)
{
  psc_stats_start(st_time_field);
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds_base = psc_mfields_get_patch(mflds, p);
    struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, HX + 3);
    if (ppsc->domain.gdims[0] == 1) {
      cuda_push_fields_H_yz(flds->p, flds);
    } else {
      assert(0);
    }
    psc_fields_put_as(flds, flds_base, HX, HX + 3);
  }
  psc_stats_stop(st_time_field);
}

// ----------------------------------------------------------------------
// psc_push_fields_cuda_step_a

static void
psc_push_fields_cuda_step_a(struct psc_push_fields *push, struct psc_mfields *mflds)
{
  step_E(mflds);
  psc_bnd_fields_fill_ghosts_a_E(push->bnd_fields, mflds);
  step_H(mflds);
  psc_bnd_fields_fill_ghosts_a_H(push->bnd_fields, mflds);
}

// ----------------------------------------------------------------------
// psc_push_fields_cuda_step_b_H

static void
psc_push_fields_cuda_step_b_H(struct psc_push_fields *push, struct psc_mfields *mflds)
{
  step_H(mflds);
}

// ----------------------------------------------------------------------
// psc_push_fields_cuda_step_b_E

static void
psc_push_fields_cuda_step_b_E(struct psc_push_fields *push, struct psc_mfields *mflds)
{
  // add and fill ghost for J
  psc_bnd_fields_add_ghosts_J(push->bnd_fields, mflds);
  psc_bnd_add_ghosts(ppsc->bnd, mflds, JXI, JXI + 3);
  psc_bnd_fill_ghosts(ppsc->bnd, mflds, JXI, JXI + 3);
  // fill ghosts for H
  psc_bnd_fields_fill_ghosts_b_H(push->bnd_fields, mflds);
  psc_bnd_fill_ghosts(ppsc->bnd, mflds, HX, HX + 3);
  
  step_E(mflds);
  
  psc_bnd_fields_fill_ghosts_b_E(push->bnd_fields, mflds);
}

// ======================================================================
// psc_push_fields: subclass "cuda"

struct psc_push_fields_ops psc_push_fields_cuda_ops = {
  .name                  = "cuda",
  .step_a                = psc_push_fields_cuda_step_a,
  .step_b_H              = psc_push_fields_cuda_step_b_H,
  .step_b_E              = psc_push_fields_cuda_step_b_E,
};
