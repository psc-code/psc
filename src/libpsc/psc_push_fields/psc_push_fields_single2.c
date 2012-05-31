
#include "psc_push_fields_private.h"
#include "psc_fields_as_single.h"

// This version "single2" avoids a lot of the communication / boundary filling.
// In it's current form, it'll only fill enough ghost points for 1st order
// particles

#include "psc_push_fields_common.c"

#include "psc_bnd.h"
#include "psc_bnd_fields.h"

static void
psc_push_fields_sub_step_a(struct psc_push_fields *push, struct psc_mfields *mflds)
{
  psc_stats_start(st_time_field);
  for (int p = 0; p < mflds->nr_patches; p++) {
    psc_push_fields_sub_push_a_E(push, psc_mfields_get_patch(mflds, p));
  }
  psc_stats_stop(st_time_field);
  
  psc_bnd_fields_fill_ghosts_a_E(push->bnd_fields, mflds);
  psc_bnd_fill_ghosts(ppsc->bnd, mflds, EX, EX + 3);
  
  psc_stats_start(st_time_field);
  for (int p = 0; p < mflds->nr_patches; p++) {
    psc_push_fields_sub_push_a_H(push, psc_mfields_get_patch(mflds, p));
  }
  psc_stats_stop(st_time_field);
  
  psc_bnd_fields_fill_ghosts_a_H(push->bnd_fields, mflds);
}

static void
psc_push_fields_sub_step_b(struct psc_push_fields *push, struct psc_mfields *mflds)
{
  psc_stats_start(st_time_field);
  for (int p = 0; p < mflds->nr_patches; p++) {
    psc_push_fields_sub_push_b_H(push, psc_mfields_get_patch(mflds, p));
  }
  psc_stats_stop(st_time_field);

  // add and fill ghost for J
  psc_bnd_fields_add_ghosts_J(push->bnd_fields, mflds);
  psc_bnd_add_ghosts(ppsc->bnd, mflds, JXI, JXI + 3);
  psc_bnd_fill_ghosts(ppsc->bnd, mflds, JXI, JXI + 3);
  // fill ghosts for H
  psc_bnd_fields_fill_ghosts_b_H(push->bnd_fields, mflds);
  psc_bnd_fill_ghosts(ppsc->bnd, mflds, HX, HX + 3);
  
  psc_stats_start(st_time_field);
  for (int p = 0; p < mflds->nr_patches; p++) {
    psc_push_fields_sub_push_b_E(push, psc_mfields_get_patch(mflds, p));
  }
  psc_stats_stop(st_time_field);
  
  psc_bnd_fields_fill_ghosts_b_E(push->bnd_fields, mflds);
  psc_bnd_fill_ghosts(ppsc->bnd, mflds, EX, EX + 3);
}

// ======================================================================
// psc_push_fields: subclass "single2"

struct psc_push_fields_ops psc_push_fields_single2_ops = {
  .name                  = "single2",
  .step_a                = psc_push_fields_sub_step_a,
  .step_b                = psc_push_fields_sub_step_b,
};
