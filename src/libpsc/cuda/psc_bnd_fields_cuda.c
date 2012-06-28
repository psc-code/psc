
#include "psc_bnd_fields_private.h"
#include "psc_fields_as_c.h"

#include "psc.h"
#include <mrc_profile.h>


EXTERN_C void cuda_conducting_wall_H_lo_hi_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_E_lo_hi_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_J_lo_hi_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_H_lo_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_H_hi_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_E_lo_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_E_hi_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_J_lo_y(int p, fields_cuda_t *pf);
EXTERN_C void cuda_conducting_wall_J_hi_y(int p, fields_cuda_t *pf);

// ----------------------------------------------------------------------
// psc_bnd_fields_cuda_fill_ghosts_E

static void
psc_bnd_fields_cuda_fill_ghosts_E(struct psc_bnd_fields *bnd, mfields_base_t *flds_base)
{
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    return;
  }

  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    mfields_cuda_t *flds_cuda = psc_mfields_get_cuda(flds_base, EX, EX + 3);
    psc_foreach_patch(ppsc, p) {
      int d = 1;
      if (ppsc->patch[p].off[d] == 0) {
	cuda_conducting_wall_E_lo_y(p, psc_mfields_get_patch_cuda(flds_cuda, p));
      }
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	cuda_conducting_wall_E_hi_y(p, psc_mfields_get_patch_cuda(flds_cuda, p));
      }
    }
    psc_mfields_put_cuda(flds_cuda, flds_base, EX, EX + 3);
    return;
  }
  assert(0);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_cuda_fill_ghosts_H

static void
psc_bnd_fields_cuda_fill_ghosts_H(struct psc_bnd_fields *bnd, mfields_base_t *flds_base)
{
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    return;
  }

  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    mfields_cuda_t *flds_cuda = psc_mfields_get_cuda(flds_base, HX, HX + 3);
    psc_foreach_patch(ppsc, p) {
      int d = 1;
      if (ppsc->patch[p].off[d] == 0) {
	cuda_conducting_wall_H_lo_y(p, psc_mfields_get_patch_cuda(flds_cuda, p));
      }
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	cuda_conducting_wall_H_hi_y(p, psc_mfields_get_patch_cuda(flds_cuda, p));
      }
    }
    psc_mfields_put_cuda(flds_cuda, flds_base, HX, HX + 3);
    return;
  }
  assert(0);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_cuda_add_ghosts_J

static void
psc_bnd_fields_cuda_add_ghosts_J(struct psc_bnd_fields *bnd, mfields_base_t *flds_base)
{
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    return;
  }

  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    mfields_cuda_t *flds_cuda = psc_mfields_get_cuda(flds_base, JXI, JXI + 3);
    psc_foreach_patch(ppsc, p) {
      int d = 1;
      if (ppsc->patch[p].off[d] == 0) {
	cuda_conducting_wall_J_lo_y(p, psc_mfields_get_patch_cuda(flds_cuda, p));
      }
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	cuda_conducting_wall_J_hi_y(p, psc_mfields_get_patch_cuda(flds_cuda, p));
      }
    }
    psc_mfields_put_cuda(flds_cuda, flds_base, JXI, JXI + 3);
    return;
  }
  assert(0);
}

// ======================================================================
// psc_bnd_fields: subclass "cuda"

struct psc_bnd_fields_ops psc_bnd_fields_cuda_ops = {
  .name                  = "cuda",
  .fill_ghosts_a_E       = psc_bnd_fields_cuda_fill_ghosts_E,
  .fill_ghosts_a_H       = psc_bnd_fields_cuda_fill_ghosts_H,
  // OPT fill_ghosts_b_E is probably not needed except for proper output
  // fill_ghosts_b_H: ?
  .fill_ghosts_b_E       = psc_bnd_fields_cuda_fill_ghosts_E,
  .fill_ghosts_b_H       = psc_bnd_fields_cuda_fill_ghosts_H,
  .add_ghosts_J          = psc_bnd_fields_cuda_add_ghosts_J,
};
