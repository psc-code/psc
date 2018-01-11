
#include "psc_bnd_fields_private.h"

#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "psc_fields_cuda.h"

#include "psc.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_bnd_fields_cuda_fill_ghosts_E

static void
psc_bnd_fields_cuda_fill_ghosts_E(struct psc_bnd_fields *bnd, struct psc_mfields *mflds_base)
{
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    return;
  }

  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(EX, EX + 3);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mf.mflds())->cmflds;
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    int d = 1;
    for (int p = 0; p < mf.n_patches(); p++) {
      if (ppsc->patch[p].off[d] == 0) {
	cuda_conducting_wall_E_lo_y(cmflds, p);
      }
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	cuda_conducting_wall_E_hi_y(cmflds, p);
      }
    }
  } else {
    assert(0);
  }
  mf.put_as(mflds_base, EX, EX + 3);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_cuda_fill_ghosts_H

static void
psc_bnd_fields_cuda_fill_ghosts_H(struct psc_bnd_fields *bnd, struct psc_mfields *mflds_base)
{
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    return;
  }

  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(HX, HX + 3);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mf.mflds())->cmflds;
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    int d = 1;
    for (int p = 0; p < mf.n_patches(); p++) {
      if (ppsc->patch[p].off[d] == 0) {
	cuda_conducting_wall_H_lo_y(cmflds, p);
      }
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	cuda_conducting_wall_H_hi_y(cmflds, p);

      }
    }
  } else {
    assert(0);
  }
  mf.put_as(mflds_base, HX, HX + 3);
}

// ----------------------------------------------------------------------
// psc_bnd_fields_cuda_add_ghosts_J

static void
psc_bnd_fields_cuda_add_ghosts_J(struct psc_bnd_fields *bnd, struct psc_mfields *mflds_base)
{
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    return;
  }

  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(JXI, JXI + 3);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mf.mflds())->cmflds;
  if (ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    int d = 1;
    for (int p = 0; p < mf.n_patches(); p++) {
      if (ppsc->patch[p].off[d] == 0) {
	cuda_conducting_wall_J_lo_y(cmflds, p);
      }
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	cuda_conducting_wall_J_hi_y(cmflds, p);
      }
    }
  } else {
    assert(0);
  }
  mf.put_as(mflds_base, JXI, JXI + 3);
}

// ======================================================================
// psc_bnd_fields: subclass "cuda"

struct psc_bnd_fields_ops_cuda : psc_bnd_fields_ops {
  psc_bnd_fields_ops_cuda() {
    name                  = "cuda";
    fill_ghosts_E         = psc_bnd_fields_cuda_fill_ghosts_E;
    fill_ghosts_H         = psc_bnd_fields_cuda_fill_ghosts_H;
    add_ghosts_J          = psc_bnd_fields_cuda_add_ghosts_J;
  }
} psc_bnd_fields_cuda_ops;
