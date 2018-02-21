
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include <mrc_io.h>

#include <stdlib.h>

// FIXME: checkpointing won't properly restore state
// FIXME: if the subclass creates objects, it'd be cleaner to have them
// be part of the subclass

// ----------------------------------------------------------------------
// fld_create
//
// FIXME, should be consolidated with psc_checks.c, and probably other places

static struct psc_mfields *
fld_create(struct psc *psc, const char *name)
{
  struct psc_mfields *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(fld, "cuda");
  psc_mfields_set_param_obj(fld, "domain", psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_set_param_int(fld, "nr_fields", 1);
  psc_mfields_setup(fld);
  psc_mfields_set_comp_name(fld, 0, name);

  return fld;
}

// ----------------------------------------------------------------------
// psc_marder_cuda_setup

static void
psc_marder_cuda_setup(struct psc_marder *marder)
{
  marder->div_e = fld_create(ppsc, "div_E");
  marder->rho = fld_create(ppsc, "rho");

  marder->bnd = psc_bnd_create(psc_marder_comm(marder));
  psc_bnd_set_name(marder->bnd, "marder_bnd");
  psc_bnd_set_type(marder->bnd, "cuda");
  psc_bnd_set_psc(marder->bnd, ppsc);
  psc_bnd_setup(marder->bnd);

  // FIXME, output_fields should be taking care of their own psc_bnd?
  marder->item_div_e = psc_output_fields_item_create(psc_comm(ppsc));
  psc_output_fields_item_set_type(marder->item_div_e, "dive_cuda");
  psc_output_fields_item_set_psc_bnd(marder->item_div_e, marder->bnd);
  psc_output_fields_item_setup(marder->item_div_e);

  marder->item_rho = psc_output_fields_item_create(psc_comm(ppsc));
  psc_output_fields_item_set_type(marder->item_rho, "rho_1st_nc_cuda");
  psc_output_fields_item_set_psc_bnd(marder->item_rho, marder->bnd);
  psc_output_fields_item_setup(marder->item_rho);

  if (marder->dump) {
    struct mrc_io *io = mrc_io_create(psc_comm(ppsc));
    mrc_io_set_type(io, "xdmf_collective");
    mrc_io_set_name(io, "mrc_io_marder");
    mrc_io_set_param_string(io, "basename", "marder");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);

    marder->io = io;
  }
}

// ----------------------------------------------------------------------
// psc_marder_cuda_destroy

static void
psc_marder_cuda_destroy(struct psc_marder *marder)
{
  psc_mfields_destroy(marder->div_e);
  psc_mfields_destroy(marder->rho);

  psc_output_fields_item_destroy(marder->item_div_e);
  psc_output_fields_item_destroy(marder->item_rho);

  psc_bnd_destroy(marder->bnd);

  if (marder->dump) {
    mrc_io_destroy(marder->io);
  }
}

// ----------------------------------------------------------------------
// psc_marder_cuda_correct
//
// Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

static void
psc_marder_cuda_correct(struct psc_marder *marder,
			struct psc_mfields *mflds_base, struct psc_mfields *mf_base)
{
  assert(ppsc->domain.gdims[0] == 1);

  const Grid_t& grid = ppsc->grid();
  // FIXME: how to choose diffusion parameter properly?
  float dx[3];
  for (int d = 0; d < 3; d++) {
    dx[d] = grid.dx[d];
  }
  float inv_sum = 0.;
  int nr_levels;
  mrc_domain_get_nr_levels(ppsc->mrc_domain, &nr_levels);
  for (int d = 0; d < 3; d++) {
      if (ppsc->domain.gdims[d] > 1) {
	inv_sum += 1. / sqr(grid.dx[d] / (1 << (nr_levels - 1)));
      }
    }
  float diffusion_max = 1. / 2. / (.5 * ppsc->dt) / inv_sum;
  float diffusion     = diffusion_max * marder->diffusion;
    
  float fac[3];
  fac[0] = 0.f;
  fac[1] = .5 * ppsc->dt * diffusion / dx[1];
  fac[2] = .5 * ppsc->dt * diffusion / dx[2];

  mfields_cuda_t mflds = mflds_base->get_as<mfields_cuda_t>(EX, EX + 3);
  mfields_cuda_t mf = mf_base->get_as<mfields_cuda_t>(0, 1);
  cuda_mfields *cmflds = mflds->cmflds;
  cuda_mfields *cmf = mf->cmflds;

  // OPT, do all patches in one kernel
  for (int p = 0; p < mf_base->nr_patches; p++) {
    int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
    int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
      if (ppsc->domain.bnd_fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
	  psc_at_boundary_lo(ppsc, p, d)) {
	l_cc[d] = -1;
	l_nc[d] = -1;
      }
      if (ppsc->domain.bnd_fld_hi[d] == BND_FLD_CONDUCTING_WALL &&
	  psc_at_boundary_hi(ppsc, p, d)) {
	r_cc[d] = -1;
	r_nc[d] = 0;
      }
    }
    
    int *ldims = ppsc->patch[p].ldims;
    
    int ly[3] = { l_nc[0], l_cc[1], l_nc[2] };
    int ry[3] = { r_nc[0] + ldims[0], r_cc[1] + ldims[1], r_nc[2] + ldims[2] };
    
    int lz[3] = { l_nc[0], l_nc[1], l_cc[2] };
    int rz[3] = { r_nc[0] + ldims[0], r_nc[1] + ldims[1], r_cc[2] + ldims[2] };
    
    cuda_marder_correct_yz(cmflds, cmf, p, fac, ly, ry, lz, rz);
  }

  mflds.put_as(mflds_base, EX, EX + 3);
  mf.put_as(mf_base, 0, 0);
}

// ======================================================================
// psc_marder: subclass "cuda"

struct psc_marder_ops_cuda : psc_marder_ops {
  psc_marder_ops_cuda() {
    name                  = "cuda";
    setup                 = psc_marder_cuda_setup;
    destroy               = psc_marder_cuda_destroy;
    correct               = psc_marder_cuda_correct;
  }
} psc_marder_cuda_ops;

