
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "fields_item.hxx"

#include <mrc_io.h>

#include <stdlib.h>

// FIXME: checkpointing won't properly restore state
// FIXME: if the subclass creates objects, it'd be cleaner to have them
// be part of the subclass

struct MarderCuda
{
  // ----------------------------------------------------------------------
  // fld_create
  //
  // FIXME, should be consolidated with psc_checks.c, and probably other places

  static struct psc_mfields *
  fld_create(struct psc *psc, const char *name)
  {
    auto mflds = PscMfields<MfieldsCuda>::create(psc_comm(psc), psc->grid(), 1, psc->ibn);
    psc_mfields_set_comp_name(mflds.mflds(), 0, name);
    return mflds.mflds();
  }

  // ----------------------------------------------------------------------
  // setup

  static void setup(struct psc_marder *marder)
  {
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
  // destroy

  static void destroy(struct psc_marder *marder)
  {
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
			  struct psc_mfields *_mflds_base, struct psc_mfields *_mf_base)
  {
    auto mflds_base = PscMfieldsBase{_mflds_base};
    auto mf_base = PscMfieldsBase{_mf_base};
    assert(mflds_base->grid().isInvar(0));

    const Grid_t& grid = ppsc->grid();
    // FIXME: how to choose diffusion parameter properly?
    float dx[3];
    for (int d = 0; d < 3; d++) {
      dx[d] = grid.domain.dx[d];
    }
    float inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!grid.isInvar(d)) {
	inv_sum += 1. / sqr(grid.domain.dx[d]);
      }
    }
    float diffusion_max = 1. / 2. / (.5 * ppsc->dt) / inv_sum;
    float diffusion     = diffusion_max * marder->diffusion;
    
    float fac[3];
    fac[0] = 0.f;
    fac[1] = .5 * ppsc->dt * diffusion / dx[1];
    fac[2] = .5 * ppsc->dt * diffusion / dx[2];

    auto& mflds = mflds_base->get_as<MfieldsCuda>(EX, EX + 3);
    auto& mf = mf_base->get_as<MfieldsCuda>(0, 1);
    cuda_mfields *cmflds = mflds.cmflds;
    cuda_mfields *cmf = mf.cmflds;

    // OPT, do all patches in one kernel
    for (int p = 0; p < mf.n_patches(); p++) {
      int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
      int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
      for (int d = 0; d < 3; d++) {
	if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
	    psc_at_boundary_lo(ppsc, p, d)) {
	  l_cc[d] = -1;
	  l_nc[d] = -1;
	}
	if (grid.bc.fld_hi[d] == BND_FLD_CONDUCTING_WALL &&
	    psc_at_boundary_hi(ppsc, p, d)) {
	  r_cc[d] = -1;
	  r_nc[d] = 0;
	}
      }
    
      const int *ldims = ppsc->grid().ldims;
    
      int ly[3] = { l_nc[0], l_cc[1], l_nc[2] };
      int ry[3] = { r_nc[0] + ldims[0], r_cc[1] + ldims[1], r_nc[2] + ldims[2] };
    
      int lz[3] = { l_nc[0], l_nc[1], l_cc[2] };
      int rz[3] = { r_nc[0] + ldims[0], r_nc[1] + ldims[1], r_cc[2] + ldims[2] };
    
      cuda_marder_correct_yz(cmflds, cmf, p, fac, ly, ry, lz, rz);
    }

    mflds_base->put_as(mflds, EX, EX + 3);
    mf_base->put_as(mf, 0, 0);
  }

  static void marder_calc_aid_fields(struct psc_marder *marder, 
				     PscMfieldsBase mflds_base, PscMparticlesBase mprts_base)
  {
    PscFieldsItemBase item_div_e(marder->item_div_e);
    PscFieldsItemBase item_rho(marder->item_rho);
    item_div_e(mflds_base, mprts_base); // FIXME, should accept NULL for particles
  
    if (marder->dump) {
      static int cnt;
      mrc_io_open(marder->io, "w", cnt, cnt);//ppsc->timestep, ppsc->timestep * ppsc->dt);
      cnt++;
      psc_mfields_write_as_mrc_fld(item_rho->mres().mflds(), marder->io);
      psc_mfields_write_as_mrc_fld(item_div_e->mres().mflds(), marder->io);
      mrc_io_close(marder->io);
    }

    item_div_e->mres()->axpy_comp(0, -1., *item_rho->mres().sub(), 0);
    // FIXME, why is this necessary?
    auto bnd = PscBndBase(marder->bnd);
    bnd.fill_ghosts(item_div_e->mres(), 0, 1);
  }

  static void run(struct psc_marder *marder, PscMfieldsBase mflds_base,
		  PscMparticlesBase mprts_base)
  {
    PscFieldsItemBase item_rho(marder->item_rho);
    PscFieldsItemBase item_div_e(marder->item_div_e);
    item_rho(mflds_base, mprts_base);
  
    // need to fill ghost cells first (should be unnecessary with only variant 1) FIXME
    auto bnd = PscBndBase(ppsc->bnd);
    bnd.fill_ghosts(mflds_base, EX, EX+3);
  
    for (int i = 0; i < marder->loop; i++) {
      marder_calc_aid_fields(marder, mflds_base, mprts_base);
      psc_marder_cuda_correct(marder, mflds_base.mflds(), item_div_e->mres().mflds());
      auto bnd = PscBndBase(ppsc->bnd);
      bnd.fill_ghosts(mflds_base, EX, EX+3);
    }
  }
};

// ======================================================================
// psc_marder: subclass "cuda"

struct psc_marder_ops_cuda : psc_marder_ops {
  psc_marder_ops_cuda() {
    name                  = "cuda";
    size                  = sizeof(MarderCuda);
    setup                 = MarderCuda::setup;
    destroy               = MarderCuda::destroy;
    run                   = MarderCuda::run;
  }
} psc_marder_cuda_ops;

