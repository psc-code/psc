
#pragma once

#include "fields_item_dive_cuda.hxx"
#include "fields_item_moments_1st_cuda.hxx"

#include <mrc_io.h>

// FIXME: checkpointing won't properly restore state
// FIXME: if the subclass creates objects, it'd be cleaner to have them
// be part of the subclass

template<typename BS>
struct MarderCuda : MarderBase
{
  using MfieldsState = MfieldsStateCuda;
  using Mfields = MfieldsCuda;
  using Mparticles = MparticlesCuda<BS>;
  using real_t = MfieldsState::real_t;
  
  MarderCuda(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : grid_{grid},
      diffusion_{diffusion},
      loop_{loop},
      dump_{dump},
      item_rho_{grid},
      item_div_e_{grid}
  {
#if 0
    bnd_ = psc_bnd_create(grid.comm());
    psc_bnd_set_name(bnd_, "marder_bnd");
    psc_bnd_set_type(bnd_, "cuda");
    psc_bnd_set_psc(bnd_, ppsc);
    psc_bnd_setup(bnd_);
#endif
    MHERE;
    return;
    assert(0);

    // FIXME, output_fields should be taking care of their own psc_bnd?
    if (dump_) {
      io_ = mrc_io_create(grid.comm());
      mrc_io_set_type(io_, "xdmf_collective");
      mrc_io_set_name(io_, "mrc_io_marder");
      mrc_io_set_param_string(io_, "basename", "marder");
      mrc_io_set_from_options(io_);
      mrc_io_setup(io_);
    }
  }

  ~MarderCuda()
  {
    //psc_bnd_destroy(bnd_);
    if (dump_) {
      mrc_io_destroy(io_);
    }
  }
  
  void calc_aid_fields(MfieldsState& mflds, Mparticles& mprts)
  {
    item_div_e_(mprts.grid(), mflds, mprts); // FIXME, should accept NULL for particles
  
    if (dump_) {
      static int cnt;
      mrc_io_open(io_, "w", cnt, cnt);//ppsc->timestep, ppsc->timestep * ppsc->dt);
      cnt++;
      // psc_mfields_write_as_mrc_fld(item_rho_.mres().mflds(), io_);
      // psc_mfields_write_as_mrc_fld(item_div_e_.mres().mflds(), io_);
      mrc_io_close(io_);
    }

    assert(0);
    //item_div_e_.mres().axpy_comp(0, -1., item_rho_.mres(), 0);
    // FIXME, why is this necessary?
    //auto bnd = PscBndBase(bnd_);
    //bnd.fill_ghosts(item_div_e_.mres(), 0, 1);
  }

  // ----------------------------------------------------------------------
  // psc_marder_cuda_correct
  //
  // Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

  void correct(MfieldsState& mflds, Mfields& mf)
  {
    assert(mflds._grid().isInvar(0));

    const Grid_t& grid = mflds._grid();
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
    float diffusion_max = 1. / 2. / (.5 * grid.dt) / inv_sum;
    float diffusion     = diffusion_max * diffusion_;
    
    float fac[3];
    fac[0] = 0.f;
    fac[1] = .5 * grid.dt * diffusion / dx[1];
    fac[2] = .5 * grid.dt * diffusion / dx[2];

    cuda_mfields *cmflds = mflds.cmflds();
    cuda_mfields *cmf = mf.cmflds();

    // OPT, do all patches in one kernel
    for (int p = 0; p < mf.n_patches(); p++) {
      int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
      int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
      for (int d = 0; d < 3; d++) {
	if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL && grid.atBoundaryLo(p, d)) {
	  l_cc[d] = -1;
	  l_nc[d] = -1;
	}
	if (grid.bc.fld_hi[d] == BND_FLD_CONDUCTING_WALL && grid.atBoundaryHi(p, d)) {
	  r_cc[d] = -1;
	  r_nc[d] = 0;
	}
      }
    
      const int *ldims = grid.ldims;
    
      int ly[3] = { l_nc[0], l_cc[1], l_nc[2] };
      int ry[3] = { r_nc[0] + ldims[0], r_cc[1] + ldims[1], r_nc[2] + ldims[2] };
    
      int lz[3] = { l_nc[0], l_nc[1], l_cc[2] };
      int rz[3] = { r_nc[0] + ldims[0], r_nc[1] + ldims[1], r_cc[2] + ldims[2] };
    
      cuda_marder_correct_yz(cmflds, cmf, p, fac, ly, ry, lz, rz);
    }
  }

  void operator()(MfieldsStateCuda& mflds, MparticlesCuda<BS>& mprts)
  {
    MHERE;
  }
  
private:
  const Grid_t& grid_;
  real_t diffusion_; //< diffusion coefficient for Marder correction
  int loop_; //< execute this many relaxation steps in a loop
  bool dump_; //< dump div_E, rho

  FieldsItemFields<Item_dive_cuda> item_div_e_;
  Moment_rho_1st_nc_cuda<MparticlesCuda<BS144>, dim_yz> item_rho_; // FIXME, hardcoded dim_yz
  mrc_io *io_; //< for debug dumping
};

