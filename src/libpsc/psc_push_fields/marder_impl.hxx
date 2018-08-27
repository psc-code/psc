
#include "marder.hxx"
#include "psc_output_fields_item.h"
#include "fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/psc_output_fields_item_moments_1st_nc.cxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"

#include <mrc_io.h>

template<typename _Mparticles, typename _MfieldsState, typename _Mfields>
struct Marder_ : MarderBase
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using fields_t = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;
  using Fields = Fields3d<fields_t>;
  using Moment_t = Moment_rho_1st_nc<Mparticles, Mfields>;

  Marder_(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : grid_{grid},
      diffusion_{diffusion},
      loop_{loop},
      dump_{dump},
      bnd_{grid, grid.ibn},
      item_rho_{grid, grid.comm()},
      item_dive_{grid, grid.comm()}
  {
    if (dump_) {
      io_ = mrc_io_create(grid.comm());
      mrc_io_set_type(io_, "xdmf_collective");
      mrc_io_set_name(io_, "mrc_io_marder");
      mrc_io_set_param_string(io_, "basename", "marder");
      mrc_io_set_from_options(io_);
      mrc_io_setup(io_);
    }
  }

  ~Marder_()
  {
    if (dump_) {
      mrc_io_destroy(io_);
    }
  }
  
  // FIXME: checkpointing won't properly restore state
  // FIXME: if the subclass creates objects, it'd be cleaner to have them
  // be part of the subclass

  // ----------------------------------------------------------------------
  // calc_aid_fields

  void calc_aid_fields(MfieldsState& mflds)
  {
    item_dive_(mflds);
	       
    if (dump_) {
      static int cnt;
      mrc_io_open(io_, "w", cnt, cnt);//ppsc->timestep, ppsc->timestep * ppsc->dt);
      cnt++;
      item_rho_.result().write_as_mrc_fld(io_, "rho", {"rho"});
      item_dive_.result().write_as_mrc_fld(io_, "dive", {"dive"});
      mrc_io_close(io_);
    }
    
    item_dive_.result().axpy_comp(0, -1., item_rho_.result(), 0);
    // FIXME, why is this necessary?
    bnd_.fill_ghosts(item_dive_.result(), 0, 1);
  }

  // ----------------------------------------------------------------------
  // correct_patch
  //
  // Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

#define define_dxdydz(dx, dy, dz)				\
  const auto& grid = flds.grid();				\
  int dx _mrc_unused = (grid.isInvar(0)) ? 0 : 1;		\
  int dy _mrc_unused = (grid.isInvar(1)) ? 0 : 1;		\
  int dz _mrc_unused = (grid.isInvar(2)) ? 0 : 1

#define psc_foreach_3d_more(psc, p, ix, iy, iz, l, r) {	\
  int __ilo[3] = { -l[0], -l[1], -l[2] };		\
  int __ihi[3] = { grid.ldims[0] + r[0],		\
		   grid.ldims[1] + r[1],		\
		   grid.ldims[2] + r[2] };		\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {	\
  for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {	\
  for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define psc_foreach_3d_more_end			\
  } } }

  void correct_patch(fields_t flds, fields_t f, int p, real_t& max_err)
  {
    Fields F(flds), FF(f);
    define_dxdydz(dx, dy, dz);

    // FIXME: how to choose diffusion parameter properly?
    //double deltax = ppsc->patch[p].dx[0];
    double deltay = grid.domain.dx[1]; // FIXME double/float
    double deltaz = grid.domain.dx[2];
    double inv_sum = 0.;
    int nr_levels;
    for (int d = 0; d < 3; d++) {
      if (!grid.isInvar(d)) {
	inv_sum += 1. / sqr(grid.domain.dx[d]);
      }
    }
    double diffusion_max = 1. / 2. / (.5 * grid.dt) / inv_sum;
    double diffusion     = diffusion_max * diffusion_;

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

#if 0
    psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
      // FIXME: F3 correct?
      F(EX, ix,iy,iz) += 
	(FF(DIVE_MARDER, ix+dx,iy,iz) - FF(DIVE_MARDER, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltax;
      F(EY, ix,iy,iz) += 
	(FF(DIVE_MARDER, ix,iy+dy,iz) - FF(DIVE_MARDER, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltay;
      F(EZ, ix,iy,iz) += 
	(FF(DIVE_MARDER, ix,iy,iz+dz) - FF(DIVE_MARDER, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltaz;
    } psc_foreach_3d_more_end;
#endif

    assert(grid.isInvar(0));

    {
      int l[3] = { l_nc[0], l_cc[1], l_nc[2] };
      int r[3] = { r_nc[0], r_cc[1], r_nc[2] };
      psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
	max_err = std::max(max_err, std::abs(FF(0, ix,iy,iz)));
	F(EY, ix,iy,iz) += 
	  (FF(0, ix,iy+dy,iz) - FF(0, ix,iy,iz))
	  * .5 *grid.dt * diffusion / deltay;
      } psc_foreach_3d_more_end;
    }

    {
      int l[3] = { l_nc[0], l_nc[1], l_cc[2] };
      int r[3] = { r_nc[0], r_nc[1], r_cc[2] };
      psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
	F(EZ, ix,iy,iz) += 
	  (FF(0, ix,iy,iz+dz) - FF(0, ix,iy,iz))
	  * .5 * grid.dt * diffusion / deltaz;
      } psc_foreach_3d_more_end;
    }
  }

#undef psc_foreach_3d_more
#undef psc_foreach_3d_more_end

  // ----------------------------------------------------------------------
  // correct

  void correct(MfieldsState& mf)
  {
    auto& mf_div_e = item_dive_.result();

    real_t max_err = 0.;
    for (int p = 0; p < mf_div_e.n_patches(); p++) {
      correct_patch(mf[p], mf_div_e[p], p, max_err);
    }

    MPI_Allreduce(MPI_IN_PLACE, &max_err, 1, Mfields_traits<Mfields>::mpi_dtype(), MPI_MAX, grid_.comm());
    mpi_printf(grid_.comm(), "marder: err %g\n", max_err);
  }

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    item_rho_.run(mprts);

    // need to fill ghost cells first (should be unnecessary with only variant 1) FIXME
    assert(0);
#if 0
    bnd_.fill_ghosts(mflds, EX, EX+3);

    for (int i = 0; i < loop_; i++) {
      calc_aid_fields(mflds);
      correct(mflds);
      bnd_.fill_ghosts(mflds, EX, EX+3);
    }
#endif
  }
  
  // ----------------------------------------------------------------------
  // run
  
  void run(MfieldsStateBase& mflds_base, MparticlesBase& mprts_base) override
  {
    auto& mflds = mflds_base.get_as<MfieldsState>(EX, EX + 3);
    auto& mprts = mprts_base.get_as<Mparticles>();

    (*this)(mflds, mprts);

    mflds_base.put_as(mflds, EX, EX + 3);
    mprts_base.put_as(mprts, MP_DONT_COPY);
  }

private:
  real_t diffusion_; //< diffusion coefficient for Marder correction
  int loop_; //< execute this many relaxation steps in a loop
  bool dump_; //< dump div_E, rho

  const Grid_t& grid_;
  Bnd_<Mfields> bnd_;
  ItemMomentLoopPatches<Moment_t> item_rho_;
  FieldsItemFields<ItemLoopPatches<Item_dive<MfieldsState, Mfields>>> item_dive_;
  mrc_io *io_; //< for debug dumping
};

