
#include "psc_checks_private.h"

#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "fields.hxx"
#include "fields_item.hxx"
#include "checks.hxx"

#include <mrc_io.h>

// FIXME, duplicated

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->grid().isInvar(0)) ? 0 : 1;      \
  int dy _mrc_unused = (ppsc->grid().isInvar(1)) ? 0 : 1;      \
  int dz _mrc_unused = (ppsc->grid().isInvar(2)) ? 0 : 1

namespace {

template<typename MP, typename MF>
struct Checks_ : ChecksParams, ChecksBase
{
  using Mfields = MF;
  using fields_t = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;
  using Fields = Fields3d<fields_t>;

  Checks_(MPI_Comm comm, const ChecksParams& params)
    : ChecksParams{params},
      comm_{comm}
  {
    rho_m = fld_create(ppsc, 1);
    rho_p = fld_create(ppsc, 1);
  }
  
  // ----------------------------------------------------------------------
  // FIXME, should be consolidated?

  static struct psc_mfields *
  fld_create(struct psc *psc, int nr_fields)
  {
    return PscMfields<Mfields>::create(psc_comm(psc), psc->grid(), nr_fields, psc->ibn).mflds();
  }

  // ----------------------------------------------------------------------
  // calc_rho

  static void
  calc_rho(struct psc *psc, struct psc_mparticles *mprts, struct psc_mfields *rho)
  {
    // FIXME, output_fields should be taking care of this?
    struct psc_bnd *bnd = psc_bnd_create(psc_comm(psc));
    psc_bnd_set_name(bnd, "psc_output_fields_bnd_calc_rho");
    psc_bnd_set_type(bnd, Mfields_traits<Mfields>::name);
    psc_bnd_set_psc(bnd, psc);
    psc_bnd_setup(bnd);

    struct psc_output_fields_item *item = psc_output_fields_item_create(psc_comm(psc));
    // makes, e.g., "rho_1st_nc_single"
    psc_output_fields_item_set_type(item, "rho_" PSC_CHECKS_ORDER "_nc_" PARTICLE_TYPE);
    psc_output_fields_item_set_psc_bnd(item, bnd);
    psc_output_fields_item_setup(item);
    PscMparticlesBase mp(mprts);
    PscFieldsItemBase _item(item);
    _item(psc->flds, mp, rho);
    psc_output_fields_item_destroy(item);

    psc_bnd_destroy(bnd);
  }

  // ======================================================================
  // psc_checks: Charge Continuity 

  // ----------------------------------------------------------------------
  // psc_calc_div_j
  //
  // FIXME, make diag_item?

  static void
  do_calc_div_j(struct psc *psc, int p, fields_t flds, fields_t div_j)
  {
    Fields F(flds), Div_J(div_j);
    define_dxdydz(dx, dy, dz);
    real_t h[3];
    for (int d = 0; d < 3; d++) {
      if (psc->grid().isInvar(d)) {
	h[d] = 0.;
      } else {
	h[d] = 1. / psc->grid().domain.dx[d];
      }
    }

    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      Div_J(0, jx,jy,jz) =
	(F(JXI, jx,jy,jz) - F(JXI, jx-dx,jy,jz)) * h[0] +
	(F(JYI, jx,jy,jz) - F(JYI, jx,jy-dy,jz)) * h[1] +
	(F(JZI, jx,jy,jz) - F(JZI, jx,jy,jz-dz)) * h[2];
    } psc_foreach_3d_end;
  }

  static void
  calc_div_j(struct psc *psc, struct psc_mfields *_mflds_base, struct psc_mfields *div_j)
  {
    auto mflds_base = PscMfieldsBase{_mflds_base};
    auto mf = mflds_base.get_as<PscMfields<Mfields>>(JXI, JXI + 3);
    PscMfields<Mfields> mf_div_j(div_j);
    psc_foreach_patch(psc, p) {
      do_calc_div_j(psc, p, mf[p], mf_div_j[p]);
    }
    mf.put_as(mflds_base, 0, 0);
  }

  // ----------------------------------------------------------------------
  // continuity

  void continuity(psc *psc, psc_mfields *rho_m, psc_mfields *rho_p)
  {
    struct psc_mfields *div_j = fld_create(psc, 1);
    psc_mfields_set_name(div_j, "div_j");
    psc_mfields_set_comp_name(div_j, 0, "div_j");
    struct psc_mfields *d_rho = fld_create(psc, 1);
    psc_mfields_set_name(d_rho, "d_rho");
    psc_mfields_set_comp_name(d_rho, 0, "d_rho");
    PscMfields<Mfields> mf_div_j(div_j), mf_d_rho(d_rho);

    mf_d_rho->axpy( 1., *PscMfields<Mfields>(rho_p).sub());
    mf_d_rho->axpy(-1., *PscMfields<Mfields>(rho_m).sub());

    calc_div_j(psc, psc->flds, div_j);
    mf_div_j->scale(psc->dt);

    double eps = continuity_threshold;
    double max_err = 0.;
    psc_foreach_patch(psc, p) {
      Fields D_rho(mf_d_rho[p]);
      Fields Div_J(mf_div_j[p]);
      psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
	double d_rho = D_rho(0, jx,jy,jz);
	double div_j = Div_J(0, jx,jy,jz);
	max_err = fmax(max_err, fabs(d_rho + div_j));
	if (fabs(d_rho + div_j) > eps) {
	  mprintf("(%d,%d,%d): %g -- %g diff %g\n", jx, jy, jz,
		  d_rho, -div_j, d_rho + div_j);
	}
      } psc_foreach_3d_end;
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (continuity_verbose || max_err >= eps) {
      mpi_printf(comm_, "continuity: max_err = %g (thres %g)\n", max_err, eps);
    }

    if (continuity_dump_always || max_err >= eps) {
      static struct mrc_io *io;
      if (!io) {
	io = mrc_io_create(psc_comm(psc));
	mrc_io_set_name(io, "mrc_io_continuity");
	mrc_io_set_param_string(io, "basename", "continuity");
	mrc_io_set_from_options(io);
	mrc_io_setup(io);
	mrc_io_view(io);
      }
      mrc_io_open(io, "w", psc->timestep, psc->timestep * psc->dt);
      psc_mfields_write_as_mrc_fld(div_j, io);
      psc_mfields_write_as_mrc_fld(d_rho, io);
      mrc_io_close(io);
    }

    assert(max_err < eps);

    psc_mfields_destroy(div_j);
    psc_mfields_destroy(d_rho);
  }

  // ----------------------------------------------------------------------
  // continuity_before_particle_push

  void continuity_before_particle_push(psc *psc) override
  {
    if (continuity_every_step < 0 || psc->timestep % continuity_every_step != 0) {
      return;
    }

    calc_rho(psc, psc->particles, rho_m);
  }

  // ----------------------------------------------------------------------
  // continuity_after_particle_push

  void continuity_after_particle_push(psc *psc) override
  {
    if (continuity_every_step < 0 || psc->timestep % continuity_every_step != 0) {
      return;
    }

    calc_rho(psc, psc->particles, rho_p);
    continuity(psc, rho_m, rho_p);
  }

  // ======================================================================
  // psc_checks: Gauss's Law

  // ----------------------------------------------------------------------
  // calc_dive

  static void
  calc_dive(struct psc *psc, struct psc_mfields *mflds, struct psc_mfields *dive)
  {
    // FIXME, output_fields should be taking care of this?
    struct psc_bnd *bnd = psc_bnd_create(psc_comm(psc));
    psc_bnd_set_name(bnd, "psc_output_fields_bnd_calc_dive");
    psc_bnd_set_type(bnd, Mfields_traits<Mfields>::name);
    psc_bnd_set_psc(bnd, psc);
    psc_bnd_setup(bnd);

    struct psc_output_fields_item *item = psc_output_fields_item_create(psc_comm(psc));
    auto s = std::string("dive_") + Mfields_traits<Mfields>::name;
    psc_output_fields_item_set_type(item, s.c_str());
    psc_output_fields_item_set_psc_bnd(item, bnd);
    psc_output_fields_item_setup(item);
    PscMparticlesBase mprts(psc->particles);
    PscFieldsItemBase _item(item);
    _item(mflds, mprts, dive); // FIXME, should accept NULL for mprts
    psc_output_fields_item_destroy(item);

    psc_bnd_destroy(bnd);
  }

  // ----------------------------------------------------------------------
  // gauss

  void gauss(psc* psc) override
  {
    if (gauss_every_step < 0 ||	psc->timestep % gauss_every_step != 0) {
      return;
    }

    struct psc_mfields *dive = fld_create(psc, 1);
    psc_mfields_set_name(dive, "div_E");
    psc_mfields_set_comp_name(dive, 0, "div_E");
    struct psc_mfields *rho = fld_create(psc, 1);
    psc_mfields_set_name(rho, "rho");
    psc_mfields_set_comp_name(rho, 0, "rho");
    PscMfields<Mfields> mf_dive(dive), mf_rho(rho);
    const auto& grid = psc->grid();
  
    calc_rho(psc, psc->particles, rho);
    calc_dive(psc, psc->flds, dive);

    double eps = gauss_threshold;
    double max_err = 0.;
    psc_foreach_patch(psc, p) {
      Fields Rho(mf_rho[p]), DivE(mf_dive[p]);

      int l[3] = {0, 0, 0}, r[3] = {0, 0, 0};
      for (int d = 0; d < 3; d++) {
	if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
	    psc_at_boundary_lo(ppsc, p, d)) {
	  l[d] = 1;
	}
      }

      psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
	if (jy < l[1] || jz < l[2] ||
	    jy >= psc->grid().ldims[1] - r[1] ||
	    jz >= psc->grid().ldims[2] - r[2]) {
	  continue;
	}
	double v_rho = Rho(0, jx,jy,jz);
	double v_dive = DivE(0, jx,jy,jz);
	max_err = fmax(max_err, fabs(v_dive - v_rho));
#if 0
	if (fabs(v_dive - v_rho) > eps) {
	  printf("(%d,%d,%d): %g -- %g diff %g\n", jx, jy, jz,
		 v_dive, v_rho, v_dive - v_rho);
	}
#endif
      } psc_foreach_3d_end;
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (gauss_verbose || max_err >= eps) {
      mpi_printf(comm_, "gauss: max_err = %g (thres %g)\n", max_err, eps);
    }

    if (gauss_dump_always || max_err >= eps) {
      static struct mrc_io *io;
      if (!io) {
	io = mrc_io_create(psc_comm(psc));
	mrc_io_set_name(io, "mrc_io_gauss");
	mrc_io_set_param_string(io, "basename", "gauss");
	mrc_io_set_from_options(io);
	mrc_io_setup(io);
	mrc_io_view(io);
      }
      mrc_io_open(io, "w", psc->timestep, psc->timestep * psc->dt);
      psc_mfields_write_as_mrc_fld(rho, io);
      psc_mfields_write_as_mrc_fld(dive, io);
      mrc_io_close(io);
    }

    assert(max_err < eps);

    psc_mfields_destroy(rho);
    psc_mfields_destroy(dive);
  }

  // state
  MPI_Comm comm_;
  psc_mfields *rho_m, *rho_p;
};

}

// ----------------------------------------------------------------------
// psc_checks_sub_read
//
// FIXME, this function exists to avoid a setup called twice error, but it's just a workaround

static void
psc_checks_sub_read(struct psc_checks *checks, struct mrc_io *io)
{
  psc_checks_read_super(checks, io);

  psc_checks_read_member_objs(checks, io);
}

