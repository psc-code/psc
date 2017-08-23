
#include "psc_checks_private.h"

#include "psc_bnd.h"
#include "psc_output_fields_item.h"

#include <mrc_io.h>

// FIXME, duplicated

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

// ----------------------------------------------------------------------
// FIXME, should be consolidated?

static struct psc_mfields *
fld_create(struct psc *psc, int nr_fields)
{
  struct psc_mfields *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(fld, FIELDS_TYPE);
  psc_mfields_set_domain(fld, psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_set_param_int(fld, "nr_fields", nr_fields);
  psc_mfields_setup(fld);

  return fld;
}

// ----------------------------------------------------------------------
// calc_rho

static void
calc_rho(struct psc *psc, struct psc_mparticles *mprts, struct psc_mfields *rho)
{
  // FIXME, output_fields should be taking care of this?
  struct psc_bnd *bnd = psc_bnd_create(psc_comm(psc));
  psc_bnd_set_name(bnd, "psc_output_fields_bnd_calc_rho");
  psc_bnd_set_type(bnd, FIELDS_TYPE);
  psc_bnd_set_psc(bnd, psc);
  psc_bnd_setup(bnd);

  struct psc_output_fields_item *item = psc_output_fields_item_create(psc_comm(psc));
  // makes, e.g., "rho_1st_nc_single"
  psc_output_fields_item_set_type(item, "rho_" PSC_CHECKS_ORDER "_nc_" PARTICLE_TYPE);
  psc_output_fields_item_set_psc_bnd(item, bnd);
  psc_output_fields_item_setup(item);
  psc_output_fields_item_run(item, psc->flds, mprts, rho);
  psc_output_fields_item_destroy(item);

  psc_bnd_destroy(bnd);
}

// ----------------------------------------------------------------------
// psc_checks_sub_setup

static void
psc_checks_sub_setup(struct psc_checks *checks)
{
  psc_mfields_set_type(checks->rho_m, FIELDS_TYPE);
  psc_mfields_set_domain(checks->rho_m, ppsc->mrc_domain);
  psc_mfields_set_param_int3(checks->rho_m, "ibn", ppsc->ibn);
  psc_mfields_set_param_int(checks->rho_m, "nr_fields", 1);
  psc_mfields_list_add(&psc_mfields_base_list, &checks->rho_m);

  psc_mfields_set_type(checks->rho_p, FIELDS_TYPE);
  psc_mfields_set_domain(checks->rho_p, ppsc->mrc_domain);
  psc_mfields_set_param_int3(checks->rho_p, "ibn", ppsc->ibn);
  psc_mfields_set_param_int(checks->rho_p, "nr_fields", 1);
  psc_mfields_list_add(&psc_mfields_base_list, &checks->rho_p);

  psc_checks_setup_member_objs(checks);
}

// ======================================================================
// psc_checks: Charge Continuity 

// ----------------------------------------------------------------------
// psc_calc_div_j
//
// FIXME, make diag_item?

static void
do_calc_div_j(struct psc *psc, int p, struct psc_fields *flds, struct psc_fields *div_j)
{
  define_dxdydz(dx, dy, dz);
  fields_real_t h[3];
  for (int d = 0; d < 3; d++) {
    if (psc->domain.gdims[d] == 1) {
      h[d] = 0.;
    } else {
      h[d] = 1. / psc->patch[p].dx[d];
    }
  }

  psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
    F3(div_j,0, jx,jy,jz) =
      (F3(flds, JXI, jx,jy,jz) - F3(flds, JXI, jx-dx,jy,jz)) * h[0] +
      (F3(flds, JYI, jx,jy,jz) - F3(flds, JYI, jx,jy-dy,jz)) * h[1] +
      (F3(flds, JZI, jx,jy,jz) - F3(flds, JZI, jx,jy,jz-dz)) * h[2];
  } psc_foreach_3d_end;
}

static void
calc_div_j(struct psc *psc, struct psc_mfields *mflds_base, struct psc_mfields *div_j)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, JXI + 3);
  psc_foreach_patch(psc, p) {
    do_calc_div_j(psc, p, psc_mfields_get_patch(mflds, p),
		  psc_mfields_get_patch(div_j, p));
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_checks_continuity

static void
psc_checks_continuity(struct psc_checks *checks, struct psc *psc,
		      struct psc_mfields *rho_m, struct psc_mfields *rho_p)
{
  struct psc_mfields *div_j = fld_create(psc, 1);
  psc_mfields_set_name(div_j, "div_j");
  psc_mfields_set_comp_name(div_j, 0, "div_j");
  struct psc_mfields *d_rho = fld_create(psc, 1);
  psc_mfields_set_name(d_rho, "d_rho");
  psc_mfields_set_comp_name(d_rho, 0, "d_rho");

  psc_mfields_axpy(d_rho,  1., rho_p);
  psc_mfields_axpy(d_rho, -1., rho_m);

  calc_div_j(psc, psc->flds, div_j);
  psc_mfields_scale(div_j, psc->dt);

  double eps = checks->continuity_threshold;
  double max_err = 0.;
  psc_foreach_patch(psc, p) {
    struct psc_fields *p_d_rho = psc_mfields_get_patch(d_rho, p);
    struct psc_fields *p_div_j = psc_mfields_get_patch(div_j, p);
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      double d_rho = F3(p_d_rho,0, jx,jy,jz);
      double div_j = F3(p_div_j,0, jx,jy,jz);
      max_err = fmax(max_err, fabs(d_rho + div_j));
      if (fabs(d_rho + div_j) > eps) {
	mprintf("(%d,%d,%d): %g -- %g diff %g\n", jx, jy, jz,
		d_rho, -div_j, d_rho + div_j);
      }
    } psc_foreach_3d_end;
  }

  // find global max
  double tmp = max_err;
  MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, psc_checks_comm(checks));

  if (checks->continuity_verbose || max_err >= eps) {
    mpi_printf(psc_checks_comm(checks),
	       "continuity: max_err = %g (thres %g)\n", max_err, eps);
  }

  if (checks->continuity_dump_always || max_err >= eps) {
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
// psc_checks_sub_continuity_before_particle_push

static void
psc_checks_sub_continuity_before_particle_push(struct psc_checks *checks, struct psc *psc)
{
  if (checks->continuity_every_step < 0 ||
      psc->timestep % checks->continuity_every_step != 0) {
    return;
  }

  calc_rho(psc, psc->particles, checks->rho_m);
}

// ----------------------------------------------------------------------
// psc_checks_sub_continuity_after_particle_push

static void
psc_checks_sub_continuity_after_particle_push(struct psc_checks *checks, struct psc *psc)
{
  if (checks->continuity_every_step < 0 ||
      psc->timestep % checks->continuity_every_step != 0) {
    return;
  }

  calc_rho(psc, psc->particles, checks->rho_p);
  psc_checks_continuity(checks, psc, checks->rho_m, checks->rho_p);
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
  psc_bnd_set_type(bnd, FIELDS_TYPE);
  psc_bnd_set_psc(bnd, psc);
  psc_bnd_setup(bnd);

  struct psc_output_fields_item *item = psc_output_fields_item_create(psc_comm(psc));
  psc_output_fields_item_set_type(item, "dive_" FIELDS_TYPE);
  psc_output_fields_item_set_psc_bnd(item, bnd);
  psc_output_fields_item_setup(item);
  psc_output_fields_item_run(item, mflds, psc->particles, dive); // FIXME, should accept NULL for mprts
  psc_output_fields_item_destroy(item);

  psc_bnd_destroy(bnd);
}

// ----------------------------------------------------------------------
// psc_checks_sub_gauss

static void
psc_checks_sub_gauss(struct psc_checks *checks, struct psc *psc)
{
  if (checks->gauss_every_step < 0 ||
      psc->timestep % checks->gauss_every_step != 0) {
    return;
  }

  struct psc_mfields *dive = fld_create(psc, 1);
  psc_mfields_set_name(dive, "div_E");
  psc_mfields_set_comp_name(dive, 0, "div_E");
  struct psc_mfields *rho = fld_create(psc, 1);
  psc_mfields_set_name(rho, "rho");
  psc_mfields_set_comp_name(rho, 0, "rho");

  calc_rho(psc, psc->particles, rho);
  calc_dive(psc, psc->flds, dive);

  double eps = checks->gauss_threshold;
  double max_err = 0.;
  psc_foreach_patch(psc, p) {
    struct psc_fields *p_rho = psc_mfields_get_patch(rho, p);
    struct psc_fields *p_dive = psc_mfields_get_patch(dive, p);

    int l[3] = {0, 0, 0}, r[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
      if (ppsc->domain.bnd_fld_lo[d] == BND_FLD_CONDUCTING_WALL && ppsc->patch[p].off[d] == 0) {
	l[d] = 1;
      }
    }

    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      if (jy < l[1] || jz < l[2] ||
	  jy >= psc->patch[p].ldims[1] - r[1] ||
	  jz >= psc->patch[p].ldims[2] - r[2]) {
	continue;
      }
      double v_rho = F3(p_rho,0, jx,jy,jz);
      double v_dive = F3(p_dive,0, jx,jy,jz);
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
  MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, psc_checks_comm(checks));

  if (checks->gauss_verbose || max_err >= eps) {
    mpi_printf(psc_checks_comm(checks), "gauss: max_err = %g (thres %g)\n", max_err, eps);
  }

  if (checks->gauss_dump_always || max_err >= eps) {
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

// ----------------------------------------------------------------------
// psc_checks_sub_ops

struct psc_checks_ops psc_checks_sub_ops = {
  // makes, e.g., "1st_single"
  .name                            = PSC_CHECKS_ORDER "_" PARTICLE_TYPE,
  .setup                           = psc_checks_sub_setup,
  .continuity_before_particle_push = psc_checks_sub_continuity_before_particle_push,
  .continuity_after_particle_push  = psc_checks_sub_continuity_after_particle_push,
  .gauss                           = psc_checks_sub_gauss,
};
