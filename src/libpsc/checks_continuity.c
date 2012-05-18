
#include "psc_testing.h"
#include "psc_bnd.h"
#include "psc_particles_c.h"
#include "psc_output_fields_item.h"

#include <math.h>

bool opt_checks_verbose = false;

typedef double creal;

static inline creal
creal_sqrt(creal x)
{
  return sqrt(x);
}

// ----------------------------------------------------------------------
// FIXME, should be consolidated?

static mfields_c_t *
fld_create(struct psc *psc, int nr_fields)
{
  mfields_c_t *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(fld, "c");
  psc_mfields_set_domain(fld, psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_set_param_int(fld, "nr_fields", nr_fields);
  psc_mfields_setup(fld);

  return fld;
}

// ----------------------------------------------------------------------
// psc_shift_particle_positions
//
// shift particle positions back / forth in time according to current
// velocity

static void
do_shift_particle_positions(particles_c_t *pp, double dt)
{
  for (int n = 0; n < pp->n_part; n++) {
    particle_c_t *part = particles_c_get_one(pp, n);
      
    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;
    part->xi += vxi * dt;
    part->yi += vyi * dt;
    part->zi += vzi * dt;
  }
}

static void
psc_shift_particle_positions(struct psc *psc, mparticles_base_t *particles_base,
			     double dt)
{
  mparticles_c_t *particles = psc_mparticles_get_c(particles_base, 0);

  psc_foreach_patch(psc, p) {
    do_shift_particle_positions(psc_mparticles_get_patch_c(particles, p), dt);
  }

  psc_mparticles_put_c(particles, particles_base, 0);
}

// ----------------------------------------------------------------------
// psc_calc_div_j

static void
do_calc_div_j(struct psc *psc, int p, fields_c_t *flds, fields_c_t *div_j)
{
  creal h[3];
  for (int d = 0; d < 3; d++) {
    if (psc->domain.gdims[d] == 1) {
      h[d] = 0.;
    } else {
      h[d] = 1. / psc->dx[d];
    }
  }

  psc_foreach_3d_g(psc, p, jx, jy, jz) {
    F3_C(div_j,0, jx,jy,jz) =
      (F3_C(flds,JXI, jx,jy,jz) - F3_C(flds,JXI, jx-1,jy,jz)) * h[0] +
      (F3_C(flds,JYI, jx,jy,jz) - F3_C(flds,JYI, jx,jy-1,jz)) * h[1] +
      (F3_C(flds,JZI, jx,jy,jz) - F3_C(flds,JZI, jx,jy,jz-1)) * h[2];
  } psc_foreach_3d_g_end;
}

static void
psc_calc_div_j(struct psc *psc, mfields_base_t *flds_base, mfields_base_t *div_j_base)
{
  mfields_c_t *flds = psc_mfields_get_c(flds_base, JXI, JXI + 3);
  mfields_c_t *div_j = psc_mfields_get_c(div_j_base, 0, 0);

  psc_foreach_patch(psc, p) {
    do_calc_div_j(psc, p, psc_mfields_get_patch_c(flds, p),
		  psc_mfields_get_patch_c(div_j, p));
  }

  psc_mfields_put_c(flds, flds_base, 0, 0);
  psc_mfields_put_c(div_j, div_j_base, 0, 1);
}

// ======================================================================
// psc_check_continuity
//
// given particles and currents at time $t^n$, checks whether discrete
// charge continuity holds

static void
psc_calc_rho(struct psc *psc, mparticles_base_t *particles, mfields_c_t *rho)
{
  mfields_c_t *dens = fld_create(psc, 3);

  struct psc_output_fields_item *item = psc_output_fields_item_create(psc_comm(psc));
  psc_output_fields_item_set_type(item, "n");
  psc_output_fields_item_set_psc_bnd(item, psc->bnd);
  psc_output_fields_item_setup(item);
  psc_output_fields_item_run(item, NULL, particles, dens);
  // rho = ni - ne
  psc_mfields_copy_comp(rho, 0, dens, 1); // FIXME, waxpy would be nicer
  psc_mfields_axpy_comp(rho, 0, -1., dens, 0);

  psc_mfields_destroy(dens);

  psc_output_fields_item_destroy(item);
}

void
psc_check_continuity(struct psc *psc, mparticles_base_t *particles,
		     mfields_base_t *flds, double eps)
{
  mfields_c_t *rho_m = fld_create(psc, 1);
  mfields_c_t *rho_p = fld_create(psc, 1);
  mfields_c_t *div_j = fld_create(psc, 1);

  psc_shift_particle_positions(psc, particles, -.5 * psc->dt);
  psc_calc_rho(psc, particles, rho_m);
  psc_shift_particle_positions(psc, particles,  1. * psc->dt);
  psc_calc_rho(psc, particles, rho_p);
  psc_shift_particle_positions(psc, particles, -.5 * psc->dt);

  psc_mfields_axpy(rho_p, -1., rho_m);
  psc_mfields_scale(rho_p, 1. / psc->dt);

  psc_bnd_fill_ghosts(psc->bnd, flds, JXI, JXI + 3);
  psc_calc_div_j(psc, flds, div_j);

  //  psc_dump_field(div_j, 0, "div_j");
  //  psc_dump_field(rho_p, 0, "dt_rho");

  // FIXME, there's gotta be a nicer way (add, norms)

  double max_err = 0.;
  psc_foreach_patch(psc, p) {
    fields_c_t *p_rho_p = psc_mfields_get_patch_c(rho_p, p);
    fields_c_t *p_div_j = psc_mfields_get_patch_c(div_j, p);
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      creal dt_rho = F3_C(p_rho_p,0, jx,jy,jz);
      creal div_j = F3_C(p_div_j,0, jx,jy,jz);
      max_err = fmax(max_err, fabs(dt_rho + div_j));
      if (fabs(dt_rho + div_j) > eps) {
	printf("(%d,%d,%d): %g -- %g diff %g\n", jx, jy, jz,
	       dt_rho, -div_j, dt_rho + div_j);
      }
    } psc_foreach_3d_end;
  }
  if (opt_checks_verbose) {
    mprintf("continuity: max_err = %g (thres %g)\n", max_err, eps);
  }

  //  psc_mfields_axpy(rho_p, +1., div_j);
  //  psc_dump_field(rho_p, 0, "cont_diff");

  assert(max_err <= eps);

  psc_mfields_destroy(rho_m);
  psc_mfields_destroy(rho_p);
  psc_mfields_destroy(div_j);
}

