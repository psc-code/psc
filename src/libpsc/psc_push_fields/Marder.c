
#include "psc.h"
#include "psc_bnd.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_double.h"
#include <math.h>
#include "psc_output_fields/common_moments.c"

enum {DIVE_MARDER, N_MARDER, NR_MARDER};

/* store parameters; not used for now
struct psc_gauss_correction_Marder_c {
  int step;
  float diffusion;
  struct psc_bnd *bnd;
};
*/

// ======================================================================
// Create/Destroy space for "aid" fields, i.e., divE and n

mfields_t *
Marder_create_aid_fields(struct psc *psc)
{
  mfields_base_t *f = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(f, FIELDS_TYPE);
  psc_mfields_set_domain(f, psc->mrc_domain);
  psc_mfields_set_param_int(f, "nr_fields", NR_MARDER);
  psc_mfields_set_param_int(f, "first_comp", DIVE_MARDER);
  psc_mfields_set_param_int3(f, "ibn", psc->ibn);
  psc_mfields_setup(f);
  return f;
}

static void
Marder_destroy_aid_fields(mfields_t *f)
{
  psc_mfields_destroy(f);
}

// ======================================================================
// Calculate divE at node centers, copied from psc_gauss_correction_item_jeh.c

#define define_dxdydz(dx, dy, dz)					\
  int dx __unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy __unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz __unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

static void
psc_push_fields_sub_push_E(struct psc_push_fields *push, struct psc_fields *flds,
			   double dt)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_real_t cnx = .5 * dt / patch->dx[0];
  fields_real_t cny = .5 * dt / patch->dx[1];
  fields_real_t cnz = .5 * dt / patch->dx[2];

  if (ppsc->domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (ppsc->domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (ppsc->domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // E-field propagation E^(n)    , H^(n), j^(n) 
  //                  -> E^(n+0.5), H^(n), j^(n)
  // Ex^{n}[-.5:+.5][-1:1][-1:1] -> Ex^{n+.5}[-.5:+.5][-1:1][-1:1]
  // using Hx^{n}[-1:1][-1.5:1.5][-1.5:1.5]
  //       jx^{n+1}[-.5:.5][-1:1][-1:1]

  psc_foreach_3d(ppsc, flds->p, ix, iy, iz, 1, 2) {
    F3(flds, EX, ix,iy,iz) +=
      cny * (F3(flds, HZ, ix,iy,iz) - F3(flds, HZ, ix,iy-1,iz)) -
      cnz * (F3(flds, HY, ix,iy,iz) - F3(flds, HY, ix,iy,iz-1)) -
      .5 * dt * F3(flds, JXI, ix,iy,iz);
    
    F3(flds, EY, ix,iy,iz) +=
      cnz * (F3(flds, HX, ix,iy,iz) - F3(flds, HX, ix,iy,iz-1)) -
      cnx * (F3(flds, HZ, ix,iy,iz) - F3(flds, HZ, ix-1,iy,iz)) -
      .5 * dt * F3(flds, JYI, ix,iy,iz);
    
    F3(flds, EZ, ix,iy,iz) +=
      cnx * (F3(flds, HY, ix,iy,iz) - F3(flds, HY, ix-1,iy,iz)) -
      cny * (F3(flds, HX, ix,iy,iz) - F3(flds, HX, ix,iy-1,iz)) -
      .5 * dt * F3(flds, JZI, ix,iy,iz);
  } foreach_3d_end;
}

static void
calc_dive_nc(struct psc_fields *flds_base, 
  struct psc_particles *prts, struct psc_fields *f, int dive_comp)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, HX + 3);
  // advance to (n+.5) dt, skipped; FIXME: should we do this?
  // psc_push_fields_sub_push_E(NULL, flds, ppsc->dt);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, dive_comp, ix,iy,iz) = 
      ((F3(flds, EX, ix,iy,iz) - F3(flds, EX, ix-dx,iy,iz)) / ppsc->patch[f->p].dx[0] +
       (F3(flds, EY, ix,iy,iz) - F3(flds, EY, ix,iy-dy,iz)) / ppsc->patch[f->p].dx[1] +
       (F3(flds, EZ, ix,iy,iz) - F3(flds, EZ, ix,iy,iz-dz)) / ppsc->patch[f->p].dx[2]);
  } foreach_3d_end;
  // back to n dt, skipped
  // psc_push_fields_sub_push_E(NULL, flds, -ppsc->dt);
  psc_fields_put_as(flds, flds_base, 0, 0);
}

// ======================================================================
// Calculate total n_nc, copied from psc_gauss_correction_moments_1st_nc.c

static void
do_n_run(int p, fields_t *pf, struct psc_particles *prts, int n_comp)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = /* particle_kind(part) */ n_comp; // calculate total n only
    DEPOSIT_TO_GRID_1ST_NC(part, pf, m,  
      ppsc->kinds[particle_kind(part)].q);
  }
}

static void
n_run(struct psc_fields *flds,
      struct psc_particles *prts_base, struct psc_fields *res, int n_comp)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, n_comp, n_comp+1);
  do_n_run(res->p, res, prts, n_comp);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

// ======================================================================

static void
Marder_calc_aid_fields(mfields_base_t *flds, mparticles_base_t *particles, mfields_t *res)
{
  for (int p = 0; p < res->nr_patches; p++) {
    calc_dive_nc(psc_mfields_get_patch(flds, p),
	     psc_mparticles_get_patch(particles, p),
	     psc_mfields_get_patch(res, p), DIVE_MARDER);
  }

  for (int p = 0; p < res->nr_patches; p++) {
    n_run(psc_mfields_get_patch(flds, p),
	     psc_mparticles_get_patch(particles, p),
	     psc_mfields_get_patch(res, p), N_MARDER);
  }
  psc_bnd_add_ghosts(ppsc->bnd, res, N_MARDER, N_MARDER+1); //+1 or not, check!!!
}

// ======================================================================
// Do the modified Marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

static void
do_Marder_correction(struct psc_fields *flds_base, struct psc_particles *prts, 
  struct psc_fields *f, int comp_dive, int comp_n)
{
  define_dxdydz(dx, dy, dz);

  // FIXME: how to choose diffusion parameter properly?
  double deltax = ppsc->patch[f->p].dx[0];
  double deltay = ppsc->patch[f->p].dx[1];
  double deltaz = ppsc->patch[f->p].dx[2];
  double inv_sum = 0.;
  int nr_levels;
  mrc_domain_get_nr_levels(ppsc->mrc_domain, &nr_levels);
  for (int d=0;d<3;d++) {
    if (ppsc->domain.gdims[d] > 1) {
      inv_sum += 1. / sqr(ppsc->patch[f->p].dx[d] / (1 << (nr_levels - 1)));
    }
  }
  double diffusion_max = 1. / 2. / (.5 * ppsc->dt) / inv_sum;
  double diffusion     = diffusion_max * .75;

  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, EX + 3);
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    // FIXME: F3 correct?
    F3(flds, EX, ix,iy,iz) += 
      (  F3(f, DIVE_MARDER, ix+dx,iy,iz) - F3(f, DIVE_MARDER, ix,iy,iz)
       - F3(f, N_MARDER,    ix+dx,iy,iz) + F3(f, N_MARDER,    ix,iy,iz)
       ) * .5 * ppsc->dt * diffusion / deltax;
    F3(flds, EY, ix,iy,iz) += 
      (  F3(f, DIVE_MARDER, ix,iy+dy,iz) - F3(f, DIVE_MARDER, ix,iy,iz)
       - F3(f, N_MARDER,    ix,iy+dy,iz) + F3(f, N_MARDER,    ix,iy,iz)
       ) * .5 * ppsc->dt * diffusion / deltay;
    F3(flds, EZ, ix,iy,iz) += 
      (  F3(f, DIVE_MARDER, ix,iy,iz+dz) - F3(f, DIVE_MARDER, ix,iy,iz)
       - F3(f, N_MARDER,    ix,iy,iz+dz) + F3(f, N_MARDER,    ix,iy,iz)
       ) * .5 * ppsc->dt * diffusion / deltaz;
  } foreach_3d_end;
  psc_fields_put_as(flds, flds_base, JXI, EX + 3);
}

void
Marder_correction_run(mfields_base_t *flds, mparticles_base_t *particles, mfields_t *res)
{
  for (int p = 0; p < res->nr_patches; p++) {
    do_Marder_correction(psc_mfields_get_patch(flds, p),
	     psc_mparticles_get_patch(particles, p),
	     psc_mfields_get_patch(res, p), DIVE_MARDER, N_MARDER);
  }
}

void
Marder_correction(mfields_base_t *flds, mparticles_base_t *particles)
{
  mfields_t *res = Marder_create_aid_fields(ppsc);
  Marder_calc_aid_fields(flds, particles, res);
  Marder_correction_run(flds, particles, res);
  Marder_destroy_aid_fields(res);
}

