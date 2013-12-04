
#include "psc.h"
#include "psc_bnd.h"
#include "psc_push_fields_private.h"
#include "psc_fields_as_single.h"
#include "psc_particles_as_single.h"
#include <math.h>
#include "../psc_output_fields/common_moments.c"

enum {
  DIVE_MARDER,
  N_MARDER,
  NR_MARDER
};

// ======================================================================
// Create/Destroy space for "aid" fields, i.e., divE and n

mfields_t *
marder_create_aid_fields(struct psc_push_fields *push, struct psc *psc)
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
marder_destroy_aid_fields(struct psc_push_fields *push, mfields_t *f)
{
  psc_mfields_destroy(f);
}

// ======================================================================
// Calculate divE at node centers, copied from psc_gauss_correction_item_jeh.c

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

static void
calc_dive_nc(struct psc_fields *flds_base, 
  struct psc_particles *prts, struct psc_fields *f, int dive_comp)
{
  define_dxdydz(dx, dy, dz);
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, HX + 3);

  // when we get called, the E fields have been advanced to (n+.5) dt
  psc_foreach_3d(ppsc, f->p, ix, iy, iz, 0, 0) {
    F3(f, dive_comp, ix,iy,iz) = 
      ((F3(flds, EX, ix,iy,iz) - F3(flds, EX, ix-dx,iy,iz)) / ppsc->patch[f->p].dx[0] +
       (F3(flds, EY, ix,iy,iz) - F3(flds, EY, ix,iy-dy,iz)) / ppsc->patch[f->p].dx[1] +
       (F3(flds, EZ, ix,iy,iz) - F3(flds, EZ, ix,iy,iz-dz)) / ppsc->patch[f->p].dx[2]);
  } foreach_3d_end;

  psc_fields_put_as(flds, flds_base, 0, 0);
}

// ======================================================================
// Calculate total n_nc, copied from psc_gauss_correction_moments_1st_nc.c

static void
do_n_run(int p, fields_t *pf, int m, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    DEPOSIT_TO_GRID_1ST_NC(part, pf, m,  
      ppsc->kinds[particle_kind(part)].q);
  }
}

static void
n_run(struct psc_fields *flds,
      struct psc_particles *prts_base, struct psc_fields *res, int m)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, m, m + 1);
  do_n_run(res->p, res, m, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

// ======================================================================

static void
marder_calc_aid_fields(struct psc_push_fields *push, 
  mfields_base_t *flds, mparticles_base_t *particles, mfields_t *res)
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
  psc_bnd_fill_ghosts(ppsc->bnd, res, DIVE_MARDER, N_MARDER+1);
}


#define psc_foreach_3d_more(psc, p, ix, iy, iz, l, r) {			\
  int __ilo[3] = { -l[0], -l[1], -l[2] };					\
  int __ihi[3] = { psc->patch[p].ldims[0] + r[0],				\
		   psc->patch[p].ldims[1] + r[1],				\
		   psc->patch[p].ldims[2] + r[2] };				\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {			\
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {			\
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define psc_foreach_3d_more_end				\
  } } }

// ======================================================================
// Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

static void
do_marder_correction(struct psc_push_fields *push, 
  struct psc_fields *flds_base, struct psc_particles *prts, 
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
  double diffusion     = diffusion_max * push->marder_diffusion;

  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, EX + 3);

  int l[3] = {0, 0, 0}, r[3] = {0, 0, 0};
  for (int d = 0; d < 3; d++) {
   if (ppsc->domain.bnd_fld_lo[d] == BND_FLD_CONDUCTING_WALL && ppsc->patch[flds->p].off[d] == 0) {
    l[d] = -2;
   }
   if (ppsc->domain.bnd_fld_hi[d] == BND_FLD_CONDUCTING_WALL && ppsc->patch[flds->p].off[d] + ppsc->patch[flds->p].ldims[d] == ppsc->domain.gdims[d]) {
    r[d] = -2;
   }
  }

  psc_foreach_3d_more(ppsc, f->p, ix, iy, iz, l, r) {
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
  } psc_foreach_3d_more_end;
  psc_fields_put_as(flds, flds_base, JXI, EX + 3);
}

void
marder_correction_run(struct psc_push_fields *push, 
  mfields_base_t *flds, mparticles_base_t *particles, mfields_t *res)
{
  for (int p = 0; p < res->nr_patches; p++) {
    do_marder_correction(push, psc_mfields_get_patch(flds, p),
	     psc_mparticles_get_patch(particles, p),
	     psc_mfields_get_patch(res, p), DIVE_MARDER, N_MARDER);
  }
}

void
marder_correction(struct psc_push_fields *push, 
  mfields_base_t *flds, mparticles_base_t *particles)
{
  if (push->marder_step < 0 || ppsc->timestep % push->marder_step != 0) 
   return;

  mfields_t *res = marder_create_aid_fields(push, ppsc);
  marder_calc_aid_fields(push, flds, particles, res);
  marder_correction_run(push, flds, particles, res);
  marder_destroy_aid_fields(push, res);
}

#undef psc_foreach_3d_more
#undef psc_foreach_3d_more_end
