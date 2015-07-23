
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_output_fields_item_private.h"

#include <math.h>

#define DEPOSIT_TO_GRID_1ST_NC(part, pf, m, val) do {			\
    particle_real_t *xi = &part->xi; /* don't shift back in time */	\
    particle_real_t u = xi[0] * dxi;					\
    particle_real_t v = xi[1] * dyi;					\
    particle_real_t w = xi[2] * dzi;					\
    int jx = particle_real_fint(u);					\
    int jy = particle_real_fint(v);					\
    int jz = particle_real_fint(w);					\
    particle_real_t h1 = u - jx;					\
    particle_real_t h2 = v - jy;					\
    particle_real_t h3 = w - jz;					\
    									\
    particle_real_t g0x = 1.f - h1;					\
    particle_real_t g0y = 1.f - h2;					\
    particle_real_t g0z = 1.f - h3;					\
    particle_real_t g1x = h1;						\
    particle_real_t g1y = h2;						\
    particle_real_t g1z = h3;						\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (ppsc->domain.gdims[0] == 1) {					\
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;				\
    }									\
    if (ppsc->domain.gdims[1] == 1) {					\
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;				\
    }									\
    if (ppsc->domain.gdims[2] == 1) {					\
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;				\
    }									\
    									\
    assert(jx >= -1 && jx < patch->ldims[0]);				\
    assert(jy >= -1 && jy < patch->ldims[1]);				\
    assert(jz >= -1 && jz < patch->ldims[2]);				\
    									\
    particle_real_t fnq = particle_wni(part) * fnqs;			\
									\
    F3(pf, m, jx    ,jy    ,jz    ) += fnq*g0x*g0y*g0z * (val);		\
    F3(pf, m, jx+jxd,jy    ,jz    ) += fnq*g1x*g0y*g0z * (val);		\
    F3(pf, m, jx    ,jy+jyd,jz    ) += fnq*g0x*g1y*g0z * (val);		\
    F3(pf, m, jx+jxd,jy+jyd,jz    ) += fnq*g1x*g1y*g0z * (val);		\
    F3(pf, m, jx    ,jy    ,jz+jzd) += fnq*g0x*g0y*g1z * (val);		\
    F3(pf, m, jx+jxd,jy    ,jz+jzd) += fnq*g1x*g0y*g1z * (val);		\
    F3(pf, m, jx    ,jy+jyd,jz+jzd) += fnq*g0x*g1y*g1z * (val);		\
    F3(pf, m, jx+jxd,jy+jyd,jz+jzd) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

// ======================================================================
// rho

static void
do_rho_run(int p, fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    DEPOSIT_TO_GRID_1ST_NC(part, pf, 0, ppsc->kinds[m].q);
  }
}

static void
rho_run_patches(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		struct psc_mparticles *mprts_base, struct psc_mfields *mres_base)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);
  struct psc_mfields *mres = psc_mfields_get_as(mres_base, FIELDS_TYPE, 0, 0);
  psc_mfields_zero_range(mres, 0, mres->nr_fields);
  for (int p = 0; p < mres->nr_patches; p++) {
    do_rho_run(p, psc_mfields_get_patch(mres, p), psc_mparticles_get_patch(mprts, p));
  }
  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
  psc_mfields_put_as(mres, mres_base, 0, 1);
}

// ======================================================================
// psc_output_fields_item: subclass "rho_1st_nc_cuda"

struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_cuda_ops = {
  .name               = "rho_1st_nc_cuda",
  .nr_comp            = 1,
  .fld_names          = { "rho_nc_cuda" }, // FIXME
  .run_patches        = rho_run_patches,
  .flags              = POFI_ADD_GHOSTS,
};

