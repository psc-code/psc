
#include "psc.h"

#include "psc_output_fields_item_private.h"
#include "psc_cuda.h"

#include <math.h>

// ======================================================================
// rho

#if 0
static void
do_rho_run(int p, struct psc_fields *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
    int m = particle_kind(part);
    particle_real_t *xi = &part->xi;
    particle_real_t v = xi[1] * dyi;
    particle_real_t w = xi[2] * dzi;
    int jy = particle_real_fint(v);
    int jz = particle_real_fint(w);
    particle_real_t h2 = v - jy;
    particle_real_t h3 = w - jz;

    particle_real_t g0y = 1.f - h2;
    particle_real_t g0z = 1.f - h3;
    particle_real_t g1y = h2;
    particle_real_t g1z = h3;

    assert(jy >= -1 && jy < patch->ldims[1]);
    assert(jz >= -1 && jz < patch->ldims[2]);

    particle_real_t fnq = particle_wni(part) * fnqs;
    particle_real_t val = ppsc->kinds[m].q;
    F3(pf, 0, 0,jy  ,jz  ) += fnq*g0y*g0z * val;
    F3(pf, 0, 0,jy+1,jz  ) += fnq*g1y*g0z * val;
    F3(pf, 0, 0,jy  ,jz+1) += fnq*g0y*g1z * val;
    F3(pf, 0, 0,jy+1,jz+1) += fnq*g1y*g1z * val;
  }
}
#endif

static void
rho_run_all(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	    struct psc_mparticles *mprts_base, struct psc_mfields *mres_base)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "cuda", 0);
  struct psc_mfields *mres = psc_mfields_get_as(mres_base, "cuda", 0, 0);
  psc_mfields_zero_range(mres, 0, mres->nr_fields);

  yz_moments_rho_1st_nc_cuda_run_patches(mprts, mres);

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
  psc_mfields_put_as(mres, mres_base, 0, 1);
}

// ======================================================================
// psc_output_fields_item: subclass "rho_1st_nc_cuda"

struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_cuda_ops = {
  .name               = "rho_1st_nc_cuda",
  .nr_comp            = 1,
  .fld_names          = { "rho_nc_cuda" }, // FIXME
  .run_all            = rho_run_all,
  .flags              = POFI_ADD_GHOSTS,
};

