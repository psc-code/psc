
#include "psc_output_fields_item_private.h"

#include "psc_fields_as_c.h"
#include "psc_particles_as_c.h"

#include <math.h>

// ======================================================================
// n_photon

static void
do_n_photon_run(int p, fields_t *pf, photons_t *photons)
{
  photon_real_t dxi = 1.f / ppsc->patch[p].dx[0], dyi = 1.f / ppsc->patch[p].dx[1], dzi = 1.f / ppsc->patch[p].dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < photons->nr; n++) {
    photon_t *p = photons_get_one(photons, n);
      
    photon_real_t u = (p->x[0] - patch->xb[0]) * dxi;
    photon_real_t v = (p->x[1] - patch->xb[1]) * dyi;
    photon_real_t w = (p->x[2] - patch->xb[2]) * dzi;
    int j1 = particle_real_nint(u);
    int j2 = particle_real_nint(v);
    int j3 = particle_real_nint(w);
    photon_real_t h1 = j1-u;
    photon_real_t h2 = j2-v;
    photon_real_t h3 = j3-w;
      
    photon_real_t gmx=.5f*(.5f+h1)*(.5f+h1);
    photon_real_t gmy=.5f*(.5f+h2)*(.5f+h2);
    photon_real_t gmz=.5f*(.5f+h3)*(.5f+h3);
    photon_real_t g0x=.75f-h1*h1;
    photon_real_t g0y=.75f-h2*h2;
    photon_real_t g0z=.75f-h3*h3;
    photon_real_t g1x=.5f*(.5f-h1)*(.5f-h1);
    photon_real_t g1y=.5f*(.5f-h2)*(.5f-h2);
    photon_real_t g1z=.5f*(.5f-h3)*(.5f-h3);
      
    if (ppsc->domain.gdims[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.;
    }
    if (ppsc->domain.gdims[1] == 1) {
      j2 = 0; gmy = 0.; g0y = 1.; g1y = 0.;
    }
    if (ppsc->domain.gdims[2] == 1) {
      j3 = 0; gmz = 0.; g0z = 1.; g1z = 0.;
    }

    int m = 0;
    photon_real_t fnq = p->wni;
    F3(pf, m, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz;
    F3(pf, m, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz;
    F3(pf, m, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz;
    F3(pf, m, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz;
    F3(pf, m, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz;
    F3(pf, m, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz;
    F3(pf, m, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz;
    F3(pf, m, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz;
    F3(pf, m, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz;
    F3(pf, m, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z;
    F3(pf, m, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z;
    F3(pf, m, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z;
    F3(pf, m, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z;
    F3(pf, m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z;
    F3(pf, m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z;
    F3(pf, m, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z;
    F3(pf, m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z;
    F3(pf, m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z;
    F3(pf, m, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z;
    F3(pf, m, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z;
    F3(pf, m, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z;
    F3(pf, m, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z;
    F3(pf, m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z;
    F3(pf, m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z;
    F3(pf, m, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z;
    F3(pf, m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z;
    F3(pf, m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z;
  }
}

static void
n_photon_run(struct psc_output_fields_item *item, struct psc_mfields *mflds,
	     struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  for (int p = 0; p < mres->nr_patches; p++) {
    struct psc_fields *res = psc_mfields_get_patch(mres, p);
    psc_fields_zero_range(res, 0, res->nr_comp);
    do_n_photon_run(p, res, &ppsc->mphotons->p[p]);
  }
}

// ======================================================================
// psc_output_fields_item: subclass "n_photon"

struct psc_output_fields_item_ops psc_output_fields_item_n_photon_ops = {
  .name               = "n_photon",
  .nr_comp            = 1,
  .fld_names          = { "n_photon" },
  .run_all            = n_photon_run,
  .flags              = POFI_ADD_GHOSTS,
};

