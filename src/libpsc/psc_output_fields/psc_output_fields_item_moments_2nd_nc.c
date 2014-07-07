
#include "psc_output_fields_item_private.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include <math.h>

// ======================================================================

typedef fields_c_real_t creal;

static void
do_n_2nd_nc_run(fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[prts->p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);
      
    creal u = part->xi * dxi;
    creal v = part->yi * dyi;
    creal w = part->zi * dzi;
    int j1 = particle_real_nint(u);
    int j2 = particle_real_nint(v);
    int j3 = particle_real_nint(w);
    creal h1 = j1-u;
    creal h2 = j2-v;
    creal h3 = j3-w;
      
    creal gmx=.5f*(.5f+h1)*(.5f+h1);
    creal gmy=.5f*(.5f+h2)*(.5f+h2);
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0x=.75f-h1*h1;
    creal g0y=.75f-h2*h2;
    creal g0z=.75f-h3*h3;
    creal g1x=.5f*(.5f-h1)*(.5f-h1);
    creal g1y=.5f*(.5f-h2)*(.5f-h2);
    creal g1z=.5f*(.5f-h3)*(.5f-h3);
      
    int j1d = 1, j2d = 1, j3d = 1;
    if (ppsc->domain.gdims[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.; j1d = 0;
    }
    if (ppsc->domain.gdims[1] == 1) {
      j2 = 0; gmy = 0.; g0y = 1.; g1y = 0.; j2d = 0;
    }
    if (ppsc->domain.gdims[2] == 1) {
      j3 = 0; gmz = 0.; g0z = 1.; g1z = 0.; j3d = 0;
    }

    creal fnq;
    int m;
    if (part->qni < 0.) {
      fnq = part->qni * part->wni * fnqs;
      m = 0;
    } else if (part->qni > 0.) {
      fnq = part->qni * part->wni * fnqs;
      m = 1;
    } else {
      fnq = part->wni * fnqs;
      m = 2;
    }
    F3(pf, m, j1-j1d,j2-j2d,j3-j3d) += fnq*gmx*gmy*gmz;
    F3(pf, m, j1    ,j2-j2d,j3-j3d) += fnq*g0x*gmy*gmz;
    F3(pf, m, j1+j1d,j2-j2d,j3-j3d) += fnq*g1x*gmy*gmz;
    F3(pf, m, j1-j1d,j2    ,j3-j3d) += fnq*gmx*g0y*gmz;
    F3(pf, m, j1    ,j2    ,j3-j3d) += fnq*g0x*g0y*gmz;
    F3(pf, m, j1+j1d,j2    ,j3-j3d) += fnq*g1x*g0y*gmz;
    F3(pf, m, j1-j1d,j2+j2d,j3-j3d) += fnq*gmx*g1y*gmz;
    F3(pf, m, j1    ,j2+j2d,j3-j3d) += fnq*g0x*g1y*gmz;
    F3(pf, m, j1+j1d,j2+j2d,j3-j3d) += fnq*g1x*g1y*gmz;
    F3(pf, m, j1-j1d,j2-j2d,j3    ) += fnq*gmx*gmy*g0z;
    F3(pf, m, j1    ,j2-j2d,j3    ) += fnq*g0x*gmy*g0z;
    F3(pf, m, j1+j1d,j2-j2d,j3    ) += fnq*g1x*gmy*g0z;
    F3(pf, m, j1-j1d,j2    ,j3    ) += fnq*gmx*g0y*g0z;
    F3(pf, m, j1    ,j2    ,j3    ) += fnq*g0x*g0y*g0z;
    F3(pf, m, j1+j1d,j2    ,j3    ) += fnq*g1x*g0y*g0z;
    F3(pf, m, j1-j1d,j2+j2d,j3    ) += fnq*gmx*g1y*g0z;
    F3(pf, m, j1    ,j2+j2d,j3    ) += fnq*g0x*g1y*g0z;
    F3(pf, m, j1+j1d,j2+j2d,j3    ) += fnq*g1x*g1y*g0z;
    F3(pf, m, j1-j1d,j2-j2d,j3+j3d) += fnq*gmx*gmy*g1z;
    F3(pf, m, j1    ,j2-j2d,j3+j3d) += fnq*g0x*gmy*g1z;
    F3(pf, m, j1+j1d,j2-j2d,j3+j3d) += fnq*g1x*gmy*g1z;
    F3(pf, m, j1-j1d,j2    ,j3+j3d) += fnq*gmx*g0y*g1z;
    F3(pf, m, j1    ,j2    ,j3+j3d) += fnq*g0x*g0y*g1z;
    F3(pf, m, j1+j1d,j2    ,j3+j3d) += fnq*g1x*g0y*g1z;
    F3(pf, m, j1-j1d,j2+j2d,j3+j3d) += fnq*gmx*g1y*g1z;
    F3(pf, m, j1    ,j2+j2d,j3+j3d) += fnq*g0x*g1y*g1z;
    F3(pf, m, j1+j1d,j2+j2d,j3+j3d) += fnq*g1x*g1y*g1z;
  }
}

static void
n_2nd_nc_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	     struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_n_2nd_nc_run(res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

// FIXME too much duplication, specialize 2d/1d

static void
do_v_2nd_nc_run(fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[prts->p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);

    creal u = part->xi * dxi;
    creal v = part->yi * dyi;
    creal w = part->zi * dzi;
    int j1 = particle_real_nint(u);
    int j2 = particle_real_nint(v);
    int j3 = particle_real_nint(w);
    creal h1 = j1-u;
    creal h2 = j2-v;
    creal h3 = j3-w;

    creal gmx=.5f*(.5f+h1)*(.5f+h1);
    creal gmy=.5f*(.5f+h2)*(.5f+h2);
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0x=.75f-h1*h1;
    creal g0y=.75f-h2*h2;
    creal g0z=.75f-h3*h3;
    creal g1x=.5f*(.5f-h1)*(.5f-h1);
    creal g1y=.5f*(.5f-h2)*(.5f-h2);
    creal g1z=.5f*(.5f-h3)*(.5f-h3);

    if (ppsc->domain.gdims[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.;
    }
    if (ppsc->domain.gdims[1] == 1) {
      j2 = 0; gmy = 0.; g0y = 1.; g1y = 0.;
    }
    if (ppsc->domain.gdims[2] == 1) {
      j3 = 0; gmz = 0.; g0z = 1.; g1z = 0.;
    }
    
    creal pxi = part->pxi;
    creal pyi = part->pyi;
    creal pzi = part->pzi;
    creal root = 1.0/sqrt(1.0+pxi*pxi+pyi*pyi+pzi*pzi);
    creal vv[3] = { pxi*root, pyi*root, pzi*root };
    creal fnq = part->wni * fnqs;
    int mm;
    if (part->qni < 0.) {
      mm = 0; // electrons
    } else if (part->qni > 0.) {
      mm = 3; // ions
    } else {
      assert(0);
    }
    for (int m = 0; m < 3; m++) {
      F3(pf, mm+m, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz * vv[m];
      F3(pf, mm+m, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz * vv[m];
      F3(pf, mm+m, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz * vv[m];
      F3(pf, mm+m, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz * vv[m];
      F3(pf, mm+m, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz * vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz * vv[m];
      F3(pf, mm+m, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz * vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz * vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz * vv[m];
      F3(pf, mm+m, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z * vv[m];
      F3(pf, mm+m, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z * vv[m];
      F3(pf, mm+m, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z * vv[m];
      F3(pf, mm+m, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z * vv[m];
      F3(pf, mm+m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z * vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z * vv[m];
      F3(pf, mm+m, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z * vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z * vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z * vv[m];
      F3(pf, mm+m, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z * vv[m];
      F3(pf, mm+m, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z * vv[m];
      F3(pf, mm+m, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z * vv[m];
      F3(pf, mm+m, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z * vv[m];
      F3(pf, mm+m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z * vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z * vv[m];
      F3(pf, mm+m, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z * vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z * vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z * vv[m];
    }  
  }
}

static void
v_2nd_nc_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	     struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_v_2nd_nc_run(res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

static void
do_vv_2nd_nc_run(fields_t *pf, struct psc_particles *prts)
{
  struct psc_patch *patch = &ppsc->patch[prts->p];
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / patch->dx[0], dyi = 1.f / patch->dx[1], dzi = 1.f / patch->dx[2];

  for (int n = 0; n < prts->n_part; n++) {
    particle_t *part = particles_get_one(prts, n);

    creal u = part->xi * dxi;
    creal v = part->yi * dyi;
    creal w = part->zi * dzi;
    int j1 = particle_real_nint(u);
    int j2 = particle_real_nint(v);
    int j3 = particle_real_nint(w);
    creal h1 = j1-u;
    creal h2 = j2-v;
    creal h3 = j3-w;

    creal gmx=.5f*(.5f+h1)*(.5f+h1);
    creal gmy=.5f*(.5f+h2)*(.5f+h2);
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0x=.75f-h1*h1;
    creal g0y=.75f-h2*h2;
    creal g0z=.75f-h3*h3;
    creal g1x=.5f*(.5f-h1)*(.5f-h1);
    creal g1y=.5f*(.5f-h2)*(.5f-h2);
    creal g1z=.5f*(.5f-h3)*(.5f-h3);

    if (ppsc->domain.gdims[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.;
    }
    if (ppsc->domain.gdims[1] == 1) {
      j2 = 0; gmy = 0.; g0y = 1.; g1y = 0.;
    }
    if (ppsc->domain.gdims[2] == 1) {
      j3 = 0; gmz = 0.; g0z = 1.; g1z = 0.;
    }
    
    creal pxi = part->pxi;
    creal pyi = part->pyi;
    creal pzi = part->pzi;
    creal root = 1.0/sqrt(1.0+pxi*pxi+pyi*pyi+pzi*pzi);
    creal vv[3] = { pxi*root, pyi*root, pzi*root };
    creal fnq = part->wni * fnqs;
    int mm;
    if (part->qni < 0.) {
      mm = 0; // electrons
    } else if (part->qni > 0.) {
      mm = 3; // ions
    } else {
      assert(0);
    }
    for (int m = 0; m < 3; m++) {
      F3(pf, mm+m, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z * vv[m]*vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z * vv[m]*vv[m];
    }  
  }
}

static void
vv_2nd_nc_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	      struct psc_particles *prts_base, struct psc_fields *res)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, PARTICLE_TYPE, 0);
  psc_fields_zero_range(res, 0, res->nr_comp);
  do_vv_2nd_nc_run(res, prts);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

// ======================================================================
// psc_output_fields_item: subclass "n_2nd_nc"

struct psc_output_fields_item_ops psc_output_fields_item_n_2nd_nc_ops = {
  .name               = "n",
  .nr_comp            = 3,
  .fld_names          = { "ne", "ni", "nn" },
  .run                = n_2nd_nc_run,
  .flags              = POFI_ADD_GHOSTS,
};

// ======================================================================
// psc_output_fields_item: subclass "v_2nd_nc"

struct psc_output_fields_item_ops psc_output_fields_item_v_2nd_nc_ops = {
  .name               = "v",
  .nr_comp            = 6,
  .fld_names          = { "vex", "vey", "vez",
			  "vix", "viy", "viz" },
  .run                = v_2nd_nc_run,
  .flags              = POFI_ADD_GHOSTS,
};

// ======================================================================
// psc_output_fields_item: subclass "vv_2nd_nc"

struct psc_output_fields_item_ops psc_output_fields_item_vv_2nd_nc_ops = {
  .name               = "vv",
  .nr_comp            = 6,
  .fld_names          = { "vexvex", "veyvey", "vezvez",
			  "vixvix", "viyviy", "vizviz" },
  .run                = vv_2nd_nc_run,
  .flags              = POFI_ADD_GHOSTS,
};

