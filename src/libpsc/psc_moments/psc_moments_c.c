
#include "psc_moments_private.h"
#include "psc_particles_as_c.h"

#include "psc_bnd.h"
#include <mrc_profile.h>
#include <math.h>

// ======================================================================

typedef fields_base_real_t creal;

static void
do_c_calc_densities(int p, fields_base_t *pf, particles_t *pp,
		    int m_NE, int m_NI, int m_NN)
{
  fields_base_zero(pf, m_NE);
  fields_base_zero(pf, m_NI);
  fields_base_zero(pf, m_NN);
  
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal dxi = 1.f / ppsc->dx[0];
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    creal u = (part->xi - patch->xb[0]) * dxi;
    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
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

    creal fnq;
    int m;
    if (part->qni < 0.) {
      fnq = part->qni * part->wni * fnqs;
      m = m_NE;
    } else if (part->qni > 0.) {
      fnq = part->qni * part->wni * fnqs;
      m = m_NI;
    } else {
      fnq = part->wni * fnqs;
      m = m_NN;
    }
    F3_BASE(pf,m, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz;
    F3_BASE(pf,m, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz;
    F3_BASE(pf,m, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz;
    F3_BASE(pf,m, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz;
    F3_BASE(pf,m, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz;
    F3_BASE(pf,m, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz;
    F3_BASE(pf,m, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz;
    F3_BASE(pf,m, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz;
    F3_BASE(pf,m, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz;
    F3_BASE(pf,m, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z;
    F3_BASE(pf,m, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z;
    F3_BASE(pf,m, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z;
    F3_BASE(pf,m, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z;
    F3_BASE(pf,m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z;
    F3_BASE(pf,m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z;
    F3_BASE(pf,m, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z;
    F3_BASE(pf,m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z;
    F3_BASE(pf,m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z;
    F3_BASE(pf,m, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z;
    F3_BASE(pf,m, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z;
    F3_BASE(pf,m, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z;
    F3_BASE(pf,m, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z;
    F3_BASE(pf,m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z;
    F3_BASE(pf,m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z;
    F3_BASE(pf,m, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z;
    F3_BASE(pf,m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z;
    F3_BASE(pf,m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z;
  }
}

static void
psc_moments_c_calc_densities(struct psc_moments *moments, mfields_base_t *flds,
			     mparticles_base_t *particles_base, mfields_base_t *res)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_densities", 1., 0, 0);
  }

  mparticles_t particles;
  psc_mparticles_get_from(&particles, particles_base);

  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_c_calc_densities(p, &res->f[p], &particles.p[p], 0, 1, 2);
  }
  prof_stop(pr);

  psc_bnd_add_ghosts(ppsc->bnd, res, 0, 3);

  psc_mparticles_put_to(&particles, particles_base); // FIXME, don't need copy-back
}

// FIXME too much duplication, specialize 2d/1d

static void
do_c_calc_v(int p, fields_base_t *pf, particles_t *pp)
{
  for (int m = 0; m < 6; m++) {
    fields_base_zero(pf, m);
  }
  
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal dxi = 1.f / ppsc->dx[0];
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    creal u = (part->xi - patch->xb[0]) * dxi;
    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
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
      F3_BASE(pf,mm+m, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z * vv[m];
      F3_BASE(pf,mm+m, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z * vv[m];
      F3_BASE(pf,mm+m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z * vv[m];
      F3_BASE(pf,mm+m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z * vv[m];
    }  
  }
}

static void
psc_moments_c_calc_v(struct psc_moments *moments, mfields_base_t *flds,
		     mparticles_base_t *particles_base, mfields_base_t *res)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_calc_v", 1., 0, 0);
  }

  mparticles_t particles;
  psc_mparticles_get_from(&particles, particles_base);

  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_c_calc_v(p, &res->f[p], &particles.p[p]);
  }
  prof_stop(pr);

  psc_bnd_add_ghosts(ppsc->bnd, res, 0, 6);
  psc_mparticles_put_to(&particles, particles_base); // FIXME, don't need copy-back
}

static void
do_c_calc_vv(int p, fields_base_t *pf, particles_t *pp)
{
  for (int m = 0; m < 6; m++) {
    fields_base_zero(pf, m);
  }
  
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal dxi = 1.f / ppsc->dx[0];
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    creal u = (part->xi - patch->xb[0]) * dxi;
    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
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
      F3_BASE(pf,mm+m, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z * vv[m]*vv[m];
      F3_BASE(pf,mm+m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z * vv[m]*vv[m];
    }  
  }
}

static void
psc_moments_c_calc_vv(struct psc_moments *moments, mfields_base_t *flds,
		      mparticles_base_t *particles_base, mfields_base_t *res)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_calc_vv", 1., 0, 0);
  }
  mparticles_t particles;
  psc_mparticles_get_from(&particles, particles_base);

  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_c_calc_vv(p, &res->f[p], &particles.p[p]);
  }
  prof_stop(pr);

  psc_bnd_add_ghosts(ppsc->bnd, res, 0, 6);
  psc_mparticles_put_to(&particles, particles_base);
}

static void
do_c_calc_photon_n(int p, fields_base_t *pf, photons_t *photons)
{
  fields_base_zero(pf, 0);
  
  creal dxi = 1.f / ppsc->dx[0];
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < photons->nr; n++) {
    photon_t *p = photons_get_one(photons, n);
      
    creal u = (p->x[0] - patch->xb[0]) * dxi;
    creal v = (p->x[1] - patch->xb[1]) * dyi;
    creal w = (p->x[2] - patch->xb[2]) * dzi;
    int j1 = particle_base_real_nint(u);
    int j2 = particle_base_real_nint(v);
    int j3 = particle_base_real_nint(w);
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

    int m = 0;
    creal fnq = p->wni;
    F3_BASE(pf,m, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz;
    F3_BASE(pf,m, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz;
    F3_BASE(pf,m, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz;
    F3_BASE(pf,m, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz;
    F3_BASE(pf,m, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz;
    F3_BASE(pf,m, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz;
    F3_BASE(pf,m, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz;
    F3_BASE(pf,m, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz;
    F3_BASE(pf,m, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz;
    F3_BASE(pf,m, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z;
    F3_BASE(pf,m, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z;
    F3_BASE(pf,m, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z;
    F3_BASE(pf,m, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z;
    F3_BASE(pf,m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z;
    F3_BASE(pf,m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z;
    F3_BASE(pf,m, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z;
    F3_BASE(pf,m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z;
    F3_BASE(pf,m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z;
    F3_BASE(pf,m, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z;
    F3_BASE(pf,m, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z;
    F3_BASE(pf,m, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z;
    F3_BASE(pf,m, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z;
    F3_BASE(pf,m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z;
    F3_BASE(pf,m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z;
    F3_BASE(pf,m, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z;
    F3_BASE(pf,m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z;
    F3_BASE(pf,m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z;
  }
}

static void
psc_moments_c_calc_photon_n(struct psc_moments *moments,
			    mphotons_t *mphotons, mfields_base_t *res)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_photon_n", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_c_calc_photon_n(p, &res->f[p], &mphotons->p[p]);
  }
  prof_stop(pr);

  psc_bnd_add_ghosts(ppsc->bnd, res, 0, 1);
}

// ======================================================================
// psc_moments: subclass "c"

struct psc_moments_ops psc_moments_c_ops = {
  .name                  = "c",
  .calc_densities        = psc_moments_c_calc_densities,
  .calc_v                = psc_moments_c_calc_v,
  .calc_vv               = psc_moments_c_calc_vv,
  .calc_photon_n         = psc_moments_c_calc_photon_n,
};
