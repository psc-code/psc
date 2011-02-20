
#include "psc.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef fields_base_real_t creal;

static inline int
nint(creal x)
{
  return (int)(x + 10.5f) - 10;
}

static void
do_c_calc_densities(fields_base_t *pf, int m_NE, int m_NI, int m_NN)
{
  fields_base_zero(pf, m_NE);
  fields_base_zero(pf, m_NI);
  fields_base_zero(pf, m_NN);
  
  creal fnqs = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
  creal dxi = 1.f / psc.dx[0];
  creal dyi = 1.f / psc.dx[1];
  creal dzi = 1.f / psc.dx[2];

  struct psc_patch *patch = &psc.patch[0];
  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *part = particles_base_get_one(&psc.pp, n);
      
    creal u = (part->xi - patch->xb[0]) * dxi;
    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
    int j1 = nint(u);
    int j2 = nint(v);
    int j3 = nint(w);
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
      
    if (psc.ihi[0] - psc.ilo[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.;
    }
    if (psc.ihi[1] - psc.ilo[1] == 1) {
      j2 = 0; gmy = 0.; g0y = 1.; g1y = 0.;
    }
    if (psc.ihi[2] - psc.ilo[2] == 1) {
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

void
c_calc_densities(fields_base_t *pf)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_densities", 1., 0, psc.pp.n_part * 12 * sizeof(creal));
  }
  prof_start(pr);
  do_c_calc_densities(pf, 0, 1, 2);
  prof_stop(pr);

  psc_add_ghosts(pf, 0, 3);
}

// FIXME too much duplication, specialize 2d/1d

static void
do_c_calc_v(fields_base_t *pf)
{
  for (int m = 0; m < 6; m++) {
    fields_base_zero(pf, m);
  }
  
  creal fnqs = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
  creal dxi = 1.f / psc.dx[0];
  creal dyi = 1.f / psc.dx[1];
  creal dzi = 1.f / psc.dx[2];

  struct psc_patch *patch = &psc.patch[0];
  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *part = particles_base_get_one(&psc.pp, n);

    creal u = (part->xi - patch->xb[0]) * dxi;
    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
    int j1 = nint(u);
    int j2 = nint(v);
    int j3 = nint(w);
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

    if (psc.ihi[0] - psc.ilo[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.;
    }
    if (psc.ihi[1] - psc.ilo[1] == 1) {
      j2 = 0; gmy = 0.; g0y = 1.; g1y = 0.;
    }
    if (psc.ihi[2] - psc.ilo[2] == 1) {
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

void
c_calc_v(fields_base_t *pf)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_calc_v", 1., 0, psc.pp.n_part * 12 * sizeof(creal));
  }
  prof_start(pr);
  do_c_calc_v(pf);
  prof_stop(pr);

  psc_add_ghosts(pf, 0, 3);
}

static void
do_c_calc_vv(fields_base_t *pf)
{
  for (int m = 0; m < 6; m++) {
    fields_base_zero(pf, m);
  }
  
  creal fnqs = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
  creal dxi = 1.f / psc.dx[0];
  creal dyi = 1.f / psc.dx[1];
  creal dzi = 1.f / psc.dx[2];

  struct psc_patch *patch = &psc.patch[0];
  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *part = particles_base_get_one(&psc.pp, n);

    creal u = (part->xi - patch->xb[0]) * dxi;
    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
    int j1 = nint(u);
    int j2 = nint(v);
    int j3 = nint(w);
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

    if (psc.ihi[0] - psc.ilo[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.;
    }
    if (psc.ihi[1] - psc.ilo[1] == 1) {
      j2 = 1; gmy = 0.; g0y = 1.; g1y = 0.;
    }
    if (psc.ihi[2] - psc.ilo[2] == 1) {
      j3 = 2; gmz = 0.; g0z = 1.; g1z = 0.;
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

void
c_calc_vv(fields_base_t *pf)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_calc_vv", 1., 0, psc.pp.n_part * 12 * sizeof(creal));
  }
  prof_start(pr);
  do_c_calc_vv(pf);
  prof_stop(pr);

  psc_add_ghosts(pf, 0, 3);
}

struct psc_moment_ops psc_moment_ops_c = {
  .name                   = "c",
  .calc_densities         = c_calc_densities,
  .calc_v                 = c_calc_v,
  .calc_vv                = c_calc_vv,
};
