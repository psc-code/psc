
#include "psc_moments_private.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include "psc_bnd.h"
#include <mrc_profile.h>
#include <math.h>

// ======================================================================

typedef fields_c_real_t creal;

static void
do_1st_calc_v(int p, fields_t *pf, particles_t *pp)
{
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal dxi = 1.f / ppsc->dx[0];
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    creal u = (part->xi - patch->xb[0]) * dxi - .5;
    creal v = (part->yi - patch->xb[1]) * dyi - .5;
    creal w = (part->zi - patch->xb[2]) * dzi - .5;
    int j1 = particle_real_fint(u);
    int j2 = particle_real_fint(v);
    int j3 = particle_real_fint(w);
    creal h1 = u-j1;
    creal h2 = v-j2;
    creal h3 = w-j3;
      
    creal g0x=1.f - h1;
    creal g0y=1.f - h2;
    creal g0z=1.f - h3;
    creal g1x=h1;
    creal g1y=h2;
    creal g1z=h3;
      
    if (ppsc->domain.gdims[0] == 1) {
      j1 = 0; g0x = 1.; g1x = 0.;
    }
    if (ppsc->domain.gdims[1] == 1) {
      j2 = 0; g0y = 1.; g1y = 0.;
    }
    if (ppsc->domain.gdims[2] == 1) {
      j3 = 0; g0z = 1.; g1z = 0.;
    }

    creal pxi = part->pxi;
    creal pyi = part->pyi;
    creal pzi = part->pzi;
    creal root = 1.0/sqrt(1.0+pxi*pxi+pyi*pyi+pzi*pzi);
    creal vv[3] = { pxi*root, pyi*root, pzi*root };
    creal fnq = part->wni * fnqs;
    int mm;
    if (part->qni < 0.) {
      mm = 0;
    } else if (part->qni > 0.) {
      mm = 3;
    } else {
      assert(0);
    }
    for (int m = 0; m < 3; m++) {
      F3(pf, mm+m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z * vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z * vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z * vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z * vv[m];
      F3(pf, mm+m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z * vv[m];
      F3(pf, mm+m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z * vv[m];
      F3(pf, mm+m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z * vv[m];
      F3(pf, mm+m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z * vv[m];
    }
  }
}

static void
do_1st_calc_densities(int p, fields_t *pf, particles_t *pp,
		    int m_NE, int m_NI, int m_NN)
{
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal dxi = 1.f / ppsc->dx[0];
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    creal u = (part->xi - patch->xb[0]) * dxi - .5;
    creal v = (part->yi - patch->xb[1]) * dyi - .5;
    creal w = (part->zi - patch->xb[2]) * dzi - .5;
    int j1 = particle_real_fint(u);
    int j2 = particle_real_fint(v);
    int j3 = particle_real_fint(w);
    creal h1 = u-j1;
    creal h2 = v-j2;
    creal h3 = w-j3;
      
    creal g0x=1.f - h1;
    creal g0y=1.f - h2;
    creal g0z=1.f - h3;
    creal g1x=h1;
    creal g1y=h2;
    creal g1z=h3;
      
    if (ppsc->domain.gdims[0] == 1) {
      j1 = 0; g0x = 1.; g1x = 0.;
    }
    if (ppsc->domain.gdims[1] == 1) {
      j2 = 0; g0y = 1.; g1y = 0.;
    }
    if (ppsc->domain.gdims[2] == 1) {
      j3 = 0; g0z = 1.; g1z = 0.;
    }

    assert(j1 >= -1 && j1 < patch->ldims[0]);
    assert(j2 >= -1 && j2 < patch->ldims[1]);
    assert(j3 >= -1 && j3 < patch->ldims[2]);
      
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
    F3(pf, m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z;
    F3(pf, m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z;
    F3(pf, m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z;
    F3(pf, m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z;
    F3(pf, m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z;
    F3(pf, m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z;
    F3(pf, m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z;
    F3(pf, m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z;
  }
}

static void
add_ghosts_reflecting_lo(mfields_c_t *res, int p, int d, int mb, int me)
{
  fields_t *pf = psc_mfields_get_patch(res, p);
  struct psc_patch *patch = ppsc->patch + p;

  if (d == 1) {
    for (int iz = 0; iz < patch->ldims[2]; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = 0; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy-1,iz);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_reflecting_hi(mfields_c_t *res, int p, int d, int mb, int me)
{
  fields_t *pf = psc_mfields_get_patch(res, p);
  struct psc_patch *patch = ppsc->patch + p;

  if (d == 1) {
    for (int iz = 0; iz < patch->ldims[2]; iz++) {
      for (int ix = 0; ix < patch->ldims[0]; ix++) {
	int iy = patch->ldims[1] - 1; {
	  for (int m = mb; m < me; m++) {
	    F3(pf, m, ix,iy,iz) += F3(pf, m, ix,iy+1,iz);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

static void
add_ghosts_boundary(mfields_c_t *res, int mb, int me)
{
  psc_foreach_patch(ppsc, p) {
    // lo
    for (int d = 0; d < 3; d++) {
      if (ppsc->patch[p].off[d] == 0) {
	if (ppsc->domain.bnd_part[d] == BND_PART_REFLECTING) {
	  add_ghosts_reflecting_lo(res, p, d, mb, me);
	}
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      if (ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d]) {
	if (ppsc->domain.bnd_fld_hi[d] == BND_PART_REFLECTING) {
	  add_ghosts_reflecting_hi(res, p, d, mb, me);
	}
      }
    }
  }
}

static void
psc_moments_1st_cc_calc_densities(struct psc_moments *moments, mfields_base_t *flds,
				  mparticles_base_t *particles_base, mfields_c_t *res)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_densities", 1., 0, 0);
  }

  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  prof_start(pr);
  psc_mfields_zero(res, 0);
  psc_mfields_zero(res, 1);
  psc_mfields_zero(res, 2);
  
  psc_foreach_patch(ppsc, p) {
    do_1st_calc_densities(p, psc_mfields_get_patch_c(res, p),
			psc_mparticles_get_patch(particles, p), 0, 1, 2);
  }
  prof_stop(pr);

  psc_mparticles_put_cf(particles, particles_base); // FIXME, don't need copy-back

  psc_bnd_add_ghosts(moments->bnd, res, 0, 3);
  add_ghosts_boundary(res, 0, 3);
}

static void
psc_moments_1st_cc_calc_v(struct psc_moments *moments, mfields_base_t *flds,
			  mparticles_base_t *particles_base, mfields_c_t *res)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_moments_v", 1., 0, 0);
  }

  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  prof_start(pr);
  for (int m = 0; m < 6; m++) {
    psc_mfields_zero(res, m); // FIXME, pass range
  }
  
  psc_foreach_patch(ppsc, p) {
    do_1st_calc_v(p, psc_mfields_get_patch_c(res, p),
		  psc_mparticles_get_patch(particles, p));
  }
  prof_stop(pr);

  psc_mparticles_put_cf(particles, particles_base); // FIXME, don't need copy-back

  psc_bnd_add_ghosts(moments->bnd, res, 0, 6);
  // FIXME, this is probably only right for densities, but needs sign
  // adjustments for velocities
  add_ghosts_boundary(res, 0, 6);
}

// ======================================================================
// psc_moments: subclass "1st_cc"

struct psc_moments_ops psc_moments_1st_cc_ops = {
  .name                  = "1st_cc",
  .calc_densities        = psc_moments_1st_cc_calc_densities,
  .calc_v                = psc_moments_1st_cc_calc_v,
};
