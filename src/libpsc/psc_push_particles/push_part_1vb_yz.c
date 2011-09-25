
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static inline void
calc_dx1(creal dx1[2], creal x[2], creal dx[2], int off[2])
{
  if (off[1] == 0) {
    dx1[0] = .5 * off[0] - x[0];
    dx1[1] = dx[1] / dx[0] * dx1[0];
  } else {
    dx1[1] = .5 * off[1] - x[1];
    dx1[0] = dx[0] / dx[1] * dx1[1];
  }
}

static inline void
curr_2d_vb_cell(fields_t *pf, int i[2], creal x[2], creal dx[2], creal fnq[2],
		creal dxt[2], int off[2])
{
  F3(JYI, 0,i[0],i[1]  ) += fnq[0] * dx[0] * (.5 - x[1] - .5 * dx[1]);
  F3(JYI, 0,i[0],i[1]+1) += fnq[0] * dx[0] * (.5 + x[1] + .5 * dx[1]);
  F3(JZI, 0,i[0],i[1]  ) += fnq[1] * dx[1] * (.5 - x[0] - .5 * dx[0]);
  F3(JZI, 0,i[0]+1,i[1]) += fnq[1] * dx[1] * (.5 + x[0] + .5 * dx[0]);
  if (dxt) {
    dxt[0] -= dx[0];
    dxt[1] -= dx[1];
    x[0] += dx[0] - off[0];
    x[1] += dx[1] - off[1];
    i[0] += off[0];
    i[1] += off[1];
  }
}

static void
do_push_part_1vb_yz(int p, fields_t *pf, particles_t *pp)
{
  creal dt = ppsc->dt;
  creal yl = .5f * dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqys = ppsc->dx[1] * fnqs / dt;
  creal fnqzs = ppsc->dx[2] * fnqs / dt;
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  fields_zero(pf, JXI);
  fields_zero(pf, JYI);
  fields_zero(pf, JZI);
  
  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // x^n, p^n -> x^(n+.5), p^n

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
    creal xm[3] = { 0.f, v, w };
    int lg2 = fint(v);
    int lg3 = fint(w);
    creal h2 = v - lg2;
    creal h3 = w - lg3;

    creal g0y = 1.f - h2;
    creal g0z = 1.f - h3;
    creal g1y = h2;
    creal g1z = h3;

    // CHARGE DENSITY AT (n+.5)*dt 

    v = (part->yi - patch->xb[1]) * dyi - .5f;
    w = (part->zi - patch->xb[2]) * dzi - .5f;
    int lh2 = fint(v);
    int lh3 = fint(w);
    h2 = v - lh2;
    h3 = w - lh3;
    creal h0y = 1.f - h2;
    creal h0z = 1.f - h3;
    creal h1y = h2;
    creal h1z = h3;

    // FIELD INTERPOLATION

#define INTERPOLATE_FIELD(m, gy, gz)					\
    (gz##0z*(gy##0y*F3(m, 0,l##gy##2  ,l##gz##3  ) +			\
	     gy##1y*F3(m, 0,l##gy##2+1,l##gz##3  )) +			\
     gz##1z*(gy##0y*F3(m, 0,l##gy##2  ,l##gz##3+1) +			\
	     gy##1y*F3(m, 0,l##gy##2+1,l##gz##3+1)))			\
      
    creal exq = INTERPOLATE_FIELD(EX, g, g);
    creal eyq = INTERPOLATE_FIELD(EY, h, g);
    creal ezq = INTERPOLATE_FIELD(EZ, g, h);

    creal hxq = INTERPOLATE_FIELD(HX, h, h);
    creal hyq = INTERPOLATE_FIELD(HY, g, h);
    creal hzq = INTERPOLATE_FIELD(HZ, h, g);

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    creal dq = dqs * part->qni / part->mni;
    creal pxm = part->pxi + dq*exq;
    creal pym = part->pyi + dq*eyq;
    creal pzm = part->pzi + dq*ezq;

    root = dq / creal_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
    creal taux = hxq*root;
    creal tauy = hyq*root;
    creal tauz = hzq*root;

    creal tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    creal pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.f*taux*tauy+2.f*tauz)*pym + 
		(2.f*taux*tauz-2.f*tauy)*pzm)*tau;
    creal pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
		(1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.f*tauy*tauz+2.f*taux)*pzm)*tau;
    creal pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
		(2.f*tauy*tauz-2.f*taux)*pym +
		(1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vxi = part->pxi * root;
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt

    v = (part->yi - patch->xb[1]) * dyi;
    w = (part->zi - patch->xb[2]) * dzi;
    int l2 = fint(v);
    int l3 = fint(w);
    h2 = v - l2;
    h3 = w - l3;

    creal fnqx = vxi * part->qni * part->wni * fnqs;
    F3(JXI, 0,l2  ,l3  ) += (1.f - h2) * (1.f - h3) * fnqx;
    F3(JXI, 0,l2+1,l3  ) += (h2      ) * (1.f - h3) * fnqx;
    F3(JXI, 0,l2  ,l3+1) += (1.f - h2) * (      h3) * fnqx;
    F3(JXI, 0,l2+1,l3+1) += (h2      ) * (      h3) * fnqx;

    // CHARGE DENSITY AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal yi = part->yi + vyi * yl;
    creal zi = part->zi + vzi * zl;

    v = (yi - patch->xb[1]) * dyi;
    w = (zi - patch->xb[2]) * dzi;
    creal xp[3] = { 0.f, v, w };
    int k2 = fint(v);
    int k3 = fint(w);
    h2 = v - k2;
    h3 = w - k3;

    // OUT OF PLANE CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt

    int i[2] = { lg2, lg3 };
    int idiff[2] = { k2 - lg2, k3 - lg3 };
    creal dx[2] = { xp[1] - xm[1], xp[2] - xm[2] };
    creal x[2] = { xm[1] - (i[0] + .5), xm[2] - (i[1] + .5) };

    creal dx1[2];
    int off[2];
    int first_dir, second_dir = -1;
    // FIXME, make sure we never div-by-zero?
    if (idiff[0] == 0 && idiff[1] == 0) {
      first_dir = -1;
    } else if (idiff[0] == 0) {
      first_dir = 1;
    } else if (idiff[1] == 0) {
      first_dir = 0;
    } else {
      dx1[0] = .5 * idiff[0] - x[0];
      dx1[1] = dx[1] / dx[0] * dx1[0];
      if (creal_abs(x[1] + dx1[1]) > .5f) {
	first_dir = 1;
      } else {
	first_dir = 0;
      }
      second_dir = 1 - first_dir;
    }

    creal fnq[2] = { part->qni * part->wni * fnqys,
		     part->qni * part->wni * fnqzs };

    if (first_dir >= 0) {
      off[1-first_dir] = 0;
      off[first_dir] = idiff[first_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(pf, i, x, dx1, fnq, dx, off);
    }

    if (second_dir >= 0) {
      off[first_dir] = 0;
      off[second_dir] = idiff[second_dir];
      calc_dx1(dx1, x, dx, off);
      curr_2d_vb_cell(pf, i, x, dx1, fnq, dx, off);
    }
    
    curr_2d_vb_cell(pf, i, x, dx, fnq, NULL, NULL);
  }
}

void
psc_push_particles_1vb_push_yz(struct psc_push_particles *push,
			       mparticles_base_t *particles_base,
			       mfields_base_t *flds_base)
{
  mfields_t flds;
  mparticles_t particles;
  psc_mfields_get_from(&flds, EX, EX + 6, flds_base);
  psc_mparticles_get_from(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("1vb_push_yz", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_push_part_1vb_yz(p, &flds.f[p], &particles.p[p]);
  }
  prof_stop(pr);

  psc_mfields_put_to(&flds, JXI, JXI + 3, flds_base);
  psc_mparticles_put_to(&particles, particles_base);
}

