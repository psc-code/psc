
#include "psc_generic_c.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_genc_push_part_y(int p, fields_t *pf, struct psc_particles *pp)
{
#define S0Y(off) s0y[off+2]
#define S1Y(off) s1y[off+2]

  creal s0y[5] = {}, s1y[5];

  creal dt = ppsc->dt;
  creal yl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqys = ppsc->patch[p].dx[1] * fnqs / dt;
  creal dyi = 1.f / ppsc->patch[p].dx[1];

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // x^n, p^n -> x^(n+.5), p^n

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->yi += vyi * yl;

    creal v = part->yi * dyi;
    int j2 = particle_real_nint(v);
    creal h2 = j2-v;

    creal gmy=.5f*(.5f+h2)*(.5f+h2);
    creal g0y=.75f-h2*h2;
    creal g1y=.5f*(.5f-h2)*(.5f-h2);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0Y(-1) = .5f*(1.5f-creal_abs(h2-1.f))*(1.5f-creal_abs(h2-1.f));
    S0Y(+0) = .75f-h2*h2;
    S0Y(+1) = .5f*(1.5f-creal_abs(h2+1.f))*(1.5f-creal_abs(h2+1.f));

    v = part->yi * dyi;
    int l2=particle_real_nint(v);
    h2=l2-v;
    creal hmy=.5f*(.5f+h2)*(.5f+h2);
    creal h0y=.75f-h2*h2;
    creal h1y=.5f*(.5f-h2)*(.5f-h2);

    // FIELD INTERPOLATION

    creal exq = 
       gmy*F3(pf, EX, 0,j2-1,0)
      +g0y*F3(pf, EX, 0,j2  ,0)
      +g1y*F3(pf, EX, 0,j2+1,0);
    creal eyq =
       gmy*F3(pf, EY, 0,l2-1,0)
      +g0y*F3(pf, EY, 0,l2  ,0)
      +g1y*F3(pf, EY, 0,l2+1,0);
    creal ezq =
       hmy*F3(pf, EZ, 0,j2-1,0)
      +h0y*F3(pf, EZ, 0,j2  ,0)
      +h1y*F3(pf, EZ, 0,j2+1,0);

    creal hxq =
       hmy*F3(pf, HX, 0,l2-1,0)
      +h0y*F3(pf, HX, 0,l2  ,0)
      +h1y*F3(pf, HX, 0,l2+1,0);
    creal hyq =
       hmy*F3(pf, HY, 0,j2-1,0)
      +h0y*F3(pf, HY, 0,j2  ,0)
      +h1y*F3(pf, HY, 0,j2+1,0);
    creal hzq =
       gmy*F3(pf, HZ, 0,l2-1,0)
      +g0y*F3(pf, HZ, 0,l2  ,0)
      +g1y*F3(pf, HZ, 0,l2+1,0);

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

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal yi = part->yi + vyi * yl;

    v = yi * dyi;
    int k2 = particle_real_nint(v);
    h2 = k2 - v;

    for (int i = -2; i <= 2; i++) {
      S1Y(i) = 0.f;
    }

    S1Y(k2-j2-1) = .5f*(1.5f-creal_abs(h2-1.f))*(1.5f-creal_abs(h2-1.f));
    S1Y(k2-j2+0) = .75f-h2*h2;
    S1Y(k2-j2+1) = .5f*(1.5f-creal_abs(h2+1.f))*(1.5f-creal_abs(h2+1.f));

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1Y(i) -= S0Y(i);
    }

    int l2min, l2max;
    
    if (k2 == j2) {
      l2min = -1; l2max = +1;
    } else if (k2 == j2 - 1) {
      l2min = -2; l2max = +1;
    } else { // (k2 == j2 + 1)
      l2min = -1; l2max = +2;
    }

    creal jxh;
    creal jyh;
    creal jzh;

    creal fnqx = vxi * part->qni * part->wni * fnqs;
    creal fnqy = part->qni * part->wni * fnqys;
    creal fnqz = vzi * part->qni * part->wni * fnqs;
    jyh = 0.f;
    for (int l2 = l2min; l2 <= l2max; l2++) {
      creal wx = S0Y(l2) + .5f * S1Y(l2);
      creal wy = S1Y(l2);
      creal wz = S0Y(l2) + .5f * S1Y(l2);
      
      jxh = fnqx*wx;
      jyh -= fnqy*wy;
      jzh = fnqz*wz;

      F3(pf, JXI, 0,j2+l2,0) += jxh;
      F3(pf, JYI, 0,j2+l2,0) += jyh;
      F3(pf, JZI, 0,j2+l2,0) += jzh;
    }
  }
}

void
psc_push_particles_generic_c_push_a_y(struct psc_push_particles *push,
				      struct psc_particles *prts,
				      struct psc_fields *flds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_y", 1., 0, 0);
  }
  
  prof_start(pr);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  do_genc_push_part_y(prts->p, flds, prts);
  prof_stop(pr);
}

