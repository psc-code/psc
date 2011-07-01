
#include "psc_generic_c.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_genc_push_part_z(int p, fields_t *pf, particles_t *pp)
{
#define S0Z(off) s0z[off+2]
#define S1Z(off) s1z[off+2]

  creal s0z[5] = {}, s1z[5];

  creal dt = ppsc->dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqzs = ppsc->dx[2] * fnqs / dt;
  creal dxi = 1.f / ppsc->dx[0];
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
    creal vyi = part->pxi * root;
    creal vzi = part->pzi * root;

    part->zi += vzi * zl;

    creal u = (part->xi - patch->xb[0]) * dxi;
    creal v = (part->yi - patch->xb[1]) * dyi;
    creal w = (part->zi - patch->xb[2]) * dzi;
    int j1 = nint(u);
    int j2 = nint(v);
    int j3 = nint(w);
    creal h1 = j1-u;
    creal h2 = j2-v;
    creal h3 = j3-w;

    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0z=.75f-h3*h3;
    creal g1z=.5f*(.5f-h3)*(.5f-h3);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0Z(-1) = .5f*(1.5f-creal_abs(h3-1.f))*(1.5f-creal_abs(h3-1.f));
    S0Z(+0) = .75f-creal_abs(h3)*creal_abs(h3);
    S0Z(+1) = .5f*(1.5f-creal_abs(h3+1.f))*(1.5f-creal_abs(h3+1.f));

    u = (part->xi - patch->xb[0]) * dxi - .5f;
    v = (part->yi - patch->xb[1]) * dyi;
    w = (part->zi - patch->xb[2]) * dzi - .5f;
    int l1=nint(u);
    int l2=nint(v);
    int l3=nint(w);
    h1=l1-u;
    h2=l2-v;
    h3=l3-w;
    creal hmz=.5f*(.5f+h3)*(.5f+h3);
    creal h0z=.75f-h3*h3;
    creal h1z=.5f*(.5f-h3)*(.5f-h3);

    // FIELD INTERPOLATION

    creal exq = 
       gmz*F3(EX, l1,j2,j3-1)
      +g0z*F3(EX, l1,j2,j3)
      +g1z*F3(EX, l1,j2,j3+1);
    creal eyq =
       gmz*F3(EY, j1,l2,j3-1)
      +g0z*F3(EY, j1,l2,j3)
      +g1z*F3(EY, j1,l2,j3+1);
    creal ezq =
       hmz*F3(EZ, j1,j2,l3-1)
      +h0z*F3(EZ, j1,j2,l3)
      +h1z*F3(EZ, j1,j2,l3+1);

    creal hxq =
       hmz*F3(HX, j1,l2,l3-1)
      +h0z*F3(HX, j1,l2,l3)
      +h1z*F3(HX, j1,l2,l3+1);
    creal hyq =
       hmz*F3(HY, l1,j2,l3-1)
      +h0z*F3(HY, l1,j2,l3)
      +h1z*F3(HY, l1,j2,l3+1);
    creal hzq =
       gmz*F3(HZ, l1,l2,j3-1)
      +g0z*F3(HZ, l1,l2,j3)
      +g1z*F3(HZ, l1,l2,j3+1);

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

    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal zi = part->zi + vzi * zl;

    w = (zi - patch->xb[2]) * dzi;
    int k3 = nint(w);
    h3 = k3 - w;

    for (int i = -2; i <= 2; i++) {
      S1Z(i) = 0.f;
    }

    S1Z(k3-j3-1) = .5f*(1.5f-creal_abs(h3-1.f))*(1.5f-creal_abs(h3-1.f));
    S1Z(k3-j3+0) = .75f-creal_abs(h3)*creal_abs(h3);
    S1Z(k3-j3+1) = .5f*(1.5f-creal_abs(h3+1.f))*(1.5f-creal_abs(h3+1.f));

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1Z(i) -= S0Z(i);
    }

    int l3min, l3max;
    
    if (k3 == j3) {
      l3min = -1; l3max = +1;
    } else if (k3 == j3 - 1) {
      l3min = -2; l3max = +1;
    } else { // (k3 == j3 + 1)
      l3min = -1; l3max = +2;
    }

    creal jxh;
    creal jyh;
    creal jzh;

    creal fnqx = vxi * part->qni * part->wni * fnqs;
    creal fnqy = vyi * part->qni * part->wni * fnqs;
    creal fnqz = part->qni * part->wni * fnqzs;
    jzh = 0.f;
    for (int l3 = l3min; l3 <= l3max; l3++) {
      creal wx = S0Z(l3) + .5f * S1Z(l3);
      creal wy = S0Z(l3) + .5f * S1Z(l3);
      creal wz = S1Z(l3);
      
      jxh = fnqx*wx;
      jyh = fnqy*wy;
      jzh -= fnqz*wz;
      
      F3(JXI, j1+l1,j2,j3+l3) += jxh;
      F3(JYI, j1+l1,j2,j3+l3) += jyh;
      F3(JZI, j1+l1,j2,j3+l3) += jzh;
    }
  }
}

void
psc_push_particles_generic_c_push_z(struct psc_push_particles *push,
				    mparticles_base_t *particles_base,
				    mfields_base_t *flds_base)
{
  mfields_t flds;
  mparticles_t particles;
  fields_get(&flds, EX, EX + 6, flds_base);
  psc_mparticles_get_from(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_z", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_genc_push_part_z(p, &flds.f[p], &particles.p[p]);
  }
  prof_stop(pr);

  fields_put(&flds, JXI, JXI + 3, flds_base);
  psc_mparticles_put_to(&particles, particles_base);
}

