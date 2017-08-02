
#include "psc_generic_c.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_genc_push_part_z(int p, fields_t *pf, particle_range_t prts)
{
#define S0Z(off) s0z[off+2]
#define S1Z(off) s1z[off+2]

  creal s0z[5] = {}, s1z[5];

  creal dt = ppsc->dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;
  creal dxi = 1.f / ppsc->patch[p].dx[0];
  creal dyi = 1.f / ppsc->patch[p].dx[1];
  creal dzi = 1.f / ppsc->patch[p].dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);

    // x^n, p^n -> x^(n+.5), p^n

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->zi += vzi * zl;

    creal u = part->xi * dxi;
    creal v = part->yi * dyi;
    creal w = part->zi * dzi;
    int j1 = particle_real_nint(u);
    int j2 = particle_real_nint(v);
    int j3 = particle_real_nint(w);
    creal h3 = j3-w;

    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0z=.75f-h3*h3;
    creal g1z=.5f*(.5f-h3)*(.5f-h3);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0Z(-1) = .5f*(1.5f-creal_abs(h3-1.f))*(1.5f-creal_abs(h3-1.f));
    S0Z(+0) = .75f-creal_abs(h3)*creal_abs(h3);
    S0Z(+1) = .5f*(1.5f-creal_abs(h3+1.f))*(1.5f-creal_abs(h3+1.f));

    u = part->xi * dxi - .5f;
    v = part->yi * dyi;
    w = part->zi * dzi - .5f;
    int l1=particle_real_nint(u);
    int l2=particle_real_nint(v);
    int l3=particle_real_nint(w);
    h3=l3-w;
    creal hmz=.5f*(.5f+h3)*(.5f+h3);
    creal h0z=.75f-h3*h3;
    creal h1z=.5f*(.5f-h3)*(.5f-h3);

    // FIELD INTERPOLATION

    creal exq = 
       gmz*F3(pf, EX, l1,j2,j3-1)
      +g0z*F3(pf, EX, l1,j2,j3)
      +g1z*F3(pf, EX, l1,j2,j3+1);
    creal eyq =
       gmz*F3(pf, EY, j1,l2,j3-1)
      +g0z*F3(pf, EY, j1,l2,j3)
      +g1z*F3(pf, EY, j1,l2,j3+1);
    creal ezq =
       hmz*F3(pf, EZ, j1,j2,l3-1)
      +h0z*F3(pf, EZ, j1,j2,l3)
      +h1z*F3(pf, EZ, j1,j2,l3+1);

    creal hxq =
       hmz*F3(pf, HX, j1,l2,l3-1)
      +h0z*F3(pf, HX, j1,l2,l3)
      +h1z*F3(pf, HX, j1,l2,l3+1);
    creal hyq =
       hmz*F3(pf, HY, l1,j2,l3-1)
      +h0z*F3(pf, HY, l1,j2,l3)
      +h1z*F3(pf, HY, l1,j2,l3+1);
    creal hzq =
       gmz*F3(pf, HZ, l1,l2,j3-1)
      +g0z*F3(pf, HZ, l1,l2,j3)
      +g1z*F3(pf, HZ, l1,l2,j3+1);

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

    w = zi * dzi;
    int k3 = particle_real_nint(w);
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
      
      F3(pf, JXI, j1+l1,j2,j3+l3) += jxh;
      F3(pf, JYI, j1+l1,j2,j3+l3) += jyh;
      F3(pf, JZI, j1+l1,j2,j3+l3) += jzh;
    }
  }
}

void
psc_push_particles_generic_c_push_mprts_z(struct psc_push_particles *push,
					  struct psc_mparticles *mprts,
					  struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_z", 1., 0, 0);
  }
  
  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);
    psc_fields_zero_range(flds, JXI, JXI + 3);
    do_genc_push_part_z(p, flds, prts);
  }
  prof_stop(pr);
}

