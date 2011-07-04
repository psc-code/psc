
#include "psc_generic_c.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_genc_push_part_xz(int p, fields_t *pf, particles_t *pp)
{
#define S0X(off) s0x[off+2]
#define S0Z(off) s0z[off+2]
#define S1X(off) s1x[off+2]
#define S1Z(off) s1z[off+2]

  creal s0x[5] = {}, s0z[5] = {}, s1x[5], s1z[5];

  creal dt = ppsc->dt;
  creal xl = .5f * dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  creal fnqxs = ppsc->dx[0] * fnqs / dt;
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
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->xi += vxi * xl;
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

    creal gmx=.5f*(.5f+h1)*(.5f+h1);
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0x=.75f-h1*h1;
    creal g0z=.75f-h3*h3;
    creal g1x=.5f*(.5f-h1)*(.5f-h1);
    creal g1z=.5f*(.5f-h3)*(.5f-h3);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0X(-1) = .5f*(1.5f-creal_abs(h1-1.f))*(1.5f-creal_abs(h1-1.f));
    S0X(+0) = .75f-creal_abs(h1)*creal_abs(h1);
    S0X(+1) = .5f*(1.5f-creal_abs(h1+1.f))*(1.5f-creal_abs(h1+1.f));
    S0Z(-1) = .5f*(1.5f-creal_abs(h3-1.f))*(1.5f-creal_abs(h3-1.f));
    S0Z(+0) = .75f-creal_abs(h3)*creal_abs(h3);
    S0Z(+1) = .5f*(1.5f-creal_abs(h3+1.f))*(1.5f-creal_abs(h3+1.f));

    u = (part->xi - patch->xb[0]) * dxi - .5f;
    v = (part->yi - patch->xb[1]) * dyi;
    w = (part->zi - patch->xb[2]) * dzi - .5f;
    int l1 = nint(u);
    int l2 = nint(v);
    int l3 = nint(w);
    h1=l1-u;
    h2=l2-v;
    h3=l3-w;
    creal hmx=.5f*(.5f+h1)*(.5f+h1);
    creal hmz=.5f*(.5f+h3)*(.5f+h3);
    creal h0x=.75f-h1*h1;
    creal h0z=.75f-h3*h3;
    creal h1x=.5f*(.5f-h1)*(.5f-h1);
    creal h1z=.5f*(.5f-h3)*(.5f-h3);

    // FIELD INTERPOLATION

    creal exq = (gmz*(hmx*F3(EX, l1-1,j2,j3-1) +
		      h0x*F3(EX, l1  ,j2,j3-1) +
		      h1x*F3(EX, l1+1,j2,j3-1)) +
		 g0z*(hmx*F3(EX, l1-1,j2,j3  ) +
		      h0x*F3(EX, l1  ,j2,j3  ) +
		      h1x*F3(EX, l1+1,j2,j3  )) +
		 g1z*(hmx*F3(EX, l1-1,j2,j3+1) +
		      h0x*F3(EX, l1  ,j2,j3+1) +
		      h1x*F3(EX, l1+1,j2,j3+1)));

    creal eyq = (gmz*(gmx*F3(EY, j1-1,l2,j3-1) +
		      g0x*F3(EY, j1  ,l2,j3-1) +
		      g1x*F3(EY, j1+1,l2,j3-1)) +
		 g0z*(gmx*F3(EY, j1-1,l2,j3  ) +
		      g0x*F3(EY, j1  ,l2,j3  ) +
		      g1x*F3(EY, j1+1,l2,j3  )) +
		 g1z*(gmx*F3(EY, j1-1,l2,j3+1) +
		      g0x*F3(EY, j1  ,l2,j3+1) +
		      g1x*F3(EY, j1+1,l2,j3+1)));

    creal ezq = (hmz*(gmx*F3(EZ, j1-1,j2,l3-1) +
		      g0x*F3(EZ, j1  ,j2,l3-1) +
		      g1x*F3(EZ, j1+1,j2,l3-1)) +
		 h0z*(gmx*F3(EZ, j1-1,j2,l3  ) +
		      g0x*F3(EZ, j1  ,j2,l3  ) +
		      g1x*F3(EZ, j1+1,j2,l3  )) +
		 h1z*(gmx*F3(EZ, j1-1,j2,l3+1) +
		      g0x*F3(EZ, j1  ,j2,l3+1) +
		      g1x*F3(EZ, j1+1,j2,l3+1)));

    creal hxq = (hmz*(gmx*F3(HX, j1-1,l2,l3-1) +
		      g0x*F3(HX, j1  ,l2,l3-1) +
		      g1x*F3(HX, j1+1,l2,l3-1)) +
		 h0z*(gmx*F3(HX, j1-1,l2,l3  ) +
		      g0x*F3(HX, j1  ,l2,l3  ) +
		      g1x*F3(HX, j1+1,l2,l3  )) +
		 h1z*(gmx*F3(HX, j1-1,l2,l3+1) +
		      g0x*F3(HX, j1  ,l2,l3+1) +
		      g1x*F3(HX, j1+1,l2,l3+1)));

    creal hyq = (hmz*(hmx*F3(HY, l1-1,j2,l3-1) +
		      h0x*F3(HY, l1  ,j2,l3-1) +
		      h1x*F3(HY, l1+1,j2,l3-1)) +
		 h0z*(hmx*F3(HY, l1-1,j2,l3  ) +
		      h0x*F3(HY, l1  ,j2,l3  ) +
		      h1x*F3(HY, l1+1,j2,l3  )) +
		 h1z*(hmx*F3(HY, l1-1,j2,l3+1) +
		      h0x*F3(HY, l1  ,j2,l3+1) +
		      h1x*F3(HY, l1+1,j2,l3+1)));

    creal hzq = (gmz*(hmx*F3(HZ, l1-1,l2,j3-1) +
		      h0x*F3(HZ, l1  ,l2,j3-1) +
		      h1x*F3(HZ, l1+1,l2,j3-1)) +
		 g0z*(hmx*F3(HZ, l1-1,l2,j3  ) +
		      h0x*F3(HZ, l1  ,l2,j3  ) +
		      h1x*F3(HZ, l1+1,l2,j3  )) +
		 g1z*(hmx*F3(HZ, l1-1,l2,j3+1) +
		      h0x*F3(HZ, l1  ,l2,j3+1) +
		      h1x*F3(HZ, l1+1,l2,j3+1)));

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

    part->xi += vxi * xl;
    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    creal xi = part->xi + vxi * xl;
    creal zi = part->zi + vzi * zl;

    u = (xi - patch->xb[0]) * dxi;
    w = (zi - patch->xb[2]) * dzi;
    int k1 = nint(u);
    int k3 = nint(w);
    h1 = k1 - u;
    h3 = k3 - w;

    for (int i = -2; i <= 2; i++) {
      S1X(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1X(k1-j1-1) = .5f*(1.5f-creal_abs(h1-1.f))*(1.5f-creal_abs(h1-1.f));
    S1X(k1-j1+0) = .75f-creal_abs(h1)*creal_abs(h1);
    S1X(k1-j1+1) = .5f*(1.5f-creal_abs(h1+1.f))*(1.5f-creal_abs(h1+1.f));
    S1Z(k3-j3-1) = .5f*(1.5f-creal_abs(h3-1.f))*(1.5f-creal_abs(h3-1.f));
    S1Z(k3-j3+0) = .75f-creal_abs(h3)*creal_abs(h3);
    S1Z(k3-j3+1) = .5f*(1.5f-creal_abs(h3+1.f))*(1.5f-creal_abs(h3+1.f));

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1X(i) -= S0X(i);
      S1Z(i) -= S0Z(i);
    }

    int l1min, l3min, l1max, l3max;
    
    if (k1 == j1) {
      l1min = -1; l1max = +1;
    } else if (k1 == j1 - 1) {
      l1min = -2; l1max = +1;
    } else { // (k1 == j1 + 1)
      l1min = -1; l1max = +2;
    }

    if (k3 == j3) {
      l3min = -1; l3max = +1;
    } else if (k3 == j3 - 1) {
      l3min = -2; l3max = +1;
    } else { // (k3 == j3 + 1)
      l3min = -1; l3max = +2;
    }

    creal jxh;
    creal jyh;
    creal jzh[5];

#define JZH(i) jzh[i+2]

    creal fnqx = part->qni * part->wni * fnqxs;
    creal fnqy = vyi * part->qni * part->wni * fnqs;
    creal fnqz = part->qni * part->wni * fnqzs;
    for (int l1 = l1min; l1 <= l1max; l1++) {
      JZH(l1) = 0.f;
    }
    for (int l3 = l3min; l3 <= l3max; l3++) {
      jxh = 0.f;
      for (int l1 = l1min; l1 <= l1max; l1++) {
	creal wx = S1X(l1) * (S0Z(l3) + .5f*S1Z(l3));
	creal wy = S0X(l1) * S0Z(l3)
	  + .5f * S1X(l1) * S0Z(l3)
	  + .5f * S0X(l1) * S1Z(l3)
	  + (1.f/3.f) * S1X(l1) * S1Z(l3);
	creal wz = S1Z(l3) * (S0X(l1) + .5f*S1X(l1));

	jxh -= fnqx*wx;
	jyh = fnqy*wy;
	JZH(l1) -= fnqz*wz;

	F3(JXI, j1+l1,j2,j3+l3) += jxh;
	F3(JYI, j1+l1,j2,j3+l3) += jyh;
	F3(JZI, j1+l1,j2,j3+l3) += JZH(l1);
      }
    }
  }
}

void
psc_push_particles_generic_c_push_xz(struct psc_push_particles *push,
				     mparticles_base_t *particles_base,
				     mfields_base_t *flds_base)
{
  mfields_t flds;
  mparticles_t particles;
  psc_mfields_get(&flds, EX, EX + 6, flds_base);
  psc_mparticles_get_from(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_xz", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_genc_push_part_xz(p, &flds.f[p], &particles.p[p]);
  }
  prof_stop(pr);

  psc_mfields_put(&flds, JXI, JXI + 3, flds_base);
  psc_mparticles_put_to(&particles, particles_base);
}

