
#include "psc_1p5_c.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
do_1p5_push_part_yz(int p, fields_t *pf, particle_range_t prts)
{
#define S0Y(off) s0y[off+2]
#define S0Z(off) s0z[off+2]
#define S1Y(off) s1y[off+2]
#define S1Z(off) s1z[off+2]

  particle_real_t s0y[5] = {}, s0z[5] = {}, s1y[5], s1z[5];

  particle_real_t dt = ppsc->dt;
  particle_real_t yl = .5f * dt;
  particle_real_t zl = .5f * dt;
  particle_real_t dqs = .5f * ppsc->coeff.eta * dt;
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t fnqys = ppsc->patch[p].dx[1] * fnqs / dt;
  particle_real_t fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;
  particle_real_t dxi = 1.f / ppsc->patch[p].dx[0];
  particle_real_t dyi = 1.f / ppsc->patch[p].dx[1];
  particle_real_t dzi = 1.f / ppsc->patch[p].dx[2];

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);

    // x^n, p^n -> x^(n+.5), p^n

    particle_real_t root = 1.f / particle_real_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    particle_real_t vxi = part->pxi * root;
    particle_real_t vyi = part->pyi * root;
    particle_real_t vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    particle_real_t u = part->xi * dxi;
    particle_real_t v = part->yi * dyi;
    particle_real_t w = part->zi * dzi;
    int j1 = particle_real_nint(u);
    int j2 = particle_real_nint(v);
    int j3 = particle_real_nint(w);
    particle_real_t h2 = j2-v;
    particle_real_t h3 = j3-w;

    particle_real_t gmy=.5f*(.5f+h2)*(.5f+h2);
    particle_real_t gmz=.5f*(.5f+h3)*(.5f+h3);
    particle_real_t g0y=.75f-h2*h2;
    particle_real_t g0z=.75f-h3*h3;
    particle_real_t g1y=.5f*(.5f-h2)*(.5f-h2);
    particle_real_t g1z=.5f*(.5f-h3)*(.5f-h3);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    S0Y(-1) = .5f*(1.5f-particle_real_abs(h2-1.f))*(1.5f-particle_real_abs(h2-1.f));
    S0Y(+0) = .75f-particle_real_abs(h2)*particle_real_abs(h2);
    S0Y(+1) = .5f*(1.5f-particle_real_abs(h2+1.f))*(1.5f-particle_real_abs(h2+1.f));
    S0Z(-1) = .5f*(1.5f-particle_real_abs(h3-1.f))*(1.5f-particle_real_abs(h3-1.f));
    S0Z(+0) = .75f-particle_real_abs(h3)*particle_real_abs(h3);
    S0Z(+1) = .5f*(1.5f-particle_real_abs(h3+1.f))*(1.5f-particle_real_abs(h3+1.f));

    u = part->xi * dxi;
    v = part->yi * dyi - .5f;
    w = part->zi * dzi - .5f;
    int l1=particle_real_nint(u);
    int l2=particle_real_nint(v);
    int l3=particle_real_nint(w);
    h2=l2-v;
    h3=l3-w;

    // 1+1/2 method from sokolov paper
    particle_real_t hmy;
    particle_real_t h0y;
    particle_real_t h1y;
    particle_real_t hmz;
    particle_real_t h0z;
    particle_real_t h1z;

    if (l2 >= v)  // h2 >= 0
      {
	hmy= h2;
	h0y= 1.f - h2;
	h1y= 0;
      }
    else  // h2 < 0
      {
	hmy= 0;
	h0y= 1.f + h2;
	h1y= -h2;
      }

    if (l3 >= w) // h3 >= 0
      {
	hmz= h3;
	h0z= 1.f - h3;
	h1z= 0;
      }
    else  // h3 < 0
      {
	hmz= 0;
	h0z= 1.f + h3;
	h1z= -h3;
      }

    // FIELD INTERPOLATION

    particle_real_t exq=gmz*(gmy*F3(pf, EX, l1,j2-1,j3-1)
		  +g0y*F3(pf, EX, l1,j2  ,j3-1)
		  +g1y*F3(pf, EX, l1,j2+1,j3-1))
             +g0z*(gmy*F3(pf, EX, l1,j2-1,j3  )
    	          +g0y*F3(pf, EX, l1,j2  ,j3  )
	          +g1y*F3(pf, EX, l1,j2+1,j3  ))
             +g1z*(gmy*F3(pf, EX, l1,j2-1,j3+1)
	          +g0y*F3(pf, EX, l1,j2  ,j3+1)
	          +g1y*F3(pf, EX, l1,j2+1,j3+1));

    particle_real_t eyq=gmz*(hmy*F3(pf, EY, j1,l2-1,j3-1)
                  +h0y*F3(pf, EY, j1,l2  ,j3-1)
                  +h1y*F3(pf, EY, j1,l2+1,j3-1))
             +g0z*(hmy*F3(pf, EY, j1,l2-1,j3  )
                  +h0y*F3(pf, EY, j1,l2  ,j3  )
                  +h1y*F3(pf, EY, j1,l2+1,j3  ))
             +g1z*(hmy*F3(pf, EY, j1,l2-1,j3+1)
                  +h0y*F3(pf, EY, j1,l2  ,j3+1)
	          +h1y*F3(pf, EY, j1,l2+1,j3+1));

    particle_real_t ezq=hmz*(gmy*F3(pf, EZ, j1,j2-1,l3-1)
		  +g0y*F3(pf, EZ, j1,j2  ,l3-1)
                  +g1y*F3(pf, EZ, j1,j2+1,l3-1))
             +h0z*(gmy*F3(pf, EZ, j1,j2-1,l3  )
                  +g0y*F3(pf, EZ, j1,j2  ,l3  )
                  +g1y*F3(pf, EZ, j1,j2+1,l3  ))
             +h1z*(gmy*F3(pf, EZ, j1,j2-1,l3+1)
                  +g0y*F3(pf, EZ, j1,j2  ,l3+1)
		  +g1y*F3(pf, EZ, j1,j2+1,l3+1));

    particle_real_t hxq=hmz*(hmy*F3(pf, HX, j1,l2-1,l3-1)
                  +h0y*F3(pf, HX, j1,l2  ,l3-1)
                  +h1y*F3(pf, HX, j1,l2+1,l3-1))
             +h0z*(hmy*F3(pf, HX, j1,l2-1,l3  )
                  +h0y*F3(pf, HX, j1,l2  ,l3  )
                  +h1y*F3(pf, HX, j1,l2+1,l3  ))
             +h1z*(hmy*F3(pf, HX, j1,l2-1,l3+1)
                  +h0y*F3(pf, HX, j1,l2  ,l3+1)
		  +h1y*F3(pf, HX, j1,l2+1,l3+1));

    particle_real_t hyq=hmz*(gmy*F3(pf, HY, l1,j2-1,l3-1)
                  +g0y*F3(pf, HY, l1,j2  ,l3-1)
                  +g1y*F3(pf, HY, l1,j2+1,l3-1))
             +h0z*(gmy*F3(pf, HY, l1,j2-1,l3  )
                  +g0y*F3(pf, HY, l1,j2  ,l3  )
                  +g1y*F3(pf, HY, l1,j2+1,l3  ))
             +h1z*(gmy*F3(pf, HY, l1,j2-1,l3+1)
                  +g0y*F3(pf, HY, l1,j2  ,l3+1)
	          +g1y*F3(pf, HY, l1,j2+1,l3+1));

    particle_real_t hzq=gmz*(hmy*F3(pf, HZ, l1,l2-1,j3-1)
		  +h0y*F3(pf, HZ, l1,l2  ,j3-1)
                  +h1y*F3(pf, HZ, l1,l2+1,j3-1))
             +g0z*(hmy*F3(pf, HZ, l1,l2-1,j3  )
                  +h0y*F3(pf, HZ, l1,l2  ,j3  )
                  +h1y*F3(pf, HZ, l1,l2+1,j3  ))
             +g1z*(hmy*F3(pf, HZ, l1,l2-1,j3+1)
                  +h0y*F3(pf, HZ, l1,l2  ,j3+1)
		  +h1y*F3(pf, HZ, l1,l2+1,j3+1));

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    particle_real_t dq = dqs * part->qni / part->mni;
    particle_real_t pxm = part->pxi + dq*exq;
    particle_real_t pym = part->pyi + dq*eyq;
    particle_real_t pzm = part->pzi + dq*ezq;

    root = dq / particle_real_sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
    particle_real_t taux = hxq*root;
    particle_real_t tauy = hyq*root;
    particle_real_t tauz = hzq*root;

    particle_real_t tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    particle_real_t pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.f*taux*tauy+2.f*tauz)*pym + 
		(2.f*taux*tauz-2.f*tauy)*pzm)*tau;
    particle_real_t pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
		(1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.f*tauy*tauz+2.f*taux)*pzm)*tau;
    particle_real_t pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
		(2.f*tauy*tauz-2.f*taux)*pym +
		(1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
    
    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1.f / particle_real_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vxi = part->pxi * root;
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    particle_real_t yi = part->yi + vyi * yl;
    particle_real_t zi = part->zi + vzi * zl;

    v = yi * dyi;
    w = zi * dzi;
    int k2 = particle_real_nint(v);
    int k3 = particle_real_nint(w);
    h2 = k2 - v;
    h3 = k3 - w;

    for (int i = -2; i <= 2; i++) {
      S1Y(i) = 0.f;
      S1Z(i) = 0.f;
    }

    S1Y(k2-j2-1) = .5f*(1.5f-particle_real_abs(h2-1.f))*(1.5f-particle_real_abs(h2-1.f));
    S1Y(k2-j2+0) = .75f-particle_real_abs(h2)*particle_real_abs(h2);
    S1Y(k2-j2+1) = .5f*(1.5f-particle_real_abs(h2+1.f))*(1.5f-particle_real_abs(h2+1.f));
    S1Z(k3-j3-1) = .5f*(1.5f-particle_real_abs(h3-1.f))*(1.5f-particle_real_abs(h3-1.f));
    S1Z(k3-j3+0) = .75f-particle_real_abs(h3)*particle_real_abs(h3);
    S1Z(k3-j3+1) = .5f*(1.5f-particle_real_abs(h3+1.f))*(1.5f-particle_real_abs(h3+1.f));

    // CURRENT DENSITY AT (n+1.0)*dt

    for (int i = -1; i <= 1; i++) {
      S1Y(i) -= S0Y(i);
      S1Z(i) -= S0Z(i);
    }

    int l2min, l3min, l2max, l3max;
    
    if (k2 == j2) {
      l2min = -1; l2max = +1;
    } else if (k2 == j2 - 1) {
      l2min = -2; l2max = +1;
    } else { // (k2 == j2 + 1)
      l2min = -1; l2max = +2;
    }

    if (k3 == j3) {
      l3min = -1; l3max = +1;
    } else if (k3 == j3 - 1) {
      l3min = -2; l3max = +1;
    } else { // (k3 == j3 + 1)
      l3min = -1; l3max = +2;
    }

    particle_real_t jxh;
    particle_real_t jyh;
    particle_real_t jzh[5];

#define JZH(i) jzh[i+2]

    particle_real_t fnqx = vxi * part->qni * part->wni * fnqs;
    particle_real_t fnqy = part->qni * part->wni * fnqys;
    particle_real_t fnqz = part->qni * part->wni * fnqzs;
    for (int l2 = l2min; l2 <= l2max; l2++) {
      JZH(l2) = 0.f;
    }
    for (int l3 = l3min; l3 <= l3max; l3++) {
      jyh = 0.f;
      for (int l2 = l2min; l2 <= l2max; l2++) {
	particle_real_t wx = S0Y(l2) * S0Z(l3)
	  + .5f * S1Y(l2) * S0Z(l3)
	  + .5f * S0Y(l2) * S1Z(l3)
	  + (1.f/3.f) * S1Y(l2) * S1Z(l3);
	particle_real_t wy = S1Y(l2) * (S0Z(l3) + .5f*S1Z(l3));
	particle_real_t wz = S1Z(l3) * (S0Y(l2) + .5f*S1Y(l2));

	jxh = fnqx*wx;
	jyh -= fnqy*wy;
	JZH(l2) -= fnqz*wz;

	F3(pf, JXI, j1,j2+l2,j3+l3) += jxh;
	F3(pf, JYI, j1,j2+l2,j3+l3) += jyh;
	F3(pf, JZI, j1,j2+l2,j3+l3) += JZH(l2);
      }
    }
  }
}

void
psc_push_particles_1p5_c_push_a_yz(struct psc_push_particles *push,
				   struct psc_particles *_prts,
				   struct psc_fields *flds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("1p5_part_yz", 1., 0, 0);
  }
  
  prof_start(pr);
  psc_fields_zero_range(flds, JXI, JXI + 3);
  particle_range_t prts = particle_range_prts(_prts);
  do_1p5_push_part_yz(_prts->p, flds, prts);
  prof_stop(pr);
}

