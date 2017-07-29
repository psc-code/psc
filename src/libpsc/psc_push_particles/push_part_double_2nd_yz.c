
#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_double.h"

#define F3_CACHE F3

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// FIXME, this really should share almost all code with its generic_c counterpart

static void
do_push_part_yz(int p, struct psc_fields *pf, struct psc_particles *pp)
{
#define S0Y(off) s0y[off+2]
#define S0Z(off) s0z[off+2]
#define S1Y(off) s1y[off+2]
#define S1Z(off) s1z[off+2]

  particle_real_t s0y[5] = {}, s0z[5] = {}, s1y[5], s1z[5];

  particle_real_t dt = ppsc->dt;
  particle_real_t dqs = .5f * ppsc->coeff.eta * dt;
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t fnqys = ppsc->patch[p].dx[1] * fnqs / dt;
  particle_real_t fnqzs = ppsc->patch[p].dx[2] * fnqs / dt;
  particle_real_t dxi = 1.f / ppsc->patch[p].dx[0];
  particle_real_t dyi = 1.f / ppsc->patch[p].dx[1];
  particle_real_t dzi = 1.f / ppsc->patch[p].dx[2];

  //int *ldims = ppsc->patch[p].ldims;

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    // x^n, p^n -> x^(n+.5), p^n

    particle_real_t root = 1.f / particle_real_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    particle_real_t vxi = part->pxi * root;
    particle_real_t vyi = part->pyi * root;
    particle_real_t vzi = part->pzi * root;

    /* assert(part->yi * dyi >= 0 && part->yi * dyi <= ldims[1]); */
    /* assert(part->zi * dzi >= 0 && part->zi * dzi <= ldims[2]); */

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
    particle_real_t hmy=.5f*(.5f+h2)*(.5f+h2);
    particle_real_t hmz=.5f*(.5f+h3)*(.5f+h3);
    particle_real_t h0y=.75f-h2*h2;
    particle_real_t h0z=.75f-h3*h3;
    particle_real_t h1y=.5f*(.5f-h2)*(.5f-h2);
    particle_real_t h1z=.5f*(.5f-h3)*(.5f-h3);

    // FIELD INTERPOLATION

    /* assert(l2 - 1 >= -1); */
    /* assert(j2 - 1 >= -1); */
    /* assert(l2 + 1 < ldims[1] + 1); */
    /* assert(j2 + 1 < ldims[1] + 2); */
    l1 = 0;
    j1 = 0;
    particle_real_t exq=gmz*(gmy*F3_CACHE(pf, EX, l1,j2-1,j3-1)
		  +g0y*F3_CACHE(pf, EX, l1,j2  ,j3-1)
		  +g1y*F3_CACHE(pf, EX, l1,j2+1,j3-1))
             +g0z*(gmy*F3_CACHE(pf, EX, l1,j2-1,j3  )
    	          +g0y*F3_CACHE(pf, EX, l1,j2  ,j3  )
	          +g1y*F3_CACHE(pf, EX, l1,j2+1,j3  ))
             +g1z*(gmy*F3_CACHE(pf, EX, l1,j2-1,j3+1)
	          +g0y*F3_CACHE(pf, EX, l1,j2  ,j3+1)
	          +g1y*F3_CACHE(pf, EX, l1,j2+1,j3+1));

    particle_real_t eyq=gmz*(hmy*F3_CACHE(pf, EY, j1,l2-1,j3-1)
                  +h0y*F3_CACHE(pf, EY, j1,l2  ,j3-1)
                  +h1y*F3_CACHE(pf, EY, j1,l2+1,j3-1))
             +g0z*(hmy*F3_CACHE(pf, EY, j1,l2-1,j3  )
                  +h0y*F3_CACHE(pf, EY, j1,l2  ,j3  )
                  +h1y*F3_CACHE(pf, EY, j1,l2+1,j3  ))
             +g1z*(hmy*F3_CACHE(pf, EY, j1,l2-1,j3+1)
                  +h0y*F3_CACHE(pf, EY, j1,l2  ,j3+1)
	          +h1y*F3_CACHE(pf, EY, j1,l2+1,j3+1));

    particle_real_t ezq=hmz*(gmy*F3_CACHE(pf, EZ, j1,j2-1,l3-1)
		  +g0y*F3_CACHE(pf, EZ, j1,j2  ,l3-1)
                  +g1y*F3_CACHE(pf, EZ, j1,j2+1,l3-1))
             +h0z*(gmy*F3_CACHE(pf, EZ, j1,j2-1,l3  )
                  +g0y*F3_CACHE(pf, EZ, j1,j2  ,l3  )
                  +g1y*F3_CACHE(pf, EZ, j1,j2+1,l3  ))
             +h1z*(gmy*F3_CACHE(pf, EZ, j1,j2-1,l3+1)
                  +g0y*F3_CACHE(pf, EZ, j1,j2  ,l3+1)
		  +g1y*F3_CACHE(pf, EZ, j1,j2+1,l3+1));

    particle_real_t hxq=hmz*(hmy*F3_CACHE(pf, HX, j1,l2-1,l3-1)
                  +h0y*F3_CACHE(pf, HX, j1,l2  ,l3-1)
                  +h1y*F3_CACHE(pf, HX, j1,l2+1,l3-1))
             +h0z*(hmy*F3_CACHE(pf, HX, j1,l2-1,l3  )
                  +h0y*F3_CACHE(pf, HX, j1,l2  ,l3  )
                  +h1y*F3_CACHE(pf, HX, j1,l2+1,l3  ))
             +h1z*(hmy*F3_CACHE(pf, HX, j1,l2-1,l3+1)
                  +h0y*F3_CACHE(pf, HX, j1,l2  ,l3+1)
		  +h1y*F3_CACHE(pf, HX, j1,l2+1,l3+1));

    particle_real_t hyq=hmz*(gmy*F3_CACHE(pf, HY, l1,j2-1,l3-1)
                  +g0y*F3_CACHE(pf, HY, l1,j2  ,l3-1)
                  +g1y*F3_CACHE(pf, HY, l1,j2+1,l3-1))
             +h0z*(gmy*F3_CACHE(pf, HY, l1,j2-1,l3  )
                  +g0y*F3_CACHE(pf, HY, l1,j2  ,l3  )
                  +g1y*F3_CACHE(pf, HY, l1,j2+1,l3  ))
             +h1z*(gmy*F3_CACHE(pf, HY, l1,j2-1,l3+1)
                  +g0y*F3_CACHE(pf, HY, l1,j2  ,l3+1)
	          +g1y*F3_CACHE(pf, HY, l1,j2+1,l3+1));

    particle_real_t hzq=gmz*(hmy*F3_CACHE(pf, HZ, l1,l2-1,j3-1)
		  +h0y*F3_CACHE(pf, HZ, l1,l2  ,j3-1)
                  +h1y*F3_CACHE(pf, HZ, l1,l2+1,j3-1))
             +g0z*(hmy*F3_CACHE(pf, HZ, l1,l2-1,j3  )
                  +h0y*F3_CACHE(pf, HZ, l1,l2  ,j3  )
                  +h1y*F3_CACHE(pf, HZ, l1,l2+1,j3  ))
             +g1z*(hmy*F3_CACHE(pf, HZ, l1,l2-1,j3+1)
                  +h0y*F3_CACHE(pf, HZ, l1,l2  ,j3+1)
		  +h1y*F3_CACHE(pf, HZ, l1,l2+1,j3+1));

     // c x^(n+.5), p^n -> x^(n+1.0), p^(n+1.0) 

    particle_real_t dq = dqs * particle_qni_div_mni(part);
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

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+.5), p^(n+1) -> x^(n+1.5), p^(n+1)

    part->yi += vyi * dt;
    part->zi += vzi * dt;

    v = part->yi * dyi;
    w = part->zi * dzi;
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

    particle_real_t fnqx = vxi * particle_qni_wni(part) * fnqs;
    particle_real_t fnqy = particle_qni_wni(part) * fnqys;
    particle_real_t fnqz = particle_qni_wni(part) * fnqzs;
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

	F3_CACHE(pf, JXI, j1,j2+l2,j3+l3) += jxh;
	F3_CACHE(pf, JYI, j1,j2+l2,j3+l3) += jyh;
	F3_CACHE(pf, JZI, j1,j2+l2,j3+l3) += JZH(l2);
      }
    }
  }
}

static struct psc_fields *
cache_fields_from_em(fields_t *pf)
{
  struct psc_fields *fld = psc_fields_create(MPI_COMM_NULL);
  psc_fields_set_type(fld, FIELDS_TYPE);
  // FIXME, can do -1 .. 2? NO!, except maybe for 1st order
  // Has to be at least -2 .. +3 because of staggering
  // FIXME, get rid of caching since it's no different from the actual
  // fields...
  psc_fields_set_param_int3(fld, "ib", pf->ib);
  psc_fields_set_param_int3(fld, "im", pf->im);
  psc_fields_set_param_int(fld, "nr_comp", 9); // JX .. HZ
  psc_fields_set_param_int(fld, "p", pf->p);
  psc_fields_setup(fld);
  for (int iz = fld->ib[2]; iz < fld->ib[2] + fld->im[2]; iz++) {
    for (int iy = fld->ib[1]; iy < fld->ib[1] + fld->im[1]; iy++) {
      F3_CACHE(fld, EX, 0,iy,iz) = F3(pf, EX, 0,iy,iz);
      F3_CACHE(fld, EY, 0,iy,iz) = F3(pf, EY, 0,iy,iz);
      F3_CACHE(fld, EZ, 0,iy,iz) = F3(pf, EZ, 0,iy,iz);
      F3_CACHE(fld, HX, 0,iy,iz) = F3(pf, HX, 0,iy,iz);
      F3_CACHE(fld, HY, 0,iy,iz) = F3(pf, HY, 0,iy,iz);
      F3_CACHE(fld, HZ, 0,iy,iz) = F3(pf, HZ, 0,iy,iz);
    }
  }
  return fld;
}

static void
cache_fields_to_j(struct psc_fields *fld, fields_t *pf)
{
  for (int iz = fld->ib[2]; iz < fld->ib[2] + fld->im[2]; iz++) {
    for (int iy = fld->ib[1]; iy < fld->ib[1] + fld->im[1]; iy++) {
      F3(pf, JXI, 0,iy,iz) += F3_CACHE(fld, JXI, 0,iy,iz);
      F3(pf, JYI, 0,iy,iz) += F3_CACHE(fld, JYI, 0,iy,iz);
      F3(pf, JZI, 0,iy,iz) += F3_CACHE(fld, JZI, 0,iy,iz);
    }
  }
}

void
psc_push_particles_2nd_double_push_a_yz(struct psc_push_particles *push,
					struct psc_particles *prts,
					struct psc_fields *flds)
{
  psc_fields_zero_range(flds, JXI, JXI + 3);
  struct psc_fields *flds_cache = cache_fields_from_em(flds);
  do_push_part_yz(prts->p, flds_cache, prts);
  cache_fields_to_j(flds_cache, flds);
  psc_fields_destroy(flds_cache);
}

