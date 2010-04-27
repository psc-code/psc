
#include "psc_generic_c.h"

#include "profile/profile.h"
#include <math.h>
#include <stdlib.h>

static inline int
nint(real x)
{
  return (int)(x + real(10.5)) - 10;
}

void
genc_push_part_yz_b()
{
  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_yz_b", 1., 0, psc.n_part * 12 * sizeof(float));
  }
  prof_start(pr);
 
  struct psc_genc *genc = psc.c_ctx;

  real dt = psc.dt;
  real yl = .5f * dt;
  real zl = .5f * dt;
  real dqs = .5f * psc.prm.eta * dt;
  real dxi = 1.f / psc.dx[0];
  real dyi = 1.f / psc.dx[1];
  real dzi = 1.f / psc.dx[2];

  for (int n = 0; n < psc.n_part; n++) {
    struct c_particle *part = &genc->part[n];

    real root = 1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    real vyi = part->pyi * root;
    real vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    real u = part->xi * dxi;
    real v = part->yi * dyi;
    real w = part->zi * dzi;
    int j1 = nint(u);
    int j2 = nint(v);
    int j3 = nint(w);
    real h1 = j1-u;
    real h2 = j2-v;
    real h3 = j3-w;

    real gmy=0.5*(0.5+h2)*(0.5+h2);
    real gmz=0.5*(0.5+h3)*(0.5+h3);
    real g0y=0.75-h2*h2;
    real g0z=0.75-h3*h3;
    real g1y=0.5*(0.5-h2)*(0.5-h2);
    real g1z=0.5*(0.5-h3)*(0.5-h3);

    u = part->xi*dxi;
    v = part->yi*dyi-0.5;
    w = part->zi*dzi-0.5;
    int l1=nint(u);
    int l2=nint(v);
    int l3=nint(w);
    h1=l1-u;
    h2=l2-v;
    h3=l3-w;
    real hmy=0.5*(0.5+h2)*(0.5+h2);
    real hmz=0.5*(0.5+h3)*(0.5+h3);
    real h0y=0.75-h2*h2;
    real h0z=0.75-h3*h3;
    real h1y=0.5*(0.5-h2)*(0.5-h2);
    real h1z=0.5*(0.5-h3)*(0.5-h3);

    // FIELD INTERPOLATION

    real exq=gmz*(gmy*F3(EX, l1,j2-1,j3-1)
		  +g0y*F3(EX, l1,j2,j3-1)
		  +g1y*F3(EX, l1,j2+1,j3-1))
      +g0z*(gmy*F3(EX, l1,j2-1,j3)
	    +g0y*F3(EX, l1,j2,j3)
	    +g1y*F3(EX, l1,j2+1,j3))
      +g1z*(gmy*F3(EX, l1,j2-1,j3+1)
	    +g0y*F3(EX, l1,j2,j3+1)
	    +g1y*F3(EX, l1,j2+1,j3+1));
    
    real eyq=gmz*(hmy*F3(EY, j1,l2-1,j3-1)
		  +h0y*F3(EY, j1,l2,j3-1)
		  +h1y*F3(EY, j1,l2+1,j3-1))
      +g0z*(hmy*F3(EY, j1,l2-1,j3)
	    +h0y*F3(EY, j1,l2,j3)
	    +h1y*F3(EY, j1,l2+1,j3))
      +g1z*(hmy*F3(EY, j1,l2-1,j3+1)
	    +h0y*F3(EY, j1,l2,j3+1)
	    +h1y*F3(EY, j1,l2+1,j3+1));

    real ezq=hmz*(gmy*F3(EZ, j1,j2-1,l3-1)
		  +g0y*F3(EZ, j1,j2,l3-1)
		  +g1y*F3(EZ, j1,j2+1,l3-1))
      +h0z*(gmy*F3(EZ, j1,j2-1,l3)
	    +g0y*F3(EZ, j1,j2,l3)
	    +g1y*F3(EZ, j1,j2+1,l3))
      +h1z*(gmy*F3(EZ, j1,j2-1,l3+1)
	    +g0y*F3(EZ, j1,j2,l3+1)
	    +g1y*F3(EZ, j1,j2+1,l3+1));

    real bxq=hmz*(hmy*F3(BX, j1,l2-1,l3-1)
		  +h0y*F3(BX, j1,l2,l3-1)
		  +h1y*F3(BX, j1,l2+1,l3-1))
      +h0z*(hmy*F3(BX, j1,l2-1,l3)
	    +h0y*F3(BX, j1,l2,l3)
	    +h1y*F3(BX, j1,l2+1,l3))
      +h1z*(hmy*F3(BX, j1,l2-1,l3+1)
	    +h0y*F3(BX, j1,l2,l3+1)
	    +h1y*F3(BX, j1,l2+1,l3+1));

    real byq=hmz*(gmy*F3(BY, l1,j2-1,l3-1)
		  +g0y*F3(BY, l1,j2,l3-1)
		  +g1y*F3(BY, l1,j2+1,l3-1))
      +h0z*(gmy*F3(BY, l1,j2-1,l3)
	    +g0y*F3(BY, l1,j2,l3)
	    +g1y*F3(BY, l1,j2+1,l3))
      +h1z*(gmy*F3(BY, l1,j2-1,l3+1)
	    +g0y*F3(BY, l1,j2,l3+1)
	    +g1y*F3(BY, l1,j2+1,l3+1));
      
    real bzq=gmz*(hmy*F3(BZ, l1,l2-1,j3-1)
		  +h0y*F3(BZ, l1,l2,j3-1)
		  +h1y*F3(BZ, l1,l2+1,j3-1))
      +g0z*(hmy*F3(BZ, l1,l2-1,j3)
	    +h0y*F3(BZ, l1,l2,j3)
	    +h1y*F3(BZ, l1,l2+1,j3))
      +g1z*(hmy*F3(BZ, l1,l2-1,j3+1)
	    +h0y*F3(BZ, l1,l2,j3+1)
	    +h1y*F3(BZ, l1,l2+1,j3+1));

     // c x^(n+0.5), p^n -> x^(n+1.0), p^(n+1.0) 

    real dq = dqs * part->qni / part->mni;
    real pxm = part->pxi + dq*exq;
    real pym = part->pyi + dq*eyq;
    real pzm = part->pzi + dq*ezq;

    root = dq / sqrtf(1.f + pxm*pxm + pym*pym + pzm*pzm);
    real taux = bxq*root;
    real tauy = byq*root;
    real tauz = bzq*root;

    real tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    real pxp = ((1.0+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.0*taux*tauy+2.0*tauz)*pym + 
		(2.0*taux*tauz-2.0*tauy)*pzm)*tau;
    real pyp = ((2.0*taux*tauy-2.0*tauz)*pxm +
		(1.0-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.0*tauy*tauz+2.0*taux)*pzm)*tau;
    real pzp = ((2.0*taux*tauz+2.0*tauy)*pxm +
		(2.0*tauy*tauz-2.0*taux)*pym +
		(1.0-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;

    part->pxi = pxp + dq * exq;
    part->pyi = pyp + dq * eyq;
    part->pzi = pzp + dq * ezq;

    root = 1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;
  }

  prof_stop(pr);
}
