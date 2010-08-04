
#include "psc_generic_c.h"

#include "util/profile.h"
#include <math.h>
#include <stdlib.h>

static void
do_genc_push_part_yz_b(psc_fields_c_t *pf, psc_particles_c_t *pp)
{
  creal dt = psc.dt;
  creal yl = .5f * dt;
  creal zl = .5f * dt;
  creal dqs = .5f * psc.coeff.eta * dt;
  creal dxi = 1.f / psc.dx[0];
  creal dyi = 1.f / psc.dx[1];
  creal dzi = 1.f / psc.dx[2];

  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_c_t *part = psc_particles_c_get_one(pp, n);

    creal root = 1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;

    creal u = part->xi * dxi;
    creal v = part->yi * dyi;
    creal w = part->zi * dzi;
    int j1 = nint(u);
    int j2 = nint(v);
    int j3 = nint(w);
    creal h1 = j1-u;
    creal h2 = j2-v;
    creal h3 = j3-w;

    creal gmy=0.5*(0.5+h2)*(0.5+h2);
    creal gmz=0.5*(0.5+h3)*(0.5+h3);
    creal g0y=0.75-h2*h2;
    creal g0z=0.75-h3*h3;
    creal g1y=0.5*(0.5-h2)*(0.5-h2);
    creal g1z=0.5*(0.5-h3)*(0.5-h3);

    u = part->xi*dxi;
    v = part->yi*dyi-0.5;
    w = part->zi*dzi-0.5;
    int l1=nint(u);
    int l2=nint(v);
    int l3=nint(w);
    h1=l1-u;
    h2=l2-v;
    h3=l3-w;
    creal hmy=0.5*(0.5+h2)*(0.5+h2);
    creal hmz=0.5*(0.5+h3)*(0.5+h3);
    creal h0y=0.75-h2*h2;
    creal h0z=0.75-h3*h3;
    creal h1y=0.5*(0.5-h2)*(0.5-h2);
    creal h1z=0.5*(0.5-h3)*(0.5-h3);

    // FIELD INTERPOLATION

    creal exq=gmz*(gmy*F3(EX, l1,j2-1,j3-1)
		  +g0y*F3(EX, l1,j2,j3-1)
		  +g1y*F3(EX, l1,j2+1,j3-1))
      +g0z*(gmy*F3(EX, l1,j2-1,j3)
	    +g0y*F3(EX, l1,j2,j3)
	    +g1y*F3(EX, l1,j2+1,j3))
      +g1z*(gmy*F3(EX, l1,j2-1,j3+1)
	    +g0y*F3(EX, l1,j2,j3+1)
	    +g1y*F3(EX, l1,j2+1,j3+1));
    
    creal eyq=gmz*(hmy*F3(EY, j1,l2-1,j3-1)
		  +h0y*F3(EY, j1,l2,j3-1)
		  +h1y*F3(EY, j1,l2+1,j3-1))
      +g0z*(hmy*F3(EY, j1,l2-1,j3)
	    +h0y*F3(EY, j1,l2,j3)
	    +h1y*F3(EY, j1,l2+1,j3))
      +g1z*(hmy*F3(EY, j1,l2-1,j3+1)
	    +h0y*F3(EY, j1,l2,j3+1)
	    +h1y*F3(EY, j1,l2+1,j3+1));

    creal ezq=hmz*(gmy*F3(EZ, j1,j2-1,l3-1)
		  +g0y*F3(EZ, j1,j2,l3-1)
		  +g1y*F3(EZ, j1,j2+1,l3-1))
      +h0z*(gmy*F3(EZ, j1,j2-1,l3)
	    +g0y*F3(EZ, j1,j2,l3)
	    +g1y*F3(EZ, j1,j2+1,l3))
      +h1z*(gmy*F3(EZ, j1,j2-1,l3+1)
	    +g0y*F3(EZ, j1,j2,l3+1)
	    +g1y*F3(EZ, j1,j2+1,l3+1));

    creal bxq=hmz*(hmy*F3(BX, j1,l2-1,l3-1)
		  +h0y*F3(BX, j1,l2,l3-1)
		  +h1y*F3(BX, j1,l2+1,l3-1))
      +h0z*(hmy*F3(BX, j1,l2-1,l3)
	    +h0y*F3(BX, j1,l2,l3)
	    +h1y*F3(BX, j1,l2+1,l3))
      +h1z*(hmy*F3(BX, j1,l2-1,l3+1)
	    +h0y*F3(BX, j1,l2,l3+1)
	    +h1y*F3(BX, j1,l2+1,l3+1));

    creal byq=hmz*(gmy*F3(BY, l1,j2-1,l3-1)
		  +g0y*F3(BY, l1,j2,l3-1)
		  +g1y*F3(BY, l1,j2+1,l3-1))
      +h0z*(gmy*F3(BY, l1,j2-1,l3)
	    +g0y*F3(BY, l1,j2,l3)
	    +g1y*F3(BY, l1,j2+1,l3))
      +h1z*(gmy*F3(BY, l1,j2-1,l3+1)
	    +g0y*F3(BY, l1,j2,l3+1)
	    +g1y*F3(BY, l1,j2+1,l3+1));
      
    creal bzq=gmz*(hmy*F3(BZ, l1,l2-1,j3-1)
		  +h0y*F3(BZ, l1,l2,j3-1)
		  +h1y*F3(BZ, l1,l2+1,j3-1))
      +g0z*(hmy*F3(BZ, l1,l2-1,j3)
	    +h0y*F3(BZ, l1,l2,j3)
	    +h1y*F3(BZ, l1,l2+1,j3))
      +g1z*(hmy*F3(BZ, l1,l2-1,j3+1)
	    +h0y*F3(BZ, l1,l2,j3+1)
	    +h1y*F3(BZ, l1,l2+1,j3+1));

     // c x^(n+0.5), p^n -> x^(n+1.0), p^(n+1.0) 

    creal dq = dqs * part->qni / part->mni;
    creal pxm = part->pxi + dq*exq;
    creal pym = part->pyi + dq*eyq;
    creal pzm = part->pzi + dq*ezq;

    root = dq / sqrtf(1.f + pxm*pxm + pym*pym + pzm*pzm);
    creal taux = bxq*root;
    creal tauy = byq*root;
    creal tauz = bzq*root;

    creal tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
    creal pxp = ((1.0+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
		(2.0*taux*tauy+2.0*tauz)*pym + 
		(2.0*taux*tauz-2.0*tauy)*pzm)*tau;
    creal pyp = ((2.0*taux*tauy-2.0*tauz)*pxm +
		(1.0-taux*taux+tauy*tauy-tauz*tauz)*pym +
		(2.0*tauy*tauz+2.0*taux)*pzm)*tau;
    creal pzp = ((2.0*taux*tauz+2.0*tauy)*pxm +
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
}

void
genc_push_part_yz_b()
{
  psc_fields_c_t pf;
  genc_fields_from_fortran(&pf);
  psc_particles_c_t pp;
  genc_particles_from_fortran(&pp);

  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_yz_b", 1., 0, psc.pp.n_part * 12 * sizeof(creal));
  }
  prof_start(pr);
  do_genc_push_part_yz_b(&pf, &pp);
  prof_stop(pr);

  genc_fields_to_fortran(&pf);
  genc_particles_to_fortran(&pp);
}
