
#include "psc_generic_c.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>

static void
do_genc_push_part_yz_b(fields_t *pf, particles_t *pp)
{
  creal dt = ppsc->dt;
  creal yl = .5f * dt;
  creal zl = .5f * dt;
  creal dqs = .5f * ppsc->coeff.eta * dt;
  creal dxi = 1.f / ppsc->dx[0];
  creal dyi = 1.f / ppsc->dx[1];
  creal dzi = 1.f / ppsc->dx[2];

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
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

    creal gmy=.5f*(.5f+h2)*(.5f+h2);
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0y=.75f-h2*h2;
    creal g0z=.75f-h3*h3;
    creal g1y=.5f*(.5f-h2)*(.5f-h2);
    creal g1z=.5f*(.5f-h3)*(.5f-h3);

    u = part->xi*dxi;
    v = part->yi*dyi-.5f;
    w = part->zi*dzi-.5f;
    int l1=nint(u);
    int l2=nint(v);
    int l3=nint(w);
    h1=l1-u;
    h2=l2-v;
    h3=l3-w;
    creal hmy=.5f*(.5f+h2)*(.5f+h2);
    creal hmz=.5f*(.5f+h3)*(.5f+h3);
    creal h0y=.75f-h2*h2;
    creal h0z=.75f-h3*h3;
    creal h1y=.5f*(.5f-h2)*(.5f-h2);
    creal h1z=.5f*(.5f-h3)*(.5f-h3);

    // FIELD INTERPOLATION

    creal exq=gmz*(gmy*F3(EX, l1,j2-1,j3-1)
		  +g0y*F3(EX, l1,j2  ,j3-1)
		  +g1y*F3(EX, l1,j2+1,j3-1))
             +g0z*(gmy*F3(EX, l1,j2-1,j3  )
    	          +g0y*F3(EX, l1,j2  ,j3  )
	          +g1y*F3(EX, l1,j2+1,j3  ))
             +g1z*(gmy*F3(EX, l1,j2-1,j3+1)
	          +g0y*F3(EX, l1,j2  ,j3+1)
	          +g1y*F3(EX, l1,j2+1,j3+1));

    creal eyq=gmz*(hmy*F3(EY, j1,l2-1,j3-1)
                  +h0y*F3(EY, j1,l2  ,j3-1)
                  +h1y*F3(EY, j1,l2+1,j3-1))
             +g0z*(hmy*F3(EY, j1,l2-1,j3  )
                  +h0y*F3(EY, j1,l2  ,j3  )
                  +h1y*F3(EY, j1,l2+1,j3  ))
             +g1z*(hmy*F3(EY, j1,l2-1,j3+1)
                  +h0y*F3(EY, j1,l2  ,j3+1)
	          +h1y*F3(EY, j1,l2+1,j3+1));

    creal ezq=hmz*(gmy*F3(EZ, j1,j2-1,l3-1)
		  +g0y*F3(EZ, j1,j2  ,l3-1)
                  +g1y*F3(EZ, j1,j2+1,l3-1))
             +h0z*(gmy*F3(EZ, j1,j2-1,l3  )
                  +g0y*F3(EZ, j1,j2  ,l3  )
                  +g1y*F3(EZ, j1,j2+1,l3  ))
             +h1z*(gmy*F3(EZ, j1,j2-1,l3+1)
                  +g0y*F3(EZ, j1,j2  ,l3+1)
		  +g1y*F3(EZ, j1,j2+1,l3+1));

    creal hxq=hmz*(hmy*F3(HX, j1,l2-1,l3-1)
                  +h0y*F3(HX, j1,l2  ,l3-1)
                  +h1y*F3(HX, j1,l2+1,l3-1))
             +h0z*(hmy*F3(HX, j1,l2-1,l3  )
                  +h0y*F3(HX, j1,l2  ,l3  )
                  +h1y*F3(HX, j1,l2+1,l3  ))
             +h1z*(hmy*F3(HX, j1,l2-1,l3+1)
                  +h0y*F3(HX, j1,l2  ,l3+1)
		  +h1y*F3(HX, j1,l2+1,l3+1));

    creal hyq=hmz*(gmy*F3(HY, l1,j2-1,l3-1)
                  +g0y*F3(HY, l1,j2  ,l3-1)
                  +g1y*F3(HY, l1,j2+1,l3-1))
             +h0z*(gmy*F3(HY, l1,j2-1,l3  )
                  +g0y*F3(HY, l1,j2  ,l3  )
                  +g1y*F3(HY, l1,j2+1,l3  ))
             +h1z*(gmy*F3(HY, l1,j2-1,l3+1)
                  +g0y*F3(HY, l1,j2  ,l3+1)
	          +g1y*F3(HY, l1,j2+1,l3+1));

    creal hzq=gmz*(hmy*F3(HZ, l1,l2-1,j3-1)
		  +h0y*F3(HZ, l1,l2  ,j3-1)
                  +h1y*F3(HZ, l1,l2+1,j3-1))
             +g0z*(hmy*F3(HZ, l1,l2-1,j3  )
                  +h0y*F3(HZ, l1,l2  ,j3  )
                  +h1y*F3(HZ, l1,l2+1,j3  ))
             +g1z*(hmy*F3(HZ, l1,l2-1,j3+1)
                  +h0y*F3(HZ, l1,l2  ,j3+1)
		  +h1y*F3(HZ, l1,l2+1,j3+1));

     // c x^(n+0.5), p^n -> x^(n+1.0), p^(n+1.0) 

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
    vyi = part->pyi * root;
    vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;
  }
}

void
psc_push_particles_generic_c_push_yz_b(struct psc_push_particles *push,
				       mparticles_base_t *particles_base,
				       mfields_base_t *flds_base)
{
  mfields_t flds;
  mparticles_t particles;
  fields_get(&flds, EX, EX + 6, flds_base);
  particles_get(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_yz_b", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_genc_push_part_yz_b(&flds.f[p], &particles.p[p]);
  }
  prof_stop(pr);

  particles_put(&particles, particles_base);
  fields_put(&flds, JXI, JXI + 3, flds_base);
}

