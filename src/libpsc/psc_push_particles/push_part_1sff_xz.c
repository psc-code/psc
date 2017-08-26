
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM DIM_XZ
#define ORDER ORDER_1ST
#define VARIANT VARIANT_SFF
#include "push_part_common.c"

static void
do_push_part(int p, fields_t flds, particle_range_t prts)
{
#include "push_part_common_vars.c"

  struct psc_patch *patch = &ppsc->patch[p];

  // FIXME, eventually no ghost points should be needed (?)
  fields_t flds_avg = fields_t_ctor((int[3]) { -1, 0, -1 },
				    (int[3]) { patch->ldims[0] + 2, 1, patch->ldims[2] + 1 },
				    6);
  
  for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {
    for (int ix = -1; ix < patch->ldims[0] + 1; ix++) {
      _F3(flds_avg, 0, ix,0,iz) = .5 * (_F3(flds, EX, ix,0,iz) + _F3(flds, EX, ix-1,0,iz));
      _F3(flds_avg, 1, ix,0,iz) = _F3(flds, EY, ix,0,iz);
      _F3(flds_avg, 2, ix,0,iz) = .5 * (_F3(flds, EZ, ix,0,iz) + _F3(flds, EZ, ix,0,iz-1));
      _F3(flds_avg, 3, ix,0,iz) = .5 * (_F3(flds, HX, ix,0,iz) + _F3(flds, HX, ix,0,iz-1));
      _F3(flds_avg, 4, ix,0,iz) = .25 * (_F3(flds, HY, ix  ,0,iz) + _F3(flds, HY, ix  ,0,iz-1) +
					_F3(flds, HY, ix-1,0,iz) + _F3(flds, HY, ix-1,0,iz-1));
      _F3(flds_avg, 5, ix,0,iz) = .5 * (_F3(flds, HZ, ix,0,iz) + _F3(flds, HZ, ix-1,0,iz));
    }
  }

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);
    particle_real_t *x = &part->xi;

    // x^n, p^n -> x^(n+.5), p^n
    particle_real_t vv[3];
    calc_v(vv, &part->pxi);
    push_x(x, vv);

    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 

    IF_DIM_X( DEPOSIT_AND_IP_COEFFS(lg1, lh1, gx, hx, 0, dxi, s0x); );
    IF_DIM_Y( DEPOSIT_AND_IP_COEFFS(lg2, lh2, gy, hy, 0, dyi, s0y); );
    IF_DIM_Z( DEPOSIT_AND_IP_COEFFS(lg3, lh3, gz, hz, 0, dzi, s0z); );

    // FIELD INTERPOLATION

    particle_real_t E[3] = { IP_FIELD(flds_avg, EX-EX, g, g, g),
			     IP_FIELD(flds_avg, EY-EX, g, g, g),
			     IP_FIELD(flds_avg, EZ-EX, g, g, g), };
    particle_real_t H[3] = { IP_FIELD(flds_avg, HX-EX, g, g, g),
			     IP_FIELD(flds_avg, HY-EX, g, g, g),
			     IP_FIELD(flds_avg, HZ-EX, g, g, g), };

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    particle_real_t dq = dqs * particle_qni_div_mni(part);
    push_p(&part->pxi, E, H, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_v(vv, &part->pxi);
    push_x(x, vv);

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)

    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)
    particle_real_t xn[3] = { x[0], x[1], x[2] };
    push_x(xn, vv);

    ZERO_S1;
    IF_DIM_X( DEPOSIT(k1, gx, 0, dxi, s1x, lg1); );
    IF_DIM_Y( DEPOSIT(k2, gy, 1, dyi, s1y, lg2); );
    IF_DIM_Z( DEPOSIT(k3, gz, 2, dzi, s1z, lg3); );

    // CURRENT DENSITY AT (n+1.0)*dt

    SUBTR_S1_S0;
    CURRENT_PREP;
    CURRENT;
  }

  fields_t_dtor(&flds_avg);
}

void
psc_push_particles_1sff_push_mprts_xz(struct psc_push_particles *push,
				      struct psc_mparticles *mprts,
				      struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register(PARTICLE_TYPE "_1sff_push_xz", 1., 0, 0);
  }

  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    fields_t flds = fields_t_mflds(mflds, p);
    fields_t_zero_range(flds, JXI, JXI + 3);
    do_push_part(p, flds, prts);
  }
  prof_stop(pr);
}

