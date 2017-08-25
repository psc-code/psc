
#include "psc_push_particles_1st.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "inc_defs.h"

#define DIM DIM_YZ
#define ORDER ORDER_1ST
#include "push_part_common.c"

#include "inc_params.c"
#include "inc_cache.c"
#include "inc_interpolate.c"
#include "inc_push.c"

static void
do_push_part_1st_yz(int p, fields_t flds, particle_range_t prts)
{
#include "push_part_common_vars.c"
  (void) dyi; // FIXME, avoid warnings
  (void) dzi;

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

    particle_real_t E[3] = { IP_FIELD(flds, EX, h, g, g),
			     IP_FIELD(flds, EY, g, h, g),
			     IP_FIELD(flds, EZ, g, g, h), };
    particle_real_t H[3] = { IP_FIELD(flds, HX, g, h, h),
			     IP_FIELD(flds, HY, h, g, h),
			     IP_FIELD(flds, HZ, h, h, g), };

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    particle_real_t dq = dqs * particle_qni_div_mni(part);
    push_p(&part->pxi, E, H, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_vxi(vv, part);
    push_xi(part, vv, .5 * dt);

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 

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
}

void
psc_push_particles_1st_push_mprts_yz(struct psc_push_particles *push,
				     struct psc_mparticles *mprts,
				     struct psc_mfields *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("1st_part_yz", 1., 0, 0);
  }

  prof_start(pr);
  params_1vb_set(ppsc, NULL, NULL);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);
    fields_t_zero_range(flds, JXI, JXI + 3);
    do_push_part_1st_yz(p, flds, prts);
  }
  prof_stop(pr);
}

