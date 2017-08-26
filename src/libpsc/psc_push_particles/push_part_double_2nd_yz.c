
#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_double.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIM DIM_YZ
#define ORDER ORDER_2ND
#include "push_part_common.c"

static void
do_push_part(int p, fields_t flds, particle_range_t prts)
{
#include "push_part_common_vars.c"

  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);
    particle_real_t *x = &part->xi;

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

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
    particle_real_t vv[3];
    calc_v(vv, &part->pxi);
    // FIXME, inelegant way of pushing full dt
    push_x(x, vv);
    push_x(x, vv);

    ZERO_S1;
    IF_DIM_X( DEPOSIT(x, k1, gx, 0, dxi, s1x, lg1); );
    IF_DIM_Y( DEPOSIT(x, k2, gy, 1, dyi, s1y, lg2); );
    IF_DIM_Z( DEPOSIT(x, k3, gz, 2, dzi, s1z, lg3); );

    // CURRENT DENSITY AT (n+1.0)*dt

    SUBTR_S1_S0;
    CURRENT_PREP;
    CURRENT;
  }
}

static fields_t
cache_fields_from_em(fields_t flds)
{
  fields_t fld = fields_t_ctor(flds.ib, flds.im, 9);
  // FIXME, can do -1 .. 2? NO!, except maybe for 1st order
  // Has to be at least -2 .. +3 because of staggering
  // FIXME, get rid of caching since it's no different from the actual
  // fields...
  for (int iz = fld.ib[2]; iz < fld.ib[2] + fld.im[2]; iz++) {
    for (int iy = fld.ib[1]; iy < fld.ib[1] + fld.im[1]; iy++) {
      _F3(fld, EX, 0,iy,iz) = _F3(flds, EX, 0,iy,iz);
      _F3(fld, EY, 0,iy,iz) = _F3(flds, EY, 0,iy,iz);
      _F3(fld, EZ, 0,iy,iz) = _F3(flds, EZ, 0,iy,iz);
      _F3(fld, HX, 0,iy,iz) = _F3(flds, HX, 0,iy,iz);
      _F3(fld, HY, 0,iy,iz) = _F3(flds, HY, 0,iy,iz);
      _F3(fld, HZ, 0,iy,iz) = _F3(flds, HZ, 0,iy,iz);
    }
  }
  return fld;
}

static void
cache_fields_to_j(fields_t fld, fields_t flds)
{
  for (int iz = fld.ib[2]; iz < fld.ib[2] + fld.im[2]; iz++) {
    for (int iy = fld.ib[1]; iy < fld.ib[1] + fld.im[1]; iy++) {
      _F3(flds, JXI, 0,iy,iz) += _F3(fld, JXI, 0,iy,iz);
      _F3(flds, JYI, 0,iy,iz) += _F3(fld, JYI, 0,iy,iz);
      _F3(flds, JZI, 0,iy,iz) += _F3(fld, JZI, 0,iy,iz);
    }
  }
}

void
psc_push_particles_2nd_double_push_mprts_yz(struct psc_push_particles *push,
					    struct psc_mparticles *mprts,
					    struct psc_mfields *mflds)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    particle_range_t prts = particle_range_mprts(mprts, p);

    // FIXME, can't we just skip this and just set j when copying back?
    fields_t_zero_range(flds, JXI, JXI + 3);
    fields_t flds_cache = cache_fields_from_em(flds);
    do_push_part(p, flds_cache, prts);
    cache_fields_to_j(flds_cache, flds);
    fields_t_dtor(&flds_cache);
  }
}

