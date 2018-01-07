
#include "inc_defs.h"

#define PRTS_STAGGERED 1

#define CACHE_EM_J 1

#define IP_DEPOSIT

#include "fields.hxx"
#include "inc_params.c"
#include "inc_push.c"
#include "inc_cache.c"
#include "inc_interpolate.c"

// ----------------------------------------------------------------------

#if ORDER == ORDER_1ST

#if DIM == DIM_XZ
#if IP_VARIANT == IP_VARIANT_SFF
#define SFX(x) x ## _1sff_xz
#define psc_push_particles_push_mprts
#define do_push_part do_push_part_1sff_xz
#define PROF_NAME "push_mprts_1sff_xz"
#else
#define SFX(x) x ## _1st_xz
#define psc_push_particles_push_mprts
#define do_push_part do_push_part_1st_xz
#define PROF_NAME "push_mprts_1st_xz"
#endif
#elif DIM == DIM_YZ
#define SFX(x) x ## _1st_yz
#define psc_push_particles_push_mprts
#define do_push_part do_push_part_1st_yz
#define PROF_NAME "push_mprts_1st_yz"
#endif

#elif ORDER == ORDER_2ND

#if DIM == DIM_Y
#define psc_push_particles_push_mprts psc_push_particles_generic_c_push_mprts_y
#define do_push_part do_push_part_genc_y
#define PROF_NAME "genc_push_mprts_y"
#elif DIM == DIM_Z
#define psc_push_particles_push_mprts psc_push_particles_generic_c_push_mprts_z
#define do_push_part do_push_part_genc_z
#define PROF_NAME "genc_push_mprts_z"
#elif DIM == DIM_XY
#define psc_push_particles_push_mprts psc_push_particles_generic_c_push_mprts_xy
#define do_push_part do_push_part_genc_xy
#define PROF_NAME "genc_push_mprts_xy"
#elif DIM == DIM_XZ
#define psc_push_particles_push_mprts psc_push_particles_generic_c_push_mprts_xz
#define do_push_part do_push_part_genc_xz
#define PROF_NAME "genc_push_mprts_xz"
#elif DIM == DIM_XYZ
#define psc_push_particles_push_mprts psc_push_particles_generic_c_push_mprts_xyz
#define do_push_part do_push_part_genc_xyz
#define PROF_NAME "genc_push_mprts_xyz"
#endif

#endif

// ----------------------------------------------------------------------
// find_l_minmax

static inline void
find_l_minmax(int *l1min, int *l1max, int k1, int lg1)
{
#if ORDER == ORDER_1ST
  if (k1 == lg1) {
    *l1min = 0; *l1max = +1;
  } else if (k1 == lg1 - 1) {
    *l1min = -1; *l1max = +1;
  } else { // (k1 == lg1 + 1)
    *l1min = 0; *l1max = +2;
  }
#elif ORDER == ORDER_2ND
  if (k1 == lg1) {
    *l1min = -1; *l1max = +1;
  } else if (k1 == lg1 - 1) {
    *l1min = -2; *l1max = +1;
  } else { // (k1 == lg1 + 1)
    *l1min = -1; *l1max = +2;
  }
#endif
}

// ======================================================================
// current

#define CURRENT_PREP_DIM(l1min, l1max, k1, lg1, fnqx, fnqxs)	\
    int l1min, l1max; find_l_minmax(&l1min, &l1max, k1, ip.lg1);	\
    particle_real_t fnqx = particle_qni_wni(part) * c_prm.fnqxs;	\

#define CURRENT_PREP							\
  IF_DIM_X( CURRENT_PREP_DIM(l1min, l1max, k1, lg1, fnqx, fnqxs); );	\
  IF_DIM_Y( CURRENT_PREP_DIM(l2min, l2max, k2, lg2, fnqy, fnqys); );	\
  IF_DIM_Z( CURRENT_PREP_DIM(l3min, l3max, k3, lg3, fnqz, fnqzs); );	\
									\
  IF_NOT_DIM_X( particle_real_t fnqxx = vv[0] * particle_qni_wni(part) * c_prm.fnqs; ); \
  IF_NOT_DIM_Y( particle_real_t fnqyy = vv[1] * particle_qni_wni(part) * c_prm.fnqs; ); \
  IF_NOT_DIM_Z( particle_real_t fnqzz = vv[2] * particle_qni_wni(part) * c_prm.fnqs; )

#define CURRENT_2ND_Y						\
  particle_real_t jyh = 0.f;					\
								\
  for (int l2 = l2min; l2 <= l2max; l2++) {			\
    particle_real_t wx = S(s0y, l2) + .5f * S(s1y, l2);		\
    particle_real_t wy = S(s1y, l2);				\
    particle_real_t wz = S(s0y, l2) + .5f * S(s1y, l2);		\
    								\
    particle_real_t jxh = fnqxx*wx;				\
    jyh -= fnqy*wy;						\
    particle_real_t jzh = fnqzz*wz;				\
    								\
    J(JXI, 0,ip.lg2+l2,0) += jxh;					\
    J(JYI, 0,ip.lg2+l2,0) += jyh;				\
    J(JZI, 0,ip.lg2+l2,0) += jzh;				\
  }

#define CURRENT_2ND_Z						\
  particle_real_t jzh = 0.f;					\
  for (int l3 = l3min; l3 <= l3max; l3++) {			\
    particle_real_t wx = S(s0z, l3) + .5f * S(s1z, l3);		\
    particle_real_t wy = S(s0z, l3) + .5f * S(s1z, l3);		\
    particle_real_t wz = S(s1z, l3);				\
    								\
    particle_real_t jxh = fnqxx*wx;				\
    particle_real_t jyh = fnqyy*wy;				\
    jzh -= fnqz*wz;						\
    								\
    J(JXI, 0,0,ip.lg3+l3) += jxh;				\
    J(JYI, 0,0,ip.lg3+l3) += jyh;				\
    J(JZI, 0,0,ip.lg3+l3) += jzh;				\
  }

#define CURRENT_2ND_XY							\
  for (int l2 = l2min; l2 <= l2max; l2++) {				\
    particle_real_t jxh = 0.f;						\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      particle_real_t wx = S(s1x, l1) * (S(s0y, l2) + .5f*S(s1y, l2));	\
      particle_real_t wz = S(s0x, l1) * S(s0y, l2)			\
	+ .5f * S(s1x, l1) * S(s0y, l2)					\
	+ .5f * S(s0x, l1) * S(s1y, l2)					\
	+ (1.f/3.f) * S(s1x, l1) * S(s1y, l2);				\
      									\
      jxh -= fnqx*wx;							\
      J(JXI, ip.lg1+l1,ip.lg2+l2,0) += jxh;				\
      J(JZI, ip.lg1+l1,ip.lg2+l2,0) += fnqzz * wz;			\
    }									\
  }									\
  for (int l1 = l1min; l1 <= l1max; l1++) {				\
    particle_real_t jyh = 0.f;						\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      particle_real_t wy = S(s1y, l2) * (S(s0x, l1) + .5f*S(s1x, l1));	\
      									\
      jyh -= fnqy*wy;							\
      J(JYI, ip.lg1+l1,ip.lg2+l2,0) += jyh;				\
    }									\
  }

#define CURRENT_XZ				\
  for (int l3 = l3min; l3 <= l3max; l3++) {	\
    particle_real_t jxh = 0.f;						\
    for (int l1 = l1min; l1 < l1max; l1++) {				\
      particle_real_t wx = S(s1x, l1) * (S(s0z, l3) + .5f*S(s1z, l3));	\
      jxh -= fnqx*wx;							\
      J(JXI, ip.lg1+l1,0,ip.lg3+l3) += jxh;				\
    }									\
  }									\
									\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      particle_real_t wy = S(s0x, l1) * S(s0z, l3)			\
	+ .5f * S(s1x, l1) * S(s0z, l3)					\
	+ .5f * S(s0x, l1) * S(s1z, l3)					\
	+ (1.f/3.f) * S(s1x, l1) * S(s1z, l3);				\
      particle_real_t jyh = fnqyy * wy;					\
      J(JYI, ip.lg1+l1,0,ip.lg3+l3) += jyh;				\
    }									\
  }									\
  for (int l1 = l1min; l1 <= l1max; l1++) {				\
    particle_real_t jzh = 0.f;						\
    for (int l3 = l3min; l3 < l3max; l3++) {				\
      particle_real_t wz = S(s1z, l3) * (S(s0x, l1) + .5f*S(s1x, l1));	\
      jzh -= fnqz*wz;							\
      J(JZI, ip.lg1+l1,0,ip.lg3+l3) += jzh;				\
    }									\
  }

#define CURRENT_1ST_YZ							\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      particle_real_t wx = S(s0y, l2) * S(s0z, l3)			\
	+ .5f * S(s1y, l2) * S(s0z, l3)					\
	+ .5f * S(s0y, l2) * S(s1z, l3)					\
	+ (1.f/3.f) * S(s1y, l2) * S(s1z, l3);				\
      particle_real_t jxh = fnqxx * wx;					\
      J(JXI, 0,ip.lg2+l2,ip.lg3+l3) += jxh;				\
    }									\
  }									\
  									\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    particle_real_t jyh = 0.f;						\
    for (int l2 = l2min; l2 < l2max; l2++) {				\
      particle_real_t wy = S(s1y, l2) * (S(s0z, l3) + .5f*S(s1z, l3));	\
      jyh -= fnqy*wy;							\
      J(JYI, 0,ip.lg2+l2,ip.lg3+l3) += jyh;				\
    }									\
  }									\
									\
  for (int l2 = l2min; l2 <= l2max; l2++) {				\
    particle_real_t jzh = 0.f;						\
    for (int l3 = l3min; l3 < l3max; l3++) {				\
      particle_real_t wz = S(s1z, l3) * (S(s0y, l2) + .5f*S(s1y, l2));	\
      jzh -= fnqz*wz;							\
      J(JZI, 0,ip.lg2+l2,ip.lg3+l3) += jzh;				\
    }									\
  }

#define JZH(i) jzh[i+2]
#define CURRENT_2ND_YZ							\
    particle_real_t jxh;						\
    particle_real_t jyh;						\
    particle_real_t jzh[5];						\
									\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      JZH(l2) = 0.f;							\
    }									\
    for (int l3 = l3min; l3 <= l3max; l3++) {				\
      jyh = 0.f;							\
      for (int l2 = l2min; l2 <= l2max; l2++) {				\
	particle_real_t wx = S(s0y, l2) * S(s0z, l3)			\
	  + .5f * S(s1y, l2) * S(s0z, l3)				\
	  + .5f * S(s0y, l2) * S(s1z, l3)				\
	+ (1.f/3.f) * S(s1y, l2) * S(s1z, l3);				\
	particle_real_t wy = S(s1y, l2) * (S(s0z, l3) + .5f*S(s1z, l3)); \
	particle_real_t wz = S(s1z, l3) * (S(s0y, l2) + .5f*S(s1y, l2)); \
									\
	jxh = fnqxx*wx;							\
	jyh -= fnqy*wy;							\
	JZH(l2) -= fnqz*wz;						\
									\
	J(JXI, 0,ip.lg2+l2,ip.lg3+l3) += jxh;				\
	J(JYI, 0,ip.lg2+l2,ip.lg3+l3) += jyh;				\
	J(JZI, 0,ip.lg2+l2,ip.lg3+l3) += JZH(l2);			\
      }									\
    }									\

#define CURRENT_2ND_XYZ							\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l2 = l2min; l2 <= l2max; l2++) {				\
      particle_real_t jxh = 0.f;					\
      for (int l1 = l1min; l1 <= l1max; l1++) {				\
	particle_real_t wx = S(s1x, l1) * (S(s0y, l2) * S(s0z, l3) +	\
					   .5f * S(s1y, l2) * S(s0z, l3) + \
					   .5f * S(s0y, l2) * S(s1z, l3) + \
					   (1.f/3.f) * S(s1y, l2) * S(s1z, l3)); \
									\
	jxh -= fnqx*wx;							\
	J(JXI, ip.lg1+l1,ip.lg2+l2,ip.lg3+l3) += jxh;			\
      }									\
    }									\
  }									\
  									\
  for (int l3 = l3min; l3 <= l3max; l3++) {				\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      particle_real_t jyh = 0.f;					\
      for (int l2 = l2min; l2 <= l2max; l2++) {				\
	particle_real_t wy = S(s1y, l2) * (S(s0x, l1) * S(s0z, l3) +	\
					   .5f * S(s1x, l1) * S(s0z, l3) + \
					   .5f * S(s0x, l1) * S(s1z, l3) + \
					   (1.f/3.f) * S(s1x, l1)*S(s1z, l3)); \
									\
	jyh -= fnqy*wy;							\
	J(JYI, ip.lg1+l1,ip.lg2+l2,ip.lg3+l3) += jyh;			\
      }									\
    }									\
  }									\
									\
  for (int l2 = l2min; l2 <= l2max; l2++) {				\
    for (int l1 = l1min; l1 <= l1max; l1++) {				\
      particle_real_t jzh = 0.f;					\
      for (int l3 = l3min; l3 <= l3max; l3++) {				\
	particle_real_t wz = S(s1z, l3) * (S(s0x, l1) * S(s0y, l2) +	\
					   .5f * S(s1x, l1) * S(s0y, l2) +\
					   .5f * S(s0x, l1) * S(s1y, l2) +\
					   (1.f/3.f) * S(s1x, l1)*S(s1y, l2)); \
									\
	jzh -= fnqz*wz;							\
	J(JZI, ip.lg1+l1,ip.lg2+l2,ip.lg3+l3) += jzh;			\
      }									\
    }									\
  }

#if DIM == DIM_Y
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_Y
#endif
#elif DIM == DIM_Z
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_Z
#endif
#elif DIM == DIM_XY
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_XY
#endif
#elif DIM == DIM_XZ
#define CURRENT CURRENT_XZ
#elif DIM == DIM_YZ

#if ORDER == ORDER_1ST
#define CURRENT CURRENT_1ST_YZ
#elif ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_YZ

#endif
#elif DIM == DIM_XYZ
#if ORDER == ORDER_2ND
#define CURRENT CURRENT_2ND_XYZ
#endif
#endif

#ifdef do_push_part

static void
do_push_part(int p, fields_t flds, particle_range_t prts)
{
  //#if (DIM & DIM_X)
  particle_real_t s0x[N_RHO] = {}, s1x[N_RHO];
  //#endif
  //#if (DIM & DIM_Y)
  particle_real_t s0y[N_RHO] = {}, s1y[N_RHO];
  //#endif
  //#if (DIM & DIM_Z)
  particle_real_t s0z[N_RHO] = {}, s1z[N_RHO];
  //#endif

  c_prm_set(ppsc);

  Fields3d<fields_t> EM(flds); // FIXME, EM and J are identical here
  Fields3d<fields_t> J(flds);

  IP_VARIANT_SFF_PREP;
  
  PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
    particle_t *part = particle_iter_deref(prt_iter);
    particle_real_t *x = &part->xi;
    particle_real_t vv[3];

#if PRTS != PRTS_STAGGERED
    // x^n, p^n -> x^(n+.5), p^n
    calc_v(vv, &part->pxi);
    push_x(x, vv, .5f * c_prm.dt);
#endif
    
    // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt 
    // FIELD INTERPOLATION

    particle_real_t xm[3];
    for (int d = 0; d < 3; d++) {
      xm[d] = x[d] * c_prm.dxi[d];
    }
    IP ip;
    INTERPOLATE_FIELDS(flds_em);

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
    particle_real_t dq = c_prm.dqs * particle_qni_div_mni(part);
    push_p(&part->pxi, ip.E, ip.H, dq);

    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0) 
    calc_v(vv, &part->pxi);

#if PRTS == PRTS_STAGGERED
    // FIXME, inelegant way of pushing full dt
    push_x(x, vv, c_prm.dt);

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    ZERO_S1;
    IF_DIM_X( DEPOSIT(x, k1, ip.gx, 0, c_prm.dxi[0], s1x, ip.lg1); );
    IF_DIM_Y( DEPOSIT(x, k2, ip.gy, 1, c_prm.dxi[1], s1y, ip.lg2); );
    IF_DIM_Z( DEPOSIT(x, k3, ip.gz, 2, c_prm.dxi[2], s1z, ip.lg3); );

#else
    push_x(x, vv, .5f * c_prm.dt);

    // x^(n+1), p^(n+1) -> x^(n+1.5f), p^(n+1)
    particle_real_t xn[3] = { x[0], x[1], x[2] };
    push_x(xn, vv, .5f * c_prm.dt);

    // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
    ZERO_S1;
    IF_DIM_X( DEPOSIT(xn, k1, gx, 0, c_prm.dxi[0], s1x, ip.lg1); );
    IF_DIM_Y( DEPOSIT(xn, k2, gy, 1, c_prm.dxi[1], s1y, ip.lg2); );
    IF_DIM_Z( DEPOSIT(xn, k3, gz, 2, c_prm.dxi[2], s1z, ip.lg3); );
#endif
    
    // CURRENT DENSITY AT (n+1.0)*dt

    SUBTR_S1_S0;
    CURRENT_PREP;
    CURRENT;
  }

  IP_VARIANT_SFF_POST;
}

#endif

#if CACHE == CACHE_EM_J

#if DIM == DIM_YZ
static fields_t
cache_fields_from_em(fields_t flds)
{
  fields_t fld = fields_t(flds.ib, flds.im, 9);
  Fields3d<fields_t> F(flds), CACHE(fld);
  // FIXME, can do -1 .. 2? NO!, except maybe for 1st order
  // Has to be at least -2 .. +3 because of staggering
  // FIXME, get rid of caching since it's no different from the actual
  // fields...
  for (int iz = fld.ib[2]; iz < fld.ib[2] + fld.im[2]; iz++) {
    for (int iy = fld.ib[1]; iy < fld.ib[1] + fld.im[1]; iy++) {
      CACHE(fld, EX, 0,iy,iz) = F(EX, 0,iy,iz);
      CACHE(fld, EY, 0,iy,iz) = F(EY, 0,iy,iz);
      CACHE(fld, EZ, 0,iy,iz) = F(EZ, 0,iy,iz);
      CACHE(fld, HX, 0,iy,iz) = F(HX, 0,iy,iz);
      CACHE(fld, HY, 0,iy,iz) = F(HY, 0,iy,iz);
      CACHE(fld, HZ, 0,iy,iz) = F(HZ, 0,iy,iz);
    }
  }
  return fld;
}

static void
cache_fields_to_j(fields_t fld, fields_t flds)
{
  Fields3d<fields_t> F(flds), CACHE(fld);
  for (int iz = fld.ib[2]; iz < fld.ib[2] + fld.im[2]; iz++) {
    for (int iy = fld.ib[1]; iy < fld.ib[1] + fld.im[1]; iy++) {
      F(JXI, 0,iy,iz) += CACHE(JXI, 0,iy,iz);
      F(JYI, 0,iy,iz) += CACHE(JYI, 0,iy,iz);
      F(JZI, 0,iy,iz) += CACHE(JZI, 0,iy,iz);
    }
  }
}
#endif

#endif

#ifdef psc_push_particles_push_mprts

// ----------------------------------------------------------------------
// psc_push_particles_push_mprts

void
#ifdef SFX
SFX(psc_push_particles_push_mprts)(struct psc_push_particles *push,
				   struct psc_mparticles *mprts,
				   struct psc_mfields *mflds_base)
#else
psc_push_particles_push_mprts(struct psc_push_particles *push,
			      struct psc_mparticles *mprts,
			      struct psc_mfields *mflds_base)
#endif
{
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 6);

  static int pr;
  if (!pr) {
    pr = prof_register(PROF_NAME, 1., 0, 0);
  }
  
  prof_start(pr);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = mf[p];
    particle_range_t prts = particle_range_mprts(mprts, p);
#if CACHE == CACHE_EM_J
    // FIXME, can't we just skip this and just set j when copying back?
    flds.zero(JXI, JXI + 3);
    fields_t flds_cache = cache_fields_from_em(flds);
    do_push_part(p, flds_cache, prts);
    cache_fields_to_j(flds_cache, flds);
    flds_cache.dtor();
#else
    flds.zero(JXI, JXI + 3);
    do_push_part(p, flds, prts);
#endif
  }
  prof_stop(pr);

  mf.put_as(mflds_base, JXI, JXI+3);
}

#endif
