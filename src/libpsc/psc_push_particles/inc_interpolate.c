
#include "psc_debug.h"

// ----------------------------------------------------------------------
// interpolation

struct ip_coeff {
#if ORDER == ORDER_1ST
  particle_real_t v0, v1;
#elif ORDER == ORDER_2ND || ORDER == ORDER_1P5
  particle_real_t vm, v0, vp, h;
#endif
};

#if ORDER == ORDER_1ST

#if DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.v0*(gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3  ) +		\
	     gx##x.v1*_F3(flds, m, l##gx##1+1,0,l##gz##3  )) +		\
   gz##z.v1*(gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3+1) +		\
	     gx##x.v1*_F3(flds, m, l##gx##1+1,0,l##gz##3+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.v0*(gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3  ) +		\
	     gy##y.v1*_F3(flds, m, 0,l##gy##2+1,l##gz##3  )) +		\
   gz##z.v1*(gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3+1) +		\
	     gy##y.v1*_F3(flds, m, 0,l##gy##2+1,l##gz##3+1)))
#define INTERPOLATE_FIELD_1ST(flds, m, gy, gz)				\
  (gz##0z*(gy##0y*F3_CACHE(flds, m, 0,l##gy[1]  ,l##gz[2]  ) +		\
	   gy##1y*F3_CACHE(flds, m, 0,l##gy[1]+1,l##gz[2]  )) +		\
   gz##1z*(gy##0y*F3_CACHE(flds, m, 0,l##gy[1]  ,l##gz[2]+1) +		\
	   gy##1y*F3_CACHE(flds, m, 0,l##gy[1]+1,l##gz[2]+1)))
#endif

#else

#if DIM == DIM_Y
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gy##y.vm*_F3(flds, m, 0,l##gy##2-1,0) +				\
   gy##y.v0*_F3(flds, m, 0,l##gy##2  ,0) +				\
   gy##y.vp*_F3(flds, m, 0,l##gy##2+1,0))
#elif DIM == DIM_Z
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*_F3(flds, m, 0,0,l##gz##3-1) +				\
   gz##z.v0*_F3(flds, m, 0,0,l##gz##3  ) +				\
   gz##z.vp*_F3(flds, m, 0,0,l##gz##3+1))
#elif DIM == DIM_XY
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,0) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,0) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,0)) +		\
   gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,0) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,0) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,0)) +		\
   gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,0) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,0) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,0)))
#elif DIM == DIM_XZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,0,l##gz##3-1) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3-1) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,0,l##gz##3-1)) +		\
   gz##z.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,0,l##gz##3  ) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3  ) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,0,l##gz##3  )) +		\
   gz##z.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,0,l##gz##3+1) +		\
	     gx##x.v0*_F3(flds, m, l##gx##1  ,0,l##gz##3+1) +		\
	     gx##x.vp*_F3(flds, m, l##gx##1+1,0,l##gz##3+1)))
#elif DIM == DIM_YZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*(gy##y.vm*_F3(flds, m, 0,l##gy##2-1,l##gz##3-1) +		\
	     gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3-1) +		\
	     gy##y.vp*_F3(flds, m, 0,l##gy##2+1,l##gz##3-1)) +		\
   gz##z.v0*(gy##y.vm*_F3(flds, m, 0,l##gy##2-1,l##gz##3  ) +		\
	     gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3  ) +		\
	     gy##y.vp*_F3(flds, m, 0,l##gy##2+1,l##gz##3  )) +		\
   gz##z.vp*(gy##y.vm*_F3(flds, m, 0,l##gy##2-1,l##gz##3+1) +		\
	     gy##y.v0*_F3(flds, m, 0,l##gy##2  ,l##gz##3+1) +		\
	     gy##y.vp*_F3(flds, m, 0,l##gy##2+1,l##gz##3+1)))
#elif DIM == DIM_XYZ
#define IP_FIELD(flds, m, gx, gy, gz)					\
  (gz##z.vm*(gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3-1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3-1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3-1)) + \
	     gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3-1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3-1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3-1)) + \
	     gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3-1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3-1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3-1))) + \
   gz##z.v0*(gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3  ) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3  ) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3  )) + \
	     gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3  ) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3  ) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3  )) + \
	     gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3  ) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3  ) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3  ))) + \
   gz##z.vp*(gy##y.vm*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2-1,l##gz##3+1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2-1,l##gz##3+1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2-1,l##gz##3+1)) + \
	     gy##y.v0*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2  ,l##gz##3+1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2  ,l##gz##3+1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2  ,l##gz##3+1)) + \
	     gy##y.vp*(gx##x.vm*_F3(flds, m, l##gx##1-1,l##gy##2+1,l##gz##3+1) + \
		       gx##x.v0*_F3(flds, m, l##gx##1  ,l##gy##2+1,l##gz##3+1) + \
		       gx##x.vp*_F3(flds, m, l##gx##1+1,l##gy##2+1,l##gz##3+1))))

#endif

#endif // ORDER


#if ORDER == ORDER_1ST

#if DIM == DIM_YZ

#if IP_VARIANT == IP_VARIANT_EC
#define INTERPOLATE_1ST(flds, E, H)					\
  do {									\
    particle_real_t g0y = 1.f - og[1];					\
    particle_real_t g0z = 1.f - og[2];					\
    particle_real_t g1y = og[1];					\
    particle_real_t g1z = og[2];					\
    									\
    E[0] = (g0z*(g0y*F3_CACHE(flds, EX, 0,lg[1]  ,lg[2]  ) +		\
		 g1y*F3_CACHE(flds, EX, 0,lg[1]+1,lg[2]  )) +		\
	    g1z*(g0y*F3_CACHE(flds, EX, 0,lg[1]  ,lg[2]+1) +		\
		 g1y*F3_CACHE(flds, EX, 0,lg[1]+1,lg[2]+1)));		\
    E[1] = (g0z*F3_CACHE(flds, EY, 0,lg[1]  ,lg[2]  ) +			\
	    g1z*F3_CACHE(flds, EY, 0,lg[1]  ,lg[2]+1));			\
    E[2] = (g0y*F3_CACHE(flds, EZ, 0,lg[1]  ,lg[2]  ) +			\
	    g1y*F3_CACHE(flds, EZ, 0,lg[1]+1,lg[2]  ));			\
									\
    H[0] = F3_CACHE(flds, HX, 0,lg[1]  ,lg[2]  );			\
    H[1] = (g0y*F3_CACHE(flds, HY, 0,lg[1]  ,lg[2]  ) +			\
	    g1y*F3_CACHE(flds, HY, 0,lg[1]+1,lg[2]  ));			\
    H[2] = (g0z*F3_CACHE(flds, HZ, 0,lg[1]  ,lg[2]  ) +			\
	    g1z*F3_CACHE(flds, HZ, 0,lg[1]  ,lg[2]+1));			\
  } while (0)
#else
#define INTERPOLATE_1ST(flds, E, H)					\
  do {									\
    particle_real_t g0y = 1.f - og[1];					\
    particle_real_t g0z = 1.f - og[2];					\
    particle_real_t g1y = og[1];					\
    particle_real_t g1z = og[2];					\
    									\
    particle_real_t h0y = 1.f - oh[1];					\
    particle_real_t h0z = 1.f - oh[2];					\
    particle_real_t h1y = oh[1];					\
    particle_real_t h1z = oh[2];					\
    									\
    E[0] = INTERPOLATE_FIELD_1ST(flds, EX, g, g);			\
    E[1] = INTERPOLATE_FIELD_1ST(flds, EY, h, g);			\
    E[2] = INTERPOLATE_FIELD_1ST(flds, EZ, g, h);			\
    									\
    H[0] = INTERPOLATE_FIELD_1ST(flds, HX, h, h);			\
    H[1] = INTERPOLATE_FIELD_1ST(flds, HY, g, h);			\
    H[2] = INTERPOLATE_FIELD_1ST(flds, HZ, h, g);			\
  } while (0)
#endif

#elif DIM == DIM_XYZ

#if IP_VARIANT == IP_VARIANT_EC
#define INTERPOLATE_1ST(flds, E, H)					\
  do {									\
    particle_real_t g0x = 1.f - og[0];					\
    particle_real_t g0y = 1.f - og[1];					\
    particle_real_t g0z = 1.f - og[2];					\
    particle_real_t g1x = og[0];					\
    particle_real_t g1y = og[1];					\
    particle_real_t g1z = og[2];					\
    									\
    E[0] = (g0z*(g0y*F3_CACHE(flds, EX, lg[0]  ,lg[1]  ,lg[2]  ) +	\
		 g1y*F3_CACHE(flds, EX, lg[0]  ,lg[1]+1,lg[2]  )) +	\
	    g1z*(g0y*F3_CACHE(flds, EX, lg[0]  ,lg[1]  ,lg[2]+1) +	\
		 g1y*F3_CACHE(flds, EX, lg[0]  ,lg[1]+1,lg[2]+1)));	\
									\
    E[1] = (g0x*(g0z*F3_CACHE(flds, EY, lg[0]  ,lg[1]  ,lg[2]  ) +	\
		 g1z*F3_CACHE(flds, EY, lg[0]  ,lg[1]  ,lg[2]+1)) +	\
	    g1x*(g0z*F3_CACHE(flds, EY, lg[0]+1,lg[1]  ,lg[2]  ) +	\
		 g1z*F3_CACHE(flds, EY, lg[0]+1,lg[1]  ,lg[2]+1)));	\
    									\
    E[2] = (g0y*(g0x*F3_CACHE(flds, EZ, lg[0]  ,lg[1]  ,lg[2]  ) +	\
		 g1x*F3_CACHE(flds, EZ, lg[0]+1,lg[1]  ,lg[2]  )) +	\
	    g1y*(g0x*F3_CACHE(flds, EZ, lg[0]  ,lg[1]+1,lg[2]  ) +	\
		 g1x*F3_CACHE(flds, EZ, lg[0]+1,lg[1]+1,lg[2]  )));	\
									\
    H[0] = (g0x*F3_CACHE(flds, HX, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	    g1x*F3_CACHE(flds, HX, lg[0]+1,lg[1]  ,lg[2]  ));		\
									\
    H[1] = (g0y*F3_CACHE(flds, HY, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	    g1y*F3_CACHE(flds, HY, lg[0]  ,lg[1]+1,lg[2]  ));		\
									\
    H[2] = (g0z*F3_CACHE(flds, HZ, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	    g1z*F3_CACHE(flds, HZ, lg[0]  ,lg[1]  ,lg[2]+1));		\
  } while (0)
#endif

#endif

#endif // ORDER == ORDER_1ST

// ======================================================================
// IP_VARIANT SFF

#if IP_VARIANT == IP_VARIANT_SFF

// FIXME, calculation of f_avg could be done at the level where we do caching, too
// (if that survives, anyway...)

#define IP_VARIANT_SFF_PREP					\
  struct psc_patch *patch = &ppsc->patch[p];			\
  								\
  /* FIXME, eventually no ghost points should be needed (?) */	\
  fields_t flds_avg = fields_t_ctor((int[3]) { -1, 0, -1 },		\
				    (int[3]) { patch->ldims[0] + 2, 1, patch->ldims[2] + 1 },\
				    6);					\
									\
  for (int iz = -1; iz < patch->ldims[2] + 1; iz++) {			\
    for (int ix = -1; ix < patch->ldims[0] + 1; ix++) {			\
      _F3(flds_avg, 0, ix,0,iz) = .5 * (_F3(flds, EX, ix,0,iz) + _F3(flds, EX, ix-1,0,iz)); \
      _F3(flds_avg, 1, ix,0,iz) = _F3(flds, EY, ix,0,iz);		\
      _F3(flds_avg, 2, ix,0,iz) = .5 * (_F3(flds, EZ, ix,0,iz) + _F3(flds, EZ, ix,0,iz-1)); \
      _F3(flds_avg, 3, ix,0,iz) = .5 * (_F3(flds, HX, ix,0,iz) + _F3(flds, HX, ix,0,iz-1)); \
      _F3(flds_avg, 4, ix,0,iz) = .25 * (_F3(flds, HY, ix  ,0,iz) + _F3(flds, HY, ix  ,0,iz-1) + \
					 _F3(flds, HY, ix-1,0,iz) + _F3(flds, HY, ix-1,0,iz-1)); \
      _F3(flds_avg, 5, ix,0,iz) = .5 * (_F3(flds, HZ, ix,0,iz) + _F3(flds, HZ, ix-1,0,iz)); \
    }									\
  }

#define IP_VARIANT_SFF_POST			\
  fields_t_dtor(&flds_avg)

#define INTERPOLATE_FIELDS						\
  /* FIXME, we don't really need h coeffs in this case, either, though	\
   * the compiler may be smart enough to figure that out */		\
  particle_real_t E[3] = { IP_FIELD(flds_avg, EX-EX, g, g, g),		\
			   IP_FIELD(flds_avg, EY-EX, g, g, g),		\
			   IP_FIELD(flds_avg, EZ-EX, g, g, g), };	\
  particle_real_t H[3] = { IP_FIELD(flds_avg, HX-EX, g, g, g),		\
			   IP_FIELD(flds_avg, HY-EX, g, g, g),		\
			   IP_FIELD(flds_avg, HZ-EX, g, g, g), }
#else

#define IP_VARIANT_SFF_PREP do {} while (0)
#define IP_VARIANT_SFF_POST do {} while (0)
#define INTERPOLATE_FIELDS						\
  particle_real_t E[3] = { IP_FIELD(flds, EX, h, g, g),			\
			   IP_FIELD(flds, EY, g, h, g),			\
			   IP_FIELD(flds, EZ, g, g, h), };		\
  particle_real_t H[3] = { IP_FIELD(flds, HX, g, h, h),			\
			   IP_FIELD(flds, HY, h, g, h),			\
			   IP_FIELD(flds, HZ, h, h, g), }

#endif

