
#include "psc_debug.h"

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

#if DIM == DIM_YZ

#define INTERPOLATE_FIELD_1ST(pf, m, gy, gz)				\
  (gz##0z*(gy##0y*F3_CACHE(pf, m, 0,l##gy[1]  ,l##gz[2]  ) +		\
	   gy##1y*F3_CACHE(pf, m, 0,l##gy[1]+1,l##gz[2]  )) +		\
   gz##1z*(gy##0y*F3_CACHE(pf, m, 0,l##gy[1]  ,l##gz[2]+1) +		\
	   gy##1y*F3_CACHE(pf, m, 0,l##gy[1]+1,l##gz[2]+1)))

#define INTERPOLATE_1ST_STD(pf, exq, eyq, ezq, hxq, hyq, hzq)		\
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
    exq = INTERPOLATE_FIELD_1ST(pf, EX, g, g);				\
    eyq = INTERPOLATE_FIELD_1ST(pf, EY, h, g);				\
    ezq = INTERPOLATE_FIELD_1ST(pf, EZ, g, h);				\
    									\
    hxq = INTERPOLATE_FIELD_1ST(pf, HX, h, h);				\
    hyq = INTERPOLATE_FIELD_1ST(pf, HY, g, h);				\
    hzq = INTERPOLATE_FIELD_1ST(pf, HZ, h, g);				\
    									\
    assert_finite(exq); assert_finite(eyq); assert_finite(ezq);		\
    assert_finite(hxq); assert_finite(hyq); assert_finite(hzq);		\
  } while (0)

#endif

#if DIM == DIM_YZ

#define INTERPOLATE_1ST_EC(pf, exq, eyq, ezq, hxq, hyq, hzq)        	\
  do {									\
    particle_real_t g0y = 1.f - og[1];					\
    particle_real_t g0z = 1.f - og[2];					\
    particle_real_t g1y = og[1];					\
    particle_real_t g1z = og[2];					\
    									\
    exq = (g0z*(g0y*F3_CACHE(pf, EX, 0,lg[1]  ,lg[2]  ) +		\
		g1y*F3_CACHE(pf, EX, 0,lg[1]+1,lg[2]  )) +		\
	   g1z*(g0y*F3_CACHE(pf, EX, 0,lg[1]  ,lg[2]+1) +		\
		g1y*F3_CACHE(pf, EX, 0,lg[1]+1,lg[2]+1)));		\
    eyq = (g0z*F3_CACHE(pf, EY, 0,lg[1]  ,lg[2]  ) +			\
	   g1z*F3_CACHE(pf, EY, 0,lg[1]  ,lg[2]+1));			\
    ezq = (g0y*F3_CACHE(pf, EZ, 0,lg[1]  ,lg[2]  ) +			\
	   g1y*F3_CACHE(pf, EZ, 0,lg[1]+1,lg[2]  ));			\
									\
    hxq = F3_CACHE(pf, HX, 0,lg[1]  ,lg[2]  );				\
    hyq = (g0y*F3_CACHE(pf, HY, 0,lg[1]  ,lg[2]  ) +			\
	   g1y*F3_CACHE(pf, HY, 0,lg[1]+1,lg[2]  ));			\
    hzq = (g0z*F3_CACHE(pf, HZ, 0,lg[1]  ,lg[2]  ) +			\
	   g1z*F3_CACHE(pf, HZ, 0,lg[1]  ,lg[2]+1));			\
    									\
    assert_finite(exq); assert_finite(eyq); assert_finite(ezq);		\
    assert_finite(hxq); assert_finite(hyq); assert_finite(hzq);		\
  } while (0)

#elif DIM == DIM_XYZ

#define INTERPOLATE_1ST_EC(pf, exq, eyq, ezq, hxq, hyq, hzq)		\
  do {									\
    particle_real_t g0x = 1.f - og[0];					\
    particle_real_t g0y = 1.f - og[1];					\
    particle_real_t g0z = 1.f - og[2];					\
    particle_real_t g1x = og[0];					\
    particle_real_t g1y = og[1];					\
    particle_real_t g1z = og[2];					\
    									\
    exq = (g0z*(g0y*F3_CACHE(pf, EX, lg[0]  ,lg[1]  ,lg[2]  ) +		\
		g1y*F3_CACHE(pf, EX, lg[0]  ,lg[1]+1,lg[2]  )) +	\
	   g1z*(g0y*F3_CACHE(pf, EX, lg[0]  ,lg[1]  ,lg[2]+1) +		\
		g1y*F3_CACHE(pf, EX, lg[0]  ,lg[1]+1,lg[2]+1)));	\
									\
    eyq = (g0x*(g0z*F3_CACHE(pf, EY, lg[0]  ,lg[1]  ,lg[2]  ) +		\
		g1z*F3_CACHE(pf, EY, lg[0]  ,lg[1]  ,lg[2]+1)) +	\
	   g1x*(g0z*F3_CACHE(pf, EY, lg[0]+1,lg[1]  ,lg[2]  ) +		\
		g1z*F3_CACHE(pf, EY, lg[0]+1,lg[1]  ,lg[2]+1)));	\
									\
    ezq = (g0y*(g0x*F3_CACHE(pf, EZ, lg[0]  ,lg[1]  ,lg[2]  ) +		\
		g1x*F3_CACHE(pf, EZ, lg[0]+1,lg[1]  ,lg[2]  )) +	\
	   g1y*(g0x*F3_CACHE(pf, EZ, lg[0]  ,lg[1]+1,lg[2]  ) +		\
		g1x*F3_CACHE(pf, EZ, lg[0]+1,lg[1]+1,lg[2]  )));	\
									\
    hxq = (g0x*F3_CACHE(pf, HX, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	   g1x*F3_CACHE(pf, HX, lg[0]+1,lg[1]  ,lg[2]  ));		\
									\
    hyq = (g0y*F3_CACHE(pf, HY, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	   g1y*F3_CACHE(pf, HY, lg[0]  ,lg[1]+1,lg[2]  ));		\
									\
    hzq = (g0z*F3_CACHE(pf, HZ, lg[0]  ,lg[1]  ,lg[2]  ) +		\
	   g1z*F3_CACHE(pf, HZ, lg[0]  ,lg[1]  ,lg[2]+1));		\
									\
    assert_finite(exq); assert_finite(eyq); assert_finite(ezq);		\
    assert_finite(hxq); assert_finite(hyq); assert_finite(hzq);		\
  } while (0)

#endif
