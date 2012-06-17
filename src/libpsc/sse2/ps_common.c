
#define SIMD_SIZE (1 << SIMD_BITS)

#if SIMD_BITS == 0
#include "simd_0.h"
#elif SIMD_BITS == 2
#include "simd_2.h"
#endif

// ======================================================================

#if SIMD_BITS == 0

#define PRT_LOAD_X(prt, xi, yi, zi, qw) do {	\
    xi = prt[0];				\
    yi = prt[1];				\
    zi = prt[2];				\
    qw = prt[3];				\
  } while (0)

#define PRT_LOAD2_X(sngl, n, xi, yi, zi, qw) do {			\
    float * restrict prt = (float *) &sngl->particles_alt[sngl->b_ids[n]]; \
    xi = prt[0];							\
    yi = prt[1];							\
    zi = prt[2];							\
    qw = prt[3];							\
  } while (0)

#define PRT_STORE_X(prt, xi, yi, zi, qw) do {		\
    prt[0] = xi;					\
    prt[1] = yi;					\
    prt[2] = zi;					\
    prt[3] = qw;					\
  } while (0)

#define PRT_LOAD_P(prt, pxi, pyi, pzi, qdm) do {	\
    pxi = prt[4];					\
    pyi = prt[5];					\
    pzi = prt[6];					\
    qdm = prt[7];					\
  } while (0)

#define PRT_LOAD2_P(sngl, n, pxi, pyi, pzi, qdm) do {		\
    float * restrict prt = (float *) &sngl->particles_alt[sngl->b_ids[n]];	\
    pxi = prt[4];						\
    pyi = prt[5];						\
    pzi = prt[6];						\
    qdm = prt[7];						\
  } while (0)

#define PRT_STORE_P(prt, pxi, pyi, pzi, qdm) do {	\
    prt[4] = pxi;					\
    prt[5] = pyi;					\
    prt[6] = pzi;					\
    prt[7] = qdm;					\
  } while (0)

#elif SIMD_BITS == 2

#define PRT_LOAD_X(prt, xi, yi, zi, qw) do {	\
    xi = v4s_load(prt);				\
    yi = v4s_load(prt + 8);			\
    zi = v4s_load(prt + 16);			\
    qw = v4s_load(prt + 24);			\
    _MM_TRANSPOSE4_PS(xi, yi, zi, qw);		\
  } while (0)

#define PRT_LOAD2_X(sngl, n, xi, yi, zi, qw) do {			\
    v4si ids = v4si_load((int *) &sngl->b_ids[n]);			\
    float * restrict prt0 = (float *) &sngl->particles_alt[v4si_extract(ids, 0)]; \
    float * restrict prt1 = (float *) &sngl->particles_alt[v4si_extract(ids, 1)]; \
    float * restrict prt2 = (float *) &sngl->particles_alt[v4si_extract(ids, 2)]; \
    float * restrict prt3 = (float *) &sngl->particles_alt[v4si_extract(ids, 3)]; \
    xi = v4s_load(prt0);						\
    yi = v4s_load(prt1);						\
    zi = v4s_load(prt2);						\
    qw = v4s_load(prt3);						\
    _MM_TRANSPOSE4_PS(xi, yi, zi, qw);					\
  } while (0)

#define PRT_STORE_X(prt, xi, yi, zi, qw) do {	\
    _MM_TRANSPOSE4_PS(xi, yi, zi, qw);		\
    v4s_store(prt     , xi);			\
    v4s_store(prt + 8 , yi);			\
    v4s_store(prt + 16, zi);			\
    v4s_store(prt + 24, qw);			\
  } while (0)

#define PRT_LOAD_P(prt, pxi, pyi, pzi, qdm) do {	\
    pxi = v4s_load(prt + 4);				\
    pyi = v4s_load(prt + 12);				\
    pzi = v4s_load(prt + 20);				\
    qdm = v4s_load(prt + 28);				\
    _MM_TRANSPOSE4_PS(pxi, pyi, pzi, qdm);		\
  } while (0)

#define PRT_LOAD2_P(sngl, n, pxi, pyi, pzi, qdm) do {		\
    v4si ids = v4si_load((int *) &sngl->b_ids[n]);			\
    float * restrict prt0 = (float *) &sngl->particles_alt[v4si_extract(ids, 0)]; \
    float * restrict prt1 = (float *) &sngl->particles_alt[v4si_extract(ids, 1)]; \
    float * restrict prt2 = (float *) &sngl->particles_alt[v4si_extract(ids, 2)]; \
    float * restrict prt3 = (float *) &sngl->particles_alt[v4si_extract(ids, 3)]; \
    pxi = v4s_load(prt0 + 4);					\
    pyi = v4s_load(prt1 + 4);					\
    pzi = v4s_load(prt2 + 4);					\
    qdm = v4s_load(prt3 + 4);					\
    _MM_TRANSPOSE4_PS(pxi, pyi, pzi, qdm);			\
  } while (0)

#define PRT_STORE_P(prt, pxi, pyi, pzi, qdm) do {	\
    _MM_TRANSPOSE4_PS(pxi, pyi, pzi, qdm);		\
    v4s_store(prt + 4  , pxi);				\
    v4s_store(prt + 12, pyi);				\
    v4s_store(prt + 20, pzi);				\
    v4s_store(prt + 28, qdm);				\
  } while (0)

#endif

// ======================================================================

#define __GET_IDX_FRAC_TRUNC(posy, lgy, ogy) do {		\
    lgy = v4s_fint(posy);					\
    ogy = posy - v4si_to_v4s(lgy);				\
  } while (0)

// OPT. try cvtps (rounding) instead

#define __GET_IDX_FRAC_ROUND(posy, lgy, ogy) do {		\
    v4s posy_floor = v4s_floor(posy);				\
    lgy = v4s_to_v4si(posy_floor);				\
    ogy = posy - posy_floor;					\
  } while (0)


#define __GET_IDX_FRAC __GET_IDX_FRAC_ROUND

#define GET_IDX_FRAC(yi, dyi, shift, lgy, ogy) do {	\
    v4s posy = yi * v4s_splat(dyi) + v4s_splat(shift);	\
    __GET_IDX_FRAC(posy, lgy, ogy);			\
  } while (0)

#define GET_POS_IDX_FRAC(yi, dyi, shift, posy, lgy, ogy) do {	\
    posy = yi * v4s_splat(dyi) + v4s_splat(shift);		\
    __GET_IDX_FRAC(posy, lgy, ogy);				\
  } while (0)

#define GET_POS_IDX(yi, dyi, shift, posy, lgy) do {	\
    posy = yi * v4s_splat(dyi) + v4s_splat(shift);	\
    v4s posy_floor = v4s_floor(posy);			\
    lgy = v4s_to_v4si(posy_floor);			\
  } while (0)

// ======================================================================

static inline v4si
calc_off(v4si ly, v4si lz, int my)
{
  // lz * my + ly
  return v4si_add(v4si_mul(lz, v4si_splat(my)), ly);
}

static inline v4si
calc_off16(v4si ly, v4si lz, int my)
{
  // 16*(lz * my + ly)
  return v4si_sll(v4si_add(v4si_mul(lz, v4si_splat(my)), ly), 4);
}

#if SIMD_BITS == 0

#define INTERPOLATE_FIELD_1ST_IP(mm, gy, gz) \
    ({									\
      v4si off;								\
      off = calc_off(l##gy##y, l##gz##z, my);				\
      fields_ip_real_t * __restrict pp = _p0 + mm * sz;			\
      (pp[off].f00  +							\
       o##gy##y * pp[off].f10 +						\
       o##gz##z * pp[off].f01 +						\
       o##gz##z * o##gy##y * pp[off].f11);				\
    })

#else

#define INTERPOLATE_FIELD_1ST_IP(mm, gy, gz) \
    ({									\
      v4si off;								\
      off = calc_off16(l##gy##y, l##gz##z, my);				\
      char * __restrict pp = (char *) (_p0 + mm * sz);			\
      v4s f00 = v4s_load((float *) (pp + v4si_extract(off, 0)));	\
      v4s f01 = v4s_load((float *) (pp + v4si_extract(off, 1)));	\
      v4s f10 = v4s_load((float *) (pp + v4si_extract(off, 2)));	\
      v4s f11 = v4s_load((float *) (pp + v4si_extract(off, 3)));	\
      _MM_TRANSPOSE4_PS(f00, f01, f10, f11);				\
      (f00  +								\
       o##gy##y * f10 +							\
       o##gz##z * f01 +							\
       o##gz##z * o##gy##y * f11);					\
    })

#endif

// ======================================================================

#define _PUSH_PXI(pxi, pyi, pzi, exq, eyq, ezq, hxq, hyq, hzq, dq) do { \
    v4s pxm = pxi + dq * exq;						\
    v4s pym = pyi + dq * eyq;						\
    v4s pzm = pzi + dq * ezq;						\
    									\
    v4s root = dq * v4s_rsqrt(vOne + sqr(pxm) + sqr(pym) + sqr(pzm));	\
    									\
    v4s taux = hxq * root;						\
    v4s tauy = hyq * root;						\
    v4s tauz = hzq * root;						\
    									\
    v4s tau = vOne / (vOne + sqr(taux) + sqr(tauy) + sqr(tauz));	\
    									\
    v4s pxp = ((vOne+taux*taux-tauy*tauy-tauz*tauz)*pxm +		\
	       (vTwo*taux*tauy+vTwo*tauz)*pym +				\
	       (vTwo*taux*tauz-vTwo*tauy)*pzm)*tau;			\
    v4s pyp = ((vTwo*taux*tauy-vTwo*tauz)*pxm +				\
	       (vOne-taux*taux+tauy*tauy-tauz*tauz)*pym +		\
	       (vTwo*tauy*tauz+vTwo*taux)*pzm)*tau;			\
    v4s pzp = ((vTwo*taux*tauz+vTwo*tauy)*pxm +				\
	       (vTwo*tauy*tauz-vTwo*taux)*pym +				\
	       (vOne-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;		\
    									\
    pxi = pxp + dq * exq;						\
    pyi = pyp + dq * eyq;						\
    pzi = pzp + dq * ezq;						\
  } while (0)

// ======================================================================

#define _GET_VYZ(vyi, vzi, pxi, pyi, pzi) do {				\
    v4s root = v4s_rsqrt(vOne +	sqr(pxi) + sqr(pyi) + sqr(pzi));	\
    vyi = pyi * root;							\
    vzi = pzi * root;							\
  } while(0)

#define _GET_VXYZ(vxi, vyi, vzi, pxi, pyi, pzi) do {			\
    v4s root = v4s_rsqrt(vOne + sqr(pxi) + sqr(pyi) + sqr(pzi));	\
    vxi = pxi * root;							\
    vyi = pyi * root;							\
    vzi = pzi * root;							\
  } while(0)

// ======================================================================

#if SIMD_BITS == 0

#define CURR_JX(p_jx, my, lfy, lfz, ofy, ofz, fnqx) do {	\
    v4s val00 = (vOne - ofy) * (vOne - ofz) * fnqx;		\
    v4s valp0 = (ofy       ) * (vOne - ofz) * fnqx;		\
    v4s val0p = (vOne - ofy) * (       ofz) * fnqx;		\
    v4s valpp = (ofy       ) * (       ofz) * fnqx;		\
    								\
    v4si off = v4si_add(v4si_mul(lfz, v4si_splat(my)), lfy);	\
    fields_ip_real_t *f = (fields_ip_real_t *) p_jx;		\
    f[off].f00 += val00;					\
    f[off].f01 += val0p;					\
    f[off].f10 += valp0;					\
    f[off].f11 += valpp;					\
  } while (0)

#elif SIMD_BITS == 2

#define CURR_JX(p_jx, my, lfy, lfz, ofy, ofz, fnqx) do {	\
    v4s val00 = (vOne - ofy) * (vOne - ofz) * fnqx;		\
    v4s valp0 = (ofy       ) * (vOne - ofz) * fnqx;		\
    v4s val0p = (vOne - ofy) * (       ofz) * fnqx;		\
    v4s valpp = (ofy       ) * (       ofz) * fnqx;		\
    								\
    v4si off = v4si_add(v4si_mul(lfz, v4si_splat(my)), lfy);	\
    _MM_TRANSPOSE4_PS(val00, val0p, valp0, valpp);		\
    p_jx[v4si_extract(off, 0)] += val00;			\
    p_jx[v4si_extract(off, 1)] += val0p;			\
    p_jx[v4si_extract(off, 2)] += valp0;			\
    p_jx[v4si_extract(off, 3)] += valpp;			\
  } while (0)

#endif

// ======================================================================

#define CALC_DX1(dx1, x, dx, off) do {					\
    v4si mask0, mask01;							\
    mask0 = v4si_cmpeq(off[0], iZero);					\
    mask01 = v4si_and(mask0, v4si_cmpeq(off[1], iZero));		\
    v4s v0, v1, vv0, vv1, o0, x0, dx_0, dx_1;				\
    o0 = v4s_blend(mask0, v4si_to_v4s(off[1]), v4si_to_v4s(off[0]));	\
    x0 = v4s_blend(mask0, x[1], x[0]);					\
    dx_0 = v4s_blend(mask0, dx[1], dx[0]);				\
    mask01 = v4si_or(mask01, v4s_cmpeq(dx_0, vZero));			\
    dx_1 = v4s_blend(mask0, dx[0], dx[1]);				\
    v0 = v4s_splat(.5f) * o0 - x0;					\
    v1 = dx_1 / dx_0 * v0;						\
    vv0 = v4s_blend(mask01, vZero, v4s_blend(mask0, v1, v0));		\
    vv1 = v4s_blend(mask01, vZero, v4s_blend(mask0, v0, v1));		\
    dx1[0] = vv0;							\
    dx1[1] = vv1;							\
  } while (0)

#if SIMD_BITS == 0

#define CURR_VB_CELL_(i, x, dx, fnq) do {				\
    v4s jy0 = fnq[0] * dx[0] * (v4s_splat(.5f) - x[1] - v4s_splat(.5f) * dx[1]); \
    v4s jyp = fnq[0] * dx[0] * (v4s_splat(.5f) + x[1] + v4s_splat(.5f) * dx[1]); \
    v4s jz0 = fnq[1] * dx[1] * (v4s_splat(.5f) - x[0] - v4s_splat(.5f) * dx[0]); \
    v4s jzp = fnq[1] * dx[1] * (v4s_splat(.5f) + x[0] + v4s_splat(.5f) * dx[0]); \
    v4si off = v4si_add(v4si_mul(i[1], v4si_splat(my)), i[0]);		\
    fields_ip_real_t *f = (fields_ip_real_t *) _p_jyz;			\
    f[off].f00 += jy0;							\
    f[off].f01 += jyp;							\
    f[off].f10 += jz0;							\
    f[off].f11 += jzp;							\
  } while (0)

#elif SIMD_BITS == 2

#define CURR_VB_CELL_(i, x, dx, fnq) do {				\
    v4s jy0 = fnq[0] * dx[0] * (v4s_splat(.5f) - x[1] - v4s_splat(.5f) * dx[1]); \
    v4s jyp = fnq[0] * dx[0] * (v4s_splat(.5f) + x[1] + v4s_splat(.5f) * dx[1]); \
    v4s jz0 = fnq[1] * dx[1] * (v4s_splat(.5f) - x[0] - v4s_splat(.5f) * dx[0]); \
    v4s jzp = fnq[1] * dx[1] * (v4s_splat(.5f) + x[0] + v4s_splat(.5f) * dx[0]); \
    v4si off = v4si_add(v4si_mul(i[1], v4si_splat(my)), i[0]);		\
    _MM_TRANSPOSE4_PS(jy0, jyp, jz0, jzp);				\
    _p_jyz[v4si_extract(off, 0)] += jy0;				\
    _p_jyz[v4si_extract(off, 1)] += jyp;				\
    _p_jyz[v4si_extract(off, 2)] += jz0;				\
    _p_jyz[v4si_extract(off, 3)] += jzp;				\
  } while (0)

#endif

#define CURR_VB_UPD(i, x, dx1, dx, off) do {				\
    dx[0] -= dx1[0];							\
    dx[1] -= dx1[1];							\
    x[0] += dx1[0] - v4si_to_v4s(off[0]);				\
    x[1] += dx1[1] - v4si_to_v4s(off[1]);				\
    i[0] = v4si_add(i[0], off[0]);					\
    i[1] = v4si_add(i[1], off[1]);					\
  } while (0)

#define CURR_JYZ_A(lgy, lgz, lfy, lfz, ogy, ogz, ym, zm, yp, zp, off0) do { \
    /* IN PLANE CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt */	\
    i[0] = lgy;								\
    i[1] = lgz;								\
    idiff[0] = v4si_sub(lfy, lgy);					\
    idiff[1] = v4si_sub(lfz, lgz);					\
    dx[0] = yp - ym;							\
    dx[1] = zp - zm;							\
    x[0] = ogy - v4s_splat(.5f);					\
    x[1] = ogz - v4s_splat(.5f);					\
    									\
    v4s x0 = v4s_xor(x[0], v4s_and(v4s_cast(idiff[0]), v4s_signmask));	\
    v4s x1 = v4s_xor(x[1], v4s_and(v4s_cast(idiff[1]), v4s_signmask));	\
    /*v4s x1 = x[1] * v4si_to_v4s(idiff[1]);*/				\
    v4si d_first = v4s_cmpge(v4s_abs(dx[1]) * (v4s_splat(.5f) - x0),	\
			     v4s_abs(dx[0]) * (v4s_splat(.5f) - x1));	\
    									\
    off0[0] = v4si_andnot(d_first, idiff[0]);				\
    off0[1] = v4si_and(d_first, idiff[1]);				\
  } while (0)

#define CURR_JYZ_(_p_jy, _p_jz, off0) do {	\
    v4s dx0[2];					\
    CALC_DX1(dx0, x, dx, off0);			\
    CURR_VB_CELL_(i, x, dx0, fnq);		\
    CURR_VB_UPD(i, x, dx0, dx, off0);		\
    						\
    v4si off1[2];				\
    off1[0] = v4si_sub(idiff[0], off0[0]);	\
    off1[1] = v4si_sub(idiff[1], off0[1]);	\
    						\
    v4s dx1[2];					\
    CALC_DX1(dx1, x, dx, off1);			\
    CURR_VB_CELL_(i, x, dx1, fnq);		\
    CURR_VB_UPD(i, x, dx1, dx, off1);		\
    						\
    CURR_VB_CELL_(i, x, dx, fnq);		\
  } while (0)

// ======================================================================


