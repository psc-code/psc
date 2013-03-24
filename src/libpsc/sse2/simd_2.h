
#include <smmintrin.h>

#define SSE2_EMU

#ifdef SSE2_EMU

typedef float v4s __attribute__ ((vector_size (16)));
typedef struct { int v[4]; } v4si;

static const v4si iZero = { .v = { 0, 0, 0, 0} };

static inline v4si
v4si_splat(int x)
{
  v4si rv = { .v = { x, x, x, x }};
  return rv;
}

static inline v4si
v4si_load(int *p)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = p[i];
  }
  return rv;
}

static inline void
v4si_store(int *p, v4si v)
{
  for (int i = 0; i < 4; i++) {
    p[i] = v.v[i];
  }
}

static inline v4si
v4si_insert(v4si v, int s, int m)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = v.v[i];
  }
  rv.v[m] = s;
  return rv;
}

static inline int
v4si_extract(v4si v, int m)
{
  return v.v[m];
}

static inline v4si
v4s_fint(v4s v)
{
  v4si rv;
  float *f = (float *) &v;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = ((int) (f[i] + 10)) - 10;
  }
  return rv;
}

static inline v4s
v4si_to_v4s(v4si v)
{
  v4s rv;
  float *f = (float *) &rv;
  for (int i = 0; i < 4; i++) {
    f[i] = v.v[i];
  }
  return rv;
}

static inline v4si
v4s_cmpge(v4s a, v4s b)
{
  v4si rv;
  float *fa = (float *) &a, *fb = (float *) &b;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (fa[i] >= fb[i] ? 0xffffffff : 0x0);
  }
  return rv;
}

static inline v4si
v4s_cmpeq(v4s a, v4s b)
{
  v4si rv;
  float *fa = (float *) &a, *fb = (float *) &b;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (fa[i] == fb[i] ? 0xffffffff : 0x0);
  }
  return rv;
}

static inline v4si
v4si_cmpeq(v4si a, v4si b)
{
  v4si rv;
  int *fa = (int *) &a, *fb = (int *) &b;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (fa[i] == fb[i] ? 0xffffffff : 0x0);
  }
  return rv;
}

static inline v4si
v4si_cmplt(v4si a, v4si b)
{
  v4si rv;
  int *fa = (int *) &a, *fb = (int *) &b;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (fa[i] < fb[i] ? 0xffffffff : 0x0);
  }
  return rv;
}

static inline v4s
v4s_blend(v4si mask, v4s a, v4s b)
{
  float rv[4];
  float *fa = (float *) &a, *fb = (float *) &b;

  for (int i = 0; i < 4; i++) {
    rv[i] = (mask.v[i] & 0x80000000) ? fa[i] : fb[i];
  }
  return (v4s) { rv[0], rv[1], rv[2], rv[3] };
}

static inline v4si
v4si_blend(v4si mask, v4si a, v4si b)
{
  v4si rv;

  for (int i = 0; i < 4; i++) {
    rv.v[i] = (mask.v[i] & 0x80000000) ? a.v[i] : b.v[i];
  }
  return rv;
}

static inline v4si
v4si_and(v4si a, v4si b)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (a.v[i] & b.v[i]);
  }
  return rv;
}

static inline v4si
v4si_or(v4si a, v4si b)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (a.v[i] | b.v[i]);
  }
  return rv;
}

static inline v4si
v4si_andnot(v4si a, v4si b)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (~a.v[i] & b.v[i]);
  }
  return rv;
}

static inline v4si
v4si_add(v4si a, v4si b)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = a.v[i] + b.v[i];
  }
  return rv;
}

static inline v4si
v4si_sub(v4si a, v4si b)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = a.v[i] - b.v[i];
  }
  return rv;
}

static inline v4si
v4si_mul(v4si a, v4si b)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = a.v[i] * b.v[i];
  }
  return rv;
}

static inline v4si
v4si_abs(v4si a)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = abs(a.v[i]);
  }
  return rv;
}

static inline v4si
v4si_sll(v4si a, int count)
{
  v4si rv;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (a.v[i] << count);
  }
  return rv;
}

static inline v4si
v4s_to_v4si(v4s a)
{
  v4si rv;
  float *fa = (float *) &a;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = (int) fa[i];
  }
  return rv;
}

static inline v4s
v4s_gather(real * __restrict pp, v4si off)
{
  float rv[4];
  for (int i = 0; i < 4; i++) {
    rv[i] = pp[off.v[i]];
  }

  return (v4s) { rv[0], rv[1], rv[2], rv[3] };
}

static inline v4s
v4s_cast(v4si a)
{
  float rv[4];
  for (int i = 0; i < 4; i++) {
    rv[i] = *((float *) &a.v[i]);
  }

  return (v4s) { rv[0], rv[1], rv[2], rv[3] };
}

static inline v4si
v4si_cast(v4s a)
{
  v4si rv;
  float *fa = (float *) &a;
  for (int i = 0; i < 4; i++) {
    rv.v[i] = *((int *) &fa[i]);
  }
  return rv;
}

#else // ----------------------------------------------------------------------

typedef float v4s __attribute__ ((vector_size (16)));
typedef int v4si __attribute__ ((vector_size (16)));

static const v4si iZero = { 0, 0, 0, 0 };

static inline v4si
v4si_splat(int x)
{
  v4si rv = { x, x, x, x };
  return rv;
}

static inline v4si
v4si_load(int *p)
{
  return (v4si) _mm_load_si128((__m128i *) p);
}

static inline void
v4si_store(int *p, v4si v)
{
  _mm_store_si128((__m128i *) p, (__m128i) v);
}

static inline int
v4si_extract(v4si v, int m)
{
  return _mm_extract_epi32((__m128i) v, m);
}

static inline v4si
v4s_fint(v4s v)
{
  return (v4si) _mm_sub_epi32(_mm_cvttps_epi32(v + v4s_splat(10.f)), (__m128i) v4si_splat(10));
}

static inline v4s
v4si_to_v4s(v4si v)
{
  return _mm_cvtepi32_ps((__m128i) v);
}

static inline v4si
v4s_cmpge(v4s a, v4s b)
{
  return (v4si) (__m128i) _mm_cmpge_ps((__m128) a, (__m128) b);
}

static inline v4si
v4s_cmpeq(v4s a, v4s b)
{
  return (v4si) (__m128i) _mm_cmpeq_ps((__m128) a, (__m128) b);
}

static inline v4si
v4si_cmpeq(v4si a, v4si b)
{
  return (v4si) _mm_cmpeq_epi32((__m128i) a, (__m128i) b);
}

static inline v4si
v4si_cmplt(v4si a, v4si b)
{
  return (v4si) _mm_cmplt_epi32((__m128i) a, (__m128i) b);
}

static inline v4s
v4s_blend(v4si mask, v4s a, v4s b)
{
  return (v4s) _mm_blendv_ps((__m128) a, (__m128) b, (__m128) (__m128i) mask);
}

static inline v4si
v4si_blend(v4si mask, v4si a, v4si b)
{
  return (v4si) _mm_blendv_epi8((__m128i) a, (__m128i) b, (__m128i) mask);
}

static inline v4si
v4si_and(v4si a, v4si b)
{
  return (v4si) _mm_and_si128((__m128i) a, (__m128i) b);
}

static inline v4si
v4si_or(v4si a, v4si b)
{
  return (v4si) _mm_or_si128((__m128i) a, (__m128i) b);
}

static inline v4si
v4si_andnot(v4si a, v4si b)
{
  return (v4si) _mm_andnot_si128((__m128i) a, (__m128i) b);
}

static inline v4si
v4si_add(v4si a, v4si b)
{
  return (v4si) _mm_add_epi32((__m128i) a, (__m128i) b);
}

static inline v4si
v4si_sub(v4si a, v4si b)
{
  return (v4si) _mm_sub_epi32((__m128i) a, (__m128i) b);
}

static inline v4si
v4si_mul(v4si a, v4si b)
{
  return (v4si) _mm_mullo_epi32((__m128i) a, (__m128i) b);
}

static inline v4si
v4si_abs(v4si a)
{
  return (v4si) _mm_abs_epi32((__m128i) a);
}

static inline v4si
v4si_sll(v4si a, int count)
{
  return (v4si) _mm_slli_epi32((__m128i) a, count);
}

static inline v4si
v4s_to_v4si(v4s a)
{
  return (v4si) _mm_cvttps_epi32((__m128) a);
}

static inline v4s
v4s_gather(real * __restrict pp, v4si _off)
{
  __m128i off = (__m128i) _off;
  __m128 res = _mm_load_ss(pp + _mm_cvtsi128_si32(off));
  res = _mm_insert_ps(res, _mm_load_ss(pp + _mm_extract_epi32(off, 1)), _MM_MK_INSERTPS_NDX(0,1,0));
  res = _mm_insert_ps(res, _mm_load_ss(pp + _mm_extract_epi32(off, 2)), _MM_MK_INSERTPS_NDX(0,2,0));
  res = _mm_insert_ps(res, _mm_load_ss(pp + _mm_extract_epi32(off, 3)), _MM_MK_INSERTPS_NDX(0,3,0));

  return (v4s) res;
}

static inline v4s
v4s_cast(v4si a)
{
  return (v4s) (__m128) (__m128i) a;
}

static inline v4si
v4si_cast(v4s a)
{
  return (v4si) (__m128i) (__m128) a;
}


#endif

static const v4s vZero = { 0.f, 0.f, 0.f, 0.f };
static const v4s vOne = { 1.f, 1.f, 1.f, 1.f };
static const v4s vTwo = { 2.f, 2.f, 2.f, 2.f };
static const v4s vThree = { 3.f, 3.f, 3.f, 3.f };
static const v4s vMp5 = { -.5f, -.5f, -.5f, -.5f };
static const v4s v4s_signmask = { -0.f, -0.f, -0.f, -0.f };

static inline v4s
v4s_splat(float x)
{
  v4s rv = { x, x, x, x };
  return rv;
}

static inline v4s
v4s_load(float *p)
{
  return (v4s) _mm_load_ps(p);
}

static inline void
v4s_prefetch(float *p)
{
  _mm_prefetch((const char *) p, _MM_HINT_T0);
}

static inline void
v4s_store(float *p, v4s v)
{
  _mm_store_ps(p, (__m128) v);
}

static inline void
v4s_stream(float *p, v4s v)
{
  _mm_stream_ps(p, (__m128) v);
}

static inline v4s
v4s_insert(v4s v, float s, int m)
{
  switch (m) {
  case 0: return (v4s) _mm_insert_ps((__m128) v, _mm_load_ss(&s), 0 << 4);
  case 1: return (v4s) _mm_insert_ps((__m128) v, _mm_load_ss(&s), 1 << 4);
  case 2: return (v4s) _mm_insert_ps((__m128) v, _mm_load_ss(&s), 2 << 4);
  case 3: return (v4s) _mm_insert_ps((__m128) v, _mm_load_ss(&s), 3 << 4);
  }
  assert(0);
}

static inline float
v4s_extract(v4s v, int m)
{
  float rv;
  _MM_EXTRACT_FLOAT(rv, (__m128) v, m);
  return rv;
}

static inline v4s
v4s_sqrt(v4s v)
{
  return _mm_sqrt_ps(v);
}

static inline v4s
v4s_rsqrt(v4s v)
{
#if 0
  v4s root = _mm_sqrt_ps(v);
  return vOne / root;
#else
  v4s x0 = _mm_rsqrt_ps(v);
  return vMp5 * x0 * (v * (x0*x0) - vThree);
#endif
}

static inline v4s
v4s_recip(v4s v)
{
#if 1
  return vOne / v;
#else
  v4s x0 = _mm_rcp_ps(v);
  return (x0 + x0) - v * (x0 * x0);
#endif
}

static inline v4s
v4s_floor(v4s a)
{
  return (v4s) _mm_floor_ps((__m128) a);
}

static inline v4s
v4s_abs(v4s x)
{
  static const v4s sign_mask = { -0.f, -0.f, -0.f, -0.f };
  return (v4s) _mm_andnot_ps((__m128) sign_mask, (__m128) x);
}

static inline v4s
v4s_xor(v4s a, v4s b)
{
  return (v4s) _mm_xor_ps((__m128) a, (__m128) b);
}

static inline v4s
v4s_and(v4s a, v4s b)
{
  return (v4s) _mm_and_ps((__m128) a, (__m128) b);
}

