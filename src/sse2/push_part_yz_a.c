#include "psc_sse2.h"

//include profiling, even though I won't use it yet
#include "profile/profile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>

void
sse2_push_part_yz_a()
{
  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz_a", 1., 0, psc.n_part * 12 * sizeof(float));
  }
  prof_start(pr);
  struct psc_sse2 *sse2 = psc.c_ctx;

  // Values that won't change from iteration to iteration

  union packed_vector dt, yl, zl, ones, half; 
  float dtfl = (float) psc.dt; // FIXME: will dt ever be small enough
                               // that is needs double.
  ones.r = _mm_set1_ps(1);
  half.r = _mm_set1_ps(.5);
  dt.r = _mm_set1_ps(dtfl);
  yl.r = _mm_mul_ps(half.r, dt.r);
  zl.r = _mm_mul_ps(half.r, dt.r);


  //  assert(psc.n_part % 4 == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += 4) {
    // FIXME : This is the wrong way to read in and pack data
    union packed_vector pxi, pyi, pzi, root, vyi, vzi, tmpx, tmpy, tmpz, yi, zi;
    for(int m = 0; m<4; m++){
      yi.v[m] = sse2->part[n+m].yi;
      zi.v[m] = sse2->part[n+m].zi;

      pxi.v[m] = sse2->part[n+m].pxi;
      pyi.v[m] = sse2->part[n+m].pyi;
      pzi.v[m] = sse2->part[n+m].pzi;
    }

    tmpx.r = _mm_mul_ps(pxi.r, pxi.r);
    tmpy.r = _mm_mul_ps(pyi.r, pyi.r);
    tmpz.r = _mm_mul_ps(pzi.r, pzi.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_add_ps(ones.r, tmpx.r);
    root.r = _mm_sqrt_ps(tmpx.r);
    
    vyi.r = _mm_div_ps(pyi.r, root.r);
    vzi.r = _mm_div_ps(pzi.r, root.r);

    tmpy.r = _mm_mul_ps(vyi.r, yl.r);
    tmpz.r = _mm_mul_ps(vzi.r, zl.r);
    
    yi.r = _mm_add_ps(yi.r, tmpy.r);
    zi.r = _mm_add_ps(zi.r, tmpz.r);
    
    // FIXME : Wrong way to unpack and return to memory
    for(int m = 0; m<4; m++){
      (sse2->part[n+m]).yi = yi.v[m];
      (sse2->part[n+m]).zi = zi.v[m];
    }
  }
  prof_stop(pr);
}
