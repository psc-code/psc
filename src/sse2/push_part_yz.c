#include "psc_sse2.h"

//include profiling, even though I won't use it yet
#include "profile/profile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>

void
sse2_push_part_yz()
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
                               // that it needs double?
  ones.r = _mm_set1_ps(1);
  half.r = _mm_set1_ps(.5);
  dt.r = _mm_set1_ps(dtfl);
  yl.r = _mm_mul_ps(half.r, dt.r);
  zl.r = _mm_mul_ps(half.r, dt.r);


  assert(psc.n_part % 4 == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += 4) {
    // FIXME : This is the wrong way to read in and pack data
    union packed_vector pxi, pyi, pzi, root, vyi, vzi, tmpx, tmpy, tmpz, yi, zi;
    for(int m = 0; m<4; m++){
    }
    // First particle
    yi.v[0] = sse2->part[n].yi;
    zi.v[0] = sse2->part[n].zi;
    
    pxi.v[0] = sse2->part[n].pxi;
    pyi.v[0] = sse2->part[n].pyi;
    pzi.v[0] = sse2->part[n].pzi;
    //Second particle
    yi.v[1] = sse2->part[n+1].yi;
    zi.v[1] = sse2->part[n+1].zi;
    
    pxi.v[1] = sse2->part[n+1].pxi;
    pyi.v[1] = sse2->part[n+1].pyi;
    pzi.v[1] = sse2->part[n+1].pzi;
    //Third Particle
    yi.v[2] = sse2->part[n+2].yi;
    zi.v[2] = sse2->part[n+2].zi;

    pxi.v[2] = sse2->part[n+2].pxi;
    pyi.v[2] = sse2->part[n+2].pyi;
    pzi.v[2] = sse2->part[n+2].pzi;
    //Fourth particle
    yi.v[3] = sse2->part[n+3].yi;
    zi.v[3] = sse2->part[n+3].zi;
    
    pxi.v[3] = sse2->part[n+3].pxi;
    pyi.v[3] = sse2->part[n+3].pyi;
    pzi.v[3] = sse2->part[n+3].pzi;
      
      
    
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
    
    (sse2->part[n]).yi = yi.v[0];
    (sse2->part[n]).zi = zi.v[0];
    (sse2->part[n+1]).yi = yi.v[1];
    (sse2->part[n+1]).zi = zi.v[1];
    (sse2->part[n+2]).yi = yi.v[2];
    (sse2->part[n+2]).zi = zi.v[2];
    (sse2->part[n+3]).yi = yi.v[3];
    (sse2->part[n+3]).zi = zi.v[3];
      
  }
  prof_stop(pr);
}
