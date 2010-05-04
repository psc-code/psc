#include "psc_sse2.h"

#include "profile/profile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <emmintrin.h>
#include <string.h>

void
sse2_push_part_yz()
{
  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz", 1., 0, psc.n_part * 12 * sizeof(float));
  }
  prof_start(pr);


//-----------------------------------------------------
// Initialization stuff (not sure what all of this is for)
  
  struct psc_sse2 *sse2 = psc.c_ctx;

  psc.p2A = 0.;
  psc.p2B = 0.;
  
  for (int m = 0; m <= JZI; m++) {
    memset(&sse2->fields[m*psc.fld_size], 0, psc.fld_size * sizeof(real));
  }  

  // Values that won't change from iteration to iteration

  pvFloat dt, yl, zl, ones, half,threefourths, onepfive,  eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs; 
  //FIXME: These are all stored as doubles in fortran!!
  float dtfl = psc.dt; 
  float etafl =  psc.prm.eta;
  float fnqsfl = sqr(psc.prm.alpha) * psc.prm.cori / etafl;
  float dxifl = 1.0f / psc.dx[0];
  float dyifl = 1.0f / psc.dx[1];
  float dzifl = 1.0f / psc.dx[2];
  float dqsfl = 0.5*etafl*dtfl;
  float fnqxsfl = psc.dx[0] * fnqsfl * psc.dt;
  float fnqysfl = psc.dx[1] * fnqsfl * psc.dt;
  float fnqzsfl = psc.dx[2] * fnqsfl * psc.dt;
  ones.r = _mm_set1_ps(1.0f);
  half.r = _mm_set1_ps(.5f);
  onepfive.r = _mm_set1_ps(1.5f);
  threefourths.r = _mm_set1_ps(.75f);
  dt.r = _mm_set1_ps(dtfl);
  eta.r = _mm_set1_ps(etafl);
  fnqs.r = _mm_set1_ps(fnqsfl);
  dqs.r = _mm_set1_ps(dqsfl);
  dxi.r = _mm_set1_ps(dxifl);
  dyi.r = _mm_set1_ps(dyifl);
  dzi.r = _mm_set1_ps(dzifl);
  yl.r = _mm_mul_ps(half.r, dt.r);
  zl.r = _mm_mul_ps(half.r, dt.r);
  fnqxs.r = _mm_set1_ps(fnqxsfl);
  fnqys.r = _mm_set1_ps(fnqysfl);
  fnqzs.r = _mm_set1_ps(fnqzsfl);

  //  assert(psc.n_part % 4 == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += 4) {
    pvFloat pxi, pyi, pzi, xi, yi, zi, qni, mni, wni;
//---------------------------------------------
// Bringing in particle specific parameters
// FIXME : This is the wrong way to read in and pack data
    // First particle
    xi.v[0] = sse2->part[n].xi;
    yi.v[0] = sse2->part[n].yi;
    zi.v[0] = sse2->part[n].zi;
    pxi.v[0] = sse2->part[n].pxi;
    pyi.v[0] = sse2->part[n].pyi;
    pzi.v[0] = sse2->part[n].pzi;
    qni.v[0] = sse2->part[n].qni;
    mni.v[0] = sse2->part[n].mni;
    wni.v[0] = sse2->part[n].wni;

    //Second particle
    xi.v[1] = sse2->part[n+1].xi;
    yi.v[1] = sse2->part[n+1].yi;
    zi.v[1] = sse2->part[n+1].zi;
    pxi.v[1] = sse2->part[n+1].pxi;
    pyi.v[1] = sse2->part[n+1].pyi;
    pzi.v[1] = sse2->part[n+1].pzi;
    qni.v[1] = sse2->part[n+1].qni;
    mni.v[1] = sse2->part[n+1].mni;
    wni.v[1] = sse2->part[n+1].wni;

    //Third Particle
    xi.v[2] = sse2->part[n+2].xi;
    yi.v[2] = sse2->part[n+2].yi;
    zi.v[2] = sse2->part[n+2].zi;
    pxi.v[2] = sse2->part[n+2].pxi;
    pyi.v[2] = sse2->part[n+2].pyi;
    pzi.v[2] = sse2->part[n+2].pzi;
    qni.v[2] = sse2->part[n+2].qni;
    mni.v[2] = sse2->part[n+2].mni;
    wni.v[2] = sse2->part[n+2].wni;

    //Fourth particle
    xi.v[3] = sse2->part[n+3].xi;
    yi.v[3] = sse2->part[n+3].yi;
    zi.v[3] = sse2->part[n+3].zi;    
    pxi.v[3] = sse2->part[n+3].pxi;
    pyi.v[3] = sse2->part[n+3].pyi;
    pzi.v[3] = sse2->part[n+3].pzi;
    qni.v[3] = sse2->part[n+3].qni;
    mni.v[3] = sse2->part[n+3].mni;
    wni.v[3] = sse2->part[n+3].wni;

    // Locals for computation      
    pvFloat vxi, vyi, vzi, tmpx, tmpy, tmpz, root;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta
    
    tmpx.r = _mm_mul_ps(pxi.r, pxi.r);
    tmpy.r = _mm_mul_ps(pyi.r, pyi.r);
    tmpz.r = _mm_mul_ps(pzi.r, pzi.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_add_ps(ones.r, tmpx.r);
    root.r = _mm_sqrt_ps(tmpx.r);
    root.r = _mm_div_ps(ones.r, root.r);

    vxi.r = _mm_mul_ps(pxi.r, root.r);
    vyi.r = _mm_mul_ps(pyi.r, root.r);
    vzi.r = _mm_mul_ps(pzi.r, root.r);

    tmpy.r = _mm_mul_ps(vyi.r, yl.r);
    tmpz.r = _mm_mul_ps(vzi.r, zl.r);
    yi.r = _mm_add_ps(yi.r, tmpy.r);
    zi.r = _mm_add_ps(zi.r, tmpz.r);


    tmpx.r = _mm_div_ps(ones.r, root.r);
    tmpx.r = _mm_sub_ps(tmpx.r, ones.r);
    tmpx.r = _mm_div_ps(tmpx.r, eta.r);
    tmpx.r = _mm_mul_ps(fnqs.r, tmpx.r);
    tmpx.r = _mm_mul_ps(mni.r, tmpx.r);
    psc.p2A += tmpx.v[0] + tmpx.v[1] + tmpx.v[2] + tmpx.v[3]; //What's this for?

// CHECKPOINT: PIC_push_part_yz.F : line 110
    
    tmpx.r = _mm_mul_ps(xi.r, dxi.r);
    tmpy.r = _mm_mul_ps(yi.r, dyi.r);
    tmpz.r = _mm_mul_ps(zi.r, dzi.r);

    pvFloat s0y[5], s0z[5], s1y[5], s1z[5];
 
// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    // Apparently this can be done in vectors. A victory!
    pvInt j1, j2, j3, l1, l2, l3;
    pvFloat j2fl, j3fl, l2fl, l3fl;
    j1.r = _mm_cvtps_epi32(tmpx.r);
    j2.r = _mm_cvtps_epi32(tmpy.r);
    j3.r = _mm_cvtps_epi32(tmpz.r);

    for(int m = 0; m < 4; m++){
      j1.v[m] = round(tmpx.v[m]);
      j2.v[m] = round(tmpy.v[m]);
      j3.v[m] = round(tmpz.v[m]);
   
    }
    
    // there must be a better way...
    j2fl.r = _mm_cvtepi32_ps(j2.r);
    j3fl.r = _mm_cvtepi32_ps(j3.r);

    tmpy.r = _mm_sub_ps(j2fl.r, tmpy.r);
    tmpz.r = _mm_sub_ps(j3fl.r, tmpz.r);


    pvFloat gmy, gmz, g0y, g0z, g1y, g1z;

    gmy.r = _mm_add_ps(half.r, tmpy.r);
    gmy.r = _mm_mul_ps(gmy.r, gmy.r);
    gmy.r = _mm_mul_ps(half.r, gmy.r);

    gmz.r = _mm_add_ps(half.r, tmpz.r);
    gmz.r = _mm_mul_ps(gmz.r, gmz.r);
    gmz.r = _mm_mul_ps(half.r, gmz.r);
    
    g0y.r = _mm_mul_ps(tmpy.r, tmpy.r);
    g0y.r = _mm_sub_ps(threefourths.r, g0y.r);

    g0z.r = _mm_mul_ps(tmpz.r, tmpz.r);
    g0z.r = _mm_sub_ps(threefourths.r, g0z.r);

    g1y.r = _mm_sub_ps(half.r, tmpy.r);
    g1y.r = _mm_mul_ps(g1y.r, g1y.r);
    g1y.r = _mm_mul_ps(half.r, g1y.r);

    g1z.r = _mm_sub_ps(half.r, tmpz.r);
    g1z.r = _mm_mul_ps(g1z.r, g1z.r);
    g1z.r = _mm_mul_ps(half.r, g1z.r);

    // indexing here departs from FORTRAN a little bit
    s0y[1].r = _mm_sub_ps(tmpy.r, ones.r);
    s0y[1].r = _mm_add_ps(onepfive.r, s0y[1].r); // h2+1.0 always <0
    s0y[1].r = _mm_mul_ps(s0y[1].r, s0y[1].r);
    s0y[1].r = _mm_mul_ps(half.r, s0y[1].r);

    s0y[2].r = _mm_mul_ps(tmpy.r, tmpy.r);
    s0y[2].r = _mm_sub_ps(threefourths.r, s0y[2].r);

    s0y[3].r = _mm_add_ps(tmpy.r, ones.r);
    s0y[3].r = _mm_sub_ps(onepfive.r, s0y[3].r); // h2+1.0 always >0
    s0y[3].r = _mm_mul_ps(s0y[3].r, s0y[3].r);
    s0y[3].r = _mm_mul_ps(half.r, s0y[3].r);

    s0z[1].r = _mm_sub_ps(tmpz.r, ones.r);
    s0z[1].r = _mm_add_ps(onepfive.r, s0z[1].r); // h3+1.0 always <0
    s0z[1].r = _mm_mul_ps(s0z[1].r, s0z[1].r);
    s0z[1].r = _mm_mul_ps(half.r, s0z[1].r);

    s0z[2].r = _mm_mul_ps(tmpz.r, tmpz.r);
    s0z[2].r = _mm_sub_ps(threefourths.r, s0z[2].r);

    s0z[3].r = _mm_add_ps(tmpz.r, ones.r);
    s0z[3].r = _mm_sub_ps(onepfive.r, s0z[3].r); // h3+1.0 always >0
    s0z[3].r = _mm_mul_ps(s0z[3].r, s0z[3].r);
    s0z[3].r = _mm_mul_ps(half.r, s0z[3].r);



// CHECKPOINT: PIC_push_part_yz.F : line 143

    tmpx.r = _mm_mul_ps(xi.r, dxi.r);
    tmpy.r = _mm_mul_ps(yi.r, dyi.r);
    tmpy.r = _mm_sub_ps(tmpy.r, half.r);
    tmpz.r = _mm_mul_ps(zi.r, dzi.r);
    tmpz.r = _mm_sub_ps(tmpz.r, half.r);

    l1.r = _mm_cvtps_epi32(tmpx.r);
    l2.r = _mm_cvtps_epi32(tmpy.r);
    l3.r = _mm_cvtps_epi32(tmpz.r);

    l2fl.r = _mm_cvtepi32_ps(l2.r);
    l3fl.r = _mm_cvtepi32_ps(l3.r);

    tmpy.r = _mm_sub_ps(l2fl.r, tmpy.r);
    tmpz.r = _mm_sub_ps(l3fl.r, tmpz.r);

    pvFloat hmy, hmz, h0y, h0z, h1y, h1z;

    hmy.r = _mm_add_ps(half.r, tmpy.r);
    hmy.r = _mm_mul_ps(hmy.r, hmy.r);
    hmy.r = _mm_mul_ps(half.r, hmy.r);

    hmz.r = _mm_add_ps(half.r, tmpz.r);
    hmz.r = _mm_mul_ps(hmz.r, hmz.r);
    hmz.r = _mm_mul_ps(half.r, hmz.r);
    
    h0y.r = _mm_mul_ps(tmpy.r, tmpy.r);
    h0y.r = _mm_sub_ps(threefourths.r, h0y.r);

    h0z.r = _mm_mul_ps(tmpz.r, tmpz.r);
    h0z.r = _mm_sub_ps(threefourths.r, h0z.r);

    h1y.r = _mm_sub_ps(half.r, tmpy.r);
    h1y.r = _mm_mul_ps(h1y.r, h1y.r);
    h1y.r = _mm_mul_ps(half.r, h1y.r);

    h1z.r = _mm_sub_ps(half.r, tmpz.r);
    h1z.r = _mm_mul_ps(h1z.r, h1z.r);
    h1z.r = _mm_mul_ps(half.r, h1z.r);

// CHECKPOINT: PIC_push_part_yz.F : line 59
// Field Interpolation

// BOTTLENECK HERE: For a given particle, none of the needed fields are stored
// next to each other in memory. Also, there's no way to tell before runtime
// exactly what fields are needed. As a consequence, I can't just load them into
// registers and shuffle. I must fetch piece by painful piece.

// A messy macro to let me get the correct data and keep notation consistent
// with the Fortran code

    pvFloat field_in, exq, eyq, ezq, bxq, byq, bzq;
     
    //exq
    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0]-1,j3.v[0]-1);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1]-1,j3.v[1]-1);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2]-1,j3.v[2]-1);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3]-1,j3.v[3]-1);
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0],j3.v[0]-1);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1],j3.v[1]-1);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2],j3.v[2]-1);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3],j3.v[3]-1);    
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0]+1,j3.v[0]-1);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1]+1,j3.v[1]-1);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2]+1,j3.v[2]-1);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3]+1,j3.v[3]-1);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    exq.r = _mm_mul_ps(gmz.r, tmpx.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0]-1,j3.v[0]);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1]-1,j3.v[1]);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2]-1,j3.v[2]);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3]-1,j3.v[3]);        
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0],j3.v[0]);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1],j3.v[1]);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2],j3.v[2]);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3],j3.v[3]);        
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0]+1,j3.v[0]);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1]+1,j3.v[1]);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2]+1,j3.v[2]);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3]+1,j3.v[3]);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(g0z.r, tmpx.r);
    exq.r = _mm_add_ps(exq.r, tmpx.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0]-1,j3.v[0]+1);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1]-1,j3.v[1]+1);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2]-1,j3.v[2]+1);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3]-1,j3.v[3]+1);        
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0],j3.v[0]+1);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1],j3.v[1]+1);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2],j3.v[2]+1);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3],j3.v[3]+1);        
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EX, l1.v[0], j2.v[0]+1,j3.v[0]+1);
    field_in.v[1] = C_FIELD(EX, l1.v[1], j2.v[1]+1,j3.v[1]+1);
    field_in.v[2] = C_FIELD(EX, l1.v[2], j2.v[2]+1,j3.v[2]+1);
    field_in.v[3] = C_FIELD(EX, l1.v[3], j2.v[3]+1,j3.v[3]+1);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(g1z.r, tmpx.r);
    exq.r = _mm_add_ps(exq.r, tmpx.r);

    //eyq
    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0]-1,j3.v[0]-1);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1]-1,j3.v[1]-1);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2]-1,j3.v[2]-1);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3]-1,j3.v[3]-1);
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0],j3.v[0]-1);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1],j3.v[1]-1);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2],j3.v[2]-1);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3],j3.v[3]-1);    
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0]+1,j3.v[0]-1);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1]+1,j3.v[1]-1);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2]+1,j3.v[2]-1);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3]+1,j3.v[3]-1);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    eyq.r = _mm_mul_ps(gmz.r, tmpx.r);
 
    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0]-1,j3.v[0]);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1]-1,j3.v[1]);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2]-1,j3.v[2]);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3]-1,j3.v[3]);        
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0],j3.v[0]);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1],j3.v[1]);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2],j3.v[2]);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3],j3.v[3]);        
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0]+1,j3.v[0]);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1]+1,j3.v[1]);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2]+1,j3.v[2]);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3]+1,j3.v[3]);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(g0z.r, tmpx.r);
    eyq.r = _mm_add_ps(eyq.r, tmpx.r);

    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0]-1,j3.v[0]+1);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1]-1,j3.v[1]+1);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2]-1,j3.v[2]+1);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3]-1,j3.v[3]+1);        
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0],j3.v[0]+1);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1],j3.v[1]+1);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2],j3.v[2]+1);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3],j3.v[3]+1);        
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EY, j1.v[0], l2.v[0]+1,j3.v[0]+1);
    field_in.v[1] = C_FIELD(EY, j1.v[1], l2.v[1]+1,j3.v[1]+1);
    field_in.v[2] = C_FIELD(EY, j1.v[2], l2.v[2]+1,j3.v[2]+1);
    field_in.v[3] = C_FIELD(EY, j1.v[3], l2.v[3]+1,j3.v[3]+1);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(g1z.r, tmpx.r);
    eyq.r = _mm_add_ps(eyq.r, tmpx.r);

    //ezq
    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0]-1,l3.v[0]-1);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1]-1,l3.v[1]-1);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2]-1,l3.v[2]-1);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3]-1,l3.v[3]-1);
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0],l3.v[0]-1);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1],l3.v[1]-1);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2],l3.v[2]-1);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3],l3.v[3]-1);    
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0]+1,l3.v[0]-1);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1]+1,l3.v[1]-1);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2]+1,l3.v[2]-1);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3]+1,l3.v[3]-1);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    ezq.r = _mm_mul_ps(hmz.r, tmpx.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0]-1,l3.v[0]);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1]-1,l3.v[1]);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2]-1,l3.v[2]);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3]-1,l3.v[3]);        
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0],l3.v[0]);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1],l3.v[1]);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2],l3.v[2]);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3],l3.v[3]);        
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0]+1,l3.v[0]);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1]+1,l3.v[1]);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2]+1,l3.v[2]);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3]+1,l3.v[3]);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(h0z.r, tmpx.r);
    ezq.r = _mm_add_ps(ezq.r, tmpx.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0]-1,l3.v[0]+1);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1]-1,l3.v[1]+1);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2]-1,l3.v[2]+1);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3]-1,l3.v[3]+1);        
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0],l3.v[0]+1);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1],l3.v[1]+1);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2],l3.v[2]+1);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3],l3.v[3]+1);        
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(EZ, j1.v[0], j2.v[0]+1,l3.v[0]+1);
    field_in.v[1] = C_FIELD(EZ, j1.v[1], j2.v[1]+1,l3.v[1]+1);
    field_in.v[2] = C_FIELD(EZ, j1.v[2], j2.v[2]+1,l3.v[2]+1);
    field_in.v[3] = C_FIELD(EZ, j1.v[3], j2.v[3]+1,l3.v[3]+1);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(h1z.r, tmpx.r);
    ezq.r = _mm_add_ps(ezq.r, tmpx.r);

    //bxq
    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0]-1,l3.v[0]-1);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1]-1,l3.v[1]-1);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2]-1,l3.v[2]-1);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3]-1,l3.v[3]-1);
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0],l3.v[0]-1);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1],l3.v[1]-1);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2],l3.v[2]-1);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3],l3.v[3]-1);    
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0]+1,l3.v[0]-1);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1]+1,l3.v[1]-1);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2]+1,l3.v[2]-1);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3]+1,l3.v[3]-1);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    bxq.r = _mm_mul_ps(hmz.r, tmpx.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0]-1,l3.v[0]);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1]-1,l3.v[1]);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2]-1,l3.v[2]);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3]-1,l3.v[3]);        
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0],l3.v[0]);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1],l3.v[1]);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2],l3.v[2]);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3],l3.v[3]);        
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0]+1,l3.v[0]);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1]+1,l3.v[1]);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2]+1,l3.v[2]);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3]+1,l3.v[3]);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(h0z.r, tmpx.r);
    bxq.r = _mm_add_ps(bxq.r, tmpx.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0]-1,l3.v[0]+1);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1]-1,l3.v[1]+1);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2]-1,l3.v[2]+1);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3]-1,l3.v[3]+1);        
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0],l3.v[0]+1);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1],l3.v[1]+1);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2],l3.v[2]+1);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3],l3.v[3]+1);        
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BX, j1.v[0], l2.v[0]+1,l3.v[0]+1);
    field_in.v[1] = C_FIELD(BX, j1.v[1], l2.v[1]+1,l3.v[1]+1);
    field_in.v[2] = C_FIELD(BX, j1.v[2], l2.v[2]+1,l3.v[2]+1);
    field_in.v[3] = C_FIELD(BX, j1.v[3], l2.v[3]+1,l3.v[3]+1);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(h1z.r, tmpx.r);
    bxq.r = _mm_add_ps(bxq.r, tmpx.r);

    //byq
    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0]-1,l3.v[0]-1);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1]-1,l3.v[1]-1);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2]-1,l3.v[2]-1);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3]-1,l3.v[3]-1);
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0],l3.v[0]-1);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1],l3.v[1]-1);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2],l3.v[2]-1);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3],l3.v[3]-1);    
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0]+1,l3.v[0]-1);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1]+1,l3.v[1]-1);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2]+1,l3.v[2]-1);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3]+1,l3.v[3]-1);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    byq.r = _mm_mul_ps(hmz.r, tmpx.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0]-1,l3.v[0]);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1]-1,l3.v[1]);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2]-1,l3.v[2]);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3]-1,l3.v[3]);        
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0],l3.v[0]);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1],l3.v[1]);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2],l3.v[2]);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3],l3.v[3]);        
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0]+1,l3.v[0]);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1]+1,l3.v[1]);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2]+1,l3.v[2]);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3]+1,l3.v[3]);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(h0z.r, tmpx.r);
    byq.r = _mm_add_ps(byq.r, tmpx.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0]-1,l3.v[0]+1);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1]-1,l3.v[1]+1);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2]-1,l3.v[2]+1);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3]-1,l3.v[3]+1);        
    tmpx.r = _mm_mul_ps(gmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0],l3.v[0]+1);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1],l3.v[1]+1);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2],l3.v[2]+1);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3],l3.v[3]+1);        
    tmpy.r = _mm_mul_ps(g0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BY, l1.v[0], j2.v[0]+1,l3.v[0]+1);
    field_in.v[1] = C_FIELD(BY, l1.v[1], j2.v[1]+1,l3.v[1]+1);
    field_in.v[2] = C_FIELD(BY, l1.v[2], j2.v[2]+1,l3.v[2]+1);
    field_in.v[3] = C_FIELD(BY, l1.v[3], j2.v[3]+1,l3.v[3]+1);        
    tmpz.r = _mm_mul_ps(g1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(h1z.r, tmpx.r);
    byq.r = _mm_add_ps(byq.r, tmpx.r);

    //bzq
    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0]-1,j3.v[0]-1);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1]-1,j3.v[1]-1);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2]-1,j3.v[2]-1);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3]-1,j3.v[3]-1);
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0],j3.v[0]-1);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1],j3.v[1]-1);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2],j3.v[2]-1);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3],j3.v[3]-1);    
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0]+1,j3.v[0]-1);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1]+1,j3.v[1]-1);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2]+1,j3.v[2]-1);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3]+1,j3.v[3]-1);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    bzq.r = _mm_mul_ps(gmz.r, tmpx.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0]-1,j3.v[0]);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1]-1,j3.v[1]);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2]-1,j3.v[2]);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3]-1,j3.v[3]);        
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0],j3.v[0]);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1],j3.v[1]);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2],j3.v[2]);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3],j3.v[3]);        
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0]+1,j3.v[0]);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1]+1,j3.v[1]);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2]+1,j3.v[2]);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3]+1,j3.v[3]);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(g0z.r, tmpx.r);
    bzq.r = _mm_add_ps(bzq.r, tmpx.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0]-1,j3.v[0]+1);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1]-1,j3.v[1]+1);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2]-1,j3.v[2]+1);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3]-1,j3.v[3]+1);        
    tmpx.r = _mm_mul_ps(hmy.r, field_in.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0],j3.v[0]+1);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1],j3.v[1]+1);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2],j3.v[2]+1);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3],j3.v[3]+1);        
    tmpy.r = _mm_mul_ps(h0y.r, field_in.r);

    field_in.v[0] = C_FIELD(BZ, l1.v[0], l2.v[0]+1,j3.v[0]+1);
    field_in.v[1] = C_FIELD(BZ, l1.v[1], l2.v[1]+1,j3.v[1]+1);
    field_in.v[2] = C_FIELD(BZ, l1.v[2], l2.v[2]+1,j3.v[2]+1);
    field_in.v[3] = C_FIELD(BZ, l1.v[3], l2.v[3]+1,j3.v[3]+1);        
    tmpz.r = _mm_mul_ps(h1y.r, field_in.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_mul_ps(g1z.r, tmpx.r);
    bzq.r = _mm_add_ps(bzq.r, tmpx.r);

// CHECKPOINT: PIC_push_part_yz.F : line 223
// Half step momentum  with E-field

    pvFloat dq;
    dq.r = _mm_div_ps(dqs.r, mni.r);
    dq.r = _mm_mul_ps(qni.r, dq.r);
    
    exq.r = _mm_mul_ps(dq.r, exq.r);
    eyq.r = _mm_mul_ps(dq.r, eyq.r);
    ezq.r = _mm_mul_ps(dq.r, ezq.r);
    
    pxi.r = _mm_add_ps(pxi.r, exq.r);
    pyi.r = _mm_add_ps(pyi.r, eyq.r);
    pzi.r = _mm_add_ps(pzi.r, ezq.r);


// CHECKPOINT: PIC_push_part_yz.F : line 228
// Rotate with B-field
    pvFloat txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;

    tmpx.r = _mm_mul_ps(pxi.r, pxi.r);
    tmpy.r = _mm_mul_ps(pyi.r, pyi.r);
    tmpz.r = _mm_mul_ps(pzi.r, pzi.r);
    root.r = _mm_add_ps(ones.r, tmpx.r);
    tmpy.r = _mm_add_ps(tmpy.r, tmpz.r);
    root.r = _mm_add_ps(root.r, tmpy.r);
    root.r = _mm_sqrt_ps(root.r);
    root.r = _mm_div_ps(dq.r, root.r);

    bxq.r = _mm_mul_ps(bxq.r, root.r);
    byq.r = _mm_mul_ps(byq.r, root.r);
    bzq.r = _mm_mul_ps(bzq.r, root.r);

    txx.r = _mm_mul_ps(bxq.r, bxq.r);
    tyy.r = _mm_mul_ps(byq.r, byq.r);
    tzz.r = _mm_mul_ps(bzq.r, bzq.r);
    t2xy.r = _mm_mul_ps(bxq.r, byq.r);
    t2xz.r = _mm_mul_ps(bxq.r, bzq.r);
    t2yz.r = _mm_mul_ps(byq.r, bzq.r);
    t2xy.r = _mm_add_ps(t2xy.r, t2xy.r);
    t2xz.r = _mm_add_ps(t2xz.r, t2xz.r);
    t2yz.r = _mm_add_ps(t2yz.r, t2yz.r);

    pvFloat tau;

    tau.r = _mm_add_ps(ones.r, txx.r);
    tmpx.r = _mm_add_ps(tyy.r, tzz.r);
    tau.r = _mm_add_ps(tau.r, tmpx.r);
    tau.r = _mm_div_ps(ones.r, tau.r); //recp is evil! evil evil evil evil!!!!
    
    bxq.r = _mm_add_ps(bxq.r, bxq.r);
    byq.r = _mm_add_ps(byq.r, byq.r);
    bzq.r = _mm_add_ps(bzq.r, bzq.r);
    
    //pxp
    tmpx.r = _mm_add_ps(ones.r, txx.r);
    tmpx.r = _mm_sub_ps(tmpx.r, tyy.r);
    tmpx.r = _mm_sub_ps(tmpx.r, tzz.r);
    tmpx.r = _mm_mul_ps(tmpx.r, pxi.r);
    
    tmpy.r = _mm_add_ps(t2xy.r, bzq.r);
    tmpy.r = _mm_mul_ps(tmpy.r, pyi.r);

    tmpz.r = _mm_sub_ps(t2xz.r, byq.r);
    tmpz.r = _mm_mul_ps(tmpz.r, pzi.r);

    pxp.r = _mm_add_ps(tmpx.r, tmpy.r);
    pxp.r = _mm_add_ps(pxp.r, tmpz.r);
    pxp.r = _mm_mul_ps(pxp.r, tau.r);
    
    //pyp
    tmpx.r = _mm_sub_ps(t2xy.r, bzq.r);
    tmpx.r = _mm_mul_ps(tmpx.r, pxi.r);

    tmpy.r = _mm_sub_ps(ones.r, txx.r);
    tmpy.r = _mm_add_ps(tmpy.r, tyy.r);
    tmpy.r = _mm_sub_ps(tmpy.r, tzz.r);
    tmpy.r = _mm_mul_ps(tmpy.r, pyi.r);
    
    tmpz.r = _mm_add_ps(t2yz.r, bxq.r);
    tmpz.r = _mm_mul_ps(tmpz.r, pzi.r);

    pyp.r = _mm_add_ps(tmpx.r, tmpy.r);
    pyp.r = _mm_add_ps(pyp.r, tmpz.r);
    pyp.r = _mm_mul_ps(pyp.r, tau.r);
   
    //pzp
    tmpx.r = _mm_add_ps(t2xz.r, byq.r);
    tmpx.r = _mm_mul_ps(tmpx.r, pxi.r);

    tmpy.r = _mm_sub_ps(t2yz.r, bxq.r);
    tmpy.r = _mm_mul_ps(tmpy.r, pyi.r);

    tmpz.r = _mm_sub_ps(ones.r, txx.r);
    tmpz.r = _mm_sub_ps(tmpz.r, tyy.r);
    tmpz.r = _mm_add_ps(tmpz.r, tzz.r);
    tmpz.r = _mm_mul_ps(tmpz.r, pzi.r);

    pzp.r = _mm_add_ps(tmpx.r, tmpy.r);
    pzp.r = _mm_add_ps(pzp.r, tmpz.r);
    pzp.r = _mm_mul_ps(pzp.r, tau.r);

// CHECKPOINT: PIC_push_part_yz.F : line 244
// Half step momentum  with E-field

    pxi.r = _mm_add_ps(pxp.r, exq.r);
    pyi.r = _mm_add_ps(pyp.r, eyq.r);
    pzi.r = _mm_add_ps(pzp.r, ezq.r);

// CHECKPOINT: PIC_push_part_yz.F : line 248
// Half step particles with new momenta
    
    tmpx.r = _mm_mul_ps(pxi.r, pxi.r);
    tmpy.r = _mm_mul_ps(pyi.r, pyi.r);
    tmpz.r = _mm_mul_ps(pzi.r, pzi.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_add_ps(ones.r, tmpx.r);
    root.r = _mm_sqrt_ps(tmpx.r);
    root.r = _mm_div_ps(ones.r, root.r);

    vxi.r = _mm_mul_ps(pxi.r, root.r);
    vyi.r = _mm_mul_ps(pyi.r, root.r);
    vzi.r = _mm_mul_ps(pzi.r, root.r);
    
    tmpy.r = _mm_mul_ps(vyi.r, yl.r);
    tmpz.r = _mm_mul_ps(vzi.r, zl.r);
    
    yi.r = _mm_add_ps(yi.r, tmpy.r);
    zi.r = _mm_add_ps(zi.r, tmpz.r);
    
    
    (sse2->part[n]).xi = xi.v[0];
    (sse2->part[n]).yi = yi.v[0];
    (sse2->part[n]).zi = zi.v[0];
    (sse2->part[n]).pxi = pxi.v[0];
    (sse2->part[n]).pyi = pyi.v[0];
    (sse2->part[n]).pzi = pzi.v[0];    

    (sse2->part[n+1]).xi = xi.v[1];
    (sse2->part[n+1]).yi = yi.v[1];
    (sse2->part[n+1]).zi = zi.v[1];
    (sse2->part[n+1]).pxi = pxi.v[1];
    (sse2->part[n+1]).pyi = pyi.v[1];
    (sse2->part[n+1]).pzi = pzi.v[1];    

    (sse2->part[n+2]).xi = xi.v[2];
    (sse2->part[n+2]).yi = yi.v[2];
    (sse2->part[n+2]).zi = zi.v[2];
    (sse2->part[n+2]).pxi = pxi.v[2];
    (sse2->part[n+2]).pyi = pyi.v[2];
    (sse2->part[n+2]).pzi = pzi.v[2];    

    (sse2->part[n+3]).xi = xi.v[3];
    (sse2->part[n+3]).yi = yi.v[3];
    (sse2->part[n+3]).zi = zi.v[3];
    (sse2->part[n+3]).pxi = pxi.v[3];
    (sse2->part[n+3]).pyi = pyi.v[3];
    (sse2->part[n+3]).pzi = pzi.v[3];    

    tmpx.r = _mm_div_ps(ones.r, root.r);
    tmpx.r = _mm_sub_ps(tmpx.r, ones.r);
    tmpx.r = _mm_div_ps(tmpx.r, eta.r);
    tmpx.r = _mm_mul_ps(fnqs.r, tmpx.r);
    tmpx.r = _mm_mul_ps(mni.r, tmpx.r);
    psc.p2B += tmpx.v[0] + tmpx.v[1] + tmpx.v[2] + tmpx.v[3]; //What's this for?

// CHECKPOINT: PIC_push_part_yz.F : line 266
// update number densities.
    tmpx.r = _mm_mul_ps(xi.r, dxi.r);
    tmpy.r = _mm_mul_ps(yi.r, dyi.r);
    tmpz.r = _mm_mul_ps(zi.r, dzi.r);

    l1.r = _mm_cvtps_epi32(tmpx.r);
    l2.r = _mm_cvtps_epi32(tmpy.r);
    l3.r = _mm_cvtps_epi32(tmpz.r);

    /* if(n == 0){ */
    /*   fprintf(stderr, "%d %d %d\n", l1.v[0], l2.v[0], l3.v[0]); */
    /* } */

    
    for(int m = 0; m < 4; m++){
      l1.v[m] = round(tmpx.v[m]);
      l2.v[m] = round(tmpy.v[m]);
      l3.v[m] = round(tmpz.v[m]);
   
    }

    /* if(n == 0){ */
    /*   fprintf(stderr, "%d %d %d\n", l1.v[0], l2.v[0], l3.v[0]); */
    /* } */

    l2fl.r = _mm_cvtepi32_ps(l2.r);
    l3fl.r = _mm_cvtepi32_ps(l3.r);

    tmpy.r = _mm_sub_ps(l2fl.r, tmpy.r);
    tmpz.r = _mm_sub_ps(l3fl.r, tmpz.r);

    gmy.r = _mm_sub_ps(tmpy.r,ones.r);
    //I'm pretty sure |h2,h3| < 0.5. If I'm wrong, I'll change this part
    gmy.r = _mm_add_ps(onepfive.r, gmy.r); // h2-1 always <0
    gmy.r = _mm_mul_ps(gmy.r, gmy.r);
    gmy.r = _mm_mul_ps(half.r, gmy.r);

    g0y.r = _mm_mul_ps(tmpy.r, tmpy.r);
    g0y.r = _mm_sub_ps(threefourths.r, g0y.r);

    g1y.r = _mm_add_ps(ones.r, tmpy.r); 
    g1y.r = _mm_sub_ps(onepfive.r, g1y.r);//h2+1 always >0
    g1y.r = _mm_mul_ps(g1y.r, g1y.r);
    g1y.r = _mm_mul_ps(half.r, g1y.r);

    gmz.r = _mm_sub_ps(tmpz.r,ones.r);
    gmz.r = _mm_add_ps(onepfive.r, gmz.r); //h3-1 always <0
    gmz.r = _mm_mul_ps(gmz.r, gmz.r);
    gmz.r = _mm_mul_ps(half.r, gmz.r);

    g0z.r = _mm_mul_ps(tmpz.r, tmpz.r);
    g0z.r = _mm_sub_ps(threefourths.r, g0z.r);

    g1z.r = _mm_add_ps(ones.r, tmpz.r);
    g1z.r = _mm_sub_ps(onepfive.r, g1z.r); //h3+1 always >0
    g1z.r = _mm_mul_ps(g1z.r, g1z.r);
    g1z.r = _mm_mul_ps(half.r, g1z.r);

// CHECKPOINT: PIC_push_part_yz.F : line 283
    //FIXME: this is, for now, a straight serial translation. I think
    //    there is a more efficient way to do this, but I might be making
    //    some changes to the the local field structure and I'll worry 
    //    about it then.

    for(int m=0; m<4; m++){ 
      float fnq;
      // This may or may not work, and may or may not help
      float *densp; 
      if (qni.v[m] < 0.0){
	densp = &(sse2->fields[NE*psc.fld_size + l1.v[m] - psc.ilo[0] + psc.ibn[0]]);
	fnq = qni.v[m] * wni.v[m] * fnqs.v[m];
      }
      else if (qni.v[m] > 0.0){
	densp = &(sse2->fields[NI*psc.fld_size + l1.v[m] - psc.ilo[0] + psc.ibn[0]]);
	fnq = qni.v[m] * wni.v[m] * fnqs.v[m];
      }
      else if (qni.v[m] == 0.0){
	densp = &(sse2->fields[NN*psc.fld_size + l1.v[m] - psc.ilo[0] + psc.ibn[0]]);
	fnq = wni.v[m] * fnqs.v[m];
      }
#define DEN_FIELD(j,k) *(densp+((j)-psc.ilo[1]+psc.ibn[1])*(psc.img[0]) + ((k) - psc.ilo[2]+psc.ibn[2])*(psc.img[0]*psc.img[1]))
      	DEN_FIELD(l2.v[m]-1, l3.v[m]-1) += fnq*gmy.v[m]*gmz.v[m];
	DEN_FIELD(l2.v[m], l3.v[m]-1) += fnq*g0y.v[m]*gmz.v[m];
	DEN_FIELD(l2.v[m]+1,l3.v[m]-1) += fnq*g1y.v[m]*gmz.v[m];
	DEN_FIELD(l2.v[m]-1, l3.v[m]) += fnq*gmy.v[m]*g0z.v[m];
	DEN_FIELD(l2.v[m], l3.v[m]) += fnq*g0y.v[m]*g0z.v[m];
	DEN_FIELD(l2.v[m]+1,l3.v[m]) += fnq*g1y.v[m]*g0z.v[m];
	DEN_FIELD(l2.v[m]-1, l3.v[m]+1) += fnq*gmy.v[m]*g1z.v[m];
	DEN_FIELD(l2.v[m], l3.v[m]+1) += fnq*g0y.v[m]*g1z.v[m];
	DEN_FIELD(l2.v[m]+1,l3.v[m]+1) += fnq*g1y.v[m]*g1z.v[m];

      /* if(n == 0){ */
      /* 	fprintf(stderr, "%g %g\n", DEN_FIELD(l2.v[m]-1, l3.v[m]-1), C_FIELD(NE,l1.v[m],l2.v[m]-1, l3.v[m]-1)); */
      /* } */
      
    }
    
    // CHECKPOINT: PIC_push_part_yz.F : line 320
    // Charge density form factor at (n+1.5)*dt
    tmpy.r = _mm_mul_ps(vyi.r, yl.r);
    tmpz.r = _mm_mul_ps(vzi.r, zl.r);
    
    yi.r = _mm_add_ps(yi.r, tmpy.r);
    zi.r = _mm_add_ps(zi.r, tmpz.r);

    pvInt k2,k3;

    pvFloat k2fl, k3fl;

    tmpy.r = _mm_mul_ps(yi.r, dyi.r);
    tmpz.r = _mm_mul_ps(zi.r, dzi.r);

    k2.r = _mm_cvtps_epi32(tmpy.r);
    k3.r = _mm_cvtps_epi32(tmpz.r);

    for(int m = 0; m < 4; m++){
      k2.v[m] = round(tmpy.v[m]);
      k3.v[m] = round(tmpz.v[m]); 
    }


    k2fl.r = _mm_cvtepi32_ps(l2.r);
    k3fl.r = _mm_cvtepi32_ps(l3.r);

    tmpy.r = _mm_sub_ps(k2fl.r, tmpy.r);
    tmpz.r = _mm_sub_ps(k3fl.r, tmpz.r);
    
    // God help me, there's some things I just can't figure out how to parallelize
    // The g-- here are just temporary variables. I can't do the assignments in parallel, 
    // but I'll be damned if I can't do the FLOPS in parallel
    gmy.r = _mm_sub_ps(tmpy.r,ones.r);
    gmy.r = _mm_add_ps(onepfive.r, gmy.r); // h2-1 always <0
    gmy.r = _mm_mul_ps(gmy.r, gmy.r);
    gmy.r = _mm_mul_ps(half.r, gmy.r);

    g0y.r = _mm_mul_ps(tmpy.r, tmpy.r);
    g0y.r = _mm_sub_ps(threefourths.r, g0y.r);

    g1y.r = _mm_add_ps(ones.r, tmpy.r); 
    g1y.r = _mm_sub_ps(onepfive.r, g1y.r);//h2+1 always >0
    g1y.r = _mm_mul_ps(g1y.r, g1y.r);
    g1y.r = _mm_mul_ps(half.r, g1y.r);

    gmz.r = _mm_sub_ps(tmpz.r,ones.r);
    gmz.r = _mm_add_ps(onepfive.r, gmz.r); //h3-1 always <0
    gmz.r = _mm_mul_ps(gmz.r, gmz.r);
    gmz.r = _mm_mul_ps(half.r, gmz.r);

    g0z.r = _mm_mul_ps(tmpz.r, tmpz.r);
    g0z.r = _mm_sub_ps(threefourths.r, g0z.r);

    g1z.r = _mm_add_ps(ones.r, tmpz.r);
    g1z.r = _mm_sub_ps(onepfive.r, g1z.r); //h3+1 always >0
    g1z.r = _mm_mul_ps(g1z.r, g1z.r);
    g1z.r = _mm_mul_ps(half.r, g1z.r);

    // All the indices below have +2 compared to the FORTRAN for C array access
    for(int p=0; p<4; p++){
      int shift = k2.v[p] - j2.v[p];
      s1y[shift + 1].v[p] = gmy.v[p];
      s1y[shift + 2].v[p] = g0y.v[p];
      s1y[shift + 3].v[p] = g1y.v[p];
      s1z[shift + 1].v[p] = gmz.v[p];
      s1z[shift + 2].v[p] = g0z.v[p];
      s1z[shift + 3].v[p] = g1z.v[p];
    }

// CHECKPOINT: PIC_push_part_yz.F : line 341
// Charge density form factor at (n+1.5)*dt

    for(int m=0; m<5; m++){
      s1y[m].r = _mm_sub_ps(s1y[m].r, s0y[m].r);
      s1z[m].r = _mm_sub_ps(s1z[m].r, s0z[m].r);
    }
    
    pvFloat fnqx, fnqy, fnqz;
    
    fnqx.r = _mm_mul_ps(wni.r, fnqs.r);
    fnqx.r = _mm_mul_ps(qni.r, fnqs.r);
    fnqx.r = _mm_mul_ps(vxi.r, fnqs.r);

    fnqy.r = _mm_mul_ps(wni.r, fnqys.r);
    fnqy.r = _mm_mul_ps(qni.r, fnqy.r);
    
    fnqz.r = _mm_mul_ps(wni.r, fnqzs.r);
    fnqz.r = _mm_mul_ps(qni.r, fnqz.r);
    
    // Again, this is extremely painful, but serial is the only way I can think of to do this
    // I have the inkling of thought on a parallel way, but it hasn't matured. If it works out, I'll implement it.
    for(int p=0; p<4; p++){
      int l2min, l2max, l3min, l3max;
      if (k2.v[p] == j2.v[p]){
    	l2min = -1;
    	l2max = 1;
      } else if(k2.v[p] == j2.v[p]-1){
    	l2min = -2;
    	l2max = 1;
      }else if(k2.v[p] == j2.v[p]+1){
    	l2min = -1;
    	l2max = 2;
      }
      if (k3.v[p] == j3.v[p]){
    	l3min = -1;
    	l3max = 1;
      } else if(k3.v[p] == j3.v[p]-1){
    	l3min = -2;
    	l3max = 1;
      } else if(k3.v[p] == j3.v[p]+1){
    	l3min = -1;
    	l3max = 2;
      }

      float jxh[5][5] = {{0.0}};
      float jyh[6][5] = {{0.0}};
      float jzh[5][6] = {{0.0}}; // Slimming these down a bit...
   
      for(int l3i=l3min+2; l3i<=(l3max+2); l3i++){
    	for(int l2i=l2min+2; l2i<=(l2max+2); l2i++){
    	  float wx = (s0y[l2i].v[p] + 0.5*s1y[l2i].v[p])*s0z[l3i].v[p]
    	    + (0.5*s0y[l2i].v[p] + (1./3.)*s1y[l2i].v[p])*s1z[l3i].v[p];
    	  float wy = s1y[l2i].v[p]*(s0z[l3i].v[p] + 0.5*s1z[l3i].v[p]);
    	  float wz = s1z[l3i].v[p]*(s0y[l2i].v[p] + 0.5*s1y[l2i].v[p]);
	  
    	  jxh[l2i][l3i] = fnqx.v[p]*wx;
    	  jyh[l2i+1][l3i] = jyh[l2i][l3i] - fnqy.v[p]*wy; // [+1][] index adjustment
    	  jzh[l2i][l3i+1] = jzh[l2i][l3i] - fnqz.v[p]*wz; //[][+1] index adjustment
	  
    	  /* if(n == 0){ */
    	  /*   fprintf(stderr, "%g %g %g \n", C_FIELD(JXI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2), C_FIELD(JYI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2), C_FIELD(JZI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2)); */
    	  /* } */

	  
    	  C_FIELD(JXI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += jxh[l2i][l3i]; //undo index adjustment of l(2,3)i
    	  C_FIELD(JYI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += jyh[l2i+1][l3i]; //undo index adjustment of l(2,3)i
    	  C_FIELD(JZI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += jzh[l2i][l3i+1]; //undo index adjustment of l(2,3)i
    	}
      }
    }   
  }
  prof_stop(pr);
}
