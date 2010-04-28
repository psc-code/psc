#include "psc_sse2.h"

#include "profile/profile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>

void
sse2_push_part_yz_b()
{
  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz_b", 1., 0, psc.n_part * 12 * sizeof(float));
  }
  prof_start(pr);


//-----------------------------------------------------
// Initialization stuff (not sure what all of this is for)
  
  struct psc_sse2 *sse2 = psc.c_ctx;

  psc.p2A = 0.;
  psc.p2B = 0.;


  // Values that won't change from iteration to iteration

  union packed_vector dt, yl, zl, ones, half,threefourths, eta, dqs,fnqs, dxi, dyi, dzi; 
  //FIXME: These are all stored as doubles in fortran!!
  float dtfl = psc.dt; 
  float etafl =  psc.prm.eta;
  float fnqsfl = sqr(psc.prm.alpha) * psc.prm.cori / etafl;
  float dxifl = 1.0f / psc.dx[0];
  float dyifl = 1.0f / psc.dx[1];
  float dzifl = 1.0f / psc.dx[2];
  float dqsfl = 0.5*etafl*dtfl;
  ones.r = _mm_set1_ps(1.0f);
  half.r = _mm_set1_ps(.5f);
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


  //  assert(psc.n_part % 4 == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += 4) {
    union packed_vector pxi, pyi, pzi, xi, yi, zi, qni, mni, wni;
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
    union packed_vector vxi, vyi, vzi, tmpx, tmpy, tmpz, root;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta
    
    tmpx.r = _mm_mul_ps(pxi.r, pxi.r);
    tmpy.r = _mm_mul_ps(pyi.r, pyi.r);
    tmpz.r = _mm_mul_ps(pzi.r, pzi.r);

    tmpx.r = _mm_add_ps(tmpx.r, tmpy.r);
    tmpx.r = _mm_add_ps(tmpx.r, tmpz.r);
    tmpx.r = _mm_add_ps(ones.r, tmpx.r);
    root.r = _mm_rsqrt_ps(tmpx.r);

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

// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    // Apparently this can be done in vectors. A victory!
      union {
      __m128i r;
      int v[4];
    } j1, j2, j3, l1, l2, l3;

    union packed_vector j2fl, j3fl, l2fl, l3fl;
    
    //Go go gadget loop unwind! Replaces u,v,w in F90 code
    j1.r = _mm_cvtps_epi32(tmpx.r);
    j2.r = _mm_cvtps_epi32(tmpy.r);
    j3.r = _mm_cvtps_epi32(tmpz.r);

    // there must be a better way...
    j2fl.r = _mm_cvtepi32_ps(j2.r);
    j3fl.r = _mm_cvtepi32_ps(j3.r);

    tmpy.r = _mm_sub_ps(j2fl.r, tmpy.r);
    tmpz.r = _mm_sub_ps(j3fl.r, tmpy.r);


    union packed_vector gmy, gmz, g0y, g0z, g1y, g1z;

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

// CHECKPOINT: PIC_push_part_yz.F : line 143

    tmpx.r = _mm_mul_ps(xi.r, dxi.r);
    tmpy.r = _mm_mul_ps(yi.r, dyi.r);
    tmpy.r = _mm_sub_ps(tmpy.r, half.r);
    tmpz.r = _mm_mul_ps(zi.r, dzi.r);
    tmpz.r = _mm_mul_ps(tmpz.r, half.r);

    l1.r = _mm_cvtps_epi32(tmpx.r);
    l2.r = _mm_cvtps_epi32(tmpy.r);
    l3.r = _mm_cvtps_epi32(tmpz.r);

    l2fl.r = _mm_cvtepi32_ps(l2.r);
    l3fl.r = _mm_cvtepi32_ps(l3.r);

    tmpy.r = _mm_sub_ps(l2fl.r, tmpy.r);
    tmpz.r = _mm_sub_ps(l3fl.r, tmpy.r);

    union packed_vector hmy, hmz, h0y, h0z, h1y, h1z;

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

#define C_FIELD(T,l,j,k) sse2->fields[(T)*psc.fld_size + ((l)-psc.ilo[0]+psc.ibn[0]) + ((j)-psc.ilo[1]+psc.ibn[1])*(psc.img[0]) + ((k) - psc.ilo[2]+psc.ibn[2])*(psc.img[0]+psc.img[1])]

    union packed_vector field_in, exq, eyq, ezq, bxq, byq, bzq;
     
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

#undef C_FIELD

    /* fprintf(stderr, "bxq: %g %g %g %g \n", bxq.v[0], bxq.v[1], bxq.v[2], bxq.v[3]); */
    /* fprintf(stderr, "byq: %g %g %g %g \n", byq.v[0], byq.v[1], byq.v[2], byq.v[3]); */
    /* fprintf(stderr, "bzq: %g %g %g %g \n", bzq.v[0], bzq.v[1], bzq.v[2], bzq.v[3]); */
    /* fprintf(stderr, "exq: %g %g %g %g \n", exq.v[0], exq.v[1], exq.v[2], exq.v[3]); */
    /* fprintf(stderr, "eyq: %g %g %g %g \n", eyq.v[0], eyq.v[1], eyq.v[2], eyq.v[3]); */
    /* fprintf(stderr, "ezq: %g %g %g %g \n", ezq.v[0], ezq.v[1], ezq.v[2], ezq.v[3]); */


// CHECKPOINT: PIC_push_part_yz.F : line 223
// Half step momentum  with E-field

    union packed_vector dq;
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
    union packed_vector txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;

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

    union packed_vector tau;

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
    root.r = _mm_rsqrt_ps(tmpx.r);

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

  }
  prof_stop(pr);
}
