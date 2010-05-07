#include "psc_sse2.h"
#include "sse2_cgen.h"
#include "profile/profile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

void
sse2_push_part_yz_a()
{
  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz", 1., 0, psc.n_part * 12 * sizeof(sse2_real));
  }
  prof_start(pr);

//-----------------------------------------------------
// Initialization stuff (not sure what all of this is for)
  
  struct psc_sse2 *sse2 = psc.c_ctx;


  pvReal dt, yl, zl, ones, half;
  //FIXME: These are all stored as doubles in fortran!!
  sse2_real dtfl = psc.dt; 
  ones.r = pv_real_SET1(1.0);
  half.r = pv_real_SET1(.5);
  dt.r = pv_real_SET1(dtfl);
  yl.r = pv_real_MUL(half.r, dt.r);
  zl.r = pv_real_MUL(half.r, dt.r);

  assert(psc.n_part % VEC_SIZE == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {

//---------------------------------------------
// Bringing in particle specific parameters
    pvReal pxi, pyi, pzi, xi, yi, zi, qni, mni, wni; 
 
    LOAD_PART(sse2,n); 

    // Locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, tmpz, root;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta
    
    tmpx.r = pv_real_MUL(pxi.r, pxi.r);
    tmpy.r = pv_real_MUL(pyi.r, pyi.r);
    tmpz.r = pv_real_MUL(pzi.r, pzi.r);

    tmpx.r = pv_real_ADD(tmpx.r, tmpy.r);
    tmpx.r = pv_real_ADD(tmpx.r, tmpz.r);
    tmpx.r = pv_real_ADD(ones.r, tmpx.r);
    root.r = pv_real_SQRT(tmpx.r);
    root.r = pv_real_DIV(ones.r, root.r);

    vxi.r = pv_real_MUL(pxi.r, root.r);
    vyi.r = pv_real_MUL(pyi.r, root.r);
    vzi.r = pv_real_MUL(pzi.r, root.r);

    tmpy.r = pv_real_MUL(vyi.r, yl.r);
    tmpz.r = pv_real_MUL(vzi.r, zl.r);
    yi.r = pv_real_ADD(yi.r, tmpy.r);
    zi.r = pv_real_ADD(zi.r, tmpz.r);

    STORE_PART_XP(sse2,n);
  }
  prof_stop(pr);
}

void
sse2_push_part_yz_b()
{
  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz", 1., 0, psc.n_part * 12 * sizeof(sse2_real));
  }
  prof_start(pr);

//-----------------------------------------------------
// Initialization stuff (not sure what all of this is for)
  
  struct psc_sse2 *sse2 = psc.c_ctx;

  psc.p2A = 0.;
  psc.p2B = 0.;
  
  // Values that won't change from iteration to iteration

  pvReal dt, yl, zl, ones, half,threefourths, onepfive,  eta, dqs,fnqs, dxi, dyi, dzi; 
  //FIXME: These are all stored as doubles in fortran!!
  sse2_real dtfl = psc.dt; 
  sse2_real etafl =  psc.prm.eta;
  sse2_real fnqsfl = sqr(psc.prm.alpha) * psc.prm.cori / etafl;
  sse2_real dxifl = 1.0 / psc.dx[0];
  sse2_real dyifl = 1.0 / psc.dx[1];
  sse2_real dzifl = 1.0 / psc.dx[2];
  sse2_real dqsfl = 0.5*etafl*dtfl;
  ones.r = pv_real_SET1(1.0);
  half.r = pv_real_SET1(.5);
  onepfive.r = pv_real_SET1(1.5);
  threefourths.r = pv_real_SET1(.75);
  dt.r = pv_real_SET1(dtfl);
  eta.r = pv_real_SET1(etafl);
  fnqs.r = pv_real_SET1(fnqsfl);
  dqs.r = pv_real_SET1(dqsfl);
  dxi.r = pv_real_SET1(dxifl);
  dyi.r = pv_real_SET1(dyifl);
  dzi.r = pv_real_SET1(dzifl);
  yl.r = pv_real_MUL(half.r, dt.r);
  zl.r = pv_real_MUL(half.r, dt.r);

  pvInt ione;
  
  ione.r = pv_int_SET1(1);


  //  assert(psc.n_part % VEC_SIZE == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {

//---------------------------------------------
// Bringing in particle specific parameters
    pvReal pxi, pyi, pzi, xi, yi, zi, qni, mni, wni; 
 
    LOAD_PART(sse2,n); 

    // Locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, tmpz, root, h2, h3;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta
    
    tmpx.r = pv_real_MUL(pxi.r, pxi.r);
    tmpy.r = pv_real_MUL(pyi.r, pyi.r);
    tmpz.r = pv_real_MUL(pzi.r, pzi.r);

    tmpx.r = pv_real_ADD(tmpx.r, tmpy.r);
    tmpx.r = pv_real_ADD(tmpx.r, tmpz.r);
    tmpx.r = pv_real_ADD(ones.r, tmpx.r);
    root.r = pv_real_SQRT(tmpx.r);
    root.r = pv_real_DIV(ones.r, root.r);

    vxi.r = pv_real_MUL(pxi.r, root.r);
    vyi.r = pv_real_MUL(pyi.r, root.r);
    vzi.r = pv_real_MUL(pzi.r, root.r);

    tmpy.r = pv_real_MUL(vyi.r, yl.r);
    tmpz.r = pv_real_MUL(vzi.r, zl.r);
    yi.r = pv_real_ADD(yi.r, tmpy.r);
    zi.r = pv_real_ADD(zi.r, tmpz.r);


    tmpx.r = pv_real_DIV(ones.r, root.r);
    tmpx.r = pv_real_SUB(tmpx.r, ones.r);
    tmpx.r = pv_real_DIV(tmpx.r, eta.r);
    tmpx.r = pv_real_MUL(fnqs.r, tmpx.r);
    tmpx.r = pv_real_MUL(mni.r, tmpx.r);
    
    for(int p = 0; p < VEC_SIZE; p++){
      psc.p2A += tmpx.v[p]; //What's this for?
    }

// CHECKPOINT: PIC_push_part_yz.F : line 110
    
    tmpx.r = pv_real_MUL(xi.r, dxi.r);
    tmpy.r = pv_real_MUL(yi.r, dyi.r);
    tmpz.r = pv_real_MUL(zi.r, dzi.r);


// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    // Apparently this can be done in vectors. A victory!
    pvInt j1, j2, j3, l1, l2, l3;
    pvReal j2fl, j3fl, l2fl, l3fl;

    //FIXME: Choose one!

    for(int m = 0; m < VEC_SIZE; m++){
      j1.v[m] = round(tmpx.v[m]);
    }

    j2.r = pv_real_to_int_CVT(tmpy.r);
    j3.r = pv_real_to_int_CVT(tmpz.r);
    
    // there must be a better way...
    j2fl.r = pv_int_to_real_CVT(j2.r);
    j3fl.r = pv_int_to_real_CVT(j3.r);

    h2.r = pv_real_SUB(j2fl.r, tmpy.r);
    h3.r = pv_real_SUB(j3fl.r, tmpz.r);

    pvInt j2pls1, j2mns1, j3pls1, j3mns1;
    j2pls1.r = pv_int_ADD(j2.r, ione.r);
    j2mns1.r = pv_int_SUB(j2.r, ione.r);
    j3pls1.r = pv_int_ADD(j3.r, ione.r);
    j3mns1.r = pv_int_SUB(j3.r, ione.r);


    pvReal gmy, gmz, gOy, gOz, gly, glz;// NB: this is gO as in octogon instead of g0, 
                                         // and gl as in library instead of g1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of g?[3], which 
                                         // seems to improve performance


    gmy.r = pv_real_ADD(half.r, h2.r);
    gmy.r = pv_real_MUL(gmy.r, gmy.r);
    gmy.r = pv_real_MUL(half.r, gmy.r);

    gmz.r = pv_real_ADD(half.r, h3.r);
    gmz.r = pv_real_MUL(gmz.r, gmz.r);
    gmz.r = pv_real_MUL(half.r, gmz.r);
    
    gOy.r = pv_real_MUL(h2.r, h2.r);
    gOy.r = pv_real_SUB(threefourths.r, gOy.r);

    gOz.r = pv_real_MUL(h3.r, h3.r);
    gOz.r = pv_real_SUB(threefourths.r, gOz.r);

    gly.r = pv_real_SUB(half.r, h2.r);
    gly.r = pv_real_MUL(gly.r, gly.r);
    gly.r = pv_real_MUL(half.r, gly.r);

    glz.r = pv_real_SUB(half.r, h3.r);
    glz.r = pv_real_MUL(glz.r, glz.r);
    glz.r = pv_real_MUL(half.r, glz.r);


// CHECKPOINT: PIC_push_part_yz.F : line 143

    tmpx.r = pv_real_MUL(xi.r, dxi.r);
    tmpy.r = pv_real_MUL(yi.r, dyi.r);
    tmpy.r = pv_real_SUB(tmpy.r, half.r);
    tmpz.r = pv_real_MUL(zi.r, dzi.r);
    tmpz.r = pv_real_SUB(tmpz.r, half.r);

    for(int m = 0; m < VEC_SIZE; m++){ // Cumbersome, but neccessary. Sometimes the intrinsic
      l1.v[m] = round(tmpx.v[m]);     // rounds .5 the wrong way. It's okay for the other
    }                                 // directions, but there's no weighting in x


    l2.r = pv_real_to_int_CVT(tmpy.r);
    l3.r = pv_real_to_int_CVT(tmpz.r);

    l2fl.r = pv_int_to_real_CVT(l2.r);
    l3fl.r = pv_int_to_real_CVT(l3.r);

    h2.r = pv_real_SUB(l2fl.r, tmpy.r);
    h3.r = pv_real_SUB(l3fl.r, tmpz.r);

    pvInt l2pls1, l2mns1, l3pls1, l3mns1;
    l2pls1.r = pv_int_ADD(l2.r, ione.r);
    l2mns1.r = pv_int_SUB(l2.r, ione.r);
    l3pls1.r = pv_int_ADD(l3.r, ione.r);
    l3mns1.r = pv_int_SUB(l3.r, ione.r);

    pvReal hmy, hmz, hOy, hOz, hly, hlz; // NB: this is hO as in octogon instead of h0, 
                                         // and hl as in library instead of h1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of h?[3], which 
                                         // seems to improve performancs

    hmy.r = pv_real_ADD(half.r, h2.r);
    hmy.r = pv_real_MUL(hmy.r, hmy.r);
    hmy.r = pv_real_MUL(half.r, hmy.r);

    hmz.r = pv_real_ADD(half.r, h3.r);
    hmz.r = pv_real_MUL(hmz.r, hmz.r);
    hmz.r = pv_real_MUL(half.r, hmz.r);
    
    hOy.r = pv_real_MUL(h2.r, h2.r);
    hOy.r = pv_real_SUB(threefourths.r, hOy.r);

    hOz.r = pv_real_MUL(h3.r, h3.r);
    hOz.r = pv_real_SUB(threefourths.r, hOz.r);

    hly.r = pv_real_SUB(half.r, h2.r);
    hly.r = pv_real_MUL(hly.r, hly.r);
    hly.r = pv_real_MUL(half.r, hly.r);

    hlz.r = pv_real_SUB(half.r, h3.r);
    hlz.r = pv_real_MUL(hlz.r, hlz.r);
    hlz.r = pv_real_MUL(half.r, hlz.r);

// CHECKPOINT: PIC_push_part_yz.F : line 59
// Field Interpolation

// BOTTLENECK HERE: For a given particle, none of the needed fields are stored
// next to each other in memory. Also, there's no way to tell before runtime
// exactly what fields are needed. As a consequence, I can't just load them into
// registers and shuffle. I must fetch piece by painful piece.

//Hopefully gcc is smart enough to unroll these for me. Otherwise, I need to figure out
// a way to generate unrolled loops
    
    pvReal exq, eyq, ezq, bxq, byq, bzq;

/*     pvReal gy[3] = {gmy, gOy, gly}; */
/*     pvReal gz[3] = {gmz, gOz, glz}; */
/*     pvReal hy[3] = {hmy, hOy, hly}; */
/*     pvReal hz[3] = {hmz, hOz, hlz}; */


    INTERP_FIELD_YZ(EX,l1,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,j1,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j1,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(BX,j1,l2,l3,h,h,bxq);
    INTERP_FIELD_YZ(BY,l1,j2,l3,h,g,byq);
    INTERP_FIELD_YZ(BZ,l1,l2,j3,g,h,bzq);

 
// CHECKPOINT: PIC_push_part_yz.F : line 223
// Half step momentum  with E-field

    pvReal dq, dqex, dqey, dqez;
    dq.r = pv_real_DIV(dqs.r, mni.r);
    dq.r = pv_real_MUL(qni.r, dq.r);
    
    dqex.r = pv_real_MUL(dq.r, exq.r);
    dqey.r = pv_real_MUL(dq.r, eyq.r);
    dqez.r = pv_real_MUL(dq.r, ezq.r);
    
    pxi.r = pv_real_ADD(pxi.r, dqex.r);
    pyi.r = pv_real_ADD(pyi.r, dqey.r);
    pzi.r = pv_real_ADD(pzi.r, dqez.r);


// CHECKPOINT: PIC_push_part_yz.F : line 228
// Rotate with B-field
    pvReal taux, tauy, tauz, txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;

    tmpx.r = pv_real_MUL(pxi.r, pxi.r);
    tmpy.r = pv_real_MUL(pyi.r, pyi.r);
    tmpz.r = pv_real_MUL(pzi.r, pzi.r);
    root.r = pv_real_ADD(ones.r, tmpx.r);
    tmpy.r = pv_real_ADD(tmpy.r, tmpz.r);
    root.r = pv_real_ADD(root.r, tmpy.r);
    root.r = pv_real_SQRT(root.r);
    root.r = pv_real_DIV(dq.r, root.r);

    taux.r = pv_real_MUL(bxq.r, root.r);
    tauy.r = pv_real_MUL(byq.r, root.r);
    tauz.r = pv_real_MUL(bzq.r, root.r);

    txx.r = pv_real_MUL(taux.r, taux.r);
    tyy.r = pv_real_MUL(tauy.r, tauy.r);
    tzz.r = pv_real_MUL(tauz.r, tauz.r);
    t2xy.r = pv_real_MUL(taux.r, tauy.r);
    t2xz.r = pv_real_MUL(taux.r, tauz.r);
    t2yz.r = pv_real_MUL(tauy.r, tauz.r);
    t2xy.r = pv_real_ADD(t2xy.r, t2xy.r);
    t2xz.r = pv_real_ADD(t2xz.r, t2xz.r);
    t2yz.r = pv_real_ADD(t2yz.r, t2yz.r);

    pvReal tau;

    tau.r = pv_real_ADD(ones.r, txx.r);
    tmpx.r = pv_real_ADD(tyy.r, tzz.r);
    tau.r = pv_real_ADD(tau.r, tmpx.r);
    tau.r = pv_real_DIV(ones.r, tau.r); //recp is evil! evil evil evil evil!!!!
    
    //Never use tau_ without a two in front
    taux.r = pv_real_ADD(taux.r, taux.r);
    tauy.r = pv_real_ADD(tauy.r, tauy.r);
    tauz.r = pv_real_ADD(tauz.r, tauz.r);
    
    //pxp
    tmpx.r = pv_real_ADD(ones.r, txx.r);
    tmpx.r = pv_real_SUB(tmpx.r, tyy.r);
    tmpx.r = pv_real_SUB(tmpx.r, tzz.r);
    tmpx.r = pv_real_MUL(tmpx.r, pxi.r);
    
    tmpy.r = pv_real_ADD(t2xy.r, tauz.r);
    tmpy.r = pv_real_MUL(tmpy.r, pyi.r);

    tmpz.r = pv_real_SUB(t2xz.r, tauy.r);
    tmpz.r = pv_real_MUL(tmpz.r, pzi.r);

    pxp.r = pv_real_ADD(tmpx.r, tmpy.r);
    pxp.r = pv_real_ADD(pxp.r, tmpz.r);
    pxp.r = pv_real_MUL(pxp.r, tau.r);
    
    //pyp
    tmpx.r = pv_real_SUB(t2xy.r, tauz.r);
    tmpx.r = pv_real_MUL(tmpx.r, pxi.r);

    tmpy.r = pv_real_SUB(ones.r, txx.r);
    tmpy.r = pv_real_ADD(tmpy.r, tyy.r);
    tmpy.r = pv_real_SUB(tmpy.r, tzz.r);
    tmpy.r = pv_real_MUL(tmpy.r, pyi.r);
    
    tmpz.r = pv_real_ADD(t2yz.r, taux.r);
    tmpz.r = pv_real_MUL(tmpz.r, pzi.r);

    pyp.r = pv_real_ADD(tmpx.r, tmpy.r);
    pyp.r = pv_real_ADD(pyp.r, tmpz.r);
    pyp.r = pv_real_MUL(pyp.r, tau.r);
   
    //pzp
    tmpx.r = pv_real_ADD(t2xz.r, tauy.r);
    tmpx.r = pv_real_MUL(tmpx.r, pxi.r);

    tmpy.r = pv_real_SUB(t2yz.r, taux.r);
    tmpy.r = pv_real_MUL(tmpy.r, pyi.r);

    tmpz.r = pv_real_SUB(ones.r, txx.r);
    tmpz.r = pv_real_SUB(tmpz.r, tyy.r);
    tmpz.r = pv_real_ADD(tmpz.r, tzz.r);
    tmpz.r = pv_real_MUL(tmpz.r, pzi.r);

    pzp.r = pv_real_ADD(tmpx.r, tmpy.r);
    pzp.r = pv_real_ADD(pzp.r, tmpz.r);
    pzp.r = pv_real_MUL(pzp.r, tau.r);

// CHECKPOINT: PIC_push_part_yz.F : line 244
// Half step momentum  with E-field

    pxi.r = pv_real_ADD(pxp.r, dqex.r);
    pyi.r = pv_real_ADD(pyp.r, dqey.r);
    pzi.r = pv_real_ADD(pzp.r, dqez.r);

// CHECKPOINT: PIC_push_part_yz.F : line 248
// Half step particles with new momenta
    
    tmpx.r = pv_real_MUL(pxi.r, pxi.r);
    tmpy.r = pv_real_MUL(pyi.r, pyi.r);
    tmpz.r = pv_real_MUL(pzi.r, pzi.r);

    tmpx.r = pv_real_ADD(tmpx.r, tmpy.r);
    tmpx.r = pv_real_ADD(tmpx.r, tmpz.r);
    tmpx.r = pv_real_ADD(ones.r, tmpx.r);
    root.r = pv_real_SQRT(tmpx.r);
    root.r = pv_real_DIV(ones.r, root.r);

    vxi.r = pv_real_MUL(pxi.r, root.r);
    vyi.r = pv_real_MUL(pyi.r, root.r);
    vzi.r = pv_real_MUL(pzi.r, root.r);
    
    tmpy.r = pv_real_MUL(vyi.r, yl.r);
    tmpz.r = pv_real_MUL(vzi.r, zl.r);
    
    yi.r = pv_real_ADD(yi.r, tmpy.r);
    zi.r = pv_real_ADD(zi.r, tmpz.r);

    STORE_PART_XP(sse2,n);

    tmpx.r = pv_real_DIV(ones.r, root.r);
    tmpx.r = pv_real_SUB(tmpx.r, ones.r);
    tmpx.r = pv_real_DIV(tmpx.r, eta.r);
    tmpx.r = pv_real_MUL(fnqs.r, tmpx.r);
    tmpx.r = pv_real_MUL(mni.r, tmpx.r);

    for(int p = 0; p < VEC_SIZE; p++){
      psc.p2B += tmpx.v[p]; //What's this for?
    }

  }
  prof_stop(pr);
}

void
sse2_push_part_yz()
{
  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz", 1., 0, psc.n_part * 12 * sizeof(sse2_real));
  }
  prof_start(pr);

//-----------------------------------------------------
// Initialization stuff (not sure what all of this is for)
  
  struct psc_sse2 *sse2 = psc.c_ctx;

  psc.p2A = 0.;
  psc.p2B = 0.;
  
  for (int m = NE; m <= JZI; m++) {
    memset(&sse2->fields[m*psc.fld_size], 0, psc.fld_size * sizeof(sse2_real));
  }  

  // Values that won't change from iteration to iteration

  pvReal dt, yl, zl, ones, half,threefourths, onepfive,  eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs; 
  //FIXME: These are all stored as doubles in fortran!!
  sse2_real dtfl = psc.dt; 
  sse2_real etafl =  psc.prm.eta;
  sse2_real fnqsfl = sqr(psc.prm.alpha) * psc.prm.cori / etafl;
  sse2_real dxifl = 1.0 / psc.dx[0];
  sse2_real dyifl = 1.0 / psc.dx[1];
  sse2_real dzifl = 1.0 / psc.dx[2];
  sse2_real dqsfl = 0.5*etafl*dtfl;
  sse2_real fnqxsfl = psc.dx[0] * fnqsfl * psc.dt;
  sse2_real fnqysfl = psc.dx[1] * fnqsfl * psc.dt;
  sse2_real fnqzsfl = psc.dx[2] * fnqsfl * psc.dt;
  ones.r = pv_real_SET1(1.0);
  half.r = pv_real_SET1(.5);
  onepfive.r = pv_real_SET1(1.5);
  threefourths.r = pv_real_SET1(.75);
  dt.r = pv_real_SET1(dtfl);
  eta.r = pv_real_SET1(etafl);
  fnqs.r = pv_real_SET1(fnqsfl);
  dqs.r = pv_real_SET1(dqsfl);
  dxi.r = pv_real_SET1(dxifl);
  dyi.r = pv_real_SET1(dyifl);
  dzi.r = pv_real_SET1(dzifl);
  yl.r = pv_real_MUL(half.r, dt.r);
  zl.r = pv_real_MUL(half.r, dt.r);
  fnqxs.r = pv_real_SET1(fnqxsfl);
  fnqys.r = pv_real_SET1(fnqysfl);
  fnqzs.r = pv_real_SET1(fnqzsfl);

  pvInt ione;
    
  ione.r = pv_int_SET1(1);

  assert(psc.n_part % VEC_SIZE == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {

//---------------------------------------------
// Bringing in particle specific parameters

    pvReal pxi, pyi, pzi, xi, yi, zi, qni, mni, wni;  
    LOAD_PART(sse2,n); 

    // Locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, tmpz, root, h2, h3;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta
    
    tmpx.r = pv_real_MUL(pxi.r, pxi.r);
    tmpy.r = pv_real_MUL(pyi.r, pyi.r);
    tmpz.r = pv_real_MUL(pzi.r, pzi.r);

    tmpx.r = pv_real_ADD(tmpx.r, tmpy.r);
    tmpx.r = pv_real_ADD(tmpx.r, tmpz.r);
    tmpx.r = pv_real_ADD(ones.r, tmpx.r);
    root.r = pv_real_SQRT(tmpx.r);
    root.r = pv_real_DIV(ones.r, root.r);

    vxi.r = pv_real_MUL(pxi.r, root.r);
    vyi.r = pv_real_MUL(pyi.r, root.r);
    vzi.r = pv_real_MUL(pzi.r, root.r);

    tmpy.r = pv_real_MUL(vyi.r, yl.r);
    tmpz.r = pv_real_MUL(vzi.r, zl.r);
    yi.r = pv_real_ADD(yi.r, tmpy.r);
    zi.r = pv_real_ADD(zi.r, tmpz.r);


    tmpx.r = pv_real_DIV(ones.r, root.r);
    tmpx.r = pv_real_SUB(tmpx.r, ones.r);
    tmpx.r = pv_real_DIV(tmpx.r, eta.r);
    tmpx.r = pv_real_MUL(fnqs.r, tmpx.r);
    tmpx.r = pv_real_MUL(mni.r, tmpx.r);
    
    for(int p = 0; p < VEC_SIZE; p++){
      psc.p2A += tmpx.v[p]; //What's this for?
    }

// CHECKPOINT: PIC_push_part_yz.F : line 110
    
    tmpx.r = pv_real_MUL(xi.r, dxi.r);
    tmpy.r = pv_real_MUL(yi.r, dyi.r);
    tmpz.r = pv_real_MUL(zi.r, dzi.r);

    pvReal s0y[5], s0z[5], s1y[5], s1z[5];
    // I'm a little scared to use memset here, though it would probably work...
    s0y[0].r = pv_real_SET1(0.0);
    s0y[4].r = pv_real_SET1(0.0);
    s0z[0].r = pv_real_SET1(0.0);
    s0z[4].r = pv_real_SET1(0.0);
    for(int m = 0; m < 5; m++){
      s1y[m].r = pv_real_SET1(0.0);
      s1z[m].r = pv_real_SET1(0.0);
    }

// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    // Apparently this can be done in vectors. A victory!
    pvInt j1, j2, j3, l1, l2, l3;
    pvReal j2fl, j3fl, l2fl, l3fl;

    //FIXME: Choose one!
    for(int m = 0; m < VEC_SIZE; m++){
      j1.v[m] = round(tmpx.v[m]);
    }

    j2.r = pv_real_to_int_CVT(tmpy.r);
    j3.r = pv_real_to_int_CVT(tmpz.r);

    // there must be a better way...
    j2fl.r = pv_int_to_real_CVT(j2.r);
    j3fl.r = pv_int_to_real_CVT(j3.r);

    h2.r = pv_real_SUB(j2fl.r, tmpy.r);
    h3.r = pv_real_SUB(j3fl.r, tmpz.r);

    pvInt j2pls1, j2mns1, j3pls1, j3mns1;
    j2pls1.r = pv_int_ADD(j2.r, ione.r);
    j2mns1.r = pv_int_SUB(j2.r, ione.r);
    j3pls1.r = pv_int_ADD(j3.r, ione.r);
    j3mns1.r = pv_int_SUB(j3.r, ione.r);


    pvReal gmy, gmz, gOy, gOz, gly, glz; // NB: this is gO as in octogon instead of g0, 
                                         // and gl as in library instead of g1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of g?[3], which 
                                         // seems to improve performance


    gmy.r = pv_real_ADD(half.r, h2.r);
    gmy.r = pv_real_MUL(gmy.r, gmy.r);
    gmy.r = pv_real_MUL(half.r, gmy.r);

    gmz.r = pv_real_ADD(half.r, h3.r);
    gmz.r = pv_real_MUL(gmz.r, gmz.r);
    gmz.r = pv_real_MUL(half.r, gmz.r);
    
    gOy.r = pv_real_MUL(h2.r, h2.r);
    gOy.r = pv_real_SUB(threefourths.r, gOy.r);

    gOz.r = pv_real_MUL(h3.r, h3.r);
    gOz.r = pv_real_SUB(threefourths.r, gOz.r);

    gly.r = pv_real_SUB(half.r, h2.r);
    gly.r = pv_real_MUL(gly.r, gly.r);
    gly.r = pv_real_MUL(half.r, gly.r);

    glz.r = pv_real_SUB(half.r, h3.r);
    glz.r = pv_real_MUL(glz.r, glz.r);
    glz.r = pv_real_MUL(half.r, glz.r);


    // indexing here departs from FORTRAN a little bit
    s0y[1].r = pv_real_SUB(h2.r, ones.r);
    s0y[1].r = pv_real_ADD(onepfive.r, s0y[1].r); // h2+1.0 always <0
    s0y[1].r = pv_real_MUL(s0y[1].r, s0y[1].r);
    s0y[1].r = pv_real_MUL(half.r, s0y[1].r);

    s0y[2].r = pv_real_MUL(h2.r, h2.r);
    s0y[2].r = pv_real_SUB(threefourths.r, s0y[2].r);

    s0y[3].r = pv_real_ADD(h2.r, ones.r);
    s0y[3].r = pv_real_SUB(onepfive.r, s0y[3].r); // h2+1.0 always >0
    s0y[3].r = pv_real_MUL(s0y[3].r, s0y[3].r);
    s0y[3].r = pv_real_MUL(half.r, s0y[3].r);

    s0z[1].r = pv_real_SUB(h3.r, ones.r);
    s0z[1].r = pv_real_ADD(onepfive.r, s0z[1].r); // h3+1.0 always <0
    s0z[1].r = pv_real_MUL(s0z[1].r, s0z[1].r);
    s0z[1].r = pv_real_MUL(half.r, s0z[1].r);

    s0z[2].r = pv_real_MUL(h3.r, h3.r);
    s0z[2].r = pv_real_SUB(threefourths.r, s0z[2].r);

    s0z[3].r = pv_real_ADD(h3.r, ones.r);
    s0z[3].r = pv_real_SUB(onepfive.r, s0z[3].r); // h3+1.0 always >0
    s0z[3].r = pv_real_MUL(s0z[3].r, s0z[3].r);
    s0z[3].r = pv_real_MUL(half.r, s0z[3].r);



// CHECKPOINT: PIC_push_part_yz.F : line 143

    tmpx.r = pv_real_MUL(xi.r, dxi.r);
    tmpy.r = pv_real_MUL(yi.r, dyi.r);
    tmpy.r = pv_real_SUB(tmpy.r, half.r);
    tmpz.r = pv_real_MUL(zi.r, dzi.r);
    tmpz.r = pv_real_SUB(tmpz.r, half.r);

    for(int m = 0; m < VEC_SIZE; m++){
      l1.v[m] = round(tmpx.v[m]);
    }


    l2.r = pv_real_to_int_CVT(tmpy.r);
    l3.r = pv_real_to_int_CVT(tmpz.r);

    l2fl.r = pv_int_to_real_CVT(l2.r);
    l3fl.r = pv_int_to_real_CVT(l3.r);

    h2.r = pv_real_SUB(l2fl.r, tmpy.r);
    h3.r = pv_real_SUB(l3fl.r, tmpz.r);

    pvInt l2pls1, l2mns1, l3pls1, l3mns1;
    l2pls1.r = pv_int_ADD(l2.r, ione.r);
    l2mns1.r = pv_int_SUB(l2.r, ione.r);
    l3pls1.r = pv_int_ADD(l3.r, ione.r);
    l3mns1.r = pv_int_SUB(l3.r, ione.r);

    pvReal hmy, hmz, hOy, hOz, hly, hlz; // NB: this is hO as in octogon instead of h0, 
                                         // and hl as in library instead of h1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of h?[3], which 
                                         // seems to improve performancs


    hmy.r = pv_real_ADD(half.r, h2.r);
    hmy.r = pv_real_MUL(hmy.r, hmy.r);
    hmy.r = pv_real_MUL(half.r, hmy.r);

    hmz.r = pv_real_ADD(half.r, h3.r);
    hmz.r = pv_real_MUL(hmz.r, hmz.r);
    hmz.r = pv_real_MUL(half.r, hmz.r);
    
    hOy.r = pv_real_MUL(h2.r, h2.r);
    hOy.r = pv_real_SUB(threefourths.r, hOy.r);

    hOz.r = pv_real_MUL(h3.r, h3.r);
    hOz.r = pv_real_SUB(threefourths.r, hOz.r);

    hly.r = pv_real_SUB(half.r, h2.r);
    hly.r = pv_real_MUL(hly.r, hly.r);
    hly.r = pv_real_MUL(half.r, hly.r);

    hlz.r = pv_real_SUB(half.r, h3.r);
    hlz.r = pv_real_MUL(hlz.r, hlz.r);
    hlz.r = pv_real_MUL(half.r, hlz.r);

// CHECKPOINT: PIC_push_part_yz.F : line 59
// Field Interpolation

// BOTTLENECK HERE: For a given particle, none of the needed fields are stored
// next to each other in memory. Also, there's no way to tell before runtime
// exactly what fields are needed. As a consequence, I can't just load them into
// registers and shuffle. I must fetch piece by painful piece.

//Hopefully gcc is smart enough to unroll these for me. Otherwise, I need to figure out
// a way to generate unrolled loops
    
    pvReal exq, eyq, ezq, bxq, byq, bzq;

/*     pvReal gy[3] = {gmy, gOy, gly}; */
/*     pvReal gz[3] = {gmz, gOz, glz}; */
/*     pvReal hy[3] = {hmy, hOy, hly}; */
/*     pvReal hz[3] = {hmz, hOz, hlz}; */

    INTERP_FIELD_YZ(EX,l1,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,j1,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j1,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(BX,j1,l2,l3,h,h,bxq);
    INTERP_FIELD_YZ(BY,l1,j2,l3,h,g,byq);
    INTERP_FIELD_YZ(BZ,l1,l2,j3,g,h,bzq);

 
// CHECKPOINT: PIC_push_part_yz.F : line 223
// Half step momentum  with E-field

    pvReal dq;
    dq.r = pv_real_DIV(dqs.r, mni.r);
    dq.r = pv_real_MUL(qni.r, dq.r);
    
    exq.r = pv_real_MUL(dq.r, exq.r);
    eyq.r = pv_real_MUL(dq.r, eyq.r);
    ezq.r = pv_real_MUL(dq.r, ezq.r);
    
    pxi.r = pv_real_ADD(pxi.r, exq.r);
    pyi.r = pv_real_ADD(pyi.r, eyq.r);
    pzi.r = pv_real_ADD(pzi.r, ezq.r);


// CHECKPOINT: PIC_push_part_yz.F : line 228
// Rotate with B-field
    pvReal txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;

    tmpx.r = pv_real_MUL(pxi.r, pxi.r);
    tmpy.r = pv_real_MUL(pyi.r, pyi.r);
    tmpz.r = pv_real_MUL(pzi.r, pzi.r);
    root.r = pv_real_ADD(ones.r, tmpx.r);
    tmpy.r = pv_real_ADD(tmpy.r, tmpz.r);
    root.r = pv_real_ADD(root.r, tmpy.r);
    root.r = pv_real_SQRT(root.r);
    root.r = pv_real_DIV(dq.r, root.r);

    bxq.r = pv_real_MUL(bxq.r, root.r);
    byq.r = pv_real_MUL(byq.r, root.r);
    bzq.r = pv_real_MUL(bzq.r, root.r);

    txx.r = pv_real_MUL(bxq.r, bxq.r);
    tyy.r = pv_real_MUL(byq.r, byq.r);
    tzz.r = pv_real_MUL(bzq.r, bzq.r);
    t2xy.r = pv_real_MUL(bxq.r, byq.r);
    t2xz.r = pv_real_MUL(bxq.r, bzq.r);
    t2yz.r = pv_real_MUL(byq.r, bzq.r);
    t2xy.r = pv_real_ADD(t2xy.r, t2xy.r);
    t2xz.r = pv_real_ADD(t2xz.r, t2xz.r);
    t2yz.r = pv_real_ADD(t2yz.r, t2yz.r);

    pvReal tau;

    tau.r = pv_real_ADD(ones.r, txx.r);
    tmpx.r = pv_real_ADD(tyy.r, tzz.r);
    tau.r = pv_real_ADD(tau.r, tmpx.r);
    tau.r = pv_real_DIV(ones.r, tau.r); //recp is evil! evil evil evil evil!!!!
    
    bxq.r = pv_real_ADD(bxq.r, bxq.r);
    byq.r = pv_real_ADD(byq.r, byq.r);
    bzq.r = pv_real_ADD(bzq.r, bzq.r);
    
    //pxp
    tmpx.r = pv_real_ADD(ones.r, txx.r);
    tmpx.r = pv_real_SUB(tmpx.r, tyy.r);
    tmpx.r = pv_real_SUB(tmpx.r, tzz.r);
    tmpx.r = pv_real_MUL(tmpx.r, pxi.r);
    
    tmpy.r = pv_real_ADD(t2xy.r, bzq.r);
    tmpy.r = pv_real_MUL(tmpy.r, pyi.r);

    tmpz.r = pv_real_SUB(t2xz.r, byq.r);
    tmpz.r = pv_real_MUL(tmpz.r, pzi.r);

    pxp.r = pv_real_ADD(tmpx.r, tmpy.r);
    pxp.r = pv_real_ADD(pxp.r, tmpz.r);
    pxp.r = pv_real_MUL(pxp.r, tau.r);
    
    //pyp
    tmpx.r = pv_real_SUB(t2xy.r, bzq.r);
    tmpx.r = pv_real_MUL(tmpx.r, pxi.r);

    tmpy.r = pv_real_SUB(ones.r, txx.r);
    tmpy.r = pv_real_ADD(tmpy.r, tyy.r);
    tmpy.r = pv_real_SUB(tmpy.r, tzz.r);
    tmpy.r = pv_real_MUL(tmpy.r, pyi.r);
    
    tmpz.r = pv_real_ADD(t2yz.r, bxq.r);
    tmpz.r = pv_real_MUL(tmpz.r, pzi.r);

    pyp.r = pv_real_ADD(tmpx.r, tmpy.r);
    pyp.r = pv_real_ADD(pyp.r, tmpz.r);
    pyp.r = pv_real_MUL(pyp.r, tau.r);
   
    //pzp
    tmpx.r = pv_real_ADD(t2xz.r, byq.r);
    tmpx.r = pv_real_MUL(tmpx.r, pxi.r);

    tmpy.r = pv_real_SUB(t2yz.r, bxq.r);
    tmpy.r = pv_real_MUL(tmpy.r, pyi.r);

    tmpz.r = pv_real_SUB(ones.r, txx.r);
    tmpz.r = pv_real_SUB(tmpz.r, tyy.r);
    tmpz.r = pv_real_ADD(tmpz.r, tzz.r);
    tmpz.r = pv_real_MUL(tmpz.r, pzi.r);

    pzp.r = pv_real_ADD(tmpx.r, tmpy.r);
    pzp.r = pv_real_ADD(pzp.r, tmpz.r);
    pzp.r = pv_real_MUL(pzp.r, tau.r);

// CHECKPOINT: PIC_push_part_yz.F : line 244
// Half step momentum  with E-field

    pxi.r = pv_real_ADD(pxp.r, exq.r);
    pyi.r = pv_real_ADD(pyp.r, eyq.r);
    pzi.r = pv_real_ADD(pzp.r, ezq.r);

// CHECKPOINT: PIC_push_part_yz.F : line 248
// Half step particles with new momenta
    
    tmpx.r = pv_real_MUL(pxi.r, pxi.r);
    tmpy.r = pv_real_MUL(pyi.r, pyi.r);
    tmpz.r = pv_real_MUL(pzi.r, pzi.r);

    tmpx.r = pv_real_ADD(tmpx.r, tmpy.r);
    tmpx.r = pv_real_ADD(tmpx.r, tmpz.r);
    tmpx.r = pv_real_ADD(ones.r, tmpx.r);
    root.r = pv_real_SQRT(tmpx.r);
    root.r = pv_real_DIV(ones.r, root.r);

    vxi.r = pv_real_MUL(pxi.r, root.r);
    vyi.r = pv_real_MUL(pyi.r, root.r);
    vzi.r = pv_real_MUL(pzi.r, root.r);

    tmpy.r = pv_real_MUL(vyi.r, yl.r);
    tmpz.r = pv_real_MUL(vzi.r, zl.r);
    
    yi.r = pv_real_ADD(yi.r, tmpy.r);
    zi.r = pv_real_ADD(zi.r, tmpz.r);

    STORE_PART_XP(sse2,n);

    tmpx.r = pv_real_DIV(ones.r, root.r);
    tmpx.r = pv_real_SUB(tmpx.r, ones.r);
    tmpx.r = pv_real_DIV(tmpx.r, eta.r);
    tmpx.r = pv_real_MUL(fnqs.r, tmpx.r);
    tmpx.r = pv_real_MUL(mni.r, tmpx.r);

    for(int p = 0; p < VEC_SIZE; p++){
      psc.p2B += tmpx.v[p]; //What's this for?
    }

// CHECKPOINT: PIC_push_part_yz.F : line 266
// update number densities.
    tmpx.r = pv_real_MUL(xi.r, dxi.r);
    tmpy.r = pv_real_MUL(yi.r, dyi.r);
    tmpz.r = pv_real_MUL(zi.r, dzi.r);

    //FIXME: Pick one!
   for(int m = 0; m < VEC_SIZE; m++){
      l1.v[m] = round(tmpx.v[m]);
    }
 
    l2.r = pv_real_to_int_CVT(tmpy.r);
    l3.r = pv_real_to_int_CVT(tmpz.r);

    l2fl.r = pv_int_to_real_CVT(l2.r);
    l3fl.r = pv_int_to_real_CVT(l3.r);

    h2.r = pv_real_SUB(l2fl.r, tmpy.r);
    h3.r = pv_real_SUB(l3fl.r, tmpz.r);


    gmy.r = pv_real_SUB(h2.r,ones.r);
    //I'm pretty sure |h2,h3| < 0.5. If I'm wrong, I'll change this part
    gmy.r = pv_real_ADD(onepfive.r, gmy.r); // h2-1 always <0
    gmy.r = pv_real_MUL(gmy.r, gmy.r);
    gmy.r = pv_real_MUL(half.r, gmy.r);

    gOy.r = pv_real_MUL(h2.r, h2.r);
    gOy.r = pv_real_SUB(threefourths.r, gOy.r);

    gly.r = pv_real_ADD(ones.r, h2.r); 
    gly.r = pv_real_SUB(onepfive.r, gly.r);//h2+1 always >0
    gly.r = pv_real_MUL(gly.r, gly.r);
    gly.r = pv_real_MUL(half.r, gly.r);

    gmz.r = pv_real_SUB(h3.r,ones.r);
    gmz.r = pv_real_ADD(onepfive.r, gmz.r); //h3-1 always <0
    gmz.r = pv_real_MUL(gmz.r, gmz.r);
    gmz.r = pv_real_MUL(half.r, gmz.r);

    gOz.r = pv_real_MUL(h3.r, h3.r);
    gOz.r = pv_real_SUB(threefourths.r, gOz.r);

    glz.r = pv_real_ADD(ones.r, h3.r);
    glz.r = pv_real_SUB(onepfive.r, glz.r); //h3+1 always >0
    glz.r = pv_real_MUL(glz.r, glz.r);
    glz.r = pv_real_MUL(half.r, glz.r);

// CHECKPOINT: PIC_push_part_yz.F : line 283
    //FIXME: this is, for now, a straight serial translation. I think
    //    there is a more efficient way to do this, but I might be making
    //    some changes to the the local field structure and I'll worry 
    //    about it then.

    for(int m=0; m < VEC_SIZE ; m++){ 
      sse2_real fnq;
      // This may or may not work, and may or may not help
      sse2_real *densp; 
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
	DEN_FIELD(l2.v[m], l3.v[m]-1) += fnq*gOy.v[m]*gmz.v[m];
	DEN_FIELD(l2.v[m]+1,l3.v[m]-1) += fnq*gly.v[m]*gmz.v[m];
	DEN_FIELD(l2.v[m]-1, l3.v[m]) += fnq*gmy.v[m]*gOz.v[m];
	DEN_FIELD(l2.v[m], l3.v[m]) += fnq*gly.v[m]*gOz.v[m];
	DEN_FIELD(l2.v[m]+1,l3.v[m]) += fnq*gly.v[m]*gOz.v[m];
	DEN_FIELD(l2.v[m]-1, l3.v[m]+1) += fnq*gmy.v[m]*glz.v[m];
	DEN_FIELD(l2.v[m], l3.v[m]+1) += fnq*gOy.v[m]*glz.v[m];
	DEN_FIELD(l2.v[m]+1,l3.v[m]+1) += fnq*gly.v[m]*glz.v[m];      
    }
    
    // CHECKPOINT: PIC_push_part_yz.F : line 320
    // Charge density form factor at (n+1.5)*dt
    tmpy.r = pv_real_MUL(vyi.r, yl.r);
    tmpz.r = pv_real_MUL(vzi.r, zl.r);
    
    yi.r = pv_real_ADD(yi.r, tmpy.r);
    zi.r = pv_real_ADD(zi.r, tmpz.r);

    pvInt k2,k3;

    pvReal k2fl, k3fl;

    tmpy.r = pv_real_MUL(yi.r, dyi.r);
    tmpz.r = pv_real_MUL(zi.r, dzi.r);

    k2.r = pv_real_to_int_CVT(tmpy.r);
    k3.r = pv_real_to_int_CVT(tmpz.r);

    k2fl.r = pv_int_to_real_CVT(k2.r);
    k3fl.r = pv_int_to_real_CVT(k3.r);

    h2.r = pv_real_SUB(k2fl.r, tmpy.r);
    h3.r = pv_real_SUB(k3fl.r, tmpz.r);


    // God help me, there's some things I just can't fgure out how to parallelize
    // The g-- here are just temporary variables. I can't do the assignments in parallel, 
    // but I'll be damned if I can't do the FLOPS in parallel

    gmy.r = pv_real_SUB(h2.r,ones.r);
    gmy.r = pv_real_ADD(onepfive.r, gmy.r); // h2-1 always <0
    gmy.r = pv_real_MUL(gmy.r, gmy.r);
    gmy.r = pv_real_MUL(half.r, gmy.r);

    gOy.r = pv_real_MUL(h2.r, h2.r);
    gOy.r = pv_real_SUB(threefourths.r, gOy.r);

    gly.r = pv_real_ADD(ones.r, h2.r); 
    gly.r = pv_real_SUB(onepfive.r, gly.r);//h2+1 always >0
    gly.r = pv_real_MUL(gly.r, gly.r);
    gly.r = pv_real_MUL(half.r, gly.r);

    gmz.r = pv_real_SUB(h3.r,ones.r);
    gmz.r = pv_real_ADD(onepfive.r, gmz.r); //h3-1 always <0
    gmz.r = pv_real_MUL(gmz.r, gmz.r);
    gmz.r = pv_real_MUL(half.r, gmz.r);

    gOz.r = pv_real_MUL(h3.r, h3.r);
    gOz.r = pv_real_SUB(threefourths.r, gOz.r);
    
    glz.r = pv_real_ADD(ones.r, h3.r);
    glz.r = pv_real_SUB(onepfive.r, glz.r); //h3+1 always >0
    glz.r = pv_real_MUL(glz.r, glz.r);
    glz.r = pv_real_MUL(half.r, glz.r);

   
    
    // All the indices below have +2 compared to the FORTRAN for C array access
    pvInt shifty, shiftz;
    shifty.r = pv_int_SUB(k2.r,j2.r);
    shiftz.r = pv_int_SUB(k3.r,j3.r);
    for(int p=0; p < VEC_SIZE; p++){
      s1y[shifty.v[p] + 1].v[p] = gmy.v[p];
      s1y[shifty.v[p] + 2].v[p] = gOy.v[p];
      s1y[shifty.v[p] + 3].v[p] = gly.v[p];
      s1z[shiftz.v[p] + 1].v[p] = gmz.v[p];
      s1z[shiftz.v[p] + 2].v[p] = gOz.v[p];
      s1z[shiftz.v[p] + 3].v[p] = glz.v[p];
    }


// CHECKPOINT: PIC_push_part_yz.F : line 341
// Charge density form factor at (n+1.5)*dt



    for(int m=0; m<5; m++){
      s1y[m].r = pv_real_SUB(s1y[m].r, s0y[m].r);
      s1z[m].r = pv_real_SUB(s1z[m].r, s0z[m].r);
    }
    
    pvReal fnqx, fnqy, fnqz;
    
    fnqx.r = pv_real_MUL(wni.r, fnqs.r);
    fnqx.r = pv_real_MUL(qni.r, fnqx.r);
    fnqx.r = pv_real_MUL(vxi.r, fnqx.r);

    fnqy.r = pv_real_MUL(wni.r, fnqys.r);
    fnqy.r = pv_real_MUL(qni.r, fnqy.r);
    
    fnqz.r = pv_real_MUL(wni.r, fnqzs.r);
    fnqz.r = pv_real_MUL(qni.r, fnqz.r);
        
    // Again, this is extremely painful, but serial is the only way I can think of to do this
    // I have the inkling of thought on a parallel way, but it hasn't matured. If it works out, I'll implement it.
    for(int p=0; p < VEC_SIZE; p++){
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

      sse2_real jzh[5] = {0.0};
   
      for(int l3i=l3min+2; l3i<=(l3max+2); l3i++){
    	sse2_real jyh = 0.0;
	for(int l2i=l2min+2; l2i<=(l2max+2); l2i++){
    	  sse2_real wx = (s0y[l2i].v[p] + 0.5*s1y[l2i].v[p])*s0z[l3i].v[p]
    	    + (0.5*s0y[l2i].v[p] + 0.3333333333*s1y[l2i].v[p])*s1z[l3i].v[p];
    	  sse2_real wy = s1y[l2i].v[p]*(s0z[l3i].v[p] + 0.5*s1z[l3i].v[p]);
    	  sse2_real wz = s1z[l3i].v[p]*(s0y[l2i].v[p] + 0.5*s1y[l2i].v[p]);
	  
    	  jyh = jyh - fnqy.v[p]*wy;
    	  jzh[l2i] = jzh[l2i] - fnqz.v[p]*wz;

    	  CF3(JXI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += fnqx.v[p]*wx; //undo index adjustment of l(2,3)i
    	  CF3(JYI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += jyh; //undo index adjustment of l(2,3)i
    	  CF3(JZI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += jzh[l2i]; //undo index adjustment of l(2,3)i
    	}
      }
    }   
  }

  prof_stop(pr);
}
