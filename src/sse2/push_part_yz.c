#include "psc_sse2.h"
#include "sse2_cgen.h"
#include "profile/profile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>



pvReal ones, half, threefourths, onepfive, third;
pvInt ione;

static void 
init_vec_numbers() {		
  ones.r = pv_real_SET1(1.0);			
  half.r = pv_real_SET1(.5);			
  onepfive.r = pv_real_SET1(1.5);		
  threefourths.r = pv_real_SET1(.75);		
  third.r = pv_real_SET1(1./3.);		
  ione.r = pv_int_SET1(1);			
}

static inline void // Root used to update p2X, so need to have it
calc_vi(struct particle_vec *p, pvReal *vxi, pvReal *vyi, pvReal *vzi, pvReal *root){
  pvReal tmpx, tmpy, tmpz, roottmp;
  
  tmpx.r = pv_real_MUL(p->pxi.r, p->pxi.r);
  tmpy.r = pv_real_MUL(p->pyi.r, p->pyi.r);
  tmpz.r = pv_real_MUL(p->pzi.r, p->pzi.r);
  
  tmpx.r = pv_real_ADD(tmpx.r, tmpy.r);
  tmpz.r = pv_real_ADD(ones.r, tmpz.r);
  tmpx.r = pv_real_ADD(tmpx.r, tmpz.r);
  root->r = pv_real_SQRT(tmpx.r);
  roottmp.r = pv_real_DIV(ones.r, root->r);  //this avoids unneccessary DIVs, which 
                                             // have a very high latency

  vxi->r = pv_real_MUL(p->pxi.r, roottmp.r);
  vyi->r = pv_real_MUL(p->pyi.r, roottmp.r);
  vzi->r = pv_real_MUL(p->pzi.r, roottmp.r);
}

static inline void
push_xi_halfdt(struct particle_vec *p, pvReal *vyi, pvReal *vzi, pvReal *yl, pvReal *zl){
  pvReal tmpy, tmpz;
  tmpy.r = pv_real_MUL(vyi->r, yl->r);
  tmpz.r = pv_real_MUL(vzi->r, zl->r);
  p->yi.r = pv_real_ADD(p->yi.r, tmpy.r);
  p->zi.r = pv_real_ADD(p->zi.r, tmpz.r);
}

// Finds index using C99 round
static inline void
find_index_Cround(pvReal *xi, pvReal *dxi, pvInt *j, pvReal *h){
  pvReal tmp;
  tmp.r = pv_real_MUL(xi->r, dxi->r);
  for(int m = 0; m < VEC_SIZE; m++){
      j->v[m] = round(tmp.v[m]);
  }
  h->r = pv_int_to_real_CVT(j->r);
  h->r = pv_real_SUB(h->r, tmp.r);
}

// Finds index using intrinsic round
static inline void
find_index_Iround(pvReal *xi, pvReal *dxi, pvInt *j, pvReal *h){
  pvReal tmp;
  tmp.r = pv_real_MUL(xi->r, dxi->r);
  j->r = pv_real_to_int_CVT(tmp.r);
  h->r = pv_int_to_real_CVT(j->r);
  h->r = pv_real_SUB(h->r, tmp.r);
}

static inline void
find_index_minus_shift(pvReal *xi, pvReal *dxi, pvInt *j, pvReal *h, pvReal *shift){
  pvReal tmp;
  tmp.r = pv_real_MUL(xi->r, dxi->r);
  tmp.r = pv_real_SUB(tmp.r, shift->r);
  j->r = pv_real_to_int_CVT(tmp.r);
  h->r = pv_int_to_real_CVT(j->r);
  h->r = pv_real_SUB(h->r, tmp.r);
}

static inline void
ip_to_grid_m(pvReal *h, pvReal *xmx){
  xmx->r = pv_real_ADD(half.r, h->r);
  xmx->r = pv_real_MUL(xmx->r, xmx->r);
  xmx->r = pv_real_MUL(half.r, xmx->r);
}

static inline void
ip_to_grid_O(pvReal *h, pvReal *xOx){
  xOx->r = pv_real_MUL(h->r, h->r);
  xOx->r = pv_real_SUB(threefourths.r, xOx->r);
}

static inline void
ip_to_grid_l(pvReal *h, pvReal *xlx){
  xlx->r = pv_real_SUB(half.r, h->r);
  xlx->r = pv_real_MUL(xlx->r, xlx->r);
  xlx->r = pv_real_MUL(half.r, xlx->r);
}

static inline void
form_factor_m(pvReal *h, pvReal *xmx){
  xmx->r = pv_real_SUB(h->r,ones.r);
  xmx->r = pv_real_ADD(onepfive.r, xmx->r); // h-1 always <0
  xmx->r = pv_real_MUL(xmx->r, xmx->r);
  xmx->r = pv_real_MUL(half.r, xmx->r);
}

static inline void
form_factor_O(pvReal *h, pvReal *xOx){
  xOx->r = pv_real_MUL(h->r, h->r);
  xOx->r = pv_real_SUB(threefourths.r, xOx->r);
}

static inline void
form_factor_l(pvReal *h, pvReal *xlx){
  xlx->r = pv_real_ADD(ones.r, h->r); 
  xlx->r = pv_real_SUB(onepfive.r, xlx->r);//h+1 always >0
  xlx->r = pv_real_MUL(xlx->r, xlx->r);
  xlx->r = pv_real_MUL(half.r, xlx->r);
}

static void
push_pi_dt(struct particle_vec *p,
	   pvReal *exq, pvReal *eyq, pvReal *ezq,
	   pvReal *bxq, pvReal *byq, pvReal *bzq,
	   pvReal *dqs){

  pvReal tmpx, tmpy, tmpz, root;
  // Half step momentum  with E-field

  pvReal dq, dqex, dqey, dqez;
  dq.r = pv_real_DIV(dqs->r, p->mni.r);
  dq.r = pv_real_MUL(p->qni.r, dq.r);
  
  dqex.r = pv_real_MUL(dq.r, exq->r);
  dqey.r = pv_real_MUL(dq.r, eyq->r);
  dqez.r = pv_real_MUL(dq.r, ezq->r);
  
  p->pxi.r = pv_real_ADD(p->pxi.r, dqex.r);
  p->pyi.r = pv_real_ADD(p->pyi.r, dqey.r);
  p->pzi.r = pv_real_ADD(p->pzi.r, dqez.r);


// CHECKPOINT: PIC_push_part_yz.F : line 228
// Rotate with B-field
  pvReal taux, tauy, tauz, txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;
  
  tmpx.r = pv_real_MUL(p->pxi.r, p->pxi.r);
  tmpy.r = pv_real_MUL(p->pyi.r, p->pyi.r);
  tmpz.r = pv_real_MUL(p->pzi.r, p->pzi.r);
  root.r = pv_real_ADD(ones.r, tmpx.r);
  tmpy.r = pv_real_ADD(tmpy.r, tmpz.r);
  root.r = pv_real_ADD(root.r, tmpy.r);
  root.r = pv_real_SQRT(root.r);
  root.r = pv_real_DIV(dq.r, root.r);

  taux.r = pv_real_MUL(bxq->r, root.r);
  tauy.r = pv_real_MUL(byq->r, root.r);
  tauz.r = pv_real_MUL(bzq->r, root.r);

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
  tmpx.r = pv_real_MUL(tmpx.r, p->pxi.r);
    
  tmpy.r = pv_real_ADD(t2xy.r, tauz.r);
  tmpy.r = pv_real_MUL(tmpy.r, p->pyi.r);

  tmpz.r = pv_real_SUB(t2xz.r, tauy.r);
  tmpz.r = pv_real_MUL(tmpz.r, p->pzi.r);

  pxp.r = pv_real_ADD(tmpx.r, tmpy.r);
  pxp.r = pv_real_ADD(pxp.r, tmpz.r);
  pxp.r = pv_real_MUL(pxp.r, tau.r);
    
  //pyp
  tmpx.r = pv_real_SUB(t2xy.r, tauz.r);
  tmpx.r = pv_real_MUL(tmpx.r, p->pxi.r);

  tmpy.r = pv_real_SUB(ones.r, txx.r);
  tmpy.r = pv_real_ADD(tmpy.r, tyy.r);
  tmpy.r = pv_real_SUB(tmpy.r, tzz.r);
  tmpy.r = pv_real_MUL(tmpy.r, p->pyi.r);
  
  tmpz.r = pv_real_ADD(t2yz.r, taux.r);
  tmpz.r = pv_real_MUL(tmpz.r, p->pzi.r);

  pyp.r = pv_real_ADD(tmpx.r, tmpy.r);
  pyp.r = pv_real_ADD(pyp.r, tmpz.r);
  pyp.r = pv_real_MUL(pyp.r, tau.r);
   
  //pzp
  tmpx.r = pv_real_ADD(t2xz.r, tauy.r);
  tmpx.r = pv_real_MUL(tmpx.r, p->pxi.r);

  tmpy.r = pv_real_SUB(t2yz.r, taux.r);
  tmpy.r = pv_real_MUL(tmpy.r, p->pyi.r);

  tmpz.r = pv_real_SUB(ones.r, txx.r);
  tmpz.r = pv_real_SUB(tmpz.r, tyy.r);
  tmpz.r = pv_real_ADD(tmpz.r, tzz.r);
  tmpz.r = pv_real_MUL(tmpz.r, p->pzi.r);

  pzp.r = pv_real_ADD(tmpx.r, tmpy.r);
  pzp.r = pv_real_ADD(pzp.r, tmpz.r);
  pzp.r = pv_real_MUL(pzp.r, tau.r);

// CHECKPOINT: PIC_push_part_yz.F : line 244
// Half step momentum  with E-field

  p->pxi.r = pv_real_ADD(pxp.r, dqex.r);
  p->pyi.r = pv_real_ADD(pyp.r, dqey.r);
  p->pzi.r = pv_real_ADD(pzp.r, dqez.r);
}


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

  // Set vector real forms of 1, 0.5, etc
  init_vec_numbers();

  pvReal dt, yl, zl;
  //FIXME: These are all stored as doubles in fortran!!
  sse2_real dtfl = psc.dt; 
  dt.r = pv_real_SET1(dtfl);
  yl.r = pv_real_MUL(half.r, dt.r);
  zl.r = pv_real_MUL(half.r, dt.r);


  assert(psc.n_part % VEC_SIZE == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {

//---------------------------------------------
// Bringing in particle specific parameters
    
    struct particle_vec p;
 
    LOAD_PART(sse2,n); 

    // Locals for computation      
    pvReal vxi, vyi, vzi, root;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);
    
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
  

  // Set vector forms of 1, 0.5, etc
  init_vec_numbers();

  // Values that won't change from iteration to iteration
  pvReal dt, yl, zl, eta, dqs,fnqs, dxi, dyi, dzi; 
  //FIXME: These are all stored as doubles in fortran!!
  sse2_real dtfl = psc.dt; 
  sse2_real etafl =  psc.prm.eta;
  sse2_real fnqsfl = sqr(psc.prm.alpha) * psc.prm.cori / etafl;
  sse2_real dxifl = 1.0 / psc.dx[0];
  sse2_real dyifl = 1.0 / psc.dx[1];
  sse2_real dzifl = 1.0 / psc.dx[2];
  sse2_real dqsfl = 0.5*etafl*dtfl;
  dt.r = pv_real_SET1(dtfl);
  eta.r = pv_real_SET1(etafl);
  fnqs.r = pv_real_SET1(fnqsfl);
  dqs.r = pv_real_SET1(dqsfl);
  dxi.r = pv_real_SET1(dxifl);
  dyi.r = pv_real_SET1(dyifl);
  dzi.r = pv_real_SET1(dzifl);
  yl.r = pv_real_MUL(half.r, dt.r);
  zl.r = pv_real_MUL(half.r, dt.r);


  //  assert(psc.n_part % VEC_SIZE == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {

//---------------------------------------------
// Bringing in particle specific parameters

    struct particle_vec p;
    //    pvReal pxi, pyi, pzi, xi, yi, zi, qni, mni, wni; 
 
    LOAD_PART(sse2,n); 

    // Locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, tmpz, root, h1, h2, h3;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta

    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    tmpx.r = pv_real_SUB(root.r, ones.r);
    tmpx.r = pv_real_DIV(tmpx.r, eta.r);
    tmpy.r = pv_real_MUL(p.mni.r, fnqs.r);
    tmpx.r = pv_real_MUL(tmpy.r, tmpx.r);
    
    for(int m = 0; m < VEC_SIZE; m++){
      psc.p2A += tmpx.v[m]; //What's this for?
    }

// CHECKPOINT: PIC_push_part_yz.F : line 110
    

// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    // Apparently this can be done in vectors. A victory!
    pvInt j1, j2, j3, l1, l2, l3;

    find_index_Cround(&(p.xi), &dxi, &j1, &h1);
    find_index_Iround(&(p.yi), &dyi, &j2, &h2);
    find_index_Iround(&(p.zi), &dzi, &j3, &h3);

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


    ip_to_grid_m(&h2, &gmy);
    ip_to_grid_m(&h3, &gmz);
     
    ip_to_grid_O(&h2, &gOy);
    ip_to_grid_O(&h3, &gOz);

    ip_to_grid_l(&h2, &gly);
    ip_to_grid_l(&h3, &glz);

// CHECKPOINT: PIC_push_part_yz.F : line 143

   
    find_index_Cround(&(p.xi), &dxi, &l1, &h1);
    find_index_minus_shift(&(p.yi), &dyi, &l2, &h2, &half);
    find_index_minus_shift(&(p.zi), &dzi, &l3, &h3, &half);

    pvInt l2pls1, l2mns1, l3pls1, l3mns1;
    l2pls1.r = pv_int_ADD(l2.r, ione.r);
    l2mns1.r = pv_int_SUB(l2.r, ione.r);
    l3pls1.r = pv_int_ADD(l3.r, ione.r);
    l3mns1.r = pv_int_SUB(l3.r, ione.r);

    pvReal hmy, hmz, hOy, hOz, hly, hlz; // NB: this is hO as in octogon instead of h0, 
                                         // and hl as in library instead of h1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of h?[3], which 
                                         // seems to improve performance


    ip_to_grid_m(&h2, &hmy);
    ip_to_grid_m(&h3, &hmz);
     
    ip_to_grid_O(&h2, &hOy);
    ip_to_grid_O(&h3, &hOz);

    ip_to_grid_l(&h2, &hly);
    ip_to_grid_l(&h3, &hlz);


// CHECKPOINT: PIC_push_part_yz.F : line 59
// Field Interpolation

// BOTTLENECK HERE: For a given particle, none of the needed fields are stored
// next to each other in memory. Also, there's no way to tell before runtime
// exactly what fields are needed. As a consequence, I can't just load them into
// registers and shuffle. I must fetch piece by painful piece.
    
    pvReal exq, eyq, ezq, bxq, byq, bzq;


    INTERP_FIELD_YZ(EX,l1,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,j1,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j1,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(BX,j1,l2,l3,h,h,bxq);
    INTERP_FIELD_YZ(BY,l1,j2,l3,h,g,byq);
    INTERP_FIELD_YZ(BZ,l1,l2,j3,g,h,bzq);

 
// CHECKPOINT: PIC_push_part_yz.F : line 223

    push_pi_dt(&p, &exq, &eyq, &ezq, &bxq, &byq, &bzq, &dqs);

// CHECKPOINT: PIC_push_part_yz.F : line 248
// Half step particles with new momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    STORE_PART_XP(sse2,n);

    tmpx.r = pv_real_SUB(root.r, ones.r);
    tmpx.r = pv_real_DIV(tmpx.r, eta.r);
    tmpy.r = pv_real_MUL(p.mni.r, fnqs.r);
    tmpx.r = pv_real_MUL(tmpy.r, tmpx.r);

    for(int m = 0; m < VEC_SIZE; m++){
      psc.p2B += tmpx.v[m]; //What's this for?
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

  //Set vector forms of numbers 1, 0.5, etc
  init_vec_numbers();

  // Values that won't change from iteration to iteration

  pvReal dt, yl, zl, eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs; 
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

  assert(psc.n_part % VEC_SIZE == 0); // Haven't implemented any padding yet
  
  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {

//---------------------------------------------
// Bringing in particle specific parameters
    
    struct particle_vec p;
    //    pvReal pxi, pyi, pzi, xi, yi, zi, qni, mni, wni;  
    LOAD_PART(sse2,n); 

    // Locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, tmpz, root, h1, h2, h3;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta

    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    tmpx.r = pv_real_SUB(root.r, ones.r);
    tmpx.r = pv_real_DIV(tmpx.r, eta.r);
    tmpy.r = pv_real_MUL(p.mni.r, fnqs.r);
    tmpx.r = pv_real_MUL(tmpy.r, tmpx.r);
    
    for(int m = 0; m < VEC_SIZE; m++){
      psc.p2A += tmpx.v[m]; //What's this for?
    }


    pvReal s0y[5], s0z[5], s1y[5], s1z[5];
    // I'm a little scared to use memset here, though it would probably work...
    for(int m = 0; m < 5; m++){
      s0y[m].r = pv_real_SET1(0.0);
      s1y[m].r = pv_real_SET1(0.0);
      s0z[m].r = pv_real_SET1(0.0);
      s1z[m].r = pv_real_SET1(0.0);
    }

// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    pvReal  H1, H2, H3;
    pvInt j1, j2, j3, l1, l2, l3;

    find_index_Cround(&(p.xi), &dxi, &j1, &H1);
    find_index_Iround(&(p.yi), &dyi, &j2, &H2);
    find_index_Iround(&(p.zi), &dzi, &j3, &H3);

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


    ip_to_grid_m(&H2, &gmy);
    ip_to_grid_m(&H3, &gmz);
     
    ip_to_grid_O(&H2, &gOy);
    ip_to_grid_O(&H3, &gOz);

    ip_to_grid_l(&H2, &gly);
    ip_to_grid_l(&H3, &glz);


// CHECKPOINT: PIC_push_part_yz.F : line 143

    find_index_Cround(&(p.xi), &dxi, &l1, &h1);
    find_index_minus_shift(&(p.yi), &dyi, &l2, &h2, &half);
    find_index_minus_shift(&(p.zi), &dzi, &l3, &h3, &half);
    
    pvInt l2pls1, l2mns1, l3pls1, l3mns1;
    l2pls1.r = pv_int_ADD(l2.r, ione.r);
    l2mns1.r = pv_int_SUB(l2.r, ione.r);
    l3pls1.r = pv_int_ADD(l3.r, ione.r);
    l3mns1.r = pv_int_SUB(l3.r, ione.r);
    
    pvReal hmy, hmz, hOy, hOz, hly, hlz; // NB: this is hO as in octogon instead of h0, 
                                         // and hl as in library instead of h1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of h?[3], which 
                                         // seems to improve performance


    ip_to_grid_m(&h2, &hmy);
    ip_to_grid_m(&h3, &hmz);
    
    ip_to_grid_O(&h2, &hOy);
    ip_to_grid_O(&h3, &hOz);
    
    ip_to_grid_l(&h2, &hly);
    ip_to_grid_l(&h3, &hlz);

// CHECKPOINT: PIC_push_part_yz.F : line 59
// Field Interpolation

// BOTTLENECK HERE: For a given particle, none of the needed fields are stored
// next to each other in memory. Also, there's no way to tell before runtime
// exactly what fields are needed. As a consequence, I can't just load them into
// registers and shuffle. I must fetch piece by painful piece.

    
    pvReal exq, eyq, ezq, bxq, byq, bzq;

    INTERP_FIELD_YZ(EX,l1,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,j1,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j1,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(BX,j1,l2,l3,h,h,bxq);
    INTERP_FIELD_YZ(BY,l1,j2,l3,h,g,byq);
    INTERP_FIELD_YZ(BZ,l1,l2,j3,g,h,bzq);


    // indexing here departs from FORTRAN a little bit
    // Simply moving this down here seems to improve performance
    // Possibly giving the processor something to do while it's waiting 
    // for fields to read in
    
    form_factor_m(&H2, &s0y[1]);
    form_factor_O(&H2, &s0y[2]);
    form_factor_l(&H2, &s0y[3]);

    form_factor_m(&H3, &s0z[1]);
    form_factor_O(&H3, &s0z[2]);
    form_factor_l(&H3, &s0z[3]);

// CHECKPOINT: PIC_push_part_yz.F : line 223

    push_pi_dt(&p, &exq, &eyq, &ezq, &bxq, &byq, &bzq, &dqs);

// CHECKPOINT: PIC_push_part_yz.F : line 248
// Half step particles with new momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    STORE_PART_XP(sse2,n);

// CHECKPOINT: PIC_push_part_yz.F : line 266
// update number densities.

#if 1
    find_index_Cround(&(p.xi), &dxi, &l1, &h1);
    find_index_Iround(&(p.yi), &dyi, &l2, &h2);
    find_index_Iround(&(p.zi), &dzi, &l3, &h3);

    form_factor_m(&h2, &gmy);
    form_factor_O(&h2, &gOy);
    form_factor_l(&h2, &gly);

    form_factor_m(&h3, &gmz);
    form_factor_O(&h3, &gOz);
    form_factor_l(&h3, &glz);

// CHECKPOINT: PIC_push_part_yz.F : line 283
    //FIXME: this is, for now, a straight serial translation. I think
    //    there is a more efficient way to do this, but I might be making
    //    some changes to the the local field structure and I'll worry 
    //    about it then.

    for(int m=0; m < VEC_SIZE ; m++){ 
      sse2_real fnq;
      // This may or may not work, and may or may not help
      sse2_real *densp; 
      if (p.qni.v[m] < 0.0){
	densp = &(sse2->fields[NE*psc.fld_size + l1.v[m] - psc.ilo[0] + psc.ibn[0]]);
	fnq = p.qni.v[m] * p.wni.v[m] * fnqs.v[m];
      }
      else if (p.qni.v[m] > 0.0){
	densp = &(sse2->fields[NI*psc.fld_size + l1.v[m] - psc.ilo[0] + psc.ibn[0]]);
	fnq = p.qni.v[m] * p.wni.v[m] * fnqs.v[m];
      }
      else if (p.qni.v[m] == 0.0){
	densp = &(sse2->fields[NN*psc.fld_size + l1.v[m] - psc.ilo[0] + psc.ibn[0]]);
	fnq = p.wni.v[m] * fnqs.v[m];
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

#endif 
    // CHECKPOINT: PIC_push_part_yz.F : line 320
    // Charge density form factor at (n+1.5)*dt

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    pvInt k2,k3;

    find_index_Iround(&(p.yi), &dyi, &k2, &h2);
    find_index_Iround(&(p.zi), &dzi, &k3, &h3);

    // God help me, there's some things I just can't fgure out how to parallelize
    // The g-- here are just temporary variables. I can't do the assignments in parallel, 
    // but I'll be damned if I can't do the FLOPS in parallel

    form_factor_m(&h2, &gmy);
    form_factor_O(&h2, &gOy);
    form_factor_l(&h2, &gly);

    form_factor_m(&h3, &gmz);
    form_factor_O(&h3, &gOz);
    form_factor_l(&h3, &glz);
       
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

    for(int m=0; m<5; m++){
      s1y[m].r = pv_real_SUB(s1y[m].r, s0y[m].r);
      s1z[m].r = pv_real_SUB(s1z[m].r, s0z[m].r);
    }

// CHECKPOINT: PIC_push_part_yz.F : line 341
// Charge density form factor at (n+1.5)*dt

    int l2min = 1;
    int l2max = 4;
    int l3min = 1;
    int l3max = 4;
    
    for(int p = 0; p < VEC_SIZE; p++){
      if(k2.v[p]==j2.v[p]){
	s1y[0].v[p] = 0.0;
	s1y[4].v[p] = 0.0;
	s0y[0].v[p] = 0.0;
	s0y[4].v[p] = 0.0;
      }
      else if (k2.v[p]==j2.v[p]-1){
	s1y[4].v[p] = 0.0;
	s0y[4].v[p] = 0.0;
	l2min = 0;
      }
      else if (k2.v[p]==j2.v[p]+1){
	s1y[0].v[p] = 0.0;
	s0y[0].v[p] = 0.0;
	l2max = 5;
      }
      if (k3.v[p]==j3.v[p]){
	s1z[0].v[p] = 0.0;
	s1z[4].v[p] = 0.0;
	s0z[0].v[p] = 0.0;
	s0z[4].v[p] = 0.0;
      }
      else if(k3.v[p]==j3.v[p]-1){
	s1z[4].v[p] = 0.0;
	s0z[4].v[p] = 0.0;
	l3min = 0;
      }
      else if(k3.v[p]==j3.v[p]+1){
	s1z[0].v[p] = 0.0;
	s0z[0].v[p] = 0.0;
	l3max = 5;
      }
    }

    pvReal fnqx, fnqy, fnqz;
    
    fnqx.r = pv_real_MUL(p.wni.r, fnqs.r);
    fnqx.r = pv_real_MUL(p.qni.r, fnqx.r);
    fnqx.r = pv_real_MUL(vxi.r, fnqx.r);

    fnqy.r = pv_real_MUL(p.wni.r, fnqys.r);
    fnqy.r = pv_real_MUL(p.qni.r, fnqy.r);
    
    fnqz.r = pv_real_MUL(p.wni.r, fnqzs.r);
    fnqz.r = pv_real_MUL(p.qni.r, fnqz.r);
        
    pvReal jzh[5];
    memset(&jzh, 0, 5 * sizeof(pvReal));
    
    for(int l3i=l3min; l3i<l3max; l3i++){
      pvReal jyh;
      jyh.r = pv_real_SET1(0.0);
      for(int l2i=l2min; l2i<l2max; l2i++){
	pvReal wx, wy, wz;
	wx.r = pv_real_MUL(half.r,s1y[l2i].r);
	wx.r = pv_real_ADD(s0y[l2i].r, wx.r);
	wx.r = pv_real_MUL(s0z[l3i].r, wx.r);
	tmpx.r = pv_real_MUL(half.r, s0y[l2i].r);
	tmpy.r = pv_real_MUL(third.r, s1y[l2i].r);
	tmpx.r = pv_real_ADD(tmpx.r, tmpy.r);
	tmpx.r = pv_real_MUL(tmpx.r, s1z[l3i].r);
	wx.r = pv_real_ADD(wx.r, tmpx.r);
	
	wy.r = pv_real_MUL(half.r, s1z[l3i].r);
	wy.r = pv_real_ADD(s0z[l3i].r, wy.r);
	wy.r = pv_real_MUL(s1y[l2i].r, wy.r);
	
	wz.r = pv_real_MUL(half.r, s1y[l2i].r);
	wz.r = pv_real_ADD(s0y[l2i].r, wz.r);
	wz.r = pv_real_MUL(s1z[l3i].r, wz.r);
	
	wx.r = pv_real_MUL(fnqx.r, wx.r);
	wy.r = pv_real_MUL(fnqy.r, wy.r);
	jyh.r = pv_real_SUB(jyh.r, wy.r);
	wz.r = pv_real_MUL(fnqz.r, wz.r);
	jzh[l2i].r = pv_real_SUB(jzh[l2i].r, wz.r);
	
	for(int p = 0; p < VEC_SIZE; p++){
	  CF3(JXI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += wx.v[p]; //undo index adjustment of l(2,3)i
	  CF3(JYI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += jyh.v[p]; //undo index adjustment of l(2,3)i
    	  CF3(JZI,j1.v[p],j2.v[p] + l2i - 2, j3.v[p] + l3i -2) += jzh[l2i].v[p]; //undo index adjustment of l(2,3)i
    	}
      }
    }   
  }
  
  prof_stop(pr);
}
