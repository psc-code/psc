#include "psc_sse2.h"
#include "simd_wrap.h"
#include "profile/profile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>


pvReal ones, half, threefourths, onepfive, third;
pvInt ione;

static void 
init_vec_numbers() {		
  ones.r = pv_set1_real(1.0);			
  half.r = pv_set1_real(.5);			
  onepfive.r = pv_set1_real(1.5);		
  threefourths.r = pv_set1_real(.75);		
  third.r = pv_set1_real(1./3.);		
  ione.r = pv_set1_int(1);			
}

#define JSX(indx2, indx3) s_jxi[((indx3) - psc.ilg[2])*psc.img[1] + ((indx2) - psc.ilg[1])]
#define JSY(indx2, indx3) s_jyi[((indx3) - psc.ilg[2])*psc.img[1] + ((indx2) - psc.ilg[1])]
#define JSZ(indx2, indx3) s_jzi[((indx3) - psc.ilg[2])*psc.img[1] + ((indx2) - psc.ilg[1])]

static inline void // Root used to update p2X, so need to have it
calc_vi(struct particle_vec *p, pvReal * restrict vxi, pvReal * restrict vyi,
	pvReal * restrict vzi, pvReal * restrict root){
  pvReal tmpx, tmpy, tmpz, roottmp;
  
  tmpx.r = pv_mul_real(p->pxi.r, p->pxi.r);
  tmpy.r = pv_mul_real(p->pyi.r, p->pyi.r);
  tmpz.r = pv_mul_real(p->pzi.r, p->pzi.r);
  
  tmpx.r = pv_add_real(tmpx.r, tmpy.r);
  tmpz.r = pv_add_real(ones.r, tmpz.r);
  tmpx.r = pv_add_real(tmpx.r, tmpz.r);
  root->r = pv_sqrt_real(tmpx.r);
  roottmp.r = pv_div_real(ones.r, root->r);
                                             
  vxi->r = pv_mul_real(p->pxi.r, roottmp.r);
  vyi->r = pv_mul_real(p->pyi.r, roottmp.r);
  vzi->r = pv_mul_real(p->pzi.r, roottmp.r);
}

static inline void
push_xi_halfdt(struct particle_vec *p, pvReal *vyi, pvReal *vzi, pvReal *yl, pvReal *zl){
  pvReal tmpy, tmpz;
  tmpy.r = pv_mul_real(vyi->r, yl->r);
  tmpz.r = pv_mul_real(vzi->r, zl->r);
  p->yi.r = pv_add_real(p->yi.r, tmpy.r);
  p->zi.r = pv_add_real(p->zi.r, tmpz.r);
}

// Finds index using C99 round
static inline void
find_index_Cround(pvReal *xi, pvReal *dxi, pvInt * restrict j, pvReal * restrict h){
  pvReal tmp;
  tmp.r = pv_mul_real(xi->r, dxi->r);
  for(int m = 0; m < VEC_SIZE; m++){
      j->v[m] = round(tmp.v[m]);
  }
  h->r = pv_cvt_int_to_real(j->r);
  h->r = pv_sub_real(h->r, tmp.r);
}

// Finds index using intrinsic round
static inline void
find_index_Iround(pvReal *xi, pvReal *dxi, pvInt * restrict j, pvReal * restrict h){
  pvReal tmp;
  tmp.r = pv_mul_real(xi->r, dxi->r);
  j->r = pv_cvt_real_to_int(tmp.r);
  h->r = pv_cvt_int_to_real(j->r);
  h->r = pv_sub_real(h->r, tmp.r);
}

static inline void
find_index_minus_shift(pvReal *xi, pvReal *dxi, pvInt * restrict j, pvReal * restrict h, pvReal *shift){
  pvReal tmp;
  tmp.r = pv_mul_real(xi->r, dxi->r);
  tmp.r = pv_sub_real(tmp.r, shift->r);
  j->r = pv_cvt_real_to_int(tmp.r);
  h->r = pv_cvt_int_to_real(j->r);
  h->r = pv_sub_real(h->r, tmp.r);
}

static inline void
ip_to_grid_m(pvReal *h, pvReal * restrict xmx){
  xmx->r = pv_add_real(half.r, h->r);
  xmx->r = pv_mul_real(xmx->r, xmx->r);
  xmx->r = pv_mul_real(half.r, xmx->r);
}

static inline void
ip_to_grid_O(pvReal *h, pvReal * restrict xOx){
  xOx->r = pv_mul_real(h->r, h->r);
  xOx->r = pv_sub_real(threefourths.r, xOx->r);
}

static inline void
ip_to_grid_l(pvReal *h, pvReal * restrict xlx){
  xlx->r = pv_sub_real(half.r, h->r);
  xlx->r = pv_mul_real(xlx->r, xlx->r);
  xlx->r = pv_mul_real(half.r, xlx->r);
}

static inline void
form_factor_m(pvReal *h, pvReal * restrict xmx){
  xmx->r = pv_sub_real(h->r,ones.r);
  xmx->r = pv_add_real(onepfive.r, xmx->r); // h-1 always <0
  xmx->r = pv_mul_real(xmx->r, xmx->r);
  xmx->r = pv_mul_real(half.r, xmx->r);
}

static inline void
form_factor_O(pvReal *h, pvReal * restrict xOx){
  xOx->r = pv_mul_real(h->r, h->r);
  xOx->r = pv_sub_real(threefourths.r, xOx->r);
}

static inline void
form_factor_l(pvReal *h, pvReal * restrict xlx){
  xlx->r = pv_add_real(ones.r, h->r); 
  xlx->r = pv_sub_real(onepfive.r, xlx->r);//h+1 always >0
  xlx->r = pv_mul_real(xlx->r, xlx->r);
  xlx->r = pv_mul_real(half.r, xlx->r);
}

static void
push_pi_dt(struct particle_vec * p,
	   pvReal * exq, pvReal * eyq,  pvReal * ezq,
	   pvReal * bxq, pvReal * byq,  pvReal * bzq,
	   pvReal * dqs){

  pvReal tmpx, tmpy, tmpz, root;
  // Half step momentum  with E-field

  pvReal dq, dqex, dqey, dqez;
  dq.r = pv_div_real(dqs->r, p->mni.r);
  dq.r = pv_mul_real(p->qni.r, dq.r);
  
  dqex.r = pv_mul_real(dq.r, exq->r);
  dqey.r = pv_mul_real(dq.r, eyq->r);
  dqez.r = pv_mul_real(dq.r, ezq->r);
  
  p->pxi.r = pv_add_real(p->pxi.r, dqex.r);
  p->pyi.r = pv_add_real(p->pyi.r, dqey.r);
  p->pzi.r = pv_add_real(p->pzi.r, dqez.r);


// CHECKPOINT: PIC_push_part_yz.F : line 228
// Rotate with B-field
  pvReal taux, tauy, tauz, txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;
  
  tmpx.r = pv_mul_real(p->pxi.r, p->pxi.r);
  tmpy.r = pv_mul_real(p->pyi.r, p->pyi.r);
  tmpz.r = pv_mul_real(p->pzi.r, p->pzi.r);
  root.r = pv_add_real(ones.r, tmpx.r);
  tmpy.r = pv_add_real(tmpy.r, tmpz.r);
  root.r = pv_add_real(root.r, tmpy.r);
  root.r = pv_sqrt_real(root.r);
  root.r = pv_div_real(dq.r, root.r);

  taux.r = pv_mul_real(bxq->r, root.r);
  tauy.r = pv_mul_real(byq->r, root.r);
  tauz.r = pv_mul_real(bzq->r, root.r);

  txx.r = pv_mul_real(taux.r, taux.r);
  tyy.r = pv_mul_real(tauy.r, tauy.r);
  tzz.r = pv_mul_real(tauz.r, tauz.r);
  t2xy.r = pv_mul_real(taux.r, tauy.r);
  t2xz.r = pv_mul_real(taux.r, tauz.r);
  t2yz.r = pv_mul_real(tauy.r, tauz.r);
  t2xy.r = pv_add_real(t2xy.r, t2xy.r);
  t2xz.r = pv_add_real(t2xz.r, t2xz.r);
  t2yz.r = pv_add_real(t2yz.r, t2yz.r);

  pvReal tau;

  tau.r = pv_add_real(ones.r, txx.r);
  tmpx.r = pv_add_real(tyy.r, tzz.r);
  tau.r = pv_add_real(tau.r, tmpx.r);
  tau.r = pv_div_real(ones.r, tau.r); //recp is evil! evil evil evil evil!!!!
    
  //Never use tau_ without a two in front
  taux.r = pv_add_real(taux.r, taux.r);
  tauy.r = pv_add_real(tauy.r, tauy.r);
  tauz.r = pv_add_real(tauz.r, tauz.r);
    
  //pxp
  tmpx.r = pv_add_real(ones.r, txx.r);
  tmpx.r = pv_sub_real(tmpx.r, tyy.r);
  tmpx.r = pv_sub_real(tmpx.r, tzz.r);
  tmpx.r = pv_mul_real(tmpx.r, p->pxi.r);
    
  tmpy.r = pv_add_real(t2xy.r, tauz.r);
  tmpy.r = pv_mul_real(tmpy.r, p->pyi.r);

  tmpz.r = pv_sub_real(t2xz.r, tauy.r);
  tmpz.r = pv_mul_real(tmpz.r, p->pzi.r);

  pxp.r = pv_add_real(tmpx.r, tmpy.r);
  pxp.r = pv_add_real(pxp.r, tmpz.r);
  pxp.r = pv_mul_real(pxp.r, tau.r);
    
  //pyp
  tmpx.r = pv_sub_real(t2xy.r, tauz.r);
  tmpx.r = pv_mul_real(tmpx.r, p->pxi.r);

  tmpy.r = pv_sub_real(ones.r, txx.r);
  tmpy.r = pv_add_real(tmpy.r, tyy.r);
  tmpy.r = pv_sub_real(tmpy.r, tzz.r);
  tmpy.r = pv_mul_real(tmpy.r, p->pyi.r);
  
  tmpz.r = pv_add_real(t2yz.r, taux.r);
  tmpz.r = pv_mul_real(tmpz.r, p->pzi.r);

  pyp.r = pv_add_real(tmpx.r, tmpy.r);
  pyp.r = pv_add_real(pyp.r, tmpz.r);
  pyp.r = pv_mul_real(pyp.r, tau.r);
   
  //pzp
  tmpx.r = pv_add_real(t2xz.r, tauy.r);
  tmpx.r = pv_mul_real(tmpx.r, p->pxi.r);

  tmpy.r = pv_sub_real(t2yz.r, taux.r);
  tmpy.r = pv_mul_real(tmpy.r, p->pyi.r);

  tmpz.r = pv_sub_real(ones.r, txx.r);
  tmpz.r = pv_sub_real(tmpz.r, tyy.r);
  tmpz.r = pv_add_real(tmpz.r, tzz.r);
  tmpz.r = pv_mul_real(tmpz.r, p->pzi.r);

  pzp.r = pv_add_real(tmpx.r, tmpy.r);
  pzp.r = pv_add_real(pzp.r, tmpz.r);
  pzp.r = pv_mul_real(pzp.r, tau.r);

// CHECKPOINT: PIC_push_part_yz.F : line 244
// Half step momentum  with E-field

  p->pxi.r = pv_add_real(pxp.r, dqex.r);
  p->pyi.r = pv_add_real(pyp.r, dqey.r);
  p->pzi.r = pv_add_real(pzp.r, dqez.r);
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
  dt.r = pv_set1_real(dtfl);
  yl.r = pv_mul_real(half.r, dt.r);
  zl.r = pv_mul_real(half.r, dt.r);

  int elements = VEC_SIZE;

  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {
    int part_left = psc.n_part - n;
    if((part_left < VEC_SIZE) && (part_left != 0)){
      elements = part_left;
    }

//---------------------------------------------
// Bringing in particle specific parameters
    
    struct particle_vec p;

    for(int m=0; m < elements; m++){
      p.xi.v[m] = sse2->part[n+m].xi;
      p.yi.v[m] = sse2->part[n+m].yi;
      p.zi.v[m] = sse2->part[n+m].zi;
      p.pxi.v[m] = sse2->part[n+m].pxi;
      p.pyi.v[m] = sse2->part[n+m].pyi;
      p.pzi.v[m] = sse2->part[n+m].pzi;
      p.qni.v[m] = sse2->part[n+m].qni;
      p.mni.v[m] = sse2->part[n+m].mni;
      p.wni.v[m] = sse2->part[n+m].wni;
    }


    // Locals for computation      
    pvReal vxi, vyi, vzi, root;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);
    
    for(int m=0; m<elements; m++){
      (sse2->part[n+m]).xi = p.xi.v[m];		
      (sse2->part[n+m]).yi = p.yi.v[m];		
      (sse2->part[n+m]).zi = p.zi.v[m];		
      (sse2->part[n+m]).pxi = p.pxi.v[m];	
      (sse2->part[n+m]).pyi = p.pyi.v[m];     
      (sse2->part[n+m]).pzi = p.pzi.v[m];	
    }

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

  //Vector versions of floating point parameters
    pvReal dt, yl, zl, eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs;
    sse2_real dxifl = 1.0 / psc.dx[0];					
    sse2_real dyifl = 1.0 / psc.dx[1];					
    sse2_real dzifl = 1.0 / psc.dx[2];					
    sse2_real fnqsfl = sqr(psc.prm.alpha) * psc.prm.cori / psc.prm.eta;
    sse2_real fnqxsfl = psc.dx[0] * fnqsfl * psc.dt;			
    sse2_real fnqysfl = psc.dx[1] * fnqsfl * psc.dt;			
    sse2_real fnqzsfl = psc.dx[2] * fnqsfl * psc.dt;			
    dt.r = pv_set1_real(psc.dt);					     
    eta.r = pv_set1_real(psc.prm.eta);					
    fnqs.r = pv_set1_real(fnqsfl);
    dqs.r = pv_set1_real(0.5*psc.prm.eta*psc.dt);			    
    dxi.r = pv_set1_real(dxifl);					
    dyi.r = pv_set1_real(dyifl);					
    dzi.r = pv_set1_real(dzifl);					
    yl.r = pv_mul_real(half.r, dt.r);					
    zl.r = pv_mul_real(half.r, dt.r);					
    fnqxs.r = pv_set1_real(fnqxsfl);				       
    fnqys.r = pv_set1_real(fnqysfl);					
    fnqzs.r = pv_set1_real(fnqzsfl);			
    
    //Vector versions of integer parameters
    pvInt ilg[3], img[3], fld_size;					
    ilg[0].r = pv_set1_int(psc.ilg[0]);					
    ilg[1].r = pv_set1_int(psc.ilg[1]);					
    ilg[2].r = pv_set1_int(psc.ilg[2]);					
    img[0].r = pv_set1_int(psc.img[0]);					
    img[1].r = pv_set1_int(psc.img[1]);					
    img[2].r = pv_set1_int(psc.img[2]);					
    fld_size.r = pv_set1_int(psc.fld_size);				
    					
    // Assign pointers to fields assuming x is uniform
    sse2_real * restrict EXpoint = &sse2->fields[EX*psc.fld_size + 1 - psc.ilg[0]]; //FIXME: this assumes xi always rounds up to 1
    sse2_real * restrict EYpoint = &sse2->fields[EY*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict EZpoint = &sse2->fields[EZ*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict BXpoint = &sse2->fields[BX*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict BYpoint = &sse2->fields[BY*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict BZpoint = &sse2->fields[BZ*psc.fld_size + 1 - psc.ilg[0]];


  int elements = VEC_SIZE;
  
  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {
    
    // This little ditty basically allows padding of the vector when it runs out of particles
    int part_left = psc.n_part - n;
    if((part_left < VEC_SIZE) && (part_left != 0)){
      elements = part_left;
    }

    pvInt itemp1, itemp2;

//---------------------------------------------
// Bringing in particle specific parameters

    struct particle_vec p;
 
    for(int m=0; m < elements; m++){
      p.xi.v[m] = sse2->part[n+m].xi;
      p.yi.v[m] = sse2->part[n+m].yi;
      p.zi.v[m] = sse2->part[n+m].zi;
      p.pxi.v[m] = sse2->part[n+m].pxi;
      p.pyi.v[m] = sse2->part[n+m].pyi;
      p.pzi.v[m] = sse2->part[n+m].pzi;
      p.qni.v[m] = sse2->part[n+m].qni;
      p.mni.v[m] = sse2->part[n+m].mni;
      p.wni.v[m] = sse2->part[n+m].wni;
    }

    // Locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, tmpz, root, h1, h2, h3;

// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta

    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    tmpx.r = pv_sub_real(root.r, ones.r);
    tmpx.r = pv_div_real(tmpx.r, eta.r);
    tmpy.r = pv_mul_real(p.mni.r, fnqs.r);
    tmpx.r = pv_mul_real(tmpy.r, tmpx.r);
    
    for(int m = 0; m < elements; m++){
      psc.p2A += tmpx.v[m]; //What's this for?
    }

// CHECKPOINT: PIC_push_part_yz.F : line 110
    

// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    // Apparently this can be done in vectors. A victory!
    pvInt j2, j3, l2, l3;

    find_index_Iround(&(p.yi), &dyi, &j2, &h2);
    find_index_Iround(&(p.zi), &dzi, &j3, &h3);

    pvInt j2pls1, j2mns1, j3pls1, j3mns1;
    j2pls1.r = pv_add_int(j2.r, ione.r);
    j2mns1.r = pv_sub_int(j2.r, ione.r);
    j3pls1.r = pv_add_int(j3.r, ione.r);
    j3mns1.r = pv_sub_int(j3.r, ione.r);


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

    find_index_minus_shift(&(p.yi), &dyi, &l2, &h2, &half);
    find_index_minus_shift(&(p.zi), &dzi, &l3, &h3, &half);

    pvInt l2pls1, l2mns1, l3pls1, l3mns1;
    l2pls1.r = pv_add_int(l2.r, ione.r);
    l2mns1.r = pv_sub_int(l2.r, ione.r);
    l3pls1.r = pv_add_int(l3.r, ione.r);
    l3mns1.r = pv_sub_int(l3.r, ione.r);

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

    INTERP_FIELD_YZ(EX,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(BX,l2,l3,h,h,bxq);
    INTERP_FIELD_YZ(BY,j2,l3,h,g,byq);
    INTERP_FIELD_YZ(BZ,l2,j3,g,h,bzq);
 
// CHECKPOINT: PIC_push_part_yz.F : line 223

    push_pi_dt(&p, &exq, &eyq, &ezq, &bxq, &byq, &bzq, &dqs);

// CHECKPOINT: PIC_push_part_yz.F : line 248
// Half step particles with new momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);


    for(int m=0; m<elements; m++){
      (sse2->part[n+m]).xi = p.xi.v[m];		
      (sse2->part[n+m]).yi = p.yi.v[m];		
      (sse2->part[n+m]).zi = p.zi.v[m];		
      (sse2->part[n+m]).pxi = p.pxi.v[m];	
      (sse2->part[n+m]).pyi = p.pyi.v[m];     
      (sse2->part[n+m]).pzi = p.pzi.v[m];	
    }


    tmpx.r = pv_sub_real(root.r, ones.r);
    tmpx.r = pv_div_real(tmpx.r, eta.r);
    tmpy.r = pv_mul_real(p.mni.r, fnqs.r);
    tmpx.r = pv_mul_real(tmpy.r, tmpx.r);

    for(int m = 0; m < elements; m++){
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

  // I don't like allocating these, but I think they could potentially be
  // too big to declare statically. If I'm wrong, having them be static will give
  // a nice performance boost
  sse2_real *s_jxi, *s_jyi, *s_jzi;
  s_jxi = calloc((psc.img[1] * psc.img[2]), sizeof(sse2_real));
  s_jyi = calloc((psc.img[1] * psc.img[2]), sizeof(sse2_real));
  s_jzi = calloc((psc.img[1] * psc.img[2]), sizeof(sse2_real));
  
  //Set vector forms of numbers 1, 0.5, etc
  init_vec_numbers();


  //Vector versions of floating point parameters
    pvReal dt, yl, zl, eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs;
    sse2_real dxifl = 1.0 / psc.dx[0];					
    sse2_real dyifl = 1.0 / psc.dx[1];					
    sse2_real dzifl = 1.0 / psc.dx[2];					
    sse2_real fnqsfl = sqr(psc.prm.alpha) * psc.prm.cori / psc.prm.eta;
    sse2_real fnqxsfl = psc.dx[0] * fnqsfl * psc.dt;			
    sse2_real fnqysfl = psc.dx[1] * fnqsfl * psc.dt;			
    sse2_real fnqzsfl = psc.dx[2] * fnqsfl * psc.dt;			
    dt.r = pv_set1_real(psc.dt);					     
    eta.r = pv_set1_real(psc.prm.eta);					
    fnqs.r = pv_set1_real(fnqsfl);
    dqs.r = pv_set1_real(0.5*psc.prm.eta*psc.dt);			    
    dxi.r = pv_set1_real(dxifl);					
    dyi.r = pv_set1_real(dyifl);					
    dzi.r = pv_set1_real(dzifl);					
    yl.r = pv_mul_real(half.r, dt.r);					
    zl.r = pv_mul_real(half.r, dt.r);					
    fnqxs.r = pv_set1_real(fnqxsfl);				       
    fnqys.r = pv_set1_real(fnqysfl);					
    fnqzs.r = pv_set1_real(fnqzsfl);			
    
    //Vector versions of integer parameters
    pvInt ilg[3], img[3], fld_size;					
    ilg[0].r = pv_set1_int(psc.ilg[0]);					
    ilg[1].r = pv_set1_int(psc.ilg[1]);					
    ilg[2].r = pv_set1_int(psc.ilg[2]);					
    img[0].r = pv_set1_int(psc.img[0]);					
    img[1].r = pv_set1_int(psc.img[1]);					
    img[2].r = pv_set1_int(psc.img[2]);					
    fld_size.r = pv_set1_int(psc.fld_size);				
    					
    //Set up some pointers for fields assuming x is uniform
    sse2_real * restrict EXpoint = &sse2->fields[EX*psc.fld_size + 1 - psc.ilg[0]]; //FIXME: this assumes xi always rounds up to 1
    sse2_real * restrict EYpoint = &sse2->fields[EY*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict EZpoint = &sse2->fields[EZ*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict BXpoint = &sse2->fields[BX*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict BYpoint = &sse2->fields[BY*psc.fld_size + 1 - psc.ilg[0]];
    sse2_real * restrict BZpoint = &sse2->fields[BZ*psc.fld_size + 1 - psc.ilg[0]];


    //    sse2_real * restrict fpoint = sse2->fields; // This keyword is awesome. 
                                                // One restrict and my laptop gets 
                                                // a 5% performance boost
    
    int elements = VEC_SIZE;

  for(int n = 0; n < psc.n_part; n += VEC_SIZE) {

    // Check if we're padding
    int part_left = psc.n_part - n;
    if((part_left < VEC_SIZE) && (part_left != 0)){
      elements = part_left;
    }


//---------------------------------------------
// Bringing in particle specific parameters
    
    struct particle_vec p;
 
    for(int m=0; m < elements; m++){
      p.xi.v[m] = sse2->part[n+m].xi;
      p.yi.v[m] = sse2->part[n+m].yi;
      p.zi.v[m] = sse2->part[n+m].zi;
      p.pxi.v[m] = sse2->part[n+m].pxi;
      p.pyi.v[m] = sse2->part[n+m].pyi;
      p.pzi.v[m] = sse2->part[n+m].pzi;
      p.qni.v[m] = sse2->part[n+m].qni;
      p.mni.v[m] = sse2->part[n+m].mni;
      p.wni.v[m] = sse2->part[n+m].wni;
    }


    // Locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, tmpz, root, h2, h3;
    pvInt itemp1, itemp2;
// CHECKPOINT: PIC_push_part_yz.F : line 104
// Start the computation
// Half step positions with current momenta

    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    tmpx.r = pv_sub_real(root.r, ones.r);
    tmpx.r = pv_div_real(tmpx.r, eta.r);
    tmpy.r = pv_mul_real(p.mni.r, fnqs.r);
    tmpx.r = pv_mul_real(tmpy.r, tmpx.r);
    
    for(int m = 0; m < elements; m++){
      psc.p2A += tmpx.v[m]; //What's this for?
    }


    pvReal s0y[5], s0z[5], s1y[5], s1z[5];
    // I'm a little scared to use memset here, though it would probably work...
    for(int m = 0; m < 5; m++){
      s0y[m].r = pv_set1_real(0.0);
      s1y[m].r = pv_set1_real(0.0);
      s0z[m].r = pv_set1_real(0.0);
      s1z[m].r = pv_set1_real(0.0);
    }

// CHECKPOINT: PIC_push_part_yz.F : line 119
// Prepare for field interpolation

    pvReal  H2, H3;
    pvInt j2, j3, l2, l3;

    find_index_Iround(&(p.yi), &dyi, &j2, &H2);
    find_index_Iround(&(p.zi), &dzi, &j3, &H3);

    pvInt j2pls1, j2mns1, j3pls1, j3mns1;
    j2pls1.r = pv_add_int(j2.r, ione.r);
    j2mns1.r = pv_sub_int(j2.r, ione.r);
    j3pls1.r = pv_add_int(j3.r, ione.r);
    j3mns1.r = pv_sub_int(j3.r, ione.r);


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

    find_index_minus_shift(&(p.yi), &dyi, &l2, &h2, &half);
    find_index_minus_shift(&(p.zi), &dzi, &l3, &h3, &half);
    
    pvInt l2pls1, l2mns1, l3pls1, l3mns1;
    l2pls1.r = pv_add_int(l2.r, ione.r);
    l2mns1.r = pv_sub_int(l2.r, ione.r);
    l3pls1.r = pv_add_int(l3.r, ione.r);
    l3mns1.r = pv_sub_int(l3.r, ione.r);
    
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

    INTERP_FIELD_YZ(EX,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(BX,l2,l3,h,h,bxq);
    INTERP_FIELD_YZ(BY,j2,l3,h,g,byq);
    INTERP_FIELD_YZ(BZ,l2,j3,g,h,bzq);


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

    for(int m=0; m < elements; m++){
      (sse2->part[n+m]).xi = p.xi.v[m];		
      (sse2->part[n+m]).yi = p.yi.v[m];		
      (sse2->part[n+m]).zi = p.zi.v[m];		
      (sse2->part[n+m]).pxi = p.pxi.v[m];	
      (sse2->part[n+m]).pyi = p.pyi.v[m];     
      (sse2->part[n+m]).pzi = p.pzi.v[m];	
    }


// CHECKPOINT: PIC_push_part_yz.F : line 266
// update number densities.

#if 1
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

    for(int m=0; m < elements ; m++){ 
      sse2_real fnq;
      // This may or may not work, and may or may not help
      sse2_real *densp; 
      if (p.qni.v[m] < 0.0){
	densp = &(sse2->fields[NE*psc.fld_size + 1 - psc.ilo[0] + psc.ibn[0]]);
	fnq = p.qni.v[m] * p.wni.v[m] * fnqs.v[m];
      }
      else if (p.qni.v[m] > 0.0){
	densp = &(sse2->fields[NI*psc.fld_size + 1 - psc.ilo[0] + psc.ibn[0]]);
	fnq = p.qni.v[m] * p.wni.v[m] * fnqs.v[m];
      }
      else if (p.qni.v[m] == 0.0){
	densp = &(sse2->fields[NN*psc.fld_size + 1 - psc.ilo[0] + psc.ibn[0]]);
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
    shifty.r = pv_sub_int(k2.r,j2.r);
    shiftz.r = pv_sub_int(k3.r,j3.r);
    for(int p=0; p < elements; p++){
      s1y[shifty.v[p] + 1].v[p] = gmy.v[p];
      s1y[shifty.v[p] + 2].v[p] = gOy.v[p];
      s1y[shifty.v[p] + 3].v[p] = gly.v[p];
      s1z[shiftz.v[p] + 1].v[p] = gmz.v[p];
      s1z[shiftz.v[p] + 2].v[p] = gOz.v[p];
      s1z[shiftz.v[p] + 3].v[p] = glz.v[p];
    }

    for(int m=0; m<5; m++){
      s1y[m].r = pv_sub_real(s1y[m].r, s0y[m].r);
      s1z[m].r = pv_sub_real(s1z[m].r, s0z[m].r);
    }

// CHECKPOINT: PIC_push_part_yz.F : line 341
// Charge density form factor at (n+1.5)*dt

    int l2min = 1;
    int l2max = 4;
    int l3min = 1;
    int l3max = 4;
    
    for(int p = 0; p < elements; p++){
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
    
    fnqx.r = pv_mul_real(p.wni.r, fnqs.r);
    fnqx.r = pv_mul_real(p.qni.r, fnqx.r);
    fnqx.r = pv_mul_real(vxi.r, fnqx.r);

    fnqy.r = pv_mul_real(p.wni.r, fnqys.r);
    fnqy.r = pv_mul_real(p.qni.r, fnqy.r);
    
    fnqz.r = pv_mul_real(p.wni.r, fnqzs.r);
    fnqz.r = pv_mul_real(p.qni.r, fnqz.r);
        
    pvReal jzh[5];
    memset(&jzh, 0, 5 * sizeof(pvReal));
    
    for(int l3i=l3min; l3i<l3max; l3i++){
      pvReal jyh;
      jyh.r = pv_set1_real(0.0);
      for(int l2i=l2min; l2i<l2max; l2i++){
	pvReal wx, wy, wz;
	wx.r = pv_mul_real(half.r,s1y[l2i].r);
	wx.r = pv_add_real(s0y[l2i].r, wx.r);
	wx.r = pv_mul_real(s0z[l3i].r, wx.r);
	tmpx.r = pv_mul_real(half.r, s0y[l2i].r);
	tmpy.r = pv_mul_real(third.r, s1y[l2i].r);
	tmpx.r = pv_add_real(tmpx.r, tmpy.r);
	tmpx.r = pv_mul_real(tmpx.r, s1z[l3i].r);
	wx.r = pv_add_real(wx.r, tmpx.r);
	
	wy.r = pv_mul_real(half.r, s1z[l3i].r);
	wy.r = pv_add_real(s0z[l3i].r, wy.r);
	wy.r = pv_mul_real(s1y[l2i].r, wy.r);
	
	wz.r = pv_mul_real(half.r, s1y[l2i].r);
	wz.r = pv_add_real(s0y[l2i].r, wz.r);
	wz.r = pv_mul_real(s1z[l3i].r, wz.r);
	
	wx.r = pv_mul_real(fnqx.r, wx.r);
	wy.r = pv_mul_real(fnqy.r, wy.r);
	jyh.r = pv_sub_real(jyh.r, wy.r);
	wz.r = pv_mul_real(fnqz.r, wz.r);
	jzh[l2i].r = pv_sub_real(jzh[l2i].r, wz.r);
	
	for(int m = 0; m < elements; m++){
	  JSX(j2.v[m] + l2i - 2, j3.v[m] + l3i - 2) += wx.v[m]; //undo index adjustment of l(2,3)i
	  JSY(j2.v[m] + l2i - 2, j3.v[m] + l3i -2) += jyh.v[m]; //undo index adjustment of l(2,3)i
    	  JSZ(j2.v[m] + l2i - 2, j3.v[m] + l3i -2) += jzh[l2i].v[m]; //undo index adjustment of l(2
    	}
      }
    }   
  }
  
  for(int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++){
    for(int iy = psc.ilg[1]; iy < psc.ihg[1]; iy++){ 
      CF3(JXI,psc.ilo[0]+1,iy,iz) = JSX(iy, iz);
      CF3(JYI,psc.ilo[0]+1,iy,iz) = JSY(iy, iz);
      CF3(JZI,psc.ilo[0]+1,iy,iz) = JSZ(iy, iz);
    }
  }
  free(s_jxi);
  free(s_jyi);
  free(s_jzi);

  prof_stop(pr);
}
