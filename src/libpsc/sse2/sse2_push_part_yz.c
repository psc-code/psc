#include "psc_sse2.h"
#include "simd_wrap.h"
#include "simd_push_common.h"
#include <mrc_profile.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

// contrary to my initial beliefs, I think the 
// push_xi function is the only one which changes
// based on the axes of the system. The rest 
// have been offloaded to simd_push_common.c

// FIXME: indices are now 0-based for the local domain --
// pretty much everything is 0-based, other than actual particle
// positions, which need subtracting patch->xb[d] still.


//---------------------------------------------
// push_xi_halfdt - advances the posistions .5dt
// using the given velocities

static inline void
push_xi_halfdt(struct particle_vec *p, pvReal *vyi, pvReal *vzi, pvReal *yl, pvReal *zl){
  pvReal tmpy, tmpz;
  tmpy.r = pv_mul_real(vyi->r, yl->r);
  tmpz.r = pv_mul_real(vzi->r, zl->r);
  p->yi.r = pv_add_real(p->yi.r, tmpy.r);
  p->zi.r = pv_add_real(p->zi.r, tmpz.r);
}


//---------------------------------------------
/// Reads in the particles
/// and half steps their positions using their 
/// current momentum

static void
do_push_part_yz_a(particles_sse2_t *pp)
{
  //-----------------------------------------------------
  // Initialization stuff 
  
  // Set vector real forms of 1, 0.5, etc
  init_vec_numbers();

  pvReal dt, yl, zl;
  sse2_real dtfl = psc.dt; 
  dt.r = pv_set1_real(dtfl);
  yl.r = pv_mul_real(half.r, dt.r);
  zl.r = pv_mul_real(half.r, dt.r);


  for(int n = 0; n < pp->n_part; n += VEC_SIZE) {

    //---------------------------------------------
    // Bringing in particle specific parameters
    
    struct particle_vec p;

    LOAD_PARTS(p);

    // Locals for computation      
    pvReal vxi, vyi, vzi, root;

    //---------------------------------------------
    // Half step positions with current momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);
    
    //---------------------------------------------
    // Store Particles


    STORE_PARTS(p);
  
  }
}

//---------------------------------------------
/// Advances the particles 
/// one full time step, including momenta, but 
/// does not update current and charge densities

static void
do_push_part_yz_b(particles_sse2_t *pp, fields_sse2_t *pf)
{
  //-----------------------------------------------------
  // Initialization stuff
  
  psc.p2A = 0.;
  psc.p2B = 0.;
  
  //---------------------------------------------
  // Set vector forms of 1, 0.5, etc
  init_vec_numbers();

  //---------------------------------------------
  // Vector versions of floating point parameters
    pvReal dt, yl, zl, eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs;
    sse2_real dxifl = 1.0 / psc.dx[0];					
    sse2_real dyifl = 1.0 / psc.dx[1];					
    sse2_real dzifl = 1.0 / psc.dx[2];					
    sse2_real fnqsfl = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
    sse2_real fnqxsfl = psc.dx[0] * fnqsfl / psc.dt;			
    sse2_real fnqysfl = psc.dx[1] * fnqsfl / psc.dt;			
    sse2_real fnqzsfl = psc.dx[2] * fnqsfl / psc.dt;			
    dt.r = pv_set1_real(psc.dt);					     
    eta.r = pv_set1_real(psc.coeff.eta);					
    fnqs.r = pv_set1_real(fnqsfl);
    dqs.r = pv_set1_real(0.5*psc.coeff.eta*psc.dt);			    
    dxi.r = pv_set1_real(dxifl);					
    dyi.r = pv_set1_real(dyifl);					
    dzi.r = pv_set1_real(dzifl);					
    yl.r = pv_mul_real(half.r, dt.r);					
    zl.r = pv_mul_real(half.r, dt.r);					
    fnqxs.r = pv_set1_real(fnqxsfl);				       
    fnqys.r = pv_set1_real(fnqysfl);					
    fnqzs.r = pv_set1_real(fnqzsfl);			
    
    //---------------------------------------------
    //Vector versions of integer parameters
    struct psc_patch *patch = &psc.patch[0];
    int sz = 1;
    for (int d = 0; d < 3; d++) {
      sz *= patch->ldims[d] + 2 * psc.ibn[d];
    }

    pvInt ilg[3], img[3], fld_size;					
    ilg[0].r = pv_set1_int(-psc.ibn[0]);					
    ilg[1].r = pv_set1_int(-psc.ibn[1]);					
    ilg[2].r = pv_set1_int(-psc.ibn[2]);					
    img[0].r = pv_set1_int(patch->ldims[0] + 2 * psc.ibn[0]);
    img[1].r = pv_set1_int(patch->ldims[1] + 2 * psc.ibn[1]);
    img[2].r = pv_set1_int(patch->ldims[2] + 2 * psc.ibn[2]);
    fld_size.r = pv_set1_int(sz);				
    		
    //---------------------------------------------			
    // Assign pointers to fields assuming x is uniform
    // FIXME : this assumes xi always rounds down to 0!
    //    sse2_real * restrict EXpoint = &pf->flds[EX*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict EYpoint = &pf->flds[EY*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict EZpoint = &pf->flds[EZ*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict HXpoint = &pf->flds[HX*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict HYpoint = &pf->flds[HY*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict HZpoint = &pf->flds[HZ*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
  
  for(int n = 0; n < pp->n_part; n += VEC_SIZE) {

    //---------------------------------------------
    // Load Particles

    struct particle_vec p;

    LOAD_PARTS(p);

    //---------------------------------------------
    // Declare locals for computation      
    pvReal vxi, vyi, vzi, tmpx, tmpy, root, h2, h3;


    //---------------------------------------------
    // Half step positions with current momenta

    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    tmpx.r = pv_sub_real(root.r, ones.r);
    tmpx.r = pv_div_real(tmpx.r, eta.r);
    tmpy.r = pv_mul_real(p.mni.r, fnqs.r);
    tmpx.r = pv_mul_real(tmpy.r, tmpx.r);
    
    for(int m = 0; m < VEC_SIZE; m++){
      psc.p2A += tmpx.v[m]; //What's this for?
    }

    //---------------------------------------------
    // Prepare for field interpolation

    // FIXME : assuming true 2D and not touching x direction

    pvInt j2, j3, l2, l3;

    find_index_Iround(&(p.yi), &dyi, &j2, &h2);
    find_index_Iround(&(p.zi), &dzi, &j3, &h3);


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


    find_index_minus_shift(&(p.yi), &dyi, &l2, &h2, &half);
    find_index_minus_shift(&(p.zi), &dzi, &l3, &h3, &half);


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

    //---------------------------------------------
    // Field Interpolation
#define print_vec( vec ) printf("%.15e %.15e\n", vec.v[0], vec.v[1]);		




    pvReal exq, eyq, ezq, hxq, hyq, hzq;

    INTERP_FIELD_YZ(EX,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(HX,l2,l3,h,h,hxq);
    INTERP_FIELD_YZ(HY,j2,l3,h,g,hyq);
    INTERP_FIELD_YZ(HZ,l2,j3,g,h,hzq);

#undef print_vec 
    //---------------------------------------------
    // Advance momenta one full dt

    push_pi_dt(&p, &exq, &eyq, &ezq, &hxq, &hyq, &hzq, &dqs);

    //---------------------------------------------
    // Half step particles with new momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);


    //---------------------------------------------
    // Store particles 

    STORE_PARTS(p);
  

    tmpx.r = pv_sub_real(root.r, ones.r);
    tmpx.r = pv_div_real(tmpx.r, eta.r);
    tmpy.r = pv_mul_real(p.mni.r, fnqs.r);
    tmpx.r = pv_mul_real(tmpy.r, tmpx.r);

    for(int m = 0; m < VEC_SIZE; m++){
      psc.p2B += tmpx.v[m]; //What's this for?
    }

  }
}

//---------------------------------------------
/// SSE2 implementation of the 
/// yz particle pushers
static void
do_push_part_yz(particles_sse2_t *pp, fields_sse2_t *pf)
{
  //-----------------------------------------------------
  // Initialization stuff (not sure what all of this is for)
  
  psc.p2A = 0.;
  psc.p2B = 0.;
  
  struct psc_patch *patch = &psc.patch[0];
  int sz = 1;
  for (int d = 0; d < 3; d++) {
    sz *= patch->ldims[d] + 2 * psc.ibn[d];
  }

  for (int m = JXI; m <= JZI; m++) {
    memset(&pf->flds[m*sz], 0, sz * sizeof(sse2_real));
  }  

  //---------------------------------------------
  // An implementation of Will's 'squished' currents
  // that excludes the x direction all together
  sse2_real * restrict s_jxi, * restrict s_jyi, * restrict s_jzi;
  int jsz = ((patch->ldims[1] + 2 * psc.ibn[1]) * 
	     (patch->ldims[2] + 2 * psc.ibn[2]));
  s_jxi = calloc(jsz, sizeof(sse2_real));
  s_jyi = calloc(jsz, sizeof(sse2_real));
  s_jzi = calloc(jsz, sizeof(sse2_real));

  // -------------------------------
  // Macros for accessing the 'squished' currents allocated above
#define JSX(indx2, indx3) s_jxi[((indx3) + psc.ibn[2])*(psc.patch[0].ldims[1] + 2*psc.ibn[1]) + ((indx2) + psc.ibn[1])]
#define JSY(indx2, indx3) s_jyi[((indx3) + psc.ibn[2])*(psc.patch[0].ldims[1] + 2*psc.ibn[1]) + ((indx2) + psc.ibn[1])]
#define JSZ(indx2, indx3) s_jzi[((indx3) + psc.ibn[2])*(psc.patch[0].ldims[1] + 2*psc.ibn[1]) + ((indx2) + psc.ibn[1])]
  
  //-----------------------------------------------------
  //Set vector forms of numbers 1, 0.5, etc
  init_vec_numbers();

  //-----------------------------------------------------
  //Vector versions of floating point parameters
    pvReal dt, yl, zl, eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs;
    sse2_real dxifl = 1.0 / psc.dx[0];					
    sse2_real dyifl = 1.0 / psc.dx[1];					
    sse2_real dzifl = 1.0 / psc.dx[2];					
    sse2_real fnqsfl = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
    sse2_real fnqxsfl = psc.dx[0] * fnqsfl / psc.dt;			
    sse2_real fnqysfl = psc.dx[1] * fnqsfl / psc.dt;			
    sse2_real fnqzsfl = psc.dx[2] * fnqsfl / psc.dt;			
    dt.r = pv_set1_real(psc.dt);					     
    eta.r = pv_set1_real(psc.coeff.eta);					
    fnqs.r = pv_set1_real(fnqsfl);
    dqs.r = pv_set1_real(0.5*psc.coeff.eta*psc.dt);			    
    dxi.r = pv_set1_real(dxifl);					
    dyi.r = pv_set1_real(dyifl);					
    dzi.r = pv_set1_real(dzifl);					
    yl.r = pv_mul_real(half.r, dt.r);					
    zl.r = pv_mul_real(half.r, dt.r);					
    fnqxs.r = pv_set1_real(fnqxsfl);				       
    fnqys.r = pv_set1_real(fnqysfl);					
    fnqzs.r = pv_set1_real(fnqzsfl);			
    
    //-----------------------------------------------------
    //Vector versions of integer parameters
    pvInt ilg[3], img[3], fld_size;					
    ilg[0].r = pv_set1_int(-psc.ibn[0]);					
    ilg[1].r = pv_set1_int(-psc.ibn[1]);					
    ilg[2].r = pv_set1_int(-psc.ibn[2]);					
    img[0].r = pv_set1_int(patch->ldims[0] + 2 * psc.ibn[0]);
    img[1].r = pv_set1_int(patch->ldims[1] + 2 * psc.ibn[1]);
    img[2].r = pv_set1_int(patch->ldims[2] + 2 * psc.ibn[2]);
    fld_size.r = pv_set1_int(sz);				
    
    //-----------------------------------------------------
    //Set up some pointers for fields assuming x is uniform
    //this assumes xi always rounds down to psc.ilo[]
    //    sse2_real * restrict EXpoint = &pf->flds[EX*psc.fld_size + psc.ilo[0] - psc.ilg[0]]; 
    //    sse2_real * restrict EYpoint = &pf->flds[EY*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict EZpoint = &pf->flds[EZ*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict HXpoint = &pf->flds[HX*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict HYpoint = &pf->flds[HY*psc.fld_size + psc.ilo[0] - psc.ilg[0]];
    //    sse2_real * restrict HZpoint = &pf->flds[HZ*psc.fld_size + psc.ilo[0] - psc.ilg[0]];

  for(int n = 0; n < pp->n_part; n += VEC_SIZE) {

    //---------------------------------------------
    // Load Particles
    
    struct particle_vec p;

    LOAD_PARTS(p);
 
    //---------------------------------------------
    // Declare locals for computation      

    pvReal vxi, vyi, vzi, tmpx, tmpy, root, h2, h3;

    //---------------------------------------------
    // Half step positions with current momenta

    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    tmpx.r = pv_sub_real(root.r, ones.r);
    tmpx.r = pv_div_real(tmpx.r, eta.r);
    tmpy.r = pv_mul_real(p.mni.r, fnqs.r);
    tmpx.r = pv_mul_real(tmpy.r, tmpx.r);
    
    for(int m = 0; m < VEC_SIZE; m++){
      psc.p2A += tmpx.v[m]; //What's this for?
    }


    //---------------------------------------------
    // Declare and zero out what will become the 
    // current density form factors
    
    pvReal s0y[5], s0z[5], s1y[5], s1z[5];
    // I'm a little scared to use memset here, though it would probably work...
    for(int m = 0; m < 5; m++){
      s0y[m].r = pv_set1_real(0.0);
      s1y[m].r = pv_set1_real(0.0);
      s0z[m].r = pv_set1_real(0.0);
      s1z[m].r = pv_set1_real(0.0);
    }

    //---------------------------------------------
    // Prepare for field interpolation

    // FIXME : assuming true 2D and not touching x direction

    pvReal  H2, H3;
    pvInt j2, j3, l2, l3;

    find_index_Iround(&(p.yi), &dyi, &j2, &H2);
    find_index_Iround(&(p.zi), &dzi, &j3, &H3);

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



    find_index_minus_shift(&(p.yi), &dyi, &l2, &h2, &half);
    find_index_minus_shift(&(p.zi), &dzi, &l3, &h3, &half);
    
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

    //---------------------------------------------
    // Field Interpolation
    
    pvReal exq, eyq, ezq, hxq, hyq, hzq;

    INTERP_FIELD_YZ(EX,j2,j3,g,g,exq);
    INTERP_FIELD_YZ(EY,l2,j3,g,h,eyq);
    INTERP_FIELD_YZ(EZ,j2,l3,h,g,ezq);
    INTERP_FIELD_YZ(HX,l2,l3,h,h,hxq);
    INTERP_FIELD_YZ(HY,j2,l3,h,g,hyq);
    INTERP_FIELD_YZ(HZ,l2,j3,g,h,hzq);

    //---------------------------------------------
    // Calculate current density form factors
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

    //---------------------------------------------
    // Advance momenta one full dt

    push_pi_dt(&p, &exq, &eyq, &ezq, &hxq, &hyq, &hzq, &dqs);

    //---------------------------------------------
    // Half step particles with new momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vyi, &vzi, &yl, &zl);

    //---------------------------------------------
    // Store particles 

    STORE_PARTS(p);

    tmpx.r = pv_sub_real(root.r, ones.r);
    tmpx.r = pv_div_real(tmpx.r, eta.r);
    tmpy.r = pv_mul_real(p.mni.r, fnqs.r);
    tmpx.r = pv_mul_real(tmpy.r, tmpx.r);

    for(int m = 0; m < VEC_SIZE; m++){
      psc.p2B += tmpx.v[m]; 
    }
    

    //---------------------------------------------
    // update charge densities removed. Now done only on outstep.
    
    //---------------------------------------------
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
    pvInt shifty, shiftz; //FIXME : need to check if this really helps
    shifty.r = pv_sub_int(k2.r,j2.r);
    shiftz.r = pv_sub_int(k3.r,j3.r);
    for(int p=0; p < VEC_SIZE; p++){
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

    int l2min = 1;
    int l2max = 4;
    int l3min = 1;
    int l3max = 4;
    
    //---------------------------------------------
    // Finding bounds for updating current densities
    // I don't like keeping these branches in, but 
    // I don't think it can be helped. The construct
    // below allows the calculations to be done in vectors
    // even if the particles have different l{2,3}{min,max}
 
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
    
    fnqx.r = pv_mul_real(p.wni.r, fnqs.r);
    fnqx.r = pv_mul_real(p.qni.r, fnqx.r);
    fnqx.r = pv_mul_real(vxi.r, fnqx.r);

    fnqy.r = pv_mul_real(p.wni.r, fnqys.r);
    fnqy.r = pv_mul_real(p.qni.r, fnqy.r);
    
    fnqz.r = pv_mul_real(p.wni.r, fnqzs.r);
    fnqz.r = pv_mul_real(p.qni.r, fnqz.r);
        
    pvReal jzh[5]; // As per Will's suggestion, using the minimal 
                   // variable here, and for jyh below, gives a nice
                   // performance boost
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
	
	
	//---------------------------------------------
	// Store into the 'squished' currents
	for(int m = 0; m < VEC_SIZE; m++){
	  JSX(j2.v[m] + l2i - 2, j3.v[m] + l3i - 2) += wx.v[m]; //undo index adjustment of l(2,3)i
	  JSY(j2.v[m] + l2i - 2, j3.v[m] + l3i -2) += jyh.v[m]; //undo index adjustment of l(2,3)i
    	  JSZ(j2.v[m] + l2i - 2, j3.v[m] + l3i -2) += jzh[l2i].v[m]; //undo index adjustment of l(2
    	}
      }
    }   
  }
  
  //---------------------------------------------
  // Store the squished currents into the global currents
  // FIXME: Assumes x-dir always rounds down to 0!
  for(int iz = -psc.ibn[2]; iz < patch->ldims[2] + psc.ibn[2]; iz++){
    for(int iy = -psc.ibn[1]; iy < patch->ldims[2] + psc.ibn[1]; iy++){ 
      F3_SSE2(pf, JXI, 0,iy,iz) = JSX(iy, iz);
      F3_SSE2(pf, JYI, 0,iy,iz) = JSY(iy, iz);
      F3_SSE2(pf, JZI, 0,iy,iz) = JSZ(iy, iz);
    }
  }
  // Always pick up your trash
  free(s_jxi);
  free(s_jyi);
  free(s_jzi);
}

#undef JSX
#undef JSY
#undef JSZ

// ======================================================================

void
psc_push_particles_sse2_push_yz_a(struct psc_push_particles *push,
				  struct psc_mparticles *particles_base,
				  struct psc_mfields *flds_base)
{
  particles_sse2_t pp;
  mparticles_sse2_get(&pp, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz_a", 1., 0, 0);
  }
  prof_start(pr);
  do_push_part_yz_a(&pp);
  prof_stop(pr);

  mparticles_sse2_put(&pp, particles_base);
}

void
psc_push_particles_sse2_push_yz_b(struct psc_push_particles *push,
				  struct psc_mparticles *particles_base,
				  struct psc_mfields *flds_base)
{
  particles_sse2_t pp;
  fields_sse2_t pf;
  mparticles_sse2_get(&pp, particles_base);
  fields_sse2_get(&pf, EX, EX + 6, flds_base);

  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz_b", 1., 0, 0);
  }
  prof_start(pr);
  do_push_part_yz_b(&pp, &pf);
  prof_stop(pr);

  mparticles_sse2_put(&pp, particles_base);
  fields_sse2_put(&pf, JXI, JXI + 3, flds_base);
}

void
psc_push_particles_sse2_push_yz(struct psc_push_particles *push,
				struct psc_mparticles *particles_base,
				struct psc_mfields *flds_base)
{
  particles_sse2_t pp;
  fields_sse2_t pf;
  mparticles_sse2_get(&pp, particles_base);
  fields_sse2_get(&pf, EX, EX + 6, flds_base);

  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_yz", 1., 0, 0);
  }
  prof_start(pr);
  do_push_part_yz(&pp, &pf);
  prof_stop(pr);

  mparticles_sse2_put(&pp, particles_base);
  fields_sse2_put(&pf, JXI, JXI + 3, flds_base);
}

/// \file sse2_push_part_yz.c SSE2 implementation of the yz particle pusher.
///
/// Includes full pusher, and part {a,b} sections for unit testing.
