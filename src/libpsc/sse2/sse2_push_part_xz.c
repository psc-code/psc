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

//---------------------------------------------
// push_xi_halfdt - advances the posistions .5dt
// using the given velocities
static inline void
push_xi_halfdt(struct particle_vec *p, pvReal *vxi, pvReal *vzi, pvReal *xl, pvReal *zl){
  pvReal tmpx, tmpz;
  tmpx.r = pv_mul_real(vxi->r, xl->r);
  tmpz.r = pv_mul_real(vzi->r, zl->r);
  p->xi.r = pv_add_real(p->xi.r, tmpx.r);
  p->zi.r = pv_add_real(p->zi.r, tmpz.r);
}



//---------------------------------------------
/// SSE2 implementation of the 
/// xz particle pusher
static void
do_push_part_xz(particles_sse2_t *pp, fields_sse2_t *pf)
{

//-----------------------------------------------------
// Initialization stuff (not sure what all of this is for)
  
  psc.p2A = 0.;
  psc.p2B = 0.;
  
  for (int m = JXI; m <= JZI; m++) {
    memset(&pf->flds[m*psc.fld_size], 0, psc.fld_size * sizeof(sse2_real));
  }  

  //---------------------------------------------
  // An implementation of Will's 'squished' currents
  // that excludes the x direction all together
  sse2_real * restrict s_jxi, * restrict s_jyi, * restrict s_jzi;
  s_jxi = calloc((psc.img[0] * psc.img[2]), sizeof(sse2_real));
  s_jyi = calloc((psc.img[0] * psc.img[2]), sizeof(sse2_real));
  s_jzi = calloc((psc.img[0] * psc.img[2]), sizeof(sse2_real));

  // -------------------------------
  // Macros for accessing the 'squished' currents allocated above
#define JSX(indx1, indx3) s_jxi[((indx3) - psc.ilg[2])*psc.img[0] + ((indx1) - psc.ilg[0])]
#define JSY(indx1, indx3) s_jyi[((indx3) - psc.ilg[2])*psc.img[0] + ((indx1) - psc.ilg[0])]
#define JSZ(indx1, indx3) s_jzi[((indx3) - psc.ilg[2])*psc.img[0] + ((indx1) - psc.ilg[0])]
  
  //-----------------------------------------------------
  //Set vector forms of numbers 1, 0.5, etc
  init_vec_numbers();

  //-----------------------------------------------------
  //Vector versions of floating point parameters
    pvReal dt, xl, zl, eta, dqs,fnqs, dxi, dyi, dzi, fnqxs, fnqys, fnqzs;
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
    xl.r = pv_mul_real(half.r, dt.r);					
    zl.r = pv_mul_real(half.r, dt.r);					
    fnqxs.r = pv_set1_real(fnqxsfl);				       
    fnqys.r = pv_set1_real(fnqysfl);					
    fnqzs.r = pv_set1_real(fnqzsfl);			
    
    //-----------------------------------------------------
    //Vector versions of integer parameters
    pvInt ilg[3], img[3], fld_size;					
    ilg[0].r = pv_set1_int(psc.ilg[0]);					
    ilg[1].r = pv_set1_int(psc.ilg[1]);					
    ilg[2].r = pv_set1_int(psc.ilg[2]);					
    img[0].r = pv_set1_int(psc.img[0]);					
    img[1].r = pv_set1_int(psc.img[1]);					
    img[2].r = pv_set1_int(psc.img[2]);					
    fld_size.r = pv_set1_int(psc.fld_size);				
    
    //-----------------------------------------------------
    //Set up some pointers for fields assuming y is uniform
    
    sse2_real * restrict EXpoint = &pf->flds[EX*psc.fld_size 
						 + (psc.ilo[1] - psc.ilg[1])*psc.img[0] ];
   
    sse2_real * restrict EYpoint = &pf->flds[EY*psc.fld_size 
						 + (psc.ilo[1] - psc.ilg[1])*psc.img[0] ];
    
    sse2_real * restrict EZpoint = &pf->flds[EZ*psc.fld_size 
						 + (psc.ilo[1] - psc.ilg[1])*psc.img[0] ];
    
    sse2_real * restrict HXpoint = &pf->flds[HX*psc.fld_size 
						 + (psc.ilo[1] - psc.ilg[1])*psc.img[0] ];
    
    sse2_real * restrict HYpoint = &pf->flds[HY*psc.fld_size 
						 + (psc.ilo[1] - psc.ilg[1])*psc.img[0] ];
    
    sse2_real * restrict HZpoint = &pf->flds[HZ*psc.fld_size 
						 + (psc.ilo[1] - psc.ilg[1])*psc.img[0] ];

  for(int n = 0; n < psc.pp.n_part; n += VEC_SIZE) {

    //---------------------------------------------
    // Load Particles
    
    struct particle_vec p;

    LOAD_PARTS(p);
 
    //---------------------------------------------
    // Declare locals for computation      

    pvReal vxi, vyi, vzi, tmpx, tmpy, root, h1, h3;

    //---------------------------------------------
    // Half step positions with current momenta

    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vxi, &vzi, &xl, &zl);

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
    
    pvReal s0x[5], s0z[5], s1x[5], s1z[5];
    // I'm a little scared to use memset here, though it would probably work...
    for(int m = 0; m < 5; m++){
      s0x[m].r = pv_set1_real(0.0);
      s1x[m].r = pv_set1_real(0.0);
      s0z[m].r = pv_set1_real(0.0);
      s1z[m].r = pv_set1_real(0.0);
    }

    //---------------------------------------------
    // Prepare for field interpolation

    pvReal  H1, H3;
    pvInt j1, j3, l1, l3;

    find_index_Iround(&(p.xi), &dxi, &j1, &H1);
    find_index_Iround(&(p.zi), &dzi, &j3, &H3);


    pvReal gmx, gmz, gOx, gOz, glx, glz;// NB: this is gO as in octogon instead of g0, 
                                         // and gl as in library instead of g1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of g?[3], which 
                                         // seems to improve performance


    ip_to_grid_m(&H1, &gmx);
    ip_to_grid_m(&H3, &gmz);
     
    ip_to_grid_O(&H1, &gOx);
    ip_to_grid_O(&H3, &gOz);

    ip_to_grid_l(&H1, &glx);
    ip_to_grid_l(&H3, &glz);



    find_index_minus_shift(&(p.xi), &dxi, &l1, &h1, &half);
    find_index_minus_shift(&(p.zi), &dzi, &l3, &h3, &half);
    
    pvReal hmx, hmz, hOx, hOz, hlx, hlz; // NB: this is hO as in octogon instead of h0, 
                                         // and hl as in library instead of h1.
                                         // This is done to let me paste in the CPP macro
                                         // instead of using an array of h?[3], which 
                                         // seems to improve performance


    ip_to_grid_m(&h1, &hmx);
    ip_to_grid_m(&h3, &hmz);
    
    ip_to_grid_O(&h1, &hOx);
    ip_to_grid_O(&h3, &hOz);
    
    ip_to_grid_l(&h1, &hlx);
    ip_to_grid_l(&h3, &hlz);


    //---------------------------------------------
    // Field Interpolation
    
    pvReal exq, eyq, ezq, hxq, hyq, hzq;

    INTERP_FIELD_XZ(EX,l1,j3,g,h,exq);
    INTERP_FIELD_XZ(EY,j1,j3,g,g,eyq);
    INTERP_FIELD_XZ(EZ,j1,l3,h,g,ezq);
    INTERP_FIELD_XZ(HX,j1,l3,h,g,hxq);
    INTERP_FIELD_XZ(HY,l1,l3,h,h,hyq);
    INTERP_FIELD_XZ(HZ,l1,j3,g,h,hzq);

    //---------------------------------------------
    // Calculate current density form factors
    // indexing here departs from FORTRAN a little bit
    // Simply moving this down here seems to improve performance
    // Possibly giving the processor something to do while it's waiting 
    // for fields to read in
    
    form_factor_m(&H1, &s0x[1]);
    form_factor_O(&H1, &s0x[2]);
    form_factor_l(&H1, &s0x[3]);

    form_factor_m(&H3, &s0z[1]);
    form_factor_O(&H3, &s0z[2]);
    form_factor_l(&H3, &s0z[3]);

    //---------------------------------------------
    // Advance momenta one full dt

    push_pi_dt(&p, &exq, &eyq, &ezq, &hxq, &hyq, &hzq, &dqs);

    //---------------------------------------------
    // Half step particles with new momenta
    
    calc_vi(&p, &vxi, &vyi, &vzi, &root);

    push_xi_halfdt(&p, &vxi, &vzi, &xl, &zl);

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
    
    push_xi_halfdt(&p, &vxi, &vzi, &xl, &zl);

    pvInt k1,k3;

    find_index_Iround(&(p.xi), &dxi, &k1, &h1);
    find_index_Iround(&(p.zi), &dzi, &k3, &h3);

    // God help me, there's some things I just can't fgure out how to parallelize
    // The g-- here are just temporary variables. I can't do the assignments in parallel, 
    // but I'll be damned if I can't do the FLOPS in parallel

    form_factor_m(&h1, &gmx);
    form_factor_O(&h1, &gOx);
    form_factor_l(&h1, &glx);

    form_factor_m(&h3, &gmz);
    form_factor_O(&h3, &gOz);
    form_factor_l(&h3, &glz);
       
    // All the indices below have +2 compared to the FORTRAN for C array access
    pvInt shiftx, shiftz; //FIXME : need to check if this really helps
    shiftx.r = pv_sub_int(k1.r,j1.r);
    shiftz.r = pv_sub_int(k3.r,j3.r);
    for(int p=0; p < VEC_SIZE; p++){
      s1x[shiftx.v[p] + 1].v[p] = gmx.v[p];
      s1x[shiftx.v[p] + 2].v[p] = gOx.v[p];
      s1x[shiftx.v[p] + 3].v[p] = glx.v[p];
      s1z[shiftz.v[p] + 1].v[p] = gmz.v[p];
      s1z[shiftz.v[p] + 2].v[p] = gOz.v[p];
      s1z[shiftz.v[p] + 3].v[p] = glz.v[p];
    }

    for(int m=0; m<5; m++){
      s1x[m].r = pv_sub_real(s1x[m].r, s0x[m].r);
      s1z[m].r = pv_sub_real(s1z[m].r, s0z[m].r);
    }

    // Macros only good for registers
#define S0X(off) s0x[(off)+2].r
#define S0Z(off) s0z[(off)+2].r
#define S1X(off) s1x[(off)+2].r
#define S1Z(off) s1z[(off)+2].r

    int l1min = -1;
    int l1max = 1;
    int l3min = -1;
    int l3max = 1;
    
    //---------------------------------------------
    // Finding bounds for updating current densities
    // I don't like keeping these branches in, but 
    // I don't think it can be helped. The construct
    // below allows the calculations to be done in vectors
    // even if the particles have different l{2,3}{min,max}
 
    for(int m = 0; m < VEC_SIZE; m++){
      if(k1.v[m]==j1.v[m]){
	s1x[0].v[m] = 0.0;
	s1x[4].v[m] = 0.0;
	s0x[0].v[m] = 0.0;
	s0x[4].v[m] = 0.0;
      }
      else if (k1.v[m]==j1.v[m]-1){
	s1x[4].v[m] = 0.0;
	s0x[4].v[m] = 0.0;
	l1min = -2;
      }
      else {// if (k1.v[p]==j1.v[p]+1)
	s1x[0].v[m] = 0.0;
	s0x[0].v[m] = 0.0;
	l1max = 2;
      }
      if (k3.v[m]==j3.v[m]){
	s1z[0].v[m] = 0.0;
	s1z[4].v[m] = 0.0;
	s0z[0].v[m] = 0.0;
	s0z[4].v[m] = 0.0;
      }
      else if(k3.v[m]==j3.v[m]-1){
	s1z[4].v[m] = 0.0;
	s0z[4].v[m] = 0.0;
	l3min = -2;
      }
      else { //(k3.v[p]==j3.v[p]+1)
	s1z[0].v[m] = 0.0;
	s0z[0].v[m] = 0.0;
	l3max = 2;
      }
    }

    pvReal fnqx, fnqy, fnqz;
    
    fnqx.r = pv_mul_real(p.wni.r, fnqxs.r);
    fnqx.r = pv_mul_real(p.qni.r, fnqx.r);

    fnqy.r = pv_mul_real(p.wni.r, fnqs.r);
    fnqy.r = pv_mul_real(p.qni.r, fnqy.r);
    fnqy.r = pv_mul_real(vyi.r, fnqy.r);

    fnqz.r = pv_mul_real(p.wni.r, fnqzs.r);
    fnqz.r = pv_mul_real(p.qni.r, fnqz.r);
        
    pvReal jzh[5]; // As per Will's suggestion, using the minimal 
                   // variable here, and for jxh below, gives a nice
                   // performance boost
    memset(&jzh, 0, 5 * sizeof(pvReal));
    
#define JZH(i) jzh[(i)+2].r

    for(int l3i=l3min; l3i<=l3max; l3i++){
      pvReal jxh;
      jxh.r = pv_set1_real(0.0);
      for(int l1i=l1min; l1i<=l1max; l1i++){
	pvReal wx, wy, wz;
	
	wx.r = pv_mul_real(half.r, S1Z(l3i));
	wx.r = pv_add_real(S0Z(l3i), wx.r);
	wx.r = pv_mul_real(S1X(l1i), wx.r);

	wy.r = pv_mul_real(half.r,S1X(l1i));
	wy.r = pv_add_real(S0X(l1i), wy.r);
	wy.r = pv_mul_real(S0Z(l3i), wy.r);
	tmpx.r = pv_mul_real(half.r, S0X(l1i));
	tmpy.r = pv_mul_real(third.r, S1X(l1i));
	tmpx.r = pv_add_real(tmpx.r, tmpy.r);
	tmpx.r = pv_mul_real(tmpx.r, S1Z(l3i));
	wy.r = pv_add_real(wy.r, tmpx.r);
		
	wz.r = pv_mul_real(half.r, S1X(l1i));
	wz.r = pv_add_real(S0X(l1i), wz.r);
	wz.r = pv_mul_real(S1Z(l3i), wz.r);
	
	wx.r = pv_mul_real(fnqx.r, wx.r);
	jxh.r = pv_sub_real(jxh.r, wx.r);
	wy.r = pv_mul_real(fnqy.r, wy.r);
	wz.r = pv_mul_real(fnqz.r, wz.r);
	JZH(l1i) = pv_sub_real(JZH(l1i), wz.r);
	
	
	//---------------------------------------------
	// Store into the 'squished' currents
	for(int m = 0; m < VEC_SIZE; m++){
	  JSX(j1.v[m] + l1i, j3.v[m] + l3i) += jxh.v[m];
	  JSY(j1.v[m] + l1i, j3.v[m] + l3i) += wy.v[m];
    	  JSZ(j1.v[m] + l1i, j3.v[m] + l3i) += jzh[l1i+2].v[m];
    	}
      }
    }   
  }
#undef JZH  
#undef S0X
#undef S0Z
#undef S1X
#undef S1Z
  //---------------------------------------------
  // Store the squished currents into the global currents
  // FIXME: Assumes x-dir always rounds down to 0!
  for(int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++){
    for(int ix = psc.ilg[0]; ix < psc.ihg[0]; ix++){ 
      F3_SSE2(pf,JXI,ix,psc.ilo[1]+0,iz) = JSX(ix, iz);
      F3_SSE2(pf,JYI,ix,psc.ilo[1]+0,iz) = JSY(ix, iz);
      F3_SSE2(pf,JZI,ix,psc.ilo[1]+0,iz) = JSZ(ix, iz);
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

void sse2_push_part_xz()
{
  particles_sse2_t pp;
  fields_sse2_t pf;
  mparticles_sse2_get(&pp);
  fields_sse2_get(&pf, EX, EX+6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("sse2_part_xz", 1., 0, psc.pp.n_part * 9 * sizeof(sse2_real));
  }
  prof_start(pr);
  do_push_part_xz(&pp, &pf);
  prof_stop(pr);
  
  mparticles_sse2_put(&pp);
  fields_sse2_put(&pf, JXI, JXI+3);
}

/// \file sse2_push_part_xz.c SSE2 implementation of the xz particle pusher.
///
/// Only the full push is included here. It has been converted from the yz pusher,
/// and errors may remain.
