#include "spu_particles.h"
#include "psc_spu.h"
#include "psc_spu_2d.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include <simdmath.h>

#ifdef __SPU__
#include <spu_mfcio.h>
#else
#include "spu_mfcio_c.h"
#endif

enum{
  ls_EX, ls_EY, ls_EZ,
  ls_HX, ls_HY, ls_HZ,
  ls_JXI, ls_JYI, ls_JZI,
  NR_LSFLDS,
};

v_real half, one, threefourths, onepfive,third,zero;

static inline void
find_index(v_real *xi, v_real *dxi, v_int *j, v_real *h)
{
  v_real tmp;
  tmp = spu_mul(*xi, *dxi);
  *h = spu_round_real(tmp);
  *j = spu_round_int(tmp);
  *h = spu_sub(*h, tmp);
} 


static inline void
find_index_minus_shift(v_real *xi, v_real *dxi, v_int *j, v_real *h, v_real *shift)
{
  v_real tmp;
  tmp = spu_msub(*xi, *dxi, *shift);
  *h = spu_round_real(tmp);
  *j = spu_round_int(tmp);
  *h = spu_sub(*h, tmp);
} 

static inline void
ip_to_grid_m(v_real *h, v_real * xmx)
{
  *xmx = spu_add(half, *h);
  *xmx = spu_mul(*xmx, *xmx);
  *xmx = spu_mul(half, *xmx);
}

static inline void
ip_to_grid_O(v_real *h, v_real *xOx)
{
  *xOx = spu_mul( *h , *h);
  *xOx = spu_sub(threefourths, *xOx);

}

static inline void
ip_to_grid_l(v_real *h, v_real * xlx)
{
  *xlx = spu_sub(half, *h);
  *xlx = spu_mul( *xlx, *xlx);
  *xlx = spu_mul(half, *xlx);
}

static inline void
form_factor_m(v_real *h, v_real * restrict xmx)
{
  *xmx = spu_sub(*h,one);
  *xmx = spu_add(onepfive,*xmx); // h-1 always <0
  *xmx = spu_mul(*xmx, *xmx);
  *xmx = spu_mul(half, *xmx);
}

static inline void
form_factor_O(v_real *h, v_real * restrict xOx)
{
  *xOx = spu_mul(*h, *h);
  *xOx = spu_sub(threefourths, *xOx);
}

static inline void
form_factor_l(v_real *h, v_real * restrict xlx)
{
  *xlx = spu_add(one, *h); 
  *xlx = spu_sub(onepfive, *xlx);//h+1 always >0
  *xlx = spu_mul(*xlx, *xlx);
  *xlx = spu_mul(half, *xlx);
}

#define F2_SPU_OFF(fldnr, jx, jy, jz)		\
      ((((((fldnr)								\
	   * psc_block.im[2] + ((jz)-psc_block.ib[2]))			\
      * psc_block.im[1] + ((jy)-psc_block.ib[1]))				\
     * psc_block.im[0] + ((jx)-psc_block.ib[0]))))


#define IP_FIELD_SPU(fldnr, jy, jz, outer_coeff, inner_coeff, field) {	\
  v_real _tmp1, _tmp2, _tmp3;						\
  v_real _in_A[2];							\
  v_real _in_B[2];							\
  int jy_0 = (int)spu_extract(jy,0);					\
  int jy_1 = (int)spu_extract(jy,1);					\
  int jz_0 = (int)spu_extract(jz,0);					\
  int jz_1 = (int)spu_extract(jz,1);					\
									\
  int off_0 = F2_SPU_OFF(fldnr, 0, jy_0 - 1, jz_0 - 1);			\
  int off_1 = F2_SPU_OFF(fldnr, 0, jy_1 - 1, jz_1 - 1);			\
									\
  int align_0 = off_0 & 1;						\
  int align_1 = off_1 & 1;						\
									\
  _in_A[0] = *(v_real *)(ls_fld + off_0 - align_0);				\
  _in_A[1] = *(v_real *)(ls_fld + off_0 + 2 - align_0);				\
  _in_B[0] = *(v_real *)(ls_fld + off_1 - align_1);				\
  _in_B[1] = *(v_real *)(ls_fld + off_1 + 2 - align_1);				\
  _tmp1 = spu_shuffle(_in_A[0], _in_B[0], fld_ip_pat[align_0][align_1]); \
  _tmp2 = spu_shuffle(_in_A[align_0], _in_B[align_1], fld_ip_pat[align_0 ^ 1][align_1 ^ 1]);\
  _tmp3 = spu_shuffle(_in_A[1], _in_B[1], fld_ip_pat[align_0][align_1]); \
  _tmp1 = spu_mul(inner_coeff##my,_tmp1);				\
  _tmp2 = spu_mul(inner_coeff##Oy,_tmp2);				\
  _tmp3 = spu_mul(inner_coeff##ly,_tmp3);				\
  _tmp1 = spu_add(_tmp1, _tmp2);					\
  _tmp1 = spu_add(_tmp1, _tmp3);					\
  field = spu_mul(outer_coeff##mz, _tmp1);				\
									\
  off_0 = F2_SPU_OFF(fldnr, 0, jy_0-1, jz_0);				\
  off_1 = F2_SPU_OFF(fldnr, 0, jy_1-1, jz_1);				\
  									\
  align_0 = off_0 & 1;							\
  align_1 = off_1 & 1;							\
									\
  _in_A[0] = *(v_real *)(ls_fld + off_0 - align_0);				\
  _in_A[1] = *(v_real *)(ls_fld + off_0 + 2 - align_0);				\
  _in_B[0] = *(v_real *)(ls_fld + off_1 - align_1);				\
  _in_B[1] = *(v_real *)(ls_fld + off_1 + 2 - align_1);				\
  _tmp1 = spu_shuffle(_in_A[0], _in_B[0], fld_ip_pat[align_0][align_1]); \
  _tmp2 = spu_shuffle(_in_A[align_0], _in_B[align_1], fld_ip_pat[align_0 ^ 1][align_1 ^ 1]);\
  _tmp3 = spu_shuffle(_in_A[1], _in_B[1], fld_ip_pat[align_0][align_1]); \
  _tmp1 = spu_mul(inner_coeff##my,_tmp1);				\
  _tmp2 = spu_mul(inner_coeff##Oy,_tmp2);				\
  _tmp3 = spu_mul(inner_coeff##ly,_tmp3);				\
  _tmp1 = spu_add(_tmp1, _tmp2);					\
  _tmp1 = spu_add(_tmp1, _tmp3);					\
  field = spu_madd(outer_coeff##Oz, _tmp1, field);			\
									\
  off_0 = F2_SPU_OFF(fldnr, 0, jy_0 - 1, jz_0+1);				\
  off_1 = F2_SPU_OFF(fldnr, 0, jy_1 - 1, jz_1+1);				\
  									\
  align_0 = off_0 & 1;							\
  align_1 = off_1 & 1;							\
									\
  _in_A[0] = *(v_real *)(ls_fld + off_0 - align_0);				\
  _in_A[1] = *(v_real *)(ls_fld + off_0 + 2 - align_0);				\
  _in_B[0] = *(v_real *)(ls_fld + off_1 - align_1);				\
  _in_B[1] = *(v_real *)(ls_fld + off_1 + 2 - align_1);				\
  _tmp1 = spu_shuffle(_in_A[0], _in_B[0], fld_ip_pat[align_0][align_1]); \
  _tmp2 = spu_shuffle(_in_A[align_0], _in_B[align_1], fld_ip_pat[align_0 ^ 1][align_1 ^ 1]); \
  _tmp3 = spu_shuffle(_in_A[1], _in_B[1], fld_ip_pat[align_0][align_1]); \
  _tmp1 = spu_mul(inner_coeff##my,_tmp1);				\
  _tmp2 = spu_mul(inner_coeff##Oy,_tmp2);				\
  _tmp3 = spu_mul(inner_coeff##ly, _tmp3);				\
  _tmp1 = spu_add(_tmp1, _tmp2);					\
  _tmp1 = spu_add(_tmp1, _tmp3);					\
  field = spu_madd(outer_coeff##lz, _tmp1, field);			\
  }

int
spu_push_part_2d(void){


#if PRINT_DEBUG

  printf("[[%#llx] start ea: %#llx end ea: %#llx \n", spu_ctx.spe_id, psc_block.part_start, psc_block.part_end);

#endif

  fields_c_real_t ls_fld[NR_LSFLDS*32*32] __attribute__((aligned(128)));

  // If we try to dma in all the field data at once, we run into the 
  // mfc 16KB per request limit. For now, we'll split it into 6 different
  // requests. Depending on whether the limit is 16KB = 1000 * 16B 
  // or 1024*16B, we could group these two by two.
  
  spu_dma_get(ls_fld, 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,EX,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_EY*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,EY,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_EZ*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,EZ,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_HX*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,HX,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_HY*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,HY,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_HZ*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,HZ,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  memset(&ls_fld[ls_JXI*32*32], 0, 3*32*32*sizeof(fields_c_real_t));

  unsigned long long cp_ea = psc_block.part_start;
  unsigned long long np_ea; 
  
  particle_cbe_t _bufferA[2] __attribute__((aligned(16))), 
    _bufferB[2] __attribute__((aligned(16))), 
    _bufferC[2] __attribute__((aligned(16)));

  buff.plb1 = &(_bufferA[0]);
  buff.plb2 = &(_bufferA[1]);
  buff.lb1 = &(_bufferB[0]); 
  buff.lb2 = &(_bufferB[1]);
  buff.sb1 = &(_bufferC[0]);
  buff.sb2 = &(_bufferC[1]);

  // Get the first two particles
  //spu_dma_get(buff.lb1, cp_ea, 2*sizeof(particle_cbe_t));

  first_preload_particle(buff.plb1, cp_ea, 2*sizeof(particle_cbe_t));

  np_ea = cp_ea + 2*sizeof(particle_cbe_t);

  // we have some stuff to do while we wait for
  // it to come in.
  
  // insert assignment, promotions, and constant 
  // calculations here.

  v_real dt, yl, zl, dyi, dzi, dqs,fnqs,fnqxs,fnqys,fnqzs;
  dt = spu_splats(spu_ctx.dt);
  half = spu_splats(0.5);
  yl = spu_mul(half, dt);
  zl = spu_mul(half, dt);
  dyi = spu_splats(1./spu_ctx.dx[1]);
  dzi = spu_splats(1./spu_ctx.dx[2]);
  one = spu_splats(1.0);
  threefourths = spu_splats(0.75);
  onepfive = spu_splats(1.5);
  third = spu_splats(1./3.);
  zero = spu_splats(0.0);
  dqs = spu_splats(0.5 * spu_ctx.eta * spu_ctx.dt);
  fnqs = spu_splats(spu_ctx.fnqs);
  fnqxs = spu_splats(spu_ctx.dx[0] * spu_ctx.fnqs / spu_ctx.dt);
  fnqys = spu_splats(spu_ctx.dx[1] * spu_ctx.fnqs / spu_ctx.dt);
  fnqzs = spu_splats(spu_ctx.dx[2] * spu_ctx.fnqs / spu_ctx.dt);
  unsigned long long end = psc_block.part_end; 

  int run = 1;
  int n = 0;

  particle_cbe_t null_part; 
    

  do {

    // rotate the buffers 
    particle_cbe_t *btmp1, *btmp2;
    btmp1 = buff.sb1;
    btmp2 = buff.sb2;
    buff.sb1 = buff.lb1;
    buff.sb2 = buff.lb2;
    buff.lb1 = buff.plb1;
    buff.lb2 = buff.plb2;
    buff.plb1 = btmp1;
    buff.plb2 = btmp2;


    // issue dma request for particle we will need 
    // next time through the loop.

    if(__builtin_expect((np_ea < end),1)) {

      loop_preload_particle(buff.plb1, np_ea, 2 * sizeof(particle_cbe_t));

    }    
    // we may need to insert some padding here, so we have to stop and check.
    // The last particle is going to be very slow (probably).
    else if(__builtin_expect(((end - cp_ea) != 2*sizeof(particle_cbe_t)),0)){ 
#if PRINT_DEBUG
      fprintf(stderr, "Detected odd particle out\n"); 
#endif
      wait_for_preload();
      null_part = *(buff.lb1);
      null_part.wni = 0.0;
      buff.lb2 = &null_part;
    } else {
      wait_for_preload();
    }


    v_real xi, yi, zi, pxi, pyi, pzi, qni, mni, wni;
  
    LOAD_PARTICLES_SPU;

    v_real vxi,vyi,vzi,root,tmpx,tmpy,tmpz; 
    
    tmpx = spu_mul(pxi,pxi);
    tmpy = spu_mul(pyi,pyi);
    tmpz = spu_mul(pzi, pzi);
    
    tmpx = spu_add(tmpx, tmpy);
    tmpz = spu_add(one, tmpz);
    tmpx = spu_add(tmpx, tmpz);
    root = spu_sqrt(tmpx);
    
    vxi = spu_div(pxi, root);
    vyi = spu_div(pyi, root);
    vzi = spu_div(pzi, root);
    
    tmpy = spu_mul(vyi, yl);
    tmpz = spu_mul(vzi, zl);
    
    yi = spu_add(yi, tmpy);
    zi = spu_add(zi, tmpz);

    v_real gmy, gmz, gOy, gOz, gly, glz, H2, H3, h2, h3;
    v_int j2, j3, l2, l3;
      
    find_index(&yi, &dyi, &j2, &H2);
    find_index(&zi, &dzi, &j3, &H3);
    
    ip_to_grid_m(&H2, &gmy);
    ip_to_grid_m(&H3, &gmz);
    
    ip_to_grid_O(&H2, &gOy);
    ip_to_grid_O(&H3, &gOz);
    
    ip_to_grid_l(&H2, &gly);
    ip_to_grid_l(&H3, &glz);
    
    find_index_minus_shift(&yi, &dyi, &l2, &h2, &half);
    find_index_minus_shift(&zi, &dzi, &l3, &h3, &half);
    
    v_real hmy, hmz, hOy, hOz, hly, hlz;
    
    ip_to_grid_m(&h2, &hmy);
    ip_to_grid_m(&h3, &hmz);
    
    ip_to_grid_O(&h2, &hOy);
    ip_to_grid_O(&h3, &hOz);
    
    ip_to_grid_l(&h2, &hly);
    ip_to_grid_l(&h3, &hlz);
    
    // Field interpolation here. urg...
    
    v_real exq, eyq, ezq, hxq, hyq, hzq;
    
    IP_FIELD_SPU(ls_EX,j2,j3,g,g,exq);
    IP_FIELD_SPU(ls_EY,l2,j3,g,h,eyq);
    IP_FIELD_SPU(ls_EZ,j2,l3,h,g,ezq);
    IP_FIELD_SPU(ls_HX,l2,l3,h,h,hxq);
    IP_FIELD_SPU(ls_HY,j2,l3,h,g,hyq);
    IP_FIELD_SPU(ls_HZ,l2,j3,g,h,hzq);
    
    v_real s0y[5], s0z[5], s1y[5], s1z[5];
    for(int mp=0; mp<5; mp++){
      s0y[mp] = spu_splats(0.0);
      s1y[mp] = spu_splats(0.0);
      s0z[mp] = spu_splats(0.0);
      s1z[mp] = spu_splats(0.0);
    }
    form_factor_m(&H2, &s0y[1]);
    form_factor_O(&H2, &s0y[2]);
    form_factor_l(&H2, &s0y[3]);
    
    form_factor_m(&H3, &s0z[1]);
    form_factor_O(&H3, &s0z[2]);
    form_factor_l(&H3, &s0z[3]);
    
    
    // It would be nice to offload the momentum to a function, but it may not
    // be necessary
    
    v_real dq;
    
    dq = spu_mul(qni, dqs);
    dq = spu_div(dq, mni);
    
    /*
      pxi = spu_madd(dq, exq, pxi);
      pyi = spu_madd(dq, eyq, pyi);
      pzi = spu_madd(dq, ezq, pzi);
      
    */
    v_real dqex, dqey, dqez;
    dqex = spu_mul(dq, exq);
    dqey = spu_mul(dq, eyq);
    dqez = spu_mul(dq, ezq);
    
    pxi = spu_add(pxi, dqex);
    pyi = spu_add(pyi, dqey);
    pzi = spu_add(pzi, dqez);
    
    v_real taux, tauy, tauz, txx, tyy, tzz, t2xy, t2xz, t2yz, pxp, pyp, pzp;
    
    tmpx = spu_mul(pxi, pxi);
    tmpy = spu_mul(pyi, pyi);
    tmpz = spu_mul(pzi, pzi);
    root = spu_add(one, tmpx);
    tmpy = spu_add(tmpy, tmpz);
    root = spu_add(root, tmpy);
    root = spu_sqrt(root);
    root = spu_div(dq, root);
    
    taux = spu_mul(hxq, root);
    tauy = spu_mul(hyq, root);
    tauz = spu_mul(hzq, root);
    
    
    txx = spu_mul(taux, taux);
    tyy = spu_mul(tauy, tauy);
    tzz = spu_mul(tauz, tauz);
    t2xy = spu_mul(taux, tauy);
    t2xz = spu_mul(taux, tauz);
    t2yz = spu_mul(tauy, tauz);
    t2xy = spu_add(t2xy, t2xy);
    t2xz = spu_add(t2xz, t2xz);
    t2yz = spu_add(t2yz, t2yz);
    
    v_real tau;
    tau = spu_add(one, txx);
    tmpx = spu_add(tyy, tzz);
    tau = spu_add(tau, tmpx);
    tau = spu_div(one, tau);
    
    taux = spu_add(taux, taux);
    tauy = spu_add(tauy, tauy);
    tauz = spu_add(tauz, tauz);
    
    //pxp
    tmpx = spu_add(one, txx);
    tmpx = spu_sub(tmpx, tyy);
    tmpx = spu_sub(tmpx, tzz);
    tmpx = spu_mul(tmpx, pxi);
    
    tmpy = spu_add(t2xy, tauz);
    tmpy = spu_mul(tmpy, pyi);
    
    tmpz = spu_sub(t2xz, tauy);
    tmpz = spu_mul(tmpz, pzi);
    
    pxp = spu_add(tmpx, tmpy);
    pxp = spu_add(pxp, tmpz);
    pxp = spu_mul(pxp, tau);
    
    //pyp
    tmpx = spu_sub(t2xy, tauz);
    tmpx = spu_mul(tmpx, pxi);
    
    tmpy = spu_sub(one, txx);
    tmpy = spu_add(tmpy, tyy);
    tmpy = spu_sub(tmpy, tzz);
    tmpy = spu_mul(tmpy, pyi);
    
    tmpz = spu_add(t2yz, taux);
    tmpz = spu_mul(tmpz, pzi);
    
    pyp = spu_add(tmpx, tmpy);
    pyp = spu_add(pyp, tmpz);
    pyp = spu_mul(pyp, tau);
    
    // pzp
    tmpx = spu_add(t2xz, tauy);
    tmpx = spu_mul(tmpx, pxi);
    
    tmpy = spu_sub(t2yz, taux);
    tmpy = spu_mul(tmpy, pyi);
    
    tmpz = spu_sub(one, txx);
    tmpz = spu_sub(tmpz, tyy);
    tmpz = spu_add(tmpz, tzz);
    tmpz = spu_mul(tmpz, pzi);
    
    pzp = spu_add(tmpx, tmpy);
    pzp = spu_add(pzp, tmpz);
    pzp = spu_mul(pzp, tau);
    
    pxi = spu_add(pxp, dqex);
    pyi = spu_add(pyp, dqey);
    pzi = spu_add(pzp, dqez);
    
    /*
      pxi = spu_madd(dq, exq, pxp);
      pyi = spu_madd(dq, eyq, pyp);
      pzi = spu_madd(dq, ezq, pzp);
    */
    // finish advancing particles
    tmpx = spu_mul(pxi,pxi);
    tmpy = spu_mul(pyi,pyi);
    tmpz = spu_mul(pzi,pzi);
    
    tmpx = spu_add(tmpx, tmpy);
    tmpz = spu_add(one, tmpz);
    tmpx = spu_add(tmpx, tmpz);
    root = spu_sqrt(tmpx);
    
    vxi = spu_div(pxi, root);
    vyi = spu_div(pyi, root);
    vzi = spu_div(pzi, root);
    
    tmpy = spu_mul(vyi, yl);
    tmpz = spu_mul(vzi, zl);
    
    yi = spu_add(yi, tmpy);
    zi = spu_add(zi, tmpz);
    

    STORE_PARTICLES_SPU;

    np_ea = cp_ea +  2 * sizeof(particle_cbe_t);
    // At this point, np_ea is one load ahead of the
    // current particle.


    if(__builtin_expect((np_ea >= end),0)) { // if we've run out of particles
      loop_store_particle(buff.lb1, cp_ea, (size_t) (end - cp_ea));
      run = 0; 
    }
    else {
      loop_store_particle(buff.lb1, cp_ea, 2 * sizeof(particle_cbe_t));
      cp_ea = np_ea;
      np_ea += 2*sizeof(particle_cbe_t);
      // cp_ea now points to the particle which will be used
      // next time in the loop. 
      // np_ea points to the one which needs to be pre-loaded.
    }

    yi = spu_add(yi, tmpy);
    zi = spu_add(zi, tmpz);

    v_int k2,k3;
    
    find_index(&yi, &dyi, &k2, &h2);
    find_index(&zi, &dzi, &k3, &h3);
    
    // God help me, there's some things I just can't fgure out how to parallelize
    // The g-- here are just temporary variables. I can't do the assignments in parallel, 
    // but I'll be damned if I can't do the FLOPS in parallel
    
    form_factor_m(&h2, &gmy);
    form_factor_O(&h2, &gOy);
    form_factor_l(&h2, &gly);
    
    form_factor_m(&h3, &gmz);
    form_factor_O(&h3, &gOz);
    form_factor_l(&h3, &glz);
    
    signed long long j2_scal[2], j3_scal[2];
    j2_scal[0] = spu_extract(j2, 0);
    j2_scal[1] = spu_extract(j2, 1);
    j3_scal[0] = spu_extract(j3, 0);
    j3_scal[1] = spu_extract(j3, 1);
    
    signed long long dfy[2];
    signed long long  dfz[2];
    dfy[0] = spu_extract(k2,0) - j2_scal[0];
    dfy[1] = spu_extract(k2,1) - j2_scal[1];
    dfz[0] = spu_extract(k3,0) - j3_scal[0];
    dfz[1] = spu_extract(k3,1) - j3_scal[1];
    
    
    for(int p=0; p < VEC_SIZE; p++){
      s1y[(int)dfy[p] + 1] = spu_sel(s1y[(int)dfy[p] + 1], gmy, (vector unsigned long long) element_assign[p]);
      s1y[(int)dfy[p] + 2] = spu_sel(s1y[(int)dfy[p] + 2], gOy, (vector unsigned long long) element_assign[p]);
      s1y[(int)dfy[p] + 3] = spu_sel(s1y[(int)dfy[p] + 3], gly, (vector unsigned long long) element_assign[p]);
      s1z[(int)dfz[p] + 1] = spu_sel(s1z[(int)dfz[p] + 1], gmz, (vector unsigned long long) element_assign[p]);
      s1z[(int)dfz[p] + 2] = spu_sel(s1z[(int)dfz[p] + 2], gOz, (vector unsigned long long) element_assign[p]);
      s1z[(int)dfz[p] + 3] = spu_sel(s1z[(int)dfz[p] + 3], glz, (vector unsigned long long) element_assign[p]);
    }
    
    for(int m=0; m<5; m++){
      s1y[m] = spu_sub(s1y[m], s0y[m]);
      s1z[m] = spu_sub(s1z[m], s0z[m]);
    }
    
    
    // I'm pretty sure this is exactly the same as the
    // branching in the serial code. Will test in 
    // sse2 code.
    
    int l2min = 1 + ((dfy[0] | dfy[1]) >> 1),
      l2max = 3 + (((dfy[0]>>1)^dfy[0]) | ((dfy[1]>>1)^dfy[1])),
      l3min = 1 + ((dfz[0] | dfz[1]) >> 1),
      l3max = 3 + (((dfz[0]>>1)^dfz[0]) | ((dfz[1]>>1)^dfz[1]));
    
    v_real fnqx, fnqy, fnqz;
    
    fnqx = spu_mul(wni, fnqs);
    fnqx = spu_mul(qni, fnqx);
    fnqx = spu_mul(vxi, fnqx);
    
    fnqy = spu_mul(wni, fnqys);
    fnqy = spu_mul(qni, fnqy);
    
    fnqz = spu_mul(wni, fnqzs);
    fnqz = spu_mul(qni, fnqz);
    
    v_real jzh[5]; // As per Will's suggestion, using the minimal 
    // variable here, and for jyh below, gives a nice
    // performance boost
    memset(jzh,0,5*sizeof(v_real));
    v_real sjx0_a, sjx0_b, sjx0_c;
    v_real sjx1_a, sjx1_b, sjx1_c;

    v_real sjy0_a, sjy0_b, sjy0_c;
    v_real sjy1_a, sjy1_b, sjy1_c;

    v_real sjz0_a, sjz0_b, sjz0_c;
    v_real sjz1_a, sjz1_b, sjz1_c;

    sjx0_c = zero;

    sjx1_c = zero;

    sjy0_c = zero;

    sjy1_c = zero;

    sjz0_c = zero;

    sjz1_c = zero;

    v_real jyh;
    for(int l3i=l3min; l3i<=l3max; l3i++){
      jyh = spu_splats(0.0);
	
      int store_off_0 = F2_SPU_OFF(ls_JXI,0,j2_scal[0] - 2, j3_scal[0] + l3i - 2);
      int store_off_1 = F2_SPU_OFF(ls_JXI,0,j2_scal[1] - 2, j3_scal[1] + l3i - 2);

      int l2i = 0;
      v_real wx, wy, wz;
      
#define CALC_J_POINT {				\
	wx = spu_mul(half,s1y[l2i]);		\
	wx = spu_add(s0y[l2i], wx);		\
	wx = spu_mul(s0z[l3i], wx);		\
	tmpx = spu_mul(half, s0y[l2i]);		\
	tmpy = spu_mul(third, s1y[l2i]);	\
	tmpx = spu_add(tmpx, tmpy);		\
	tmpx = spu_mul(tmpx, s1z[l3i]);		\
	wx = spu_add(wx, tmpx);			\
						\
	wy = spu_mul(half, s1z[l3i]);		\
	wy = spu_add(s0z[l3i], wy);		\
	wy = spu_mul(s1y[l2i], wy);		\
						\
	wz = spu_mul(half, s1y[l2i]);		\
	wz = spu_add(s0y[l2i], wz);		\
	wz = spu_mul(s1z[l3i], wz);		\
						\
	wx = spu_mul(fnqx, wx);			\
	wy = spu_mul(fnqy, wy);			\
	jyh = spu_sub(jyh, wy);			\
	wz = spu_mul(fnqz, wz);			\
	jzh[l2i] = spu_sub(jzh[l2i], wz);	\
	}



      CALC_J_POINT;
      
      sjx0_a = spu_sel(zero, wx, (vector unsigned long long) element_assign[0]);
      sjy0_a = spu_sel(zero, jyh, (vector unsigned long long) element_assign[0]);
      sjz0_a = spu_sel(zero, jzh[l2i], (vector unsigned long long) element_assign[0]);
      
      sjx1_a = spu_shuffle(wx, zero, uphi_pat);
      sjy1_a = spu_shuffle(jyh, zero, uphi_pat);
      sjz1_a = spu_shuffle(jzh[l2i], zero, uphi_pat);
      
      l2i = 1;
      
      CALC_J_POINT; 
      
      sjx0_a = spu_shuffle(sjx0_a, wx, uplo_pat);
      sjy0_a = spu_shuffle(sjy0_a, jyh, uplo_pat);
      sjz0_a = spu_shuffle(sjz0_a, jzh[l2i], uplo_pat);
      
      sjx1_a = spu_sel(sjx1_a, wx, (vector unsigned long long) element_assign[1]);
      sjy1_a = spu_sel(sjy1_a, jyh, (vector unsigned long long) element_assign[1]);
      sjz1_a = spu_sel(sjz1_a, jzh[l2i], (vector unsigned long long) element_assign[1]);
      
      l2i = 2;
      
      CALC_J_POINT;
      
      sjx0_b = spu_sel(zero, wx, (vector unsigned long long) element_assign[0]);
      sjy0_b = spu_sel(zero, jyh, (vector unsigned long long) element_assign[0]);
      sjz0_b = spu_sel(zero, jzh[l2i], (vector unsigned long long) element_assign[0]);
      
      sjx1_b = spu_shuffle(wx, zero, uphi_pat);
      sjy1_b = spu_shuffle(jyh, zero, uphi_pat);
      sjz1_b = spu_shuffle(jzh[l2i], zero, uphi_pat);
      
      l2i = 3;
      
      CALC_J_POINT; 
      
      sjx0_b = spu_shuffle(sjx0_b, wx, uplo_pat);
      sjy0_b = spu_shuffle(sjy0_b, jyh, uplo_pat);
      sjz0_b = spu_shuffle(sjz0_b, jzh[l2i], uplo_pat);
      
      sjx1_b = spu_sel(sjx1_b, wx, (vector unsigned long long) element_assign[1]);
      sjy1_b = spu_sel(sjy1_b, jyh, (vector unsigned long long) element_assign[1]);
      sjz1_b = spu_sel(sjz1_b, jzh[l2i], (vector unsigned long long) element_assign[1]);
      
      
      l2i = 4;
      
      CALC_J_POINT;
      
      sjx0_c = spu_sel(zero, wx, (vector unsigned long long) element_assign[0]);
      sjy0_c = spu_sel(zero, jyh, (vector unsigned long long) element_assign[0]);
      sjz0_c = spu_sel(zero, jzh[l2i], (vector unsigned long long) element_assign[0]);
      
      sjx1_c = spu_shuffle(wx, zero, uphi_pat);
      sjy1_c = spu_shuffle(jyh, zero, uphi_pat);
      sjz1_c = spu_shuffle(jzh[l2i], zero, uphi_pat);


      if((store_off_0 & 1) == 0) {
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[0] - 2, j3_scal[0] - 2 ))) += sjx0_a;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[0] + VEC_SIZE - 2, j3_scal[0] - 2 ))) += sjx0_b;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[0] + 2 * VEC_SIZE - 2, j3_scal[0] - 2 ))) += sjx0_c;

	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[0] - 2, j3_scal[0] - 2 ))) += sjy0_a;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[0] + VEC_SIZE - 2, j3_scal[0] - 2 ))) += sjy0_b;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[0] + 2 * VEC_SIZE - 2, j3_scal[0] - 2 ))) += sjy0_c;

	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[0] - 2, j3_scal[0] - 2 ))) += sjz0_a;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[0] + VEC_SIZE - 2, j3_scal[0] - 2 ))) += sjz0_b;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[0] + 2 * VEC_SIZE - 2, j3_scal[0] - 2 ))) += sjz0_c;
    } else { // if we're unaligned
	v_real newx1, newx2, newx3;
	v_real newy1, newy2, newy3;
	v_real newz1, newz2, newz3;
	
	newx1 = spu_shuffle(zero, sjx0_a, uplo_pat);
	newx2 = spu_shuffle(sjx0_a, sjx0_b, fld_ip_pat[1][0]);
	newx3 = spu_shuffle(sjx0_b, sjx0_c, fld_ip_pat[1][0]);
	
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[0] - 2 - 1, j3_scal[0] - 2 ))) += newx1;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[0] + VEC_SIZE - 2 - 1, j3_scal[0] - 2 ))) += newx2;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[0] + 2 * VEC_SIZE - 2 - 1, j3_scal[0] - 2 ))) += newx3;

	newy1 = spu_shuffle(zero, sjy0_a, uplo_pat);
	newy2 = spu_shuffle(sjy0_a, sjy0_b, fld_ip_pat[1][0]);
	newy3 = spu_shuffle(sjy0_b, sjy0_c, fld_ip_pat[1][0]);
	
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[0] - 2 - 1, j3_scal[0] - 2 ))) += newy1;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[0] + VEC_SIZE - 2 - 1, j3_scal[0] - 2 ))) += newy2;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[0] + 2 * VEC_SIZE - 2 - 1, j3_scal[0] - 2 ))) += newy3;

	newz1 = spu_shuffle(zero, sjz0_a, uplo_pat);
	newz2 = spu_shuffle(sjz0_a, sjz0_b, fld_ip_pat[1][0]);
	newz3 = spu_shuffle(sjz0_b, sjz0_c, fld_ip_pat[1][0]);
	
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[0] - 2 - 1, j3_scal[0] - 2 ))) += newz1;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[0] + VEC_SIZE - 2 - 1, j3_scal[0] - 2 ))) += newz2;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[0] + 2 * VEC_SIZE - 2 - 1, j3_scal[0] - 2 ))) += newz3;
      }

    if((store_off_1 & 1) == 0) {
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[1] - 2, j3_scal[1] - 2 ))) += sjx1_a;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[1] + VEC_SIZE - 2, j3_scal[1] - 2 ))) += sjx1_b;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[1] + 2 * VEC_SIZE - 2, j3_scal[1] - 2 ))) += sjx1_c;

	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[1] - 2, j3_scal[1] - 2 ))) += sjy1_a;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[1] + VEC_SIZE - 2, j3_scal[1] - 2 ))) += sjy1_b;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[1] + 2 * VEC_SIZE - 2, j3_scal[1] - 2 ))) += sjy1_c;

	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[1] - 2, j3_scal[1] - 2 ))) += sjz1_a;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[1] + VEC_SIZE - 2, j3_scal[1] - 2 ))) += sjz1_b;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[1] + 2 * VEC_SIZE - 2, j3_scal[1] - 2 ))) += sjz1_c;
    }  else { // if we're unaligned
	v_real newx1, newx2, newx3;
	v_real newy1, newy2, newy3;
	v_real newz1, newz2, newz3;
	
	newx1 = spu_shuffle(zero, sjx1_a, uplo_pat);
	newx2 = spu_shuffle(sjx1_a, sjx1_b, fld_ip_pat[1][0]);
	newx3 = spu_shuffle(sjx1_b, sjx1_c, fld_ip_pat[1][0]);
	
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[1] - 2 - 1, j3_scal[1] - 2 ))) += newx1;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[1] + VEC_SIZE - 2 - 1, j3_scal[1] - 2 ))) += newx2;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JXI,0,j2_scal[1] + 2 * VEC_SIZE - 2 - 1, j3_scal[1] - 2 ))) += newx3;

	newy1 = spu_shuffle(zero, sjy1_a, uplo_pat);
	newy2 = spu_shuffle(sjy1_a, sjy1_b, fld_ip_pat[1][0]);
	newy3 = spu_shuffle(sjy1_b, sjy1_c, fld_ip_pat[1][0]);
	
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[1] - 2 - 1, j3_scal[1] - 2 ))) += newy1;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[1] + VEC_SIZE - 2 - 1, j3_scal[1] - 2 ))) += newy2;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JYI,0,j2_scal[1] + 2 * VEC_SIZE - 2 - 1, j3_scal[1] - 2 ))) += newy3;

	newz1 = spu_shuffle(zero, sjz1_a, uplo_pat);
	newz2 = spu_shuffle(sjz1_a, sjz1_b, fld_ip_pat[1][0]);
	newz3 = spu_shuffle(sjz1_b, sjz1_c, fld_ip_pat[1][0]);
	
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[1] - 2 - 1, j3_scal[1] - 2 ))) += newz1;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[1] + VEC_SIZE - 2 - 1, j3_scal[1] - 2 ))) += newz2;
	*((v_real *)(ls_fld + 
		     F2_SPU_OFF(ls_JZI,0,j2_scal[1] + 2 * VEC_SIZE - 2 - 1, j3_scal[1] - 2 ))) += newz3;
      }

      
#if 0		
      for(int l2i=l2min; l2i<=l2max; l2i++){
	CALC_J_POINT;
	int m = 0; 
	cbe_real wx_s = spu_extract(wx,m);
	ls_fld[F2_SPU_OFF(ls_JXI,0,j2_scal[m] + l2i - 2, j3_scal[m] + l3i - 2 )] += wx_s;
	cbe_real jyh_s = spu_extract(jyh,m);
	ls_fld[F2_SPU_OFF(ls_JYI,0,j2_scal[m] + l2i - 2, j3_scal[m] + l3i - 2 )] += jyh_s;
	cbe_real jzh_s = spu_extract(jzh[l2i],m);
	ls_fld[F2_SPU_OFF(ls_JZI,0,j2_scal[m] + l2i - 2, j3_scal[m] + l3i - 2 )] += jzh_s;
	m = 1; 
	wx_s = spu_extract(wx,m);
	ls_fld[F2_SPU_OFF(ls_JXI,0,j2_scal[m] + l2i - 2, j3_scal[m] + l3i - 2 )] += wx_s;
	jyh_s = spu_extract(jyh,m);
	ls_fld[F2_SPU_OFF(ls_JYI,0,j2_scal[m] + l2i - 2, j3_scal[m] + l3i - 2 )] += jyh_s;
	jzh_s = spu_extract(jzh[l2i],m);
	ls_fld[F2_SPU_OFF(ls_JZI,0,j2_scal[m] + l2i - 2, j3_scal[m] + l3i - 2 )] += jzh_s;
      }	
#endif

    }  
    
    n += 2; 
    
  } while(__builtin_expect((run),1));


  spu_dma_put(&ls_fld[ls_JXI*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,JXI,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_put(&ls_fld[ls_JYI*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,JYI,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_put(&ls_fld[ls_JZI*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F2_OFF_BLOCK(&psc_block,JZI,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      psc_block.ib[2])),
	      32*32*sizeof(fields_c_real_t));
  

  end_wait_particles_stored();
#if PRINT_DEBUG
  fprintf(stderr, "[[%#llx] ran %d particles\n", spu_ctx.spe_id, n);
#endif
  return 0;

}
