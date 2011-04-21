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

// Things are indexed a little differently on the
// local store, so we want our own enum for each
// compute kernel.
enum{
  ls_EX, ls_EY, ls_EZ,
  ls_HX, ls_HY, ls_HZ,
  ls_JXI, ls_JYI, ls_JZI,
  NR_LSFLDS,
};


// Some floating points constants which are needed
// in vector form. 
v_real half, one, two,  threefourths, onepfive,third,zero;

// A bunch of standard functions, set out here as static
// inlines for ease of reading the actual source. 
// These are essential slightly modified versions of the
// SSE2 functions.



#define find_index(dir, xi, dxi, j, h) {		\
    v_real tmp;					\
    tmp = spu_sub( xi, xb##dir);		\
    tmp = spu_mul( tmp, dxi);			\
    h = spu_round_real(tmp);			\
    j = spu_round_int(tmp);			\
    h = spu_sub(h, tmp);			\
} 


#define find_index_minus_shift(dir, xi,dxi,j,h,shift){	\
    v_real tmp;						\
    tmp = spu_sub( xi, xb##dir);			\
    tmp = spu_msub(tmp, dxi, shift);			\
    h = spu_round_real(tmp);				\
    j = spu_round_int(tmp);				\
    h = spu_sub(h, tmp);				\
} 


#define ip_to_grid_m(h,xmx){			\
    xmx = spu_add(half, h);			\
    xmx = spu_mul(xmx, xmx);			\
    xmx = spu_mul(half, xmx);			\
}

#define ip_to_grid_O(h, xOx) {			\
    xOx = spu_mul( h , h);			\
    xOx = spu_sub(threefourths, xOx);		\
}

#define ip_to_grid_l(h, xlx) {			\
    xlx = spu_sub(half, h);			\
    xlx = spu_mul( xlx, xlx);			\
    xlx = spu_mul(half, xlx);			\
}

#define form_factor_m(h, xmx) {			\
    xmx = spu_sub(h,one);			\
    xmx = spu_add(onepfive,xmx);					\
    xmx = spu_mul(xmx, xmx);						\
    xmx = spu_mul(half, xmx);						\
}


#define form_factor_O(h, xOx) {			\
    xOx = spu_mul(h, h);			\
    xOx = spu_sub(threefourths, xOx);		\
}

#define form_factor_l(h, xlx) {			\
    xlx = spu_add(one, h);			\
    xlx = spu_sub(onepfive, xlx);		\
    xlx = spu_mul(xlx, xlx);			\
    xlx = spu_mul(half, xlx);			\
}


// For dma purposes, we need C fields access macro here on the spu, but I don't
// want to contaminate the spu code, so it's just copied here. 
#define F3_OFF_C(pf, fldnr, jx,jy,jz)					\
  ((((((fldnr)								\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))



// FIXME Assumes xy plane!!
// We need a macro to access the offsets in the local store. 
// The fields in the local store will only be 2d, and right now I'm going to 
// force them to be 2d in the xy plane, life is too messy otherwise.
#define F2_SPU_OFF(fldnr, jx, jy)					\
  (((((fldnr)								\
      * psc_block.im[1] + ((jy)-psc_block.ib[1]))			\
     * psc_block.im[0] + ((jx)-psc_block.ib[0]))))


// A very ugly macro. This does the field interpolation. As far as I can tell,
// it does it in a very effecient manner. Could be wrong about that to.
// This macro could use a careful look in asmvis. If precision switching
// is going to work, we're going to need a single precision version of this too.

#define IP_FIELD_SPU(fldnr, jx, jy, outer_coeff, inner_coeff, field) {	\
  v_real _tmp1, _tmp2, _tmp3;						\
  v_real _in_A[2];							\
  v_real _in_B[2];							\
  int jx_0 = (int)spu_extract(jx,0);					\
  int jx_1 = (int)spu_extract(jx,1);					\
  int jy_0 = (int)spu_extract(jy,0);					\
  int jy_1 = (int)spu_extract(jy,1);					\
									\
  int off_0 = F2_SPU_OFF(fldnr, jx_0 - 1, jy_0 - 1);			\
  int off_1 = F2_SPU_OFF(fldnr, jx_1 - 1, jy_1 - 1);			\
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
  _tmp1 = spu_mul(inner_coeff##mx,_tmp1);				\
  _tmp2 = spu_mul(inner_coeff##Ox,_tmp2);				\
  _tmp3 = spu_mul(inner_coeff##lx,_tmp3);				\
  _tmp1 = spu_add(_tmp1, _tmp2);					\
  _tmp1 = spu_add(_tmp1, _tmp3);					\
  field = spu_mul(outer_coeff##my, _tmp1);				\
									\
  off_0 = F2_SPU_OFF(fldnr, jx_0-1, jy_0);				\
  off_1 = F2_SPU_OFF(fldnr, jx_1-1, jy_1);				\
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
  _tmp1 = spu_mul(inner_coeff##mx,_tmp1);				\
  _tmp2 = spu_mul(inner_coeff##Ox,_tmp2);				\
  _tmp3 = spu_mul(inner_coeff##lx,_tmp3);				\
  _tmp1 = spu_add(_tmp1, _tmp2);					\
  _tmp1 = spu_add(_tmp1, _tmp3);					\
  field = spu_madd(outer_coeff##Oy, _tmp1, field);			\
									\
  off_0 = F2_SPU_OFF(fldnr, jx_0 - 1, jy_0+1);				\
  off_1 = F2_SPU_OFF(fldnr, jx_1 - 1, jy_1+1);				\
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
  _tmp1 = spu_mul(inner_coeff##mx,_tmp1);				\
  _tmp2 = spu_mul(inner_coeff##Ox,_tmp2);				\
  _tmp3 = spu_mul(inner_coeff##lx, _tmp3);				\
  _tmp1 = spu_add(_tmp1, _tmp2);					\
  _tmp1 = spu_add(_tmp1, _tmp3);					\
  field = spu_madd(outer_coeff##ly, _tmp1, field);			\
  }

// The actual compute kernel function.
int
spu_push_part_2d(void){


#if PRINT_DEBUG

  printf("[[%#llx] start ea: %#llx end ea: %#llx \n", spu_ctx.spe_id, psc_block.part_start, psc_block.part_end);

#endif

  // Local store for the fields. At this point, contains E,H, and J.
  fields_c_real_t ls_fld[NR_LSFLDS*32*32] __attribute__((aligned(128)));

  // If we try to dma in all the field data at once, we run into the 
  // mfc 16KB per request limit. For now, we'll split it into 6 different
  // requests. Depending on whether the limit is 16KB = 1000 * 16B 
  // or 1024*16B, we could group these two by two.
  
  // There might be a faster way to do this using dma lists. One would
  // think we're stalling pretty bad here, but until I get a good idea 
  // of the bandwidth and latency of the mfc I can't really say for sure.

  // This could be pretty easily generalized for precision switching, but 
  // right now it won't work with single precision. 

  // FIXME: This assumes the index in the symmetry direction is 0.

  spu_dma_get(ls_fld, 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,EX,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_EY*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,EY,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_EZ*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,EZ,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_HX*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,HX,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_HY*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,HY,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_get(&ls_fld[ls_HZ*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,HZ,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));


  // We need to make sure the currents are zeroed out before we start work.
  memset(&ls_fld[ls_JXI*32*32], 0, 3*32*32*sizeof(fields_c_real_t));

  // These two pointers reference the current particle being 
  // worked on by the spu (cp_ea) and the next particle which needs
  // to be preloaded (np_ea). The addresses are in main memory.
  unsigned long long cp_ea = psc_block.part_start;
  unsigned long long np_ea; 
  
  // We're triple buffering the particles here: One preloading, 
  // one being worked on, and one storing. (when I say one, I actually 
  // mean one vector set of particles. Two for double, four for single.)
  //
  // Kai feels we might ineffeciently using the mfc by only loading one set
  // of particles at a time, feeling instead that each buffer should be some
  // larger group of particles to cut down on the latency. 
  //
  // He may be correct, however I am of the opinion that the code is computationally
  // bound, so we're not really going to be speeding things up, and life will just become
  // more difficult. It's really going to take more testing, because I really can't
  // say for now. For now let's assume Kai is correct (the more likely case) and I need
  // to work on the buffering system a bit. 
  particle_cbe_t _bufferA[2] __attribute__((aligned(16))), 
    _bufferB[2] __attribute__((aligned(16))), 
    _bufferC[2] __attribute__((aligned(16)));

  buff.plb1 = &(_bufferA[0]);
  buff.plb2 = &(_bufferA[1]);
  buff.lb1 = &(_bufferB[0]); 
  buff.lb2 = &(_bufferB[1]);
  buff.sb1 = &(_bufferC[0]);
  buff.sb2 = &(_bufferC[1]);

  // Get the first two particles, we need to preload them
  // before we actually start work on the loop.
  // Doing it requires a slightly different function call than 
  // we'll use inside the loop. I'm not really happy with having 
  // this off-loaded to a function in spu_dma.c, but this is a pretty 
  // small complaint. 
  first_preload_particle(buff.plb1, cp_ea, 2*sizeof(particle_cbe_t));

  np_ea = cp_ea + 2*sizeof(particle_cbe_t);

  // we have some stuff to do while we wait for
  // it to come in.
  
  // insert assignment, promotions, and constant 
  // calculations here.
  // When we change precision, this area will need to be modified. 
  v_real dt, yl, xl, dyi, dxi, dqs,fnqs,fnqxs,fnqys,fnqzs;
  v_real xbx, xby, xbz;
  dt = spu_splats(spu_ctx.dt);
  half = spu_splats(0.5);
  xl = spu_mul(half, dt);
  yl = spu_mul(half, dt);
  dxi = spu_splats(1./spu_ctx.dx[0]);
  dyi = spu_splats(1./spu_ctx.dx[1]);
  one = spu_splats(1.0);
  two = spu_splats(2.0);
  threefourths = spu_splats(0.75);
  onepfive = spu_splats(1.5);
  third = spu_splats(1./3.);
  zero = spu_splats(0.0);
  dqs = spu_splats(0.5 * spu_ctx.eta * spu_ctx.dt);
  fnqs = spu_splats(spu_ctx.fnqs);
  fnqxs = spu_splats(spu_ctx.dx[0] * spu_ctx.fnqs / spu_ctx.dt);
  fnqys = spu_splats(spu_ctx.dx[1] * spu_ctx.fnqs / spu_ctx.dt);
  fnqzs = spu_splats(spu_ctx.dx[2] * spu_ctx.fnqs / spu_ctx.dt);
  xbx = spu_splats(psc_block.xb[0]);
  xby = spu_splats(psc_block.xb[1]);
  xbz = spu_splats(psc_block.xb[2]);
  unsigned long long end = psc_block.part_end; 

  int run = 1;
  int n = 0;

  // We might need an empty particle to pad out the last 
  // time through the loop if we have an odd number of particles. 
  particle_cbe_t null_part; 
    

  // Time to actually do the loop. This is a while loop, instead of
  // a for. It will terminate when the variable 'run', defined above, 
  // is set to '0'.
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
      // Final branch: If this is the last particle, 
      // and we don't need to add any padding. 
      wait_for_preload();
    }


    // The actual computational loop. This is pretty straight
    // forward, but I need to test if using the atomic intrinsics
    // really helps. It's really a question as to whether I'm better at 
    // writing assembly than the compiler.
    // Will test soon.

    v_real xi, yi, zi, pxi, pyi, pzi, qni, mni, wni;
  
    // The particle loading/storing macro was ported directly
    // from the SSE2 implementation. We should really take a look
    // at the assembly and see if it could be optimized more for 
    // the cell.
    LOAD_PARTICLES_SPU;

    v_real vxi,vyi,vzi,root,tmpx,tmpy,tmpz; 
    
    // Regular part a, as expected.
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
    
    tmpx = spu_mul(vxi, xl);
    tmpy = spu_mul(vyi, yl);
    
    xi = spu_add(xi, tmpx);
    yi = spu_add(yi, tmpy);

    // A little part b loving. 

    v_real gmx, gmy, gOx, gOy, glx, gly, H1, H2, h1, h2;
    v_int j1, j2, l1, l2;
      
    find_index(x,xi, dxi, j1, H1);
    find_index(y,yi, dyi, j2, H2);
    
    ip_to_grid_m(H1, gmx);
    ip_to_grid_m(H2, gmy);
    
    ip_to_grid_O(H1, gOx);
    ip_to_grid_O(H2, gOy);
    
    ip_to_grid_l(H1, glx);
    ip_to_grid_l(H2, gly);
    
    find_index_minus_shift(x,xi, dxi, l1, h1, half);
    find_index_minus_shift(y,yi, dyi, l2, h2, half);
    
    v_real hmx, hmy, hOx, hOy, hlx, hly;

    ip_to_grid_m(h1, hmx);
    ip_to_grid_m(h2, hmy);
    
    ip_to_grid_O(h1, hOx);
    ip_to_grid_O(h2, hOy);
    
    ip_to_grid_l(h1, hlx);
    ip_to_grid_l(h2, hly);
    
    // Field interpolation here. urg...
    
    v_real exq, eyq, ezq, hxq, hyq, hzq;
    
    IP_FIELD_SPU(ls_EX,l1,j2,g,h,exq);
    IP_FIELD_SPU(ls_EY,j1,l2,h,g,eyq);
    IP_FIELD_SPU(ls_EZ,j1,j2,g,g,ezq);
    IP_FIELD_SPU(ls_HX,j1,l2,h,g,hxq);
    IP_FIELD_SPU(ls_HY,l1,j2,g,h,hyq);
    IP_FIELD_SPU(ls_HZ,l1,l2,h,h,hzq);
    
    
    // These arrays may be sort of a problem.
    // It all sorts of depends on how smart the compiler is. 
    // Allow me to explain in some depth. Arrays are groups
    // of sequential memory, and indexing the array involves 
    // adding the index to the address of the first element. 
    // The problem is, because registers don't have an address, 
    // you can't really create an array of registers. So, these variables
    // are probably always on the stack, and any operations we do are 
    // using them off the stack, which is slower than just keeping them 
    // in registers. The compiler, however, may be smart enough to unroll 
    // the loops and use them as registers. Again, I'll have to look
    // at the assembly and figure out exactly what's happening. At the 
    // very least, we could probably manually unroll the loop for the
    // s0* and make sure they stay in registers. Unfortunately, I'm not really
    // sure it will be possible to do that with the s1*.
    v_real s0x[5], s0y[5], s1x[5], s1y[5];
    for(int mp=0; mp<5; mp++){
      s0x[mp] = spu_splats(0.0);
      s1x[mp] = spu_splats(0.0);
      s0y[mp] = spu_splats(0.0);
      s1y[mp] = spu_splats(0.0);
    }
    form_factor_m(H1, s0x[1]);
    form_factor_O(H1, s0x[2]);
    form_factor_l(H1, s0x[3]);
    
    form_factor_m(H2, s0y[1]);
    form_factor_O(H2, s0y[2]);
    form_factor_l(H2, s0y[3]);
    
    // FIXME: I'm missing a step here. I need to look into it. 
    // should just be some sort of reduction. 
    
    // It would be nice to offload the momentum to a function, but it may not
    // be necessary
    
    v_real dq, pxp, pyp, pzp, taux, tauy, tauz, tau;
    
    // My version is faster than the fast-math 
    // version ( ie mine = .86 * fast-math)
    // so even though this is harder to read, 
    // I'm going to use it instead.

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
    
    v_real txx, tyy, tzz, t2xy, t2xz, t2yz;
    
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
    
    tmpx = spu_mul(vxi, xl);
    tmpy = spu_mul(vyi, yl);
    
    xi = spu_add(xi, tmpx);
    yi = spu_add(yi, tmpy);

    STORE_PARTICLES_SPU;

    np_ea = cp_ea +  2 * sizeof(particle_cbe_t);
    // At this point, np_ea is one load ahead of the
    // current particle.


    // Need to issue the particle dma out requests. 
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

    // The final piece, current calculations. 
    xi = spu_add(xi, tmpx);
    yi = spu_add(yi, tmpy);

    v_int k1,k2;
    
    find_index(x,xi, dxi, k1, h1);
    find_index(y,yi, dyi, k2, h2);
    
    // God help me, there's some things I just can't fgure out how to parallelize
    // The g-- here are just temporary variables. I can't do the assignments in parallel, 
    // but I'll be damned if I can't do the FLOPS in parallel
    
    form_factor_m(h1, gmx);
    form_factor_O(h1, gOx);
    form_factor_l(h1, glx);
    
    form_factor_m(h2, gmy);
    form_factor_O(h2, gOy);
    form_factor_l(h2, gly);
    

    // This next part is kind of a pain. The two particles can be in 
    // different cells and move different directions, so we need to 
    // be able to get the cell index for each one. This prompts these
    // spu_extracts into scalar variables. Extracts are slow, as each one 
    // is a combination of some number of rotates and what not. Sadly, 
    // I can't really think of any alteranative. 
    signed long long j1_scal[2], j2_scal[2];
    j1_scal[0] = spu_extract(j1, 0);
    j1_scal[1] = spu_extract(j1, 1);
    j2_scal[0] = spu_extract(j2, 0);
    j2_scal[1] = spu_extract(j2, 1);
    
    signed long long dfx[2];
    signed long long  dfy[2];
    dfx[0] = spu_extract(k1,0) - j1_scal[0];
    dfx[1] = spu_extract(k1,1) - j1_scal[1];
    dfy[0] = spu_extract(k2,0) - j2_scal[0];
    dfy[1] = spu_extract(k2,1) - j2_scal[1];
    
    
    // This is the main section which is preventing me from 
    // eliminating the arrays. It may be a major bottleneck, 
    // depending on how the compiler has implemented it. 
    for(int p=0; p < VEC_SIZE; p++){
      s1x[(int)dfx[p] + 1] = spu_sel(s1x[(int)dfx[p] + 1], gmx, (vector unsigned long long) element_assign[p]);
      s1x[(int)dfx[p] + 2] = spu_sel(s1x[(int)dfx[p] + 2], gOx, (vector unsigned long long) element_assign[p]);
      s1x[(int)dfx[p] + 3] = spu_sel(s1x[(int)dfx[p] + 3], glx, (vector unsigned long long) element_assign[p]);
      s1y[(int)dfy[p] + 1] = spu_sel(s1y[(int)dfy[p] + 1], gmy, (vector unsigned long long) element_assign[p]);
      s1y[(int)dfy[p] + 2] = spu_sel(s1y[(int)dfy[p] + 2], gOy, (vector unsigned long long) element_assign[p]);
      s1y[(int)dfy[p] + 3] = spu_sel(s1y[(int)dfy[p] + 3], gly, (vector unsigned long long) element_assign[p]);
    }
    
    for(int m=0; m<5; m++){
      s1x[m] = spu_sub(s1x[m], s0x[m]);
      s1y[m] = spu_sub(s1y[m], s0y[m]);
    }
    
    
    // This section appears to effectively eliminate
    // the branching section in the original code. 
    // The l2min/max part is not, at this moment, needed
    // as I have unrolled the inner loop in the current assign.    
    int l1min = 1 + ((dfx[0] | dfx[1]) >> 1),
      l1max = 3 + (((dfx[0]>>1)^dfx[0]) | ((dfx[1]>>1)^dfx[1])),
      l2min = 1 + ((dfy[0] | dfy[1]) >> 1),
      l2max = 3 + (((dfy[0]>>1)^dfy[0]) | ((dfy[1]>>1)^dfy[1]));
    
    v_real fnqx, fnqy, fnqz;
    
    fnqx = spu_mul(wni, fnqxs);
    fnqx = spu_mul(qni, fnqx);
    
    fnqy = spu_mul(wni, fnqys);
    fnqy = spu_mul(qni, fnqy);
    
    fnqz = spu_mul(wni, fnqs);
    fnqz = spu_mul(qni, fnqz);
    fnqz = spu_mul(vzi, fnqz);

    
    // This is an ugly section, and it 
    // doesn't translate well to single 
    // precision. The += to the local 
    // store are a very bad idea. They stall
    // the processor rather significantly. As
    // such, we need to do some acrobatics. 
    // First, we're going to shuffle each current element
    // for each particle into some registers lined up on the
    // fast running index. 
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

    v_real jyh[5]; // As per Will's suggestion, using the minimal 
    // variable here, and for jyh below, gives a nice
    // performance boost
    memset(jyh,0,5*sizeof(v_real));


    for(int l2i=l2min; l2i<=l2max; l2i++){
      v_real jxh = spu_splats(0.0);
	
      // Let's find the offsets in the local store to which 
      // we will store the new currents. 
      long int store_off_0 = F2_SPU_OFF(ls_JXI,j1_scal[0] - 2, j2_scal[0] + l2i - 2);
      long int store_off_1 = F2_SPU_OFF(ls_JXI,j1_scal[1] - 2, j2_scal[1] + l2i - 2);


      v_real wx, wy, wz;
      
#define CALC_J_X {				\
	wx = spu_mul(half, s1y[l2i]);		\
	wx = spu_add(s0y[l2i], wx);		\
	wx = spu_mul(s1x[l1i], wx);		\
	wx = spu_mul(fnqx, wx);			\
	jxh = spu_sub(jxh, wx);			\
      }
      
#define CALC_J_Y {				\
	wy = spu_mul(half, s1x[l1i]);		\
	wy = spu_add(s0x[l1i], wy);		\
	wy = spu_mul(s1y[l2i], wy);		\
	wy = spu_mul(fnqy, wy);			\
	jyh[l1i] = spu_sub(jyh[l1i],wy);	\
      }
      
#define CALC_J_Z {					\
	wz = spu_mul(half,s1x[l1i]);			\
	wz = spu_add(s0x[l1i], wz);			\
	wz = spu_mul(s0y[l2i], wz);			\
	tmpx = spu_mul(half, s0x[l1i]);			\
	tmpy = spu_mul(third, s1x[l1i]);		\
	tmpx = spu_add(tmpx, tmpy);			\
	tmpx = spu_mul(tmpx, s1y[l2i]);			\
	wz = spu_add(wz, tmpx);				\
	wz = spu_mul(fnqz, wz);				\
      }

      // Now we need to preload in the currents in
      // local store that will be affected by each particle. 
      v_real l0x_a, l0x_b, l0x_c;
      v_real l0y_a, l0y_b, l0y_c;
      v_real l0z_a, l0z_b, l0z_c;
      v_real l1x_a, l1x_b, l1x_c;
      v_real l1y_a, l1y_b, l1y_c;
      v_real l1z_a, l1z_b, l1z_c;


      // need to load these up here and hope the compiler
      // is smart enough to realize we're preloading for when
      // we get to the bottom.
      
      l0x_a = *((v_real *)(ls_fld + (store_off_0 & ~1)));
      
      l0x_b = *((v_real *)(ls_fld + (store_off_0 & ~1) + 2));

      l0x_c = *((v_real *)(ls_fld + (store_off_0 & ~1) + 4));

      l0y_a = *((v_real *)(ls_fld + ((store_off_0 + 32*32) & ~1)));
      
      l0y_b = *((v_real *)(ls_fld + ((store_off_0 + 32*32) & ~1) + 2));

      l0y_c = *((v_real *)(ls_fld + ((store_off_0 + 32*32) & ~1) + 4));

      l0z_a = *((v_real *)(ls_fld + ((store_off_0 + 2*32*32) & ~1)));
      
      l0z_b = *((v_real *)(ls_fld + ((store_off_0 + 2*32*32) & ~1) + 2));

      l0z_c = *((v_real *)(ls_fld + ((store_off_0 + 2*32*32) & ~1) + 4));

      
      // now for second particle

      l1x_a = *((v_real *)(ls_fld + (store_off_1 & ~1)));
      
      l1x_b = *((v_real *)(ls_fld + (store_off_1 & ~1) + 2));

      l1x_c = *((v_real *)(ls_fld + (store_off_1 & ~1) + 4));

      l1y_a = *((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1)));
      
      l1y_b = *((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 2));

      l1y_c = *((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 4));

      l1z_a = *((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1)));
      
      l1z_b = *((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 2));

      l1z_c = *((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 4));

      

#if 1
      // Now for the unrolled loop. We calculate the 
      // value of each particles current at each point, 
      // and then shuffle into the appropriate place in the 
      // registers we allocated before. Note, we use spu_sel whenever
      // possible. A shuffle is 4 cycles, a select (ie pass through)
      // is only 2. 
      int l1i = 0;
      CALC_J_X;
      CALC_J_Y;
      CALC_J_Z;

      sjx0_a = spu_sel(zero, jxh, (vector unsigned long long) element_assign[0]);
      sjy0_a = spu_sel(zero, jyh[l1i], (vector unsigned long long) element_assign[0]);
      sjz0_a = spu_sel(zero, wz, (vector unsigned long long) element_assign[0]);
      
      sjx1_a = spu_shuffle(jxh, zero, uphi_pat);
      sjy1_a = spu_shuffle(jyh[l1i], zero, uphi_pat);
      sjz1_a = spu_shuffle(wz, zero, uphi_pat);
      
      l1i = 1;
      
      CALC_J_X;
      CALC_J_Y;
      CALC_J_Z;
      
      sjx0_a = spu_shuffle(sjx0_a, jxh, uplo_pat);
      sjy0_a = spu_shuffle(sjy0_a, jyh[l1i], uplo_pat);
      sjz0_a = spu_shuffle(sjz0_a, wz, uplo_pat);
      
      sjx1_a = spu_sel(sjx1_a, jxh, (vector unsigned long long) element_assign[1]);
      sjy1_a = spu_sel(sjy1_a, jyh[l1i], (vector unsigned long long) element_assign[1]);
      sjz1_a = spu_sel(sjz1_a, wz, (vector unsigned long long) element_assign[1]);
      
      l1i = 2;
      
      CALC_J_X;
      CALC_J_Y;
      CALC_J_Z;
            
      sjx0_b = spu_sel(zero, jxh, (vector unsigned long long) element_assign[0]);
      sjy0_b = spu_sel(zero, jyh[l1i], (vector unsigned long long) element_assign[0]);
      sjz0_b = spu_sel(zero, wz, (vector unsigned long long) element_assign[0]);
      
      sjx1_b = spu_shuffle(jxh, zero, uphi_pat);
      sjy1_b = spu_shuffle(jyh[l1i], zero, uphi_pat);
      sjz1_b = spu_shuffle(wz, zero, uphi_pat);
      
      l1i = 3;
      
      CALC_J_X;
      CALC_J_Y;
      CALC_J_Z;
      
      sjx0_b = spu_shuffle(sjx0_b, jxh, uplo_pat);
      sjy0_b = spu_shuffle(sjy0_b, jyh[l1i], uplo_pat);
      sjz0_b = spu_shuffle(sjz0_b, wz, uplo_pat);
      
      sjx1_b = spu_sel(sjx1_b, jxh, (vector unsigned long long) element_assign[1]);
      sjy1_b = spu_sel(sjy1_b, jyh[l1i], (vector unsigned long long) element_assign[1]);
      sjz1_b = spu_sel(sjz1_b, wz, (vector unsigned long long) element_assign[1]);
      
      
      l1i = 4;
      
      CALC_J_X;
      CALC_J_Y;
      CALC_J_Z;
      
      sjx0_c = spu_sel(zero, jxh, (vector unsigned long long) element_assign[0]);
      sjy0_c = spu_sel(zero, jyh[l1i], (vector unsigned long long) element_assign[0]);
      sjz0_c = spu_sel(zero, wz, (vector unsigned long long) element_assign[0]);
      
      sjx1_c = spu_shuffle(jxh, zero, uphi_pat);
      sjy1_c = spu_shuffle(jyh[l1i], zero, uphi_pat);
      sjz1_c = spu_shuffle(wz, zero, uphi_pat);


      // I call this pain, and it is. Basically, we need 
      // to account for the currents from the two particles
      // overlapping each other. There are six cases we need to consider,
      // and this variable will let us differentiate them. This 
      // variable is calculated up here so the branches can be predicted. 
      // I'm not sure if I specifically have to issue some instruction
      // for the spe to do the branch prediction. 
      int pain = store_off_1/2 - store_off_0/2;

      // First, if particle is unaligned, we need to 
      // shift the new currents one to the right in our 
      // store registers. This is not easy, as the currents
      // to be stored span three registers. With a smart use of shuffles, 
      // however, it's not too expensive. 
      if((store_off_0 & 1) == 0) { // we're aligned.

	l0x_a += sjx0_a;
	l0x_b += sjx0_b;
	l0x_c += sjx0_c;

	l0y_a += sjy0_a;
	l0y_b += sjy0_b;
	l0y_c += sjy0_c;

	l0z_a += sjz0_a;
	l0z_b += sjz0_b;
	l0z_c += sjz0_c;
	
      } 
    
      else { // if we're unaligned
	v_real newx1, newx2, newx3;
	v_real newy1, newy2, newy3;
	v_real newz1, newz2, newz3;
	
	newx1 = spu_shuffle(zero, sjx0_a, uplo_pat);
	newx2 = spu_shuffle(sjx0_a, sjx0_b, fld_ip_pat[1][0]);
	newx3 = spu_shuffle(sjx0_b, sjx0_c, fld_ip_pat[1][0]);
	
	l0x_a += newx1;
	l0x_b += newx2;
	l0x_c += newx3;

	newy1 = spu_shuffle(zero, sjy0_a, uplo_pat);
	newy2 = spu_shuffle(sjy0_a, sjy0_b, fld_ip_pat[1][0]);
	newy3 = spu_shuffle(sjy0_b, sjy0_c, fld_ip_pat[1][0]);
	
	l0y_a += newy1;
	l0y_b += newy2;
	l0y_c += newy3;

	newz1 = spu_shuffle(zero, sjz0_a, uplo_pat);
	newz2 = spu_shuffle(sjz0_a, sjz0_b, fld_ip_pat[1][0]);
	newz3 = spu_shuffle(sjz0_b, sjz0_c, fld_ip_pat[1][0]);
	
	l0z_a += newz1;
	l0z_b += newz2;
	l0z_c += newz3;
      }
      

      // If the second particle is also unaligned, we need
      // to do the same shifting buisness. 
      if((store_off_1 & 1) == 1) { // if we're unaligned
	v_real newx1, newx2, newx3;
	v_real newy1, newy2, newy3;
	v_real newz1, newz2, newz3;
	
	newx1 = spu_shuffle(zero, sjx1_a, uplo_pat);
	newx2 = spu_shuffle(sjx1_a, sjx1_b, fld_ip_pat[1][0]);
	newx3 = spu_shuffle(sjx1_b, sjx1_c, fld_ip_pat[1][0]);
	
	sjx1_a = newx1;
	sjx1_b = newx2;
	sjx1_c = newx3;

	newy1 = spu_shuffle(zero, sjy1_a, uplo_pat);
	newy2 = spu_shuffle(sjy1_a, sjy1_b, fld_ip_pat[1][0]);
	newy3 = spu_shuffle(sjy1_b, sjy1_c, fld_ip_pat[1][0]);
	
	sjy1_a = newy1;
	sjy1_b = newy2;
	sjy1_c = newy3;

	newz1 = spu_shuffle(zero, sjz1_a, uplo_pat);
	newz2 = spu_shuffle(sjz1_a, sjz1_b, fld_ip_pat[1][0]);
	newz3 = spu_shuffle(sjz1_b, sjz1_c, fld_ip_pat[1][0]);
	
	sjz1_a = newz1;
	sjz1_b = newz2;
	sjz1_c = newz3;
      }
      
      // Now we need to account for the six
      // cases of relative position between the two 
      // particles. 

      if( (pain <= -3) || (pain >= 3)) {
	// There is no overlap between the particles, 
	// so we can store them seperately. 
	l1x_a += sjx1_a;
	l1x_b += sjx1_b;
	l1x_c += sjx1_c;

	l1y_a += sjy1_a;
	l1y_b += sjy1_b;
	l1y_c += sjy1_c;

	l1z_a += sjz1_a;
	l1z_b += sjz1_b;
	l1z_c += sjz1_c;

	*((v_real *)(ls_fld + (store_off_1 & ~1))) = l1x_a;
	
	*((v_real *)(ls_fld + (store_off_1 & ~1) + 2)) = l1x_b;
	
	*((v_real *)(ls_fld + (store_off_1 & ~1) + 4)) = l1x_c;
	
	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1))) = l1y_a;
	
	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 2)) = l1y_b;
	
	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 4)) = l1y_c;
	
	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1))) = l1z_a;
	
	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 2)) = l1z_b;
	
	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 4)) = l1z_c;
	
	
      }
      else if (pain == -2) {
	// The particles overlap at one point, as shown below. 
	// | 0a | 0b | 0c |    |
	// |    |    | 1a | 1b | 1c |
	l1x_a += sjx1_a;
	l1x_b += sjx1_b;

	l1y_a += sjy1_a;
	l1y_b += sjy1_b;

	l1z_a += sjz1_a;
	l1z_b += sjz1_b;
	
	*((v_real *)(ls_fld + (store_off_1 & ~1))) = l1x_a;
	
	*((v_real *)(ls_fld + (store_off_1 & ~1) + 2)) = l1x_b;

	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1))) = l1y_a;
	
	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 2)) = l1y_b;

	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1))) = l1z_a;
	
	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 2)) = l1z_b;

	l0x_a += sjx1_c;
	l0y_a += sjy1_c;
	l0z_a += sjz1_c;

      } 
      else if (pain == 2) {
	// Overlap at one point:
	// |    |    | 0a | 0b | 0c | 
	// | 1a | 1b | 1c |    |    |
	l1x_b += sjx1_b;
	l1x_c += sjx1_c;

	l1y_b += sjy1_b;
	l1y_c += sjy1_c;

	l1z_b += sjz1_b;
	l1z_c += sjz1_c;
	
	*((v_real *)(ls_fld + (store_off_1 & ~1) + 2)) = l1x_b;
	
	*((v_real *)(ls_fld + (store_off_1 & ~1) + 4)) = l1x_c;

	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 2)) = l1y_b;
	
	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 4)) = l1y_c;

	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 2)) = l1z_b;
	
	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 4)) = l1z_c;

	l0x_c += sjx1_a;
	l0y_c += sjy1_a;
	l0z_c += sjz1_a;

      } 
      else if (pain == -1) {
	// Overlap at two points:
	// | 0a | 0b | 0c |    |
	// |    | 1a | 1b | 1c |

	l1x_a += sjx1_a;

	l1y_a += sjy1_a;

	l1z_a += sjz1_a;

	*((v_real *)(ls_fld + (store_off_1 & ~1))) = l1x_a;

	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1))) = l1y_a;
      
	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1))) = l1z_a;

	l0x_a += sjx1_b;
	l0x_b += sjx1_c;
	l0y_a += sjy1_b;
	l0y_b += sjy1_c;
	l0z_a += sjz1_b;
	l0z_b += sjz1_c;

      } 
      else if (pain == 1) {
	// Overlap at two points:
	// |    | 0a | 0b | 0c |
	// | 1a | 1b | 1c |    |
	l1x_c += sjx1_c;

	l1y_c += sjy1_c;

	l1z_c += sjz1_c;

	*((v_real *)(ls_fld + (store_off_1 & ~1) + 4 )) = l1x_c;

	*((v_real *)(ls_fld + ((store_off_1 + 32*32) & ~1) + 4 )) = l1y_c;
      
	*((v_real *)(ls_fld + ((store_off_1 + 2*32*32) & ~1) + 4)) = l1z_c;

	l0x_b += sjx1_a;
	l0x_c += sjx1_b;
	l0y_b += sjy1_a;
	l0y_c += sjy1_b;
	l0z_b += sjz1_a;
	l0z_c += sjz1_b;
	
      }
      else if (pain == 0) {
	// Overlap at all three points
	// | 0a | 0b | 0c |
	// | 1a | 1b | 1c |
	l0x_a += sjx1_a;
	l0x_b += sjx1_b;
	l0x_c += sjx1_c;
	l0y_a += sjy1_a;
	l0y_b += sjy1_b;
	l0y_c += sjy1_c;
	l0z_a += sjz1_a;
	l0z_b += sjz1_b;
	l0z_c += sjz1_c;
      }

      // Time to store and hope for the best.


      *((v_real *)(ls_fld + (store_off_0 & ~1))) = l0x_a ;
      
      *((v_real *)(ls_fld + (store_off_0 & ~1) + 2)) = l0x_b; 

      *((v_real *)(ls_fld + (store_off_0 & ~1) + 4)) = l0x_c;

      *((v_real *)(ls_fld + ((store_off_0 + 32*32) & ~1))) = l0y_a;
      
      *((v_real *)(ls_fld + ((store_off_0 + 32*32) & ~1) + 2)) = l0y_b;

      *((v_real *)(ls_fld + ((store_off_0 + 32*32) & ~1) + 4)) = l0y_c;

      *((v_real *)(ls_fld + ((store_off_0 + 2*32*32) & ~1))) = l0z_a;
      
      *((v_real *)(ls_fld + ((store_off_0 + 2*32*32) & ~1) + 2)) = l0z_b;

      *((v_real *)(ls_fld + ((store_off_0 + 2*32*32) & ~1) + 4)) = l0z_c;




#endif

    }  
    
    // This just counts how many particles we've run.
    // It's only here for debugging purposes. 
    n += 2; 
    
  } while(__builtin_expect((run),1));


  // Now it's time to write out the currents

  // FIXME: This assumes the index of the plane in 0 in the symmetry direction.
  spu_dma_put(&ls_fld[ls_JXI*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,JXI,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_put(&ls_fld[ls_JYI*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,JYI,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));

  spu_dma_put(&ls_fld[ls_JZI*32*32], 
	      (psc_block.wb_flds + 
	       sizeof(fields_c_real_t) * F3_OFF_C(&psc_block,JZI,
						      psc_block.ib[0],
						      psc_block.ib[1],
						      0)),
	      32*32*sizeof(fields_c_real_t));
  

  // Then we just wait to make sure the particles have 
  // finished storing. 
  end_wait_particles_stored();
#if PRINT_DEBUG
  fprintf(stderr, "[[%#llx] ran %d particles\n", spu_ctx.spe_id, n);
#endif
  return 0;

}
