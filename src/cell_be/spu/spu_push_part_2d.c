#include "spu_particles.h"
#include "psc_spu.h"
#include "psc_spu_2d.h"
#include <stddef.h>
#include <stdio.h>

#include <simdmath.h>

#ifdef __SPU__
#include <spu_mfcio.h>
#else
#include "spu_mfcio_c.h"
#endif

enum{
  ls_EX, ls_EY, ls_EZ,
  ls_HX, ls_HY, ls_HZ,
  NR_LSFLDS,
};


int
spu_push_part_2d(void){


#if PRINT_DEBUG
  printf("[[%#llx] start ea: %#llx end ea: %#llx \n", spu_ctx.spe_id, psc_block.part_start, psc_block.part_end);
#endif

  fields_c_real_t ls_fld[6*32*32] __attribute__((aligned(128)));

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

  v_real dt, yl, zl, half, one;
  dt = spu_splats(spu_ctx.dt);
  half = spu_splats(0.5);
  yl = spu_mul(half, dt);
  zl = spu_mul(half, dt);
  one = spu_splats(1.0);
  

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
  n += 2; 

  } while(__builtin_expect((run),1));

  end_wait_particles_stored();
#if PRINT_DEBUG
  fprintf(stderr, "[[%#llx] ran %d particles\n", spu_ctx.spe_id, n);
#endif
  return 0;

}
