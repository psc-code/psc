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


int
spu_push_part_2d(void){

  printf("[[%#llx] start ea: %#llx end ea: %#llx \n", spu_ctx.spe_id, psc_block.part_start, psc_block.part_end);

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

    mfc_get(buff.lb1, cp_ea, 2*sizeof(particle_cbe_t), 
    	  tag_pget, 0, 0);
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
  
  mfc_write_tag_mask(1 << tag_pget);

  unsigned int mask = mfc_read_tag_status_any();

  //  fprintf(stderr, "mask %d \n", mask);
  particle_cbe_t null_part = *(buff.lb1);
  null_part.wni = 0.0;


  unsigned long long end = psc_block.part_end; 

  int run = 1;
  int n = 0;

  do {

    // issue dma request for particle we will need 
    // next time through the loop.
    if(__builtin_expect((np_ea < end),1)) {
      //      fprintf(stderr,"Preloading particle\n");
      mfc_get(buff.plb1, np_ea, 2 * sizeof(particle_cbe_t), 
	      tag_pget, 0, 0);
      mfc_write_tag_mask(1 << tag_pget);
      mask = mfc_read_tag_status_any();

    }    
    // we may need to insert some padding here, so we have to stop and check.
    // The last particle is going to be very slow (probably).
    else if(__builtin_expect(((end - cp_ea) != sizeof(particle_cbe_t)),0)){ 
      buff.lb2 = &null_part;
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
    // At this point, np_ea is one ahead of the
    // current particle.

    
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
    
    if(__builtin_expect((np_ea >= end),0)) { // if we've run out of particles
      fprintf(stderr, "Storing particle\n");
      mfc_put(buff.sb1, cp_ea, (size_t) (end - cp_ea),
	      tag_pput, 0, 0);
      mfc_write_tag_mask(1 << tag_pput);
      mask = mfc_read_tag_status_any();
      run = 0; 
    }
    else {
      mfc_put(buff.sb1, cp_ea, 2 * sizeof(particle_cbe_t),
	      tag_pput, 0, 0);

      cp_ea = np_ea;
      np_ea += 2*sizeof(particle_cbe_t);
      mfc_write_tag_mask(1 << tag_pput);
      mask = mfc_read_tag_status_any();

      // cp_ea now points to the particle which will be used
      // next time in the loop. 
      // np_ea points to the one which needs to be pre-loaded.
    }
  n += 2; 

  } while(__builtin_expect((run),1));

  fprintf(stderr, "[[%#llx] ran %d particles\n", spu_ctx.spe_id, n);
  return 0;

}
