
#include <stdio.h>
#include <assert.h>

#include "psc_spu.h"

#include "psc_spu_2d.h"
// In case I need cell emulation

#ifdef __SPU__
#include <spu_mfcio.h>
#else
#include "spu_mfcio_c.h"
#endif

psc_cell_ctx_t spu_ctx;
psc_cell_block_t psc_block;


// Much of this is taken from Kai's openggcm spu implementation

// Error handler. Want these passed out the 
// PPE to handle. 
static void 
spu_error(void)
{
  spu_write_out_mbox(SPU_ERROR);
}


static int
spu_run_job(unsigned long long ea)
{
  spu_dma_get(&psc_block, ea, sizeof(psc_block));
  switch(psc_block.job) { 
  case SPU_PART:
    return spu_push_part_2d();
    break;
  case SPU_FIELD_A:
    return spu_push_field_a_nopml();
    break;
  case SPU_FIELD_B:
    return spu_push_field_b_nopml();
    break;
  default:
#ifndef NDEBUG
    fprintf(stderr, "spu: unknown job %d\n", psc_block.job);
#endif
    spu_error();
    return -1;
  }
} 


int
spu_main(unsigned long long spe_id, unsigned long long spu_comm_ea, 
	 unsigned long long env)
{
  unsigned int msg_in, msg_out; 

#if PRINT_DEBUG
    fprintf(stderr,"spu main [%#llx]\n", spe_id);
#endif
    
    spu_dma_get(&spu_ctx, env, sizeof(spu_ctx));
#if PRINT_DEBUG
    fprintf(stderr,"context %#llx, dt: %g\n", env, spu_ctx.dt);
#endif
    spu_ctx.spe_id = spe_id; 

    int rc; 
    msg_out = SPE_READY; 
    spu_write_out_mbox(msg_out);
    
    for(;;){
#if PRINT_DEBUG
      fprintf(stderr, "Sitting in Loop\n");
#endif
      msg_in = spu_read_in_mbox();
      
      switch(msg_in) {
      case SPU_QUIT:
#if PRINT_DEBUG
	fprintf(stderr, "Got msg SPU_QUIT\n");
#endif
	return 0;
	
      case SPU_RUNJOB:
#if PRINT_DEBUG
	fprintf(stderr, "Got msg SPU_RUNJOB\n");
#endif
	rc = spu_run_job(spu_comm_ea);
	assert(rc == 0); 
	msg_out = SPE_IDLE;
	spu_write_out_mbox(msg_out); 
	break;
	  
      case SPE_CLEAR:
#if PRINT_DEBUG
	fprintf(stderr, "Got msg SPU_ClEAR\n");
#endif
	msg_out = SPE_CLEAR;
	spu_write_out_mbox(msg_out);
	break;
	
      default:
#if PRINT_DEBUG
	fprintf(stderr, "spu: unknown msg %#x\n", msg_in);
#endif
	spu_error();
	return -1;
      }
    }
    
    return 0;
}
