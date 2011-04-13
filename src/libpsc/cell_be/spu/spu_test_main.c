#include <stdio.h>
#include <assert.h>

#include "psc_spu.h"
#include "spu_test_func.h"
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
  case SPU_HELLO:
    return spu_psc_hello();
    break;
  case SPU_BYE:
    return spu_psc_goodbye();
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

#ifndef NDEBUG
    printf("spu main [%#llx]\n", spe_id);
#endif
    //    fprintf(stderr, "spu main [%#llx] ea %p env %p \n", spe_id, spu_comm_ea, env);
    //    fprintf(stderr, "spu main [%#llx] spu_ctx %p size %d \n", spe_id, &spu_ctx, sizeof(spu_ctx));
    spu_dma_get(&spu_ctx, env, sizeof(spu_ctx));
    spu_ctx.spe_id = spe_id; 

    int rc; 
    msg_out = SPE_IDLE; 
    spu_write_out_mbox(msg_out);
    
    for(;;){
      msg_in = spu_read_in_mbox();
      
      switch(msg_in) {
      case SPU_QUIT:
	return 0;
	
      case SPU_RUNJOB:
	rc = spu_run_job(spu_comm_ea);
	assert(rc == 0); 
	msg_out = SPE_IDLE;
	spu_write_out_mbox(msg_out); 
	break;
	
      default:
#ifndef NDEBUG
	fprintf(stderr, "spu: unknown msg %#x\n", msg_in);
#endif
	spu_error();
	return -1;
      }
      
      // exit successfully
      //    spu_write_out_mbox(msg);
    }
    
    return 0;
}
