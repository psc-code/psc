#include <stdio.h>
#include <assert.h>

// In case I need cell emulation

#ifdef __SPU__
#include <spu_mfcio.h>
#else
#include "spu_mfcio_c.h"
#endif

// Much of this is taken from Kai's openggcm spu implementation

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

// Error handler. Want these passed out the 
// PPE to handle. 
static void 
spu_error(void)
{
  spu_write_out_mbox(SPU_ERROR);
}



int
spu_main(unsigned long long spe_id, unsigned long long spu_comm_ea, 
	 unsigned long long env)
{
  unsigned int msg_in, msg_out; 

#ifndef NDEBUG
    printf("spu main [%#llx]\n", spe_id);
#endif
    spu_dma_get(&psc_env, env, sizeof(psc_env));
    psc_env.spe_id = spe_id; 
  // Using loop lets us run switches 
  // inside an infinite loop.
  int rc; 
  for(;;){
    msg_in = spu_read_in_mbox();

    switch(msg_in) {
    case SPU_QUIT:
      return 0;

    case SPU_RUNJOB:
      rc = spu_run_job(spu_comm_ea);
      assert(rc == 0); 
      msg_out = SPU_IDLE;
      rc = spu_write_out_mbox(msg_out); 
      break;

    default:
#ifndef NDEBUG
      fprintf(stderr, "spu: unknown msg %#x\n", msg);
#endif
      spu_error();
      return -1;
    }

    // exit successfully
    //    spu_write_out_mbox(msg);
  }

  return 0;
}
