#include <stdlib.h>
#include <assert.h>
#include <libspe2.h>
#include <sys/time.h>
#include "psc.h"
#include "psc_ppu.h"


#include <pthread.h>
#include <string.h>

#if CELLEMU
spe_program_handle_t test_handle = spu_main;
#endif

//extern spe_program_handle_t test_handle;
extern spe_program_handle_t spu_2d_handle; 

static pthread_t thread_id[NR_SPE];

struct psc_spu_ops spu_ctl; 


///////
/// Computes the global parameters for the task context and registers their existence
/// in the task description
///
/// Call after the task has been fully described but before task creation.

/*
void 
psc_alf_init_push_task_context(push_task_context_t * tcs, alf_task_desc_handle_t * task_desc, fields_cbe_t *pf)
{
  // Calculate the values, and stick them into the struct
  tcs->dxi = 1.0 / psc.dx[0];
  tcs->dyi = 1.0 / psc.dx[1];
  tcs->dzi = 1.0 / psc.dx[2];
  tcs->dt = psc.dt;
  tcs->xl = 0.5 * psc.dt;
  tcs->yl = 0.5 * psc.dt;
  tcs->zl = 0.5 * psc.dt;
  tcs->eta = psc.coeff.eta;
  tcs->fnqs = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
  tcs->fnqxs = psc.dx[0] * tcs->fnqs / psc.dt;
  tcs->fnqys = psc.dx[1] * tcs->fnqs / psc.dt;
  tcs->fnqzs = psc.dx[2] * tcs->fnqs / psc.dt;
  tcs->dqs = 0.5*psc.coeff.eta*psc.dt;
  tcs->ilg[0] = psc.ilg[0];
  tcs->ilg[1] = psc.ilg[1];
  tcs->ilg[2] = psc.ilg[2];
  tcs->img[0] = psc.img[0];
  tcs->img[1] = psc.img[1];
  tcs->img[2] = psc.img[2];
  tcs->p_fields = (unsigned long long) pf->flds;
  // Register the number and type of the parameters in the handle
  int ierr;
  ierr = alf_task_desc_ctx_entry_add(*task_desc, ALF_DATA_BYTE, sizeof(push_task_context_t) ); ACE;
}

*/



static void 
global_ctx_create()
{
  int rc; 
  void *m; 
  rc = posix_memalign(&m, 128, sizeof(psc_cell_ctx_t));
  assert(rc == 0);
  spu_ctl.global_ctx = (psc_cell_ctx_t *) m;
  psc_cell_ctx_t * global_ctx = spu_ctl.global_ctx; 
  global_ctx->spe_id = 0;
  global_ctx->dx[0] = psc.dx[0];
  global_ctx->dx[1] = psc.dx[1];
  global_ctx->dx[2] = psc.dx[2];
  global_ctx->dt = psc.dt;
  global_ctx->eta = psc.coeff.eta;
  global_ctx->fnqs = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
}


static void *
spe_thread_function(void *data)
{
  int i = (unsigned long) data;
  int rc; 
  unsigned int entry = SPE_DEFAULT_ENTRY;
  
  //  fprintf(stderr, "block pointer %p\n", spe_blocks[i]);
  do {
    rc = spe_context_run(spu_ctl.spe_id[i], &entry, 0,
			 spu_ctl.spe_blocks[i], spu_ctl.global_ctx, NULL);
  } while (rc > 0);

  pthread_exit(NULL);
}

void 
psc_init_spes(void)
{
  
  assert(!spu_ctl.spes_inited);

  //  spu_ctl.spu_test = test_handle; 
  spu_ctl.spu_2d = spu_2d_handle; 
  //  assert(sizeof(psc_cell_ctx_t) % 16 == 0);
  int rc; 
  spe_program_handle_t spu_prog;

  global_ctx_create();

  spu_prog = spu_ctl.spu_2d;  

  for (int i = 0; i < NR_SPE; i++){
    void *m;
    rc = posix_memalign(&m, 128, sizeof(psc_cell_block_t)); 
    assert(rc == 0);
    spu_ctl.spe_blocks[i] = (psc_cell_block_t *) m;
    assert(spu_ctl.spe_blocks[i] != NULL);
    spu_ctl.spe_id[i] = spe_context_create(0,NULL);
    spe_program_load(spu_ctl.spe_id[i],&spu_prog);
    rc = pthread_create(&thread_id[i], NULL, spe_thread_function,
			(void *)(unsigned long) i); 
    assert(rc == 0);
  }
  
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double start =  tv.tv_sec;

  int ready = 0; 
  while(ready < NR_SPE){
    for(int spe=0; spe < NR_SPE; spe++){
      unsigned int msg;
      spe_out_mbox_read(spu_ctl.spe_id[spe], &msg,1);
      if(msg == SPE_READY){
	spu_ctl.spe_state[spe] = SPE_IDLE;
	ready++;
      }
      assert(msg != SPU_ERROR);
      gettimeofday(&tv, NULL);
      if((tv.tv_sec - start) > 30) {
	fprintf(stderr, "Did not obtian ready from all spes after 30s. Exiting\n");
	assert(0);
      }
    }
  }
  
  spu_ctl.active_spes = 0;
  
  spu_ctl.spes_inited = 1; 
}

void
psc_kill_spes(void)
{
  
  assert(spu_ctl.spes_inited);
  
  unsigned int msg; 
  for(int i = 0; i<NR_SPE; i++){
    msg = SPU_QUIT;
    spe_in_mbox_write(spu_ctl.spe_id[i], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
    free(spu_ctl.spe_blocks[i]);
    spu_ctl.spe_blocks[i] = NULL;
  }
}
  


int
get_spe(void)
{
  assert(spu_ctl.spes_inited);
  
  int spe; 
  
  // Look for first idle spe
  for(spe = 0; spe <= NR_SPE; spe++){
    if (spu_ctl.spe_state[spe] == SPE_IDLE)
      break;
  }
  
  // This assert checks that we never call this function
  // when all the SPEs are being used ( I think)
  assert(spe < NR_SPE);
  
  spu_ctl.spe_state[spe] = SPE_RUN;
  spu_ctl.active_spes++;
  
  assert(spu_ctl.active_spes <= NR_SPE);

  return spe;
}


void
put_spe(int spe)
{
  spu_ctl.spe_state[spe] = SPE_IDLE;
  spu_ctl.active_spes--;
}

void
update_spes_status(void)
{
  for(int spe=0; spe < NR_SPE; spe++){
    unsigned int msg;
    spe_out_mbox_read(spu_ctl.spe_id[spe], &msg,1);
    assert(msg != SPU_ERROR);
    if((msg == SPE_IDLE) && (spu_ctl.spe_state[spe] == SPE_RUN)){
#if PRINT_DEBUG
      fprintf(stderr, "[ppe] Got SPE_IDLE\n");
#endif
      put_spe(spe);
      unsigned int msg_out = SPE_CLEAR;
      spe_in_mbox_write(spu_ctl.spe_id[spe], &msg_out, 1, SPE_MBOX_ANY_BLOCKING);
      // could eliminate this by checking the return value 
      while(msg != SPE_CLEAR){
	spe_out_mbox_read(spu_ctl.spe_id[spe], &msg,1);
      }
    }
  }
}



static void 
cbe_create(void)
{

  // Because I don't want to take the performance
  // it of storing the cached fields to the full domain
  // between each mode, I'm going to force a consistency 
  // check. If cbe_create is being called (you're using
  // the cbe particle pusher) you damn well better be using
  // the other cell mods too. Right now, they're just going to 
  // wrap the C versions (excepting the sort), but they 
  // will be implemented eventually. 

  spu_ctl.spes_inited = 0;
  spu_ctl.blocks_inited = 0;
  spu_ctl.particles_coarse_sorted = 0; 

  // To make sure every thing is layed out as planned,
  // we'll set anything to be allocated in the control 
  // struct to null:
  spu_ctl.layout = NULL; 
  spu_ctl.block_list = NULL;
  spu_ctl.cnts = NULL; 
  spu_ctl.global_ctx = NULL; 

}


static void
cbe_destroy(void)
{
  
  psc_kill_spes(); 

  cbe_blocks_destroy(); 
  
  free(spu_ctl.global_ctx);
  spu_ctl.global_ctx = NULL;
  
  free(spu_ctl.cnts);
  spu_ctl.cnts = NULL; 

  free(spu_ctl.layout);
  spu_ctl.layout = NULL; 
   
}


struct psc_ops psc_ops_cbe = {
  .name                   = "cbe",
  .create                 = cbe_create,
  .destroy                = cbe_destroy, 
  .push_part_yz           = cbe_push_part_2d,
};



