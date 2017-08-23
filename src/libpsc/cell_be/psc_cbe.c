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

// Handles to the spu code in different dimensions. Accessed through the
// spu_psc_cfs.o library. 
//extern spe_program_handle_t spu_1d; ///< Handle for 1-dimensional spe executable
static spe_program_handle_t spu_2d; ///< Handle for 2-dimensional spe executable
//extern spe_program_handle_t spu_3d; ///< Handle for 3-dimensional spe executable

// Some static variables dealing with spe control which are only needed
// in this file (which is why they are static, duh)
static spe_context_ptr_t spe_id[NR_SPE]; ///< Address of each spe 
static pthread_t thread_id[NR_SPE]; ///< The thread corresponding to each spe
static int spe_state[NR_SPE]; ///< Which spes are doing what
static int active_spes; ///< Number of spes currently running a job
static psc_cell_block_t *spe_blocks[NR_SPE]; ///< Each spe's job description load/store cache
static psc_cell_ctx_t *global_ctx; ///< Params and such needed on each spu for entire run.

// Apparently the standard states that variable declared 'static'
// are initialized to 0. I'm not realy happy relying on this behaviour, 
// but at the moment I think is the best way to do it and get rid of
// the cbe_create function.
static int spes_inited; ///< Has the task been loaded onto the spes


static void 
global_ctx_create()
{
  int rc; 
  void *m; 
  rc = posix_memalign(&m, 128, sizeof(psc_cell_ctx_t));
  assert(rc == 0);
  global_ctx = (psc_cell_ctx_t *) m;
  global_ctx->spe_id = 0;
  global_ctx->dx[0] = psc.dx[0];
  global_ctx->dx[1] = psc.dx[1];
  global_ctx->dx[2] = psc.dx[2];
  global_ctx->dt = psc.dt;
  global_ctx->eta = psc.coeff.eta;
  global_ctx->fnqs = sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta;
}


// I'm not entirely sure what the point of only passing this
// function a void pointer is...
// (which is odd, because I wrote it)
static void *
spe_thread_function(void *data)
{
  int i = (unsigned long) data;
  int rc; 
  unsigned int entry = SPE_DEFAULT_ENTRY;
  
  //  fprintf(stderr, "block pointer %p\n", spe_blocks[i]);
  do {
    rc = spe_context_run(spe_id[i], &entry, 0,
			 spe_blocks[i], global_ctx, NULL);
  } while (rc > 0);

  pthread_exit(NULL);
}

void 
psc_init_spes(void)
{
  

  if(!spes_inited){

    //  spu_test = test_handle; 
    spu_2d = spu_2d_handle; 
    //  assert(sizeof(psc_cell_ctx_t) % 16 == 0);
    int rc; 
    spe_program_handle_t spu_prog;
    
    global_ctx_create();
    
    spu_prog = spu_2d;  
    
    for (int i = 0; i < NR_SPE; i++){
      void *m;
      rc = posix_memalign(&m, 128, sizeof(psc_cell_block_t)); 
      assert(rc == 0);
      spe_blocks[i] = (psc_cell_block_t *) m;
      assert(spe_blocks[i] != NULL);
      spe_id[i] = spe_context_create(0,NULL);
      spe_program_load(spe_id[i],&spu_prog);
      rc = pthread_create(&thread_id[i], NULL, spe_thread_function,
			(void *)(unsigned long) i); 
      assert(rc == 0);
    }
  
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double start =  tv.tv_sec;


    // FIXME: There may be something wrong with this construct, and it
    // feels kind of clumsy. Better way to check that spe's actually loaded
    // the program?
    //
    // This may be a good place to use some of the more sophisticated spe
    // synchroniztion or error handeling function.
    int ready = 0; 
    while(ready < NR_SPE){
      for(int spe=0; spe < NR_SPE; spe++){
	if (spe_out_mbox_status(spe_id[spe]) > 0) {
	  unsigned int msg;
	  int nmesg = spe_out_mbox_read(spe_id[spe], &msg,1);
	  while(nmesg > 0) {
	    if(msg == SPE_READY){
	      spe_state[spe] = SPE_IDLE;
	      ready++;
	    }
	    assert(msg != SPU_ERROR);
	    nmesg = spe_out_mbox_read(spe_id[spe], &msg,1);
	  }
	}
	gettimeofday(&tv, NULL);
	if((tv.tv_sec - start) > 30) {
	  fprintf(stderr, "Did not obtian ready from all spes after 30s. Exiting\n");
	  assert(0);
	}
      }
    }
  
    active_spes = 0;
  
    spes_inited = 1; 
  }
}

void
psc_kill_spes(void)
{
 
  assert(active_spes == 0);
  if(spes_inited) {
    // Shutdown the spe threads
    unsigned int msg; 
    for(int i = 0; i<NR_SPE; i++){
      msg = SPU_QUIT;
      int nmesg = spe_in_mbox_write(spe_id[i], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
      assert(nmesg == 1);
      free(spe_blocks[i]);
      spe_blocks[i] = NULL;
    }
    // free the global ctx
    free(global_ctx);
    global_ctx = NULL;
    spes_inited = 0;
  }
}
  


static int
get_spe(void)
{
  assert(spes_inited);
  
  int spe; 
  
  // Look for first idle spe
  for(spe = 0; spe <= NR_SPE; spe++){
    if (spe_state[spe] == SPE_IDLE)
      break;
  }
  
  // This assert checks that we never call this function
  // when all the SPEs are being used ( I think)
  // seems a bit... paranoid.
  assert(spe < NR_SPE);
  
  spe_state[spe] = SPE_RUN;
  active_spes++;
  
  assert(active_spes <= NR_SPE);

  return spe;
}


static void
put_spe(int spe)
{
  spe_state[spe] = SPE_IDLE;
  active_spes--;
}

void
update_spes_status(void)
{
  for(int spe=0; spe < NR_SPE; spe++){
    unsigned int msg;
    int nmesg; 
    nmesg = spe_out_mbox_read(spe_id[spe], &msg,1);
    assert(msg != SPU_ERROR);
    while(nmesg > 0){
      if (msg == SPE_IDLE) {
#if PRINT_DEBUG
	fprintf(stderr, "[ppe] Got SPE_IDLE\n");
#endif
	put_spe(spe);
	nmesg = spe_out_mbox_read(spe_id[spe], &msg,1);
	// could eliminate this by checking the return value 
      }
      assert(msg != SPU_ERROR);
    }
  }
}


void cell_run_patch(int p,struct psc_fields *pf, particles_t *pp, int job) 
{

  while(active_spes == NR_SPE) {
    update_spes_status();
  }
  
  int spe = get_spe();
  
  spe_blocks[spe]->job = job;
  for(int i = 0; i < 3; i++){
    spe_blocks[spe]->ib[i] = pf->ib[i];
    spe_blocks[spe]->im[i] = pf->im[i];
    spe_blocks[spe]->xb[i] = psc.patch[p].xb[i];
  }
  spe_blocks[spe]->wb_flds = (unsigned long long) pf->flds;
  spe_blocks[spe]->part_start = (unsigned long long) pp->particles;
  spe_blocks[spe]->part_end = (unsigned long long) (pp->particles + pp->n_part);

#ifdef DEBUG_PATCH_SCHED
  fprintf(stderr, "queing patch %d\n", p);
  fprintf(stderr, "patch characteristics:\n");
  fprintf(stderr, "Job ID: %d\n", job);
  fprintf(stderr, "ib: [%d %d %d]\n",spe_blocks[spe]->ib[0],spe_blocks[spe]->ib[1],spe_blocks[spe]->ib[2]);
  fprintf(stderr, "im: [%d %d %d]\n",spe_blocks[spe]->im[0],spe_blocks[spe]->im[1],spe_blocks[spe]->im[2]);
  fprintf(stderr, "xb: [%g %g %g]\n",spe_blocks[spe]->xb[0],spe_blocks[spe]->xb[1],spe_blocks[spe]->xb[2]);
  fprintf(stderr, "part_start: %lld part_end: %lld\n",spe_blocks[spe]->part_start,spe_blocks[spe]->part_end);
#endif

  unsigned int msg = SPU_RUNJOB;
  int nmesg =  spe_in_mbox_write(spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
  assert(nmesg == 1);

}

void wait_all_spe(void)
{
  while(active_spes != 0) {
    update_spes_status();
  }
}



