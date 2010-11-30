#include <stdlib.h>
#include <assert.h>
#include "psc.h"
#include "psc_cbe.h"

#ifdef CELLEMU
#include <libspe2_c.h>
#else
#include <libspe2.h>
#endif 

#include <pthread.h>
#include <string.h>

#ifdef CELLEMU
#define NR_SPE (1)
#else
#define NR_SPE (8)
#endif 



enum { 
  SPE_IDLE,
  SPE_RUN,
};

struct psc_spu_ops spu_progs = {
  .test = test_handle,
};

psc_cell_ctx_t global_ctx __attribute__((aligned(128)));
psc_cell_block_t *spe_blocks[NR_SPE];
block_node_t * staged_blocks; 
block_node_t * finished_blocks; 

static int spes_inited;
static int active_spes;
static int spe_state[NR_SPE];
static spe_context_ptr_t spe_id[NR_SPE];
static pthread_t thread_id[NR_SPE];


///////
/// Computes the global parameters for the task context and registers their existence
/// in the task description
///
/// Call after the task has been fully described but before task creation.
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

void
wb_current_cache_init(spu_curr_cache_t * cache, push_wb_context_t * blk)
{
   cache->lg[0] = blk->wb_lg[0];
  cache->lg[1] = blk->wb_lg[1];
  cache->lg[2] = blk->wb_lg[2];
  cache->hg[0] = blk->wb_hg[0];
  cache->hg[1] = blk->wb_hg[1];
  cache->hg[2] = blk->wb_hg[2];
  
  int fld_size = (blk->wb_hg[0] - blk->wb_lg[0] + 1)
    * (blk->wb_hg[1] - blk->wb_lg[1] + 1)
    * (blk->wb_hg[2] - blk->wb_lg[2] + 1);

  void *m;
  int ierr = posix_memalign(&m, 128, 4*fld_size*sizeof(fields_cbe_real_t));
  assert(ierr == 0);
  memset(m,0,4*fld_size*sizeof(fields_cbe_real_t));
  cache->flds = (fields_cbe_real_t *)m;
  blk->p_cache = (unsigned long long) cache->flds;
}

void
wb_current_cache_store(fields_cbe_t *pf, spu_curr_cache_t * cache)
{
#define JC_OFF(jx,jy,jz)			\
  (((((jz) - lg[2])				\
     *img[1] + ((jy) - lg[1]))			\
    *img[0] + ((jx) - lg[0]))			\
   *4)

  int lg[3];
  int hg[3];
  lg[0] = cache->lg[0];
  lg[1] = cache->lg[1];
  lg[2] = cache->lg[2];
  hg[0] = cache->hg[0];
  hg[1] = cache->hg[1];
  hg[2] = cache->hg[2];

  int img[3] = {hg[0] - lg[0] + 1,
		hg[1] - lg[1] + 1,
		hg[2] - lg[2] + 1};

  for(int jz = lg[2]; jz <= hg[2]; jz++){
    for(int jy = lg[1]; jy <= hg[1]; jy++){
      for(int jx = lg[0]; jx <= hg[0]; jx++){
      	F3_CBE(pf, JXI,jx, jy, jz) += cache->flds[0 + JC_OFF(jx,jy,jz)];
	F3_CBE(pf, JYI,jx, jy, jz) += cache->flds[1 + JC_OFF(jx,jy,jz)];
	F3_CBE(pf, JZI,jx, jy, jz) += cache->flds[2 + JC_OFF(jx,jy,jz)];
      }
    }
  }
}
#undef JC_OFF

static void *
spe_thread_function(void *data)
{
  int i = (unsigned long) data;
  int rc; 
  unsigned int entry = SPE_DEFAULT_ENTRY;

  do {
    rc = spe_context_run(spe_id[i], &entry, 0,
			 spe_blocks[i], &global_ctx, NULL);
  } while (rc > 0);

  pthread_exit(NULL);
}

static void 
cbe_create(void)
{
  assert(sizeof(psc_cell_env_t) % 16 == 0);
  int rc; 
  spu_program_handle_t spu_prog;

  spu_prog = spu_progs.test;
  
  for (int i = 0; i < NR_SPE; i++){
    void *m;
    rc = posix_memalign(&m, 128, sizeof(psc_cell_block_t)); 
    assert(rc == 0);
    spe_blocks[i] = (psc_cell_block_t *) m;
    spe_id[i] = spe_context_create(0,NULL);
    spe_program_load(spe_id[i],&spu_prog);
    rc = pthread_create(&thread_id[i], NULL, spe_thread_function,
			(void *)(unsigned long) i); 
    assert(rc == 0);
    spe_state[i] = SPE_IDLE;
  }
  active_spes = 0;
  
  spes_inited = 1; 
  }
}


int
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
  assert(spe < NR_SPE);
  
  spe_state[spe] = SPE_RUN;
  active_spes++;
  
  assert(active_spes <= NR_SPE);

  return spe;
}

void
put_spe(int spe)
{
  spe_state[spe] = SPE_IDLE;
  active_spes--;
}

block_node_t *
block_ll_pop(block_node_t ** head) 
{
  block_node_t * old_head = *head;
  *head = *head->next;
  return old_head;
}

void
block_ll_push(block_node_t * new, block_node_t ** head)
{
  new->next = *head;
  *head = new; 
}

block_node_t *
block_create(int * blo, int * bhi)
{
  block_node_t * node = malloc(sizeof(block_node_t));
  node->block = calloc(1,sizeof(psc_cell_block_t));
  assert(node->block != NULL);
  node->next = NULL;
  return node;
}
 
void
block_ll_create(void)
{
  for(int i = 0; i < NR_SPES * 2; i++){
    block_ll_push(block_create(1,1), &staged_blocks);
  }
}
 
/// \FIXME Do the particle module destroy functions every get called?
/// If not, we're going to be leaking an ALF environment...

static void
cbe_destroy(void)
{
  if( alf_env_shared != ALF_NULL_HANDLE ){
    psc_alf_env_shared_destroy(6000);
  }
}


struct psc_ops psc_ops_cbe = {
  .name                   = "cbe",
  .create                 = cbe_create,
  .destroy                = cbe_destroy,
  .push_part_xy           = cbe_push_part_2d,
};
