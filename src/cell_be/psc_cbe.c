#include <stdlib.h>
#include <alf.h>
#include <assert.h>
#include "psc.h"
#include "psc_cbe.h"
#include "psc_alf.h"

alf_handle_t alf_env_shared = ALF_NULL_HANDLE;

///////
/// Set up a shared ALF environment.
/// 
/// Queries the system for number of available spus.
/// Pass max_nodes = 0 to use all available spus.
/// If max_nodes is greater than than the number of 
/// available nodes, all available nodes are used.
/// If max_nodes is less than the number of available
/// nodes, only max_nodes are used. On certain systems, 
/// for example Blade Centers with dual CBE configurations, 
/// it may be useful to restrict the number of nodes used 
/// to the ALF environment to only on chip (8) spus, if
/// memory bandwidth is a major limitation.
void
psc_alf_env_shared_create(unsigned int max_nodes)
{

  int ierr;
  unsigned int avail_nodes;
  unsigned int nodes;

  
  // Create alf environment
  // Right now, path to shared accelerator code libraries is gathered from 
  // environment variable ALF_LIBRARY_PATH
  ierr = alf_init_shared(NULL, &alf_env_shared); ACE; 
  
  //Find out how many spu's we have available.
  ierr = alf_query_system_info(alf_env_shared, ALF_QUERY_NUM_ACCEL, 0 ,&avail_nodes); ACE; 
  assert(avail_nodes > 0);
  if( max_nodes == 0) {
    nodes = avail_nodes;
  } else if ( max_nodes < avail_nodes ) {
    nodes = max_nodes;
  } else {
    nodes = avail_nodes;
  }

  // Tell alf how many spus to use.
  ierr = alf_num_instances_set(alf_env_shared, nodes); ACE;
}

///////
/// Shut down shared ALF environment.
///
/// Will wait at most 'wait_time' milliseconds for
/// any running tasks to finish before forcing shutdown
/// and setting alf_env_shared to ALF_NULL_HANDLE.
/// Passing 0 forces immediate shutdown. 
/// Passing <0 will wait forever (*not* recommended).
void
psc_alf_env_shared_destroy(int wait_time)
{
    int ierr = alf_exit(alf_env_shared, ALF_EXIT_POLICY_WAIT, wait_time);ACE;
    alf_env_shared = ALF_NULL_HANDLE;
}

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

static void 
cbe_create(void)
{
  if( alf_env_shared == ALF_NULL_HANDLE ){
    psc_alf_env_shared_create(0);
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
  .push_part_yz           = cbe_push_part_yz,
};
