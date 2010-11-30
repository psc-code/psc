#ifndef PSC_CBE_H
#define PSC_CBE_H

#include "psc.h"
#include "psc_particles_cbe.h"
#include "psc_fields_cbe.h"

enum kern {
  SPU_HELLO,
  SPU_BYE,
  NR_KERN,
};

// The SPU can only load off 16byte boundaries, and as a 
// consequence each particle need to occupy some multiple
// of 16B. There might be another way around this, but
// for now this is a hacky workaround.

  
/// Parameters which are the same across all blocks and compute
/// kernels. 
typedef struct _psc_cell_ctx
{
  const char hello[8];
  const char bye[8];
  unsigned long long spe_id; 

  /*
  cbe_real dxi; ///< 1/dx
  cbe_real dyi; ///< 1/dy
  cbe_real dzi; ///< 1/dz
  cbe_real dt; ///< timestep size
  cbe_real xl; ///< 0.5dt
  cbe_real yl; ///< 0.5dt
  cbe_real zl; ///< 0.5dt
  cbe_real eta; ///< psc.coeff.eta
  cbe_real fnqs; ///< sqr(psc.coeff.alpha) * psc.coeff.cori / psc.coeff.eta
  cbe_real fnqxs; ///< psc.dx[0] * fnqsfl / psc.dt;
  cbe_real fnqys; ///< psc.dx[1] * fnqsfl / psc.dt;
  cbe_real fnqzs; ///< psc.dx[2] * fnqsfl / psc.dt;
  cbe_real dqs; ///< 0.5*psc.coeff.eta*psc.dt
  int ilg[3]; 
  int img[3];
  unsigned long long p_fields; ///< address of field array
  */
} psc_cell_ctx_t;

/// Parameters specific to each work block
typedef struct _psc_cell_block
{
  enum kern job;
  /*
  unsigned long long p_start; ///< starting address of WB particles. 
  unsigned long long p_pad; ///< address of WB padding particles. 
  unsigned long long p_cache; ///< address of WB current cache
  unsigned int lo_npart; ///< number of particles proccesed by the work block
  unsigned int maxpart; ///< upper limit of particles on spu at one time
  int wb_lg[3]; ///< lower work block index in xyz w/ ghosts
  int wb_hg[3]; ///< upper work block index in xyz w/ ghosts
  */
} psc_cell_block_t;

typedef struct
{
  fields_cbe_real_t * flds;
  int lg[3];
  int hg[3];
} spu_curr_cache_t;


struct psc_spu_ops {
  void (*spu_test)(void);
  void (*spu_1d)(void);
  void (*spu_2d)(void);
  void (*spu_3d)(void);
}
void cbe_push_part_yz(void);
void psc_alf_env_shared_create(unsigned int max_nodes);
void psc_alf_env_shared_destroy(int wait_time);
void wb_current_cache_init(spu_curr_cache_t * cache, push_wb_context_t * blk);
void wb_current_cache_store(fields_cbe_t * pf, spu_curr_cache_t * cache);

#endif
