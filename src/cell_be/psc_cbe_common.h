#ifndef PSC_CBE_COMMON_H
#define PSC_CBE_COMMON_H

#ifdef CELLEMU
#define NR_SPE (1)
#else
#define NR_SPE (8)
#endif 


enum kern {
  SPU_HELLO,
  SPU_BYE,
  NR_KERN,
};

enum { 
  SPU_QUIT,
  SPU_ERROR,
  SPU_RUNJOB,
  SPE_IDLE,
  SPE_RUN,
};


// The SPU can only load off 16byte boundaries, and as a 
// consequence each particle need to occupy some multiple
// of 16B. There might be another way around this, but
// for now this is a hacky workaround.

  
/// Parameters which are the same across all blocks and compute
/// kernels. 
typedef struct _psc_cell_ctx
{
  char hello[8];
  char bye[8];
  unsigned long long spe_id; 
  unsigned long long padding; 
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
  unsigned long long job;
  unsigned long long padding; 
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

#endif
