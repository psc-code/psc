#ifndef PSC_CBE_H
#define PSC_CBE_H

// Order of includes is important!
// There are preprocessor sections in 
// psc_cbe_common.h which depend on 
// which include files preceded it. 
// I have a feeling this is very bad form.
#include "psc.h"
#include "psc_particles_as_cbe.h"
#include "psc_fields_as_c.h"

#include "psc_cbe_common.h"

#ifdef CELLEMU
#include "libspe2_c.h"
#else
#include <libspe2.h>
#endif 

#ifdef CELLEMU
#define NR_SPE (1)
#else
#define NR_SPE (8)
#endif 


struct cbe_block_layout {
  int block_size[3]; ///< The dimensions of the sub-domains (blocks) to be handed to spes.
  int block_grid[3]; ///< The number of subdomains in each direction 
  int nblocks; ///< Number of sub-domains. 
};

struct psc_spu_ops {
  // SPE Stuff
  spe_program_handle_t spu_test; ///< A wee test executable 
  spe_program_handle_t spu_1d; ///< Handle for 1-dimensional spe executable
  spe_program_handle_t spu_2d; ///< Handle for 2-dimensional spe executable
  spe_program_handle_t spu_3d; ///< Handle for 3-dimensional spe executable
  spe_context_ptr_t spe_id[NR_SPE]; ///< Address of each spe 
  int spes_inited; ///< Has the task been loaded onto the spes
  volatile int active_spes; ///< Number of spes currently running a job
  int spe_state[NR_SPE]; ///< Which spes are doing what
  pthread_t thread_id[NR_SPE]; ///< The thread corresponding to each spe
  psc_cell_block_t *spe_blocks[NR_SPE]; ///< Each spe's job description load/store cache

  // Work Block Stuff
  struct cbe_block_layout *layout; ///< Stuff related to the domain decomposition
  psc_cell_block_t ** block_list; ///< sub-domains to be handed off to the spes
  int blocks_inited; ///< Have the sub-domains been allocated.
  unsigned int * cnts; ///< The number of particles in each block. Easier to just save this here. 
  int particles_coarse_sorted; ///< Have the particles been sorted by work-block

  psc_cell_ctx_t *global_ctx; ///< Params and such needed on each spu for entire run.
};

int spu_main(unsigned long long spe_id, unsigned long long spu_comm_ea,
	     unsigned long long env);

void cbe_push_part_2d(void);

extern struct psc_spu_ops spu_ctl; ///< Information related to domains and spu programs

// Spe handeling functions from psc_cbe.c
void psc_init_spes(void);
void psc_kill_spes(void);
int get_spe(void);
void put_spe(int spe);
void update_spes_status(void);

// Block handeling functions from cbe_blocks.
void cbe_setup_layout(void);
void cbe_blocks_create(void);
void cbe_blocks_destroy(void);
void cbe_assign_parts_to_blocks(particles_cbe_t * pp);
void cbe_field_blocks_get(fields_c_t *pf, int mb, int me);
void cbe_field_blocks_get(fields_c_t *pf, int mb, int me);
void cbe_currents_put(fields_c_t *pf);


#endif
