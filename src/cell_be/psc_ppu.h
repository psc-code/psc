#ifndef PSC_CBE_H
#define PSC_CBE_H

#include "psc.h"
#include "psc_particles_cbe.h"
#include "psc_fields_cbe.h"

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


// The SPU can only load off 16byte boundaries, and as a 
// consequence each particle need to occupy some multiple
// of 16B. There might be another way around this, but
// for now this is a hacky workaround.

struct psc_spu_ops {
  spe_program_handle_t spu_test; ///< A wee test executable 
  spe_program_handle_t spu_1d; ///< Handle for 1-dimensional spe executable
  spe_program_handle_t spu_2d; ///< Handle for 2-dimensional spe executable
  spe_program_handle_t spu_3d; ///< Handle for 3-dimensional spe executable
  int block_size[3]; ///< The dimensions of the sub-domains (blocks) to be handed to spes.
  int block_grid[3]; ///< The number of subdomains in each direction 
  int nblocks; ///< Number of sub-domains. 
  int blocks_ready; ///< A flag to check if the blocks are ready to be passed to the spes (ie sorted, domain decomposed, etc)
  int * cnts; ///< The number of particles in each block. Easier to just save this here. 
};

int spu_main(unsigned long long spe_id, unsigned long long spu_comm_ea,
	     unsigned long long env);

void cbe_push_part_2d(void);

extern psc_cell_block_t *spe_blocks[NR_SPE]; ///< Block Cache for each spe
extern psc_cell_block_t **block_list;  ///< List of all the work blocks
extern int active_spes; ///< spes currently running a job.
extern spe_context_ptr_t spe_id[NR_SPE]; ///< The address of each spe
extern struct psc_spu_ops spu_ctl; ///< Information related to domains and spu programs

int get_spe(void);
void update_idle_spes(void);
// Deprecated alf stuff
/*
void cbe_push_part_yz(void);
void psc_alf_env_shared_create(unsigned int max_nodes);
void psc_alf_env_shared_destroy(int wait_time);
void wb_current_cache_init(spu_curr_cache_t * cache, push_wb_context_t * blk);
void wb_current_cache_store(fields_cbe_t * pf, spu_curr_cache_t * cache);
*/
#endif
