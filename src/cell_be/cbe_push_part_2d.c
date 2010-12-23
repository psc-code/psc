#include <string.h>
#include "psc_ppu.h"
#include "util/profile.h"
#include "psc_particles_cbe.h"

#ifdef CELLEMU
#include "libspe2_c.h"
#else
#include <libspe2.h>
#endif

static void
do_cbe_push_part_2d(particles_cbe_t *pp, fields_c_t *pf)
 {
  
  unsigned int msg; 
  int state;
  psc_cell_block_t ** active_blk; 

  cbe_assign_parts_to_blocks(pp);

  active_blk = spu_ctl.block_list; 

  while( (*active_blk != NULL) || (spu_ctl.active_spes != 0)){
    while(spu_ctl.active_spes < NR_SPE && *active_blk != NULL){
      int spe = get_spe();
      (*active_blk)->job = SPU_PART;
      memcpy(spu_ctl.spe_blocks[spe], *active_blk, sizeof(psc_cell_block_t));
      msg = SPU_RUNJOB;
      spe_in_mbox_write(spu_ctl.spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
#if PRINT_DEBUG
      fprintf(stderr, "[ppe queing %p] start %#llx end %#llx\n", spu_ctl.spe_id[spe],
	      spu_ctl.spe_blocks[spe]->part_start, spu_ctl.spe_blocks[spe]->part_end);
#endif
      active_blk++;
    }
    update_spes_status();
    fflush(stderr);
    fflush(stdout);
  }

  spu_ctl.particles_coarse_sorted = 0; 
}  

void
cbe_push_part_2d(void)
{
  fields_c_t pf;
  particles_cbe_t pp;

  particles_cbe_get(&pp);
  fields_c_get(&pf, EX, EX+6);

  if(!spu_ctl.spes_inited)
    psc_init_spes();

  cbe_field_blocks_get(&pf,EX,EX+6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_part_yz", 1., 0, psc.pp.n_part * 9 * sizeof(cbe_real));
  }
  prof_start(pr);
  do_cbe_push_part_2d(&pp, &pf);
  prof_stop(pr);

  cbe_currents_put(&pf);

  particles_cbe_put(&pp);
  fields_c_put(&pf,JXI,JXI+3);
}
