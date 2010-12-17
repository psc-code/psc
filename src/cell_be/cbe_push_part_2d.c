#include <string.h>
#include "psc_ppu.h"
#include "util/profile.h"
#include "psc_particles_cbe.h"

#ifdef CELLEMU
#include "libspe2_c.h"
const unsigned int all_spes = 1;
#else
#include <libspe2.h>
const unsigned int all_spes = 255;
#endif

static void
do_cbe_push_part_2d(particles_cbe_t *pp, fields_cbe_t *pf)
 {
  
  unsigned int msg; 
  int state;
  psc_cell_block_t ** active_blk; 

  assert(spu_ctl.cnts);

  unsigned int *cnts = spu_ctl.cnts;

  particle_cbe_t *fp = pp->particles; 

  psc_cell_block_t ** curr = block_list;
  (*curr)->part_start =  fp;
  (*curr)->part_end = (fp + cnts[1]);
  curr++;
  
  fprintf(stderr, "nblocks %i \n", spu_ctl.nblocks);

  for(int i = 1; i < spu_ctl.nblocks; i++){
    fprintf(stderr, "[%d] ctns[i] %d cnts[i+1] %d\n", i, cnts[i], cnts[i+1]);
    (*curr)->part_start = (fp + cnts[i]);
    (*curr)->part_end = (fp + cnts[i+1]);
    fprintf(stderr, "[%d] start %p end %p\n", i, (*curr)->part_start, (*curr)->part_end);
    curr++;
  }
    
    spu_ctl.blocks_ready = 1;


  active_blk = block_list; 

  while( (*active_blk != NULL) || (active_spes != 0)){
    while(active_spes < NR_SPE && *active_blk != NULL){
      int spe = get_spe();
      (*active_blk)->job = SPU_PART;
      printf(" Running block %p \n", *active_blk);
      memcpy(spe_blocks[spe], *active_blk, sizeof(psc_cell_block_t));
      fprintf(stderr, "[ppe queing %#llx] start %p end %p\n", spe_id[spe],
	      spe_blocks[spe]->part_start, spe_blocks[spe]->part_end);
      msg = SPU_RUNJOB;
      spe_in_mbox_write(spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
      fprintf(stderr, "mass of first p %g \n", 
	      ((particle_cbe_t *)spe_blocks[spe]->part_start)->mni);
      active_blk++;
    }
    //    fprintf(stderr, " active spes: %d \n", active_spes);
    update_idle_spes();
    //    fprintf(stderr, " active spes: %d \n", active_spes);
    fflush(stderr);
    fflush(stdout);
  }
}  

void
cbe_push_part_2d(void)
{
  fields_cbe_t pf;
  particles_cbe_t pp;

  particles_cbe_get(&pp);
  fields_cbe_get(&pf, EX, EX+6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("cbe_part_yz", 1., 0, psc.pp.n_part * 9 * sizeof(cbe_real));
  }
  prof_start(pr);
  do_cbe_push_part_2d(&pp, &pf);
  prof_stop(pr);

  particles_cbe_put(&pp);
  fields_cbe_put(&pf,JXI,JXI+3);
}
