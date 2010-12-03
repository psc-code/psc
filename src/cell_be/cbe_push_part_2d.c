#include <string.h>
#include "psc_ppu.h"

#ifdef CELLEMU
#include "libspe2_c.h"
const unsigned int all_spes = 1;
#else
#include <libspe2.h>
const unsigned int all_spes = 255;
#endif

void
cbe_push_part_2d(void)
{
  
  unsigned int msg; 
  int state;
  psc_cell_block_t ** active_blk; 

  active_blk = block_list; 

  while( (*active_blk != NULL) || (active_spes != 0)){
    while(active_spes < NR_SPE && *active_blk != NULL){
      HERE;
      int spe = get_spe();
      (*active_blk)->job = SPU_HELLO;
      memcpy(spe_blocks[spe], *active_blk, sizeof(psc_cell_block_t));
      msg = SPU_RUNJOB;
      spe_in_mbox_write(spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
      active_blk++;
    }
    
    update_idle_spes();
    
  }
  
  active_blk = block_list; 
  
  while( (*active_blk != NULL) || (active_spes != 0)){
    while(active_spes < NR_SPE && *active_blk != NULL){
      HERE;
      int spe = get_spe();
      (*active_blk)->job = SPU_BYE;
      memcpy(spe_blocks[spe], *active_blk, sizeof(psc_cell_block_t));
      msg = SPU_RUNJOB;
      spe_in_mbox_write(spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
      active_blk++;
    }
    
    update_idle_spes();    
  }  
}
