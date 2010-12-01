#include <string.h>
#include "psc_ppu.h"

#ifdef CELLEMU
#include "libspe2_c.h"
#else
#include <libspe2.h>
#endif

void
cbe_push_part_2d(void)
{
  
  unsigned int msg; 
  psc_cell_block_t ** active; 

  active = block_list; 

  while( (*active != NULL) || (active_spes != 0)){
    for(int spe = 0; spe < NR_SPE; spe++){
      spe_out_mbox_read(spe_id[spe], &msg,1);
      //      printf("Looking at spe %d \n", spe);
      if((msg == SPE_IDLE) && (active != NULL)){
	(*active)->job = SPU_HELLO;
	memcpy(spe_blocks[spe], *active, sizeof(psc_cell_block_t));
	msg = SPU_RUNJOB;
	spe_in_mbox_write(spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
	active++;
      }
    }
  }
  
  active = block_list; 

  while( (*active != NULL) || (active_spes != 0)){
    for(int spe = 0; spe < NR_SPE; spe++){
      spe_out_mbox_read(spe_id[spe], &msg,1);
      if((msg == SPE_IDLE) && (active != NULL)){
	(*active)->job = SPU_BYE;
	memcpy(spe_blocks[spe], *active, sizeof(psc_cell_block_t));
	msg = SPU_RUNJOB;
	spe_in_mbox_write(spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
	active++;
      }
    }
  }


}
