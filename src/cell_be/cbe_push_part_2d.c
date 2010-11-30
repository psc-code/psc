#include <string.h>

void
cbe_push_part_2d(void)
{
  block_ll_create();
  
  unsigned int msg; 
  while( (staged_blocks->next != NULL) && (active_spes != 0)){
    for(int spe = 0; spe < NR_SPES; spe++){
      spe_out_mbox_read(spe_id[spe], &msg,1);
      if((msg == SPU_IDLE) && (staged_blocks->next != NULL)){
	block_node_t * todo = block_ll_pop(&staged_blocks);
	todo->block->job = SPU_HELLO;
	memcpy(spe_blocks[i], todo->block, sizeof(psc_cell_block_t));
	msg = SPU_RUNJOB;
	spe_in_mbox_write(spe_id[spe], &msg, 1, SPE_MBOX_ANY_NONBLOCKING);
	block_ll_push(todo, &finished_blocks);
      }
    }
  }


}
