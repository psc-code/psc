#include "psc.h"
#include "psc_ppu.h"



void cbe_setup_layout(void)
{
  assert(!spu_ctl.layout);
  spu_ctl.layout = calloc(1, sizeof(struct cbe_block_layout));
  
  assert(psc.img[0] == 7);
  assert(psc.img[1] == 110);
  assert(psc.img[2] == 110);

  struct cbe_block_layout *layout = spu_ctl.layout;  

  layout.nblocks = 16;

  layout.block_size[0] = 1;
  layout.block_size[1] = 26;
  layout.block_size[2] = 26;

  layout.block_grid[0] = 1;
  layout.block_grid[1] = 4;
  layout.block_grid[2] = 4;
}


void 
cbe_blocks_create(void)
{

  assert(!spu_ctl.blocks_inited);
  
  if(!spu_ctl.layout) 
    cbe_setup_layout();

  int nblocks = spu_ctl.nblocks; 

  spu_ctl.block_list = calloc(nblocks+1, sizeof(psc_cell_block_t*));

  psc_cell_block_t **curr_block = spu_ctl.block_list;

  for(int i = 0; i < nblocks; i++){
    *curr_block = calloc(1,sizeof(psc_cell_block_t));
    assert(*curr_block != NULL);
    curr_block++;
  }
  
  *curr_block = NULL;

  spu_ctl.blocks_inited = 1; 

}

void 
cbe_blocks_destroy(void)
{

  assert(spu_ctl.blocks_inited); 

  psc_cell_block_t ** active; 
  active = spu_ctl.block_list; 

  while(*active != NULL){
    free(*active);
    *active == NULL; 
    active++;
  }
}

void 
cbe_assign_parts_to_blocks(particles_cbe_t * pp)
{
  if(! spu_ctl.layout)
    cbe_setup_layout();
  
  if(! spu_ctl.blocks_inited)
    cbe_blocks_create();

  if(! spu_ctl.particles_coarse_sorted)
    psc_sort();

  unsigned int *cnts = spu_ctl.cnts;

  particle_cbe_t *fp = pp->particles; 

  psc_cell_block_t ** curr = spu_ctl.block_list;

  fprintf(stderr, "nblocks %d fp: %p \n", spu_ctl.nblocks, fp);

  (*curr)->part_start =  fp;
  (*curr)->part_end = (fp + cnts[0]);
  fprintf(stderr, "[0] cnts[0] %d\n", cnts[0]);
  fprintf(stderr, "[0] start %p end %p\n", 0, (*curr)->part_start, (*curr)->part_end);
  curr++;
  
  for(int i = 1; i < spu_ctl.nblocks; i++){
    fprintf(stderr, "[%d] ctns[%d - 1] %d cnts[%d] %d\n", i, i, cnts[i-1], i, cnts[i]);
    (*curr)->part_start = (fp + cnts[i-1]);
    (*curr)->part_end = (fp + cnts[i]);
    fprintf(stderr, "[%d] start %p end %p\n", i, (*curr)->part_start, (*curr)->part_end);
    curr++;
  }
  
}

