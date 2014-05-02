
#include "mrc_block_factory_private.h"

#include <mrc_domain.h> // will become just mrc_domain.h
#include <mrc_params.h>
#include <mrc_io.h>
#include <string.h>
#include <mrc_bits.h>

#define mrc_domain_mb(domain) mrc_to_subobj(domain, struct mrc_domain_mb)
#define pfb
#define pfr return(0)
#define CE assert(ierr == 0)

// ======================================================================
// mrc_block_factory_run
//
// Takes a multiblock domain and subdivides it into blocks based
// on the factory mapping, then creates the block mapping
// members of the domain and assigns them.

static int
check_blocks(struct mrc_domain_mb *mb)
{
  pfb;
  for (int n = 0; n < mb->nr_blocks; n++) {
    for (int f = 0; f < NR_FACES; f++) {
      struct MB_block *block = &mb->mb_blocks[n];
      struct MB_face *face = &block->faces[f];
      if (face->block < 0) {
	continue;
      }
      assert(face->block < mb->nr_blocks);
      assert(face->face  < NR_FACES);
      struct MB_block *nblock = &mb->mb_blocks[face->block];
      struct MB_face *nface = &nblock->faces[face->face];
      assert(nface->block == n);
      assert(nface->face  == f);
      int dir = face2dir(f);
      assert(face->map[dir] == 0);
      for (int i = 0; i < 3; i++) {
	if (i == dir)
	  continue;

	int ndir = map2dir(face->map[i]);
	assert(map2dir(nface->map[ndir]) == i);
      }
    }
  }
  pfr;
}


void
mrc_block_factory_run(struct mrc_block_factory *fac, struct mrc_domain *domain)
{
  assert(strcmp(mrc_domain_type(domain), "mb") == 0);
  struct mrc_block_factory_ops *ops = mrc_block_factory_ops(fac);
  assert(ops && ops->run);
  ops->run(fac, domain);
  
  // Some housekeeping that used to be done in mb_set_params
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  for (int n = 0; n < mb->nr_blocks; n++) {
    mb->mb_blocks[n].nr_block = n;
    for (int f = 0; f < NR_FACES; f++) {
      int *map = mb->mb_blocks[n].faces[f].map;
      // Not entirely sure what this does, but it was in the original
      if (map[0]+map[1]+map[2] == 0)
	mb->mb_blocks[n].faces[f].block = -1;
    }
  }
  int ierr = check_blocks(mb); CE;
  for (int d = 0; d < 3; d++) {
    // this may not be needed now...
    mb->ppb[d] = MAX(1, mb->ppb[d]);
  }

}

// ----------------------------------------------------------------------
// mrc_block_factory_init

static void
mrc_block_factory_init(void)
{
  mrc_class_register_subclass(&mrc_class_mrc_block_factory, &mrc_block_factory_simple2d);
  mrc_class_register_subclass(&mrc_class_mrc_block_factory, &mrc_block_factory_simple3d);
  mrc_class_register_subclass(&mrc_class_mrc_block_factory, &mrc_block_factory_cylindrical);
  mrc_class_register_subclass(&mrc_class_mrc_block_factory, &mrc_block_factory_half_cylinder);
}

// ----------------------------------------------------------------------
// mrc_block_factory class

struct mrc_class_mrc_block_factory mrc_class_mrc_block_factory = {
  .name             = "mrc_block_factory",
  .size             = sizeof(struct mrc_block_factory),
  .init             = mrc_block_factory_init,
};
