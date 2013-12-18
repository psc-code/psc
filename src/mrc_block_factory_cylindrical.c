
#include "mrc_block_factory_private.h"
#include <mrc_domain.h> // will become just mrc_domain.h

#include <stdlib.h>

#define mrc_domain_mb(domain) mrc_to_subobj(domain, struct mrc_domain_mb)
// ======================================================================
// mrc_block_factory_cylindrical
//
// Constrcut the mapping needed for a cylindrical 2d trafo.
// 
// Diagram ahoy!
// 2 blocks
// Numbers map to each other, 'o' indicates an outer boundary
//
// BLOCK 0   BLOCK 1
//
// ---1---  ---3---
// |     |  |     |
// 2     o  o     2
// |     |  |     |
// ---3---  ---1---
//
//

// ----------------------------------------------------------------------
// mrc_block_factory_cylindrical_run

static void 
mrc_block_factory_cylindrical_run(struct mrc_block_factory *fac, struct mrc_domain *domain)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  
  int nr_blocks = 2;
  mb->nr_blocks = nr_blocks;
  struct MB_block *blocks = (struct MB_block *) calloc(nr_blocks, sizeof(struct MB_block));
  
  mb->mb_blocks = blocks;
  
  blocks[0] = (struct MB_block) {
    .mx = { gdims[0], gdims[1]/2, 1 },
    .faces = {
      [FACE_LEFT  ] = { 1, FACE_LEFT  , .map = { 0, MB_Y, MB_Z },
			.btype = BTYPE_SP    },
      [FACE_RIGHT ] = { .btype = BTYPE_OUTER },
      [FACE_BOTTOM] = { 1, FACE_TOP   , .map = { MB_X, 0, MB_Z } },
      [FACE_TOP   ] = { 1, FACE_BOTTOM, .map = { MB_X, 0, MB_Z } },
    },
  };
    
  blocks[1] = (struct MB_block) {
    .mx = { gdims[0], gdims[1]/2, 1 },
    .faces = {
      [FACE_LEFT  ] = { 0, FACE_LEFT  , .map = { 0, MB_Y, MB_Z },
			.btype = BTYPE_SP    },
      [FACE_RIGHT ] = { .btype = BTYPE_OUTER },
      [FACE_BOTTOM] = { 0, FACE_TOP   , .map = { MB_X, 0, MB_Z } },
      [FACE_TOP   ] = { 0, FACE_BOTTOM, .map = { MB_X, 0, MB_Z } },
    },
  };
}

// ----------------------------------------------------------------------
// mrc_block_factory subclass "cylindrical"

struct mrc_block_factory_ops mrc_block_factory_cylindrical = {
  .name        = "cylindrical",
  .run         = mrc_block_factory_cylindrical_run,
};

