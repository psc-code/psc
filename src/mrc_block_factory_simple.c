
#include "mrc_block_factory_private.h"
#include <mrc_domain.h> // will become just mrc_domain.h

#include <stdlib.h>

#define mrc_domain_mb(domain) mrc_to_subobj(domain, struct mrc_domain_mb)
// ======================================================================
// mrc_block_factory_simple[2d, 3d]
//
// A simple 2D and 3D domains with one block.

// ----------------------------------------------------------------------
// mrc_block_factory_simple2d_run

static void 
mrc_block_factory_simple2d_run(struct mrc_block_factory *fac, struct mrc_domain *domain)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  
  int nr_blocks = 1;
  mb->nr_blocks = nr_blocks;
  struct MB_block *blocks = (struct MB_block *) calloc(nr_blocks, sizeof(struct MB_block));
  
  mb->mb_blocks = blocks;

  assert(gdims[2] == 1);
  blocks[0] = (struct MB_block) {
    .mx = { gdims[0], gdims[1], 1 }, // Dimensions of the block
    .faces = { // Default boundary mapping for the faces
      [FACE_LEFT  ] = { .btype = BTYPE_OUTER },
      [FACE_RIGHT ] = { .btype = BTYPE_OUTER },
      [FACE_BOTTOM] = { .btype = BTYPE_OUTER },
      [FACE_TOP   ] = { .btype = BTYPE_OUTER },
    },
  };
  
  // check if any of our directions is periodic. If so, assign 
  // the appropriate mapping.

  if (mb->bc[0] == BC_PERIODIC) {
    blocks[0].faces[FACE_LEFT] = 
      (struct MB_face) { 0, FACE_RIGHT , .map = { 0, MB_Y, MB_Z } };
    blocks[0].faces[FACE_RIGHT] = 
      (struct MB_face) { 0, FACE_LEFT  , .map = { 0, MB_Y, MB_Z } };
  }

  if (mb->bc[1] == BC_PERIODIC) {
    blocks[0].faces[FACE_BOTTOM] = 
      (struct MB_face) { 0, FACE_TOP   , .map = { MB_X, 0, MB_Z } };
    blocks[0].faces[FACE_TOP] = 
      (struct MB_face) { 0, FACE_BOTTOM, .map = { MB_X, 0, MB_Z } };
  }

}

// ----------------------------------------------------------------------
// mrc_block_factory subclass "simple2d"

struct mrc_block_factory_ops mrc_block_factory_simple2d = {
  .name        = "simple2d",
  .run         = mrc_block_factory_simple2d_run,
};


// ----------------------------------------------------------------------
// mrc_block_factory_simple3d_run

static void 
mrc_block_factory_simple3d_run(struct mrc_block_factory *fac, struct mrc_domain *domain)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  
  int nr_blocks = 1;
  mb->nr_blocks = nr_blocks;
  struct MB_block *blocks = (struct MB_block *) calloc(nr_blocks, sizeof(struct MB_block));
  
  mb->mb_blocks = blocks;

   blocks[0] = (struct MB_block){
     .mx = { gdims[0], gdims[1], gdims[2] },
     .faces = {
       [FACE_LEFT  ] = { .btype = BTYPE_OUTER },
       [FACE_RIGHT ] = { .btype = BTYPE_OUTER },
       [FACE_BOTTOM] = { .btype = BTYPE_OUTER },
       [FACE_TOP   ] = { .btype = BTYPE_OUTER },
       [FACE_FRONT ] = { .btype = BTYPE_OUTER },
       [FACE_BACK  ] = { .btype = BTYPE_OUTER },
     },
   };

  if (mb->bc[0] == BC_PERIODIC) {
    blocks[0].faces[FACE_LEFT] = 
      (struct MB_face) { 0, FACE_RIGHT , .map = { 0, MB_Y, MB_Z } };
    blocks[0].faces[FACE_RIGHT] = 
      (struct MB_face) { 0, FACE_LEFT  , .map = { 0, MB_Y, MB_Z } };
  }
  if (mb->bc[1] == BC_PERIODIC) {
    blocks[0].faces[FACE_BOTTOM] = 
      (struct MB_face) { 0, FACE_TOP   , .map = { MB_X, 0, MB_Z } };
    blocks[0].faces[FACE_TOP] = 
      (struct MB_face) { 0, FACE_BOTTOM, .map = { MB_X, 0, MB_Z } };
  }
  if (mb->bc[2] == BC_PERIODIC) {
    blocks[0].faces[FACE_FRONT] = 
      (struct MB_face) { 0, FACE_BACK  , .map = { MB_X, MB_Y, 0 } };
    blocks[0].faces[FACE_BACK] = 
      (struct MB_face) { 0, FACE_FRONT , .map = { MB_X, MB_Y, 0 } };
  }
}

// ----------------------------------------------------------------------
// mrc_block_factory subclass "simple3d"

struct mrc_block_factory_ops mrc_block_factory_simple3d = {
  .name        = "simple3d",
  .run         = mrc_block_factory_simple3d_run,
};


