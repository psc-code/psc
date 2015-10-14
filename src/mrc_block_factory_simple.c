
#include "mrc_block_factory_private.h"
#include <mrc_domain.h> // will become just mrc_domain.h
#include <mrc_crds_gen.h>

#include <stdlib.h>

#define mrc_domain_mb(domain) mrc_to_subobj(domain, struct mrc_domain_mb)
// ======================================================================
// mrc_block_factory_simple[2d, 3d]
//
// A simple 2D and 3D domains with one block.

// ----------------------------------------------------------------------
// mrc_block_factory_simple2d_run

struct mrc_bf_simple {
  double xb, xe;
  double yb, ye;
  double zb, ze;
  int dims[3];
  struct mrc_crds_gen *coord_gen[3];
};

static void
_bf_simple_create(struct mrc_block_factory *fac)
{
  struct mrc_bf_simple *sub = mrc_to_subobj(fac, struct mrc_bf_simple);
  mrc_crds_gen_set_name(sub->coord_gen[0], "coord_gen_x");
  mrc_crds_gen_set_name(sub->coord_gen[1], "coord_gen_y");
  mrc_crds_gen_set_name(sub->coord_gen[2], "coord_gen_z");
}

static void
_mrc_block_factory_read(struct mrc_block_factory *fac, struct mrc_io *io)
{
  mrc_block_factory_read_member_objs(fac, io);
}


#define VAR(x) (void *)offsetof(struct mrc_bf_simple, x)
static struct param mrc_bf_simple_param_descr[] = {
  { "xb"              , VAR(xb)             , PARAM_DOUBLE(0)         },
  { "xe"              , VAR(xe)             , PARAM_DOUBLE(1.0  )     },
  { "yb"              , VAR(yb)             , PARAM_DOUBLE(0)         },
  { "ye"              , VAR(ye)             , PARAM_DOUBLE(1.0  )     },
  { "zb"              , VAR(zb)             , PARAM_DOUBLE(0)         },
  { "ze"              , VAR(ze)             , PARAM_DOUBLE(1.0  )     },
  // Fixme: I'd rather use an int3 here, but apparently the fancy namex
  // appending only works on the command line
  { "mx"                , VAR(dims[0])          , PARAM_INT(16 )     },
  { "my"                , VAR(dims[1])          , PARAM_INT(16 )     },
  { "mz"                , VAR(dims[2])          , PARAM_INT(16 )     },
  { "coord_gen_x"       , VAR(coord_gen[0])    , MRC_VAR_OBJ(mrc_crds_gen)},
  { "coord_gen_y"       , VAR(coord_gen[1])    , MRC_VAR_OBJ(mrc_crds_gen)},
  { "coord_gen_z"       , VAR(coord_gen[2])    , MRC_VAR_OBJ(mrc_crds_gen)},
  {}
};
#undef VAR


// FIXME: This is kind of redundent, could be handled by the 3D version just fine.
static void 
mrc_block_factory_simple2d_run(struct mrc_block_factory *fac, struct mrc_domain *domain)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  struct mrc_bf_simple *sub = mrc_to_subobj(fac, struct mrc_bf_simple);

  int nr_blocks = 1;
  mb->nr_blocks = nr_blocks;
  struct MB_block *blocks = (struct MB_block *) calloc(nr_blocks, sizeof(struct MB_block));
  
  mb->mb_blocks = blocks;

  blocks[0] = (struct MB_block) {
    .mx = { sub->dims[0], sub->dims[1], 1 }, // Dimensions of the block
    .faces = { // Default boundary mapping for the faces
      [FACE_LEFT  ] = { .btype = BTYPE_OUTER },
      [FACE_RIGHT ] = { .btype = BTYPE_OUTER },
      [FACE_BOTTOM] = { .btype = BTYPE_OUTER },
      [FACE_TOP   ] = { .btype = BTYPE_OUTER },
    },
    .coord_gen = {sub->coord_gen[0] , sub->coord_gen[1], sub->coord_gen[2]},
    .xl = { sub->xb, sub->yb, sub->zb },
    .xh = { sub->xe, sub->ye, sub->ze },   
  };
  
  // check if any of our directions is periodic. If so, assign 
  // the appropriate mapping.

  if (domain->bc[0] == BC_PERIODIC) {
    blocks[0].faces[FACE_LEFT] = 
      (struct MB_face) { 0, FACE_RIGHT , .map = { 0, MB_Y, MB_Z } };
    blocks[0].faces[FACE_RIGHT] = 
      (struct MB_face) { 0, FACE_LEFT  , .map = { 0, MB_Y, MB_Z } };
  }

  if (domain->bc[1] == BC_PERIODIC) {
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
  .size        = sizeof(struct mrc_bf_simple),
  .create      = _bf_simple_create,
  .read        = _mrc_block_factory_read,
  .param_descr = mrc_bf_simple_param_descr,
  .run         = mrc_block_factory_simple2d_run,
};


// ----------------------------------------------------------------------
// mrc_block_factory_simple3d_run

static void 
mrc_block_factory_simple3d_run(struct mrc_block_factory *fac, struct mrc_domain *domain)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  struct mrc_bf_simple *sub = mrc_to_subobj(fac, struct mrc_bf_simple);
  
  int nr_blocks = 1;
  mb->nr_blocks = nr_blocks;
  struct MB_block *blocks = (struct MB_block *) calloc(nr_blocks, sizeof(struct MB_block));
  
  mb->mb_blocks = blocks;

   blocks[0] = (struct MB_block){
     .mx = { sub->dims[0], sub->dims[1], sub->dims[2] }, // Dimensions of the block
     .faces = {
       [FACE_LEFT  ] = { .btype = BTYPE_OUTER },
       [FACE_RIGHT ] = { .btype = BTYPE_OUTER },
       [FACE_BOTTOM] = { .btype = BTYPE_OUTER },
       [FACE_TOP   ] = { .btype = BTYPE_OUTER },
       [FACE_FRONT ] = { .btype = BTYPE_OUTER },
       [FACE_BACK  ] = { .btype = BTYPE_OUTER },
     },
     .coord_gen = {sub->coord_gen[0] , sub->coord_gen[1], sub->coord_gen[2]},
     .xl = { sub->xb, sub->yb, sub->zb },
     .xh = { sub->xe, sub->ye, sub->ze } ,   
   };

  if (domain->bc[0] == BC_PERIODIC) {
    blocks[0].faces[FACE_LEFT] = 
      (struct MB_face) { 0, FACE_RIGHT , .map = { 0, MB_Y, MB_Z } };
    blocks[0].faces[FACE_RIGHT] = 
      (struct MB_face) { 0, FACE_LEFT  , .map = { 0, MB_Y, MB_Z } };
  }
  if (domain->bc[1] == BC_PERIODIC) {
    blocks[0].faces[FACE_BOTTOM] = 
      (struct MB_face) { 0, FACE_TOP   , .map = { MB_X, 0, MB_Z } };
    blocks[0].faces[FACE_TOP] = 
      (struct MB_face) { 0, FACE_BOTTOM, .map = { MB_X, 0, MB_Z } };
  }
  if (domain->bc[2] == BC_PERIODIC) {
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
  .size        = sizeof(struct mrc_bf_simple),
  .create      = _bf_simple_create,
  .read             = _mrc_block_factory_read,
  .param_descr = mrc_bf_simple_param_descr,
  .run         = mrc_block_factory_simple3d_run,
};


