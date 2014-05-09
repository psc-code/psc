
#include "mrc_block_factory_private.h"

#include <mrc_domain.h> // will become just mrc_domain.h
#include <mrc_crds_gen.h>

#include <stdlib.h>
#include <math.h>

#define mrc_domain_mb(domain) mrc_to_subobj(domain, struct mrc_domain_mb)
// ======================================================================
// mrc_block_factory_half_cylinder
//
// Constrcut the mapping needed for a cylinder with even symmetry
// 
// Diagram ahoy!
// 1 block
// Face 2 maps to itself. The y direction is periodic. Right face is outer bound. 
//
// BLOCK 0 
//
// ---1---
// |     |
// 2     o
// |     |
// ---1---
//
//

// ----------------------------------------------------------------------
// mrc_block_factory_hc_run
// This function commmented to the point of absurdity because it is used as a tutorial
// function in the documentation.

struct mrc_bf_hc {
  double rb, re; // Bounds in the radial (coord_gen_r) direction
  double zb, ze; // Bounds in the Z (coord_gen_z) direction
  int dims[2]; // FIXME: just 2D for now
  struct mrc_crds_gen *coord_gen[3]; // Coordinate generation objects 
                                     //for the 3 domain coordinates.
};

#define VAR(x) (void *)offsetof(struct mrc_bf_hc, x)
static struct param mrc_bf_hc_param_descr[] = {
  { "rb"              , VAR(rb)             , PARAM_DOUBLE(1e-15)     },
  { "re"              , VAR(re)             , PARAM_DOUBLE(1.0  )     },
  { "zb"              , VAR(zb)             , PARAM_DOUBLE(0)         },
  { "ze"              , VAR(ze)             , PARAM_DOUBLE(1.0  )     },
  { "mr"                , VAR(dims[0])          , PARAM_INT(16)          },
  { "mth"               , VAR(dims[1])          , PARAM_INT(16)          },
  { "coord_gen_r"       , VAR(coord_gen[0])    , MRC_VAR_OBJ(mrc_crds_gen)},
  { "coord_gen_th"      , VAR(coord_gen[1])    , MRC_VAR_OBJ(mrc_crds_gen)},
  { "coord_gen_z"       , VAR(coord_gen[2])    , MRC_VAR_OBJ(mrc_crds_gen)},
  {}
};
#undef VAR
// FIXME: Just a 2D plane for now..
//   { "mz"                , VAR(dims[2])          , PARAM_INT(1)           },

static void
_bf_hc_create(struct mrc_block_factory *fac)
{
  struct mrc_bf_hc *sub = mrc_to_subobj(fac, struct mrc_bf_hc);
  mrc_crds_gen_set_name(sub->coord_gen[0], "coord_gen_r");
  mrc_crds_gen_set_name(sub->coord_gen[1], "coord_gen_th");
  mrc_crds_gen_set_name(sub->coord_gen[2], "coord_gen_z");
}

// Unfortunately, this has to be here. The subclass create *cannot* be called on read
// since member objects don't get created.
static void
_mrc_block_factory_read(struct mrc_block_factory *fac, struct mrc_io *io)
{
  mrc_block_factory_read_member_objs(fac, io);
}


static void 
mrc_block_factory_hc_run(struct mrc_block_factory *fac, struct mrc_domain *domain)
{

  // Extract the domain sub-context
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  struct mrc_bf_hc *sub = mrc_to_subobj(fac, struct mrc_bf_hc);
  // Set the number of blocks, allocate space, and assign elements in the domain context.
  
  int nr_blocks = 1;
  mb->nr_blocks = nr_blocks;
  struct MB_block *blocks = (struct MB_block *) calloc(nr_blocks, sizeof(struct MB_block));
  mb->mb_blocks = blocks;
  
  // Construct the block mappings.

  // Block Structs have the following definition
  /////////////////
  //  struct MB_block {
  //    int mx[3]; ------------------------ number of points in each dim
  //    struct MB_face faces[NR_FACES]; --- Face mappings (6 faces total)
  //    int nr_block; --------------------- The number of this block
  //    char (*coord_names)[3]; // What coordinates the 0,1,2 dims are (used for crds_gen_parsing)
  //    double xl[3]; // Lower bounds of this block in coordinate space
  //   double xh[3]; // Upper bounds of this block in coordinate space
  //  };
  ///////

  // Each face mapping has the following structure
  ///////////////
  //  struct MB_face {
  //    int block; --------------- Which block do we exchange points for?
  //    int face; ---------------- Which face on that block has my points?
  //    int map[3]; Which dimension [x, y, z] on that face maps to each of my dimensions [0,1,2]?
  //    int btype;  What boundary condition should I use when exchanging data?
  //  };
  ////////////////////

  double pi = M_PI;
  // First block
  
  blocks[0] = (struct MB_block) {
    .nr_block = 0,
    .mx = { sub->dims[0], sub->dims[1], 1 },
    .faces = {
      // The Left face lies at the coordinate singularity in a circle, and has a special
      // boundary type.
      [FACE_LEFT  ] = { 0, FACE_LEFT  , .map = { 0, MB_Y, MB_Z },
			.btype = BTYPE_SP    },
      // The right face lies at the outer edge of cylindrical and doesn't exchange points with 
      // anyone
      [FACE_RIGHT ] = { .btype = BTYPE_OUTER },
      // The default boundary type is BTYPE_NONE, which indicates a simple ghost exchange
      [FACE_BOTTOM] = { 0, FACE_TOP   , .map = { MB_X, 0, MB_Z } },
      [FACE_TOP   ] = { 0, FACE_BOTTOM, .map = { MB_X, 0, MB_Z } },
    },
    // Assign which of the coordinate generators should be used for each local axis.
    .coord_gen = {sub->coord_gen[0] , sub->coord_gen[1], sub->coord_gen[2]},
    // Assign the coordinate boundaries for each local axis.
    .xl = { sub->rb,  0, 0 },
    .xh = { sub->re, pi, sub->ze },
  };
    
}


// ----------------------------------------------------------------------
// mrc_block_factory subclass "half_cylinder"

struct mrc_block_factory_ops mrc_block_factory_half_cylinder = {
  .name        = "half_cylinder",
  .size        = sizeof(struct mrc_bf_hc),
  .create      = _bf_hc_create,
  .read        = _mrc_block_factory_read,
  .param_descr = mrc_bf_hc_param_descr,
  .run         = mrc_block_factory_hc_run,
};

