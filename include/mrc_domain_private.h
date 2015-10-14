
#ifndef MRC_DOMAIN_PRIVATE_H
#define MRC_DOMAIN_PRIVATE_H

#include <mrc_domain.h>

#include <bitfield3d.h>
#include <assert.h>

// ======================================================================
// mrc_domain

struct mrc_domain {
  struct mrc_obj obj;
  int rank, size;
  struct mrc_crds *crds;
  struct mrc_ddc *ddc;
  
  int bc[3];
};

struct mrc_domain_ops {
  MRC_SUBCLASS_OPS(struct mrc_domain);
  struct mrc_patch *(*get_patches)(struct mrc_domain *domain, int *nr_patches);
  void (*get_global_dims)(struct mrc_domain *domain, int *dims);
  void (*get_nr_procs)(struct mrc_domain *domain, int *nr_procs);
  int  (*get_neighbor_rank)(struct mrc_domain *, int shift[3]);
  void (*get_neighbor_rank_patch)(struct mrc_domain *domain, int p, int dir[3],
				  int *nei_rank, int *nei_patch);
  void (*get_nr_global_patches)(struct mrc_domain *domain, int *nr_global_patches);
  void (*get_global_patch_info)(struct mrc_domain *domain, int gp,
				struct mrc_patch_info *info);
  void (*get_local_patch_info)(struct mrc_domain *domain, int p,
			       struct mrc_patch_info *info);
  void (*get_level_idx3_patch_info)(struct mrc_domain *domain, int level, int idx[3],
				    struct mrc_patch_info *info);
  void (*get_nr_levels)(struct mrc_domain *domain, int *p_nr_levels);
  void (*plot)(struct mrc_domain *domain);
  struct mrc_ddc *(*create_ddc)(struct mrc_domain *);
  int* (*get_offset)(struct mrc_domain* domain);
  void (*add_patch)(struct mrc_domain *domain, int l, int idx3[3]);
};

// ======================================================================
// mrc_domain_simple

struct mrc_domain_simple {
  int gdims[3];
  struct mrc_patch patch;
  int nr_procs[3];

  int proc[3];
};

extern struct mrc_domain_ops mrc_domain_simple_ops;

// ======================================================================
// mrc_sfc

enum {
  CURVE_BYDIM,
  CURVE_MORTON,
  CURVE_HILBERT,
};

struct mrc_sfc {
  int curve_type; //< type of space filling curve

  // for sfc simple
  int np[3];

  // for sfc morton, hilbert
  int nbits[3];
  int nbits_max;

  // for sfc hilbert
  int hilbert_nr_dims;
  int hilbert_dim[3];
  int *hilbert_map_h_to_p;
  int *hilbert_map_p_to_h;
};

void sfc_setup(struct mrc_sfc *sfc, int *np);
void sfc_destroy(struct mrc_sfc *sfc);
int  sfc_idx3_to_idx(struct mrc_sfc *sfc, const int p[3]);
void sfc_idx_to_idx3(struct mrc_sfc *sfc, int idx, int p[3]);

// ======================================================================
// mrc_domain_multi

#include <bintree.h>

struct mrc_domain_multi {
  int nr_global_patches;	//Number of global patches
  int gdims[3];
  int nr_patches;
  int gpatch_off; //< global patch # on this proc is gpatch_off..gpatch_off+nr_patches
  struct mrc_patch *patches;
  int np[3]; //< # of patches per direction
  int *gpatch_off_all;
  int *ldims[3]; //< patch dims for all patches by direction
  int *off[3]; //< offsets for all patches by direction
  struct mrc_sfc sfc;

  // map
  struct bintree g_patches;	//Provides a mapping gpatch -> gpatchinfo / patches
  int* gp;	//Maps [0..nr_gpatches] -> gpatch

  bool have_activepatches; //Otherwise, we're assuming all patches are active
  struct bitfield3d* p_activepatches;	//Only used as a parameter. Will be invalid after setup()
  struct bitfield3d activepatches; 	//Index of which patches are active
};

extern struct mrc_domain_ops mrc_domain_multi_ops;

// ======================================================================
// mrc_domain_amr

struct mrc_amr_patch {
  int l;
  int idx3[3];
  list_t entry;
};

struct mrc_domain_amr {
  int nr_global_patches;	//< Number of global patches
  int nr_levels;                //< Number of AMR levels
  int gdims[3];                 //< Patch dimensions == global dims on level 0
  int nr_patches;
  int gpatch_off; //< global patch # on this proc is gpatch_off..gpatch_off+nr_patches
  struct mrc_patch *patches;
  int *gpatch_off_all;
  struct mrc_sfc sfc;

  list_t global_patch_list;     //< List of patches added, only before setup()

  // map
  struct mrc_amr_level_sfc_idx {
    int l, sfc_idx;
  } *map_gpatch_to_sfc;	//Maps [0..nr_gpatches] -> l, sfc_idx
};

#define mrc_domain_amr(domain) mrc_to_subobj(domain, struct mrc_domain_amr)

extern struct mrc_domain_ops mrc_domain_amr_ops;

// ======================================================================
// mrc_domain_mb (multi-block)

enum {
  MB_NONE,
  MB_X, 
  MB_Y,
  MB_Z,
  MB_R  = 4,
  MB_XR = MB_X | MB_R,
  MB_YR = MB_Y | MB_R,
  MB_ZR = MB_Z | MB_R,
};

  
// These are too clever for their own good: Look at value of face enums
// to see why they works.

// ----------------------------------------------------------------------
// return dir (0,1,2) perpendicular to given face

static inline int
face2dir(int face)
{
  return face / 2;
}

// ----------------------------------------------------------------------
// return 1 for boundary at positive side, 0 at zero side

static inline int
face2bnd(int face)
{
  return face % 2;
}

static inline int
map2dir(int map)
{
  return (map & 3) - 1;
}

// ----------------------------------------------------------------------


struct mrc_domain_mb { 
  // unique block related bits
  struct mrc_block_factory * _block_factory;
  int ppb[3]; // how to subdivide blocks into patches
  int mb_ppb_total; // total patches per block
  int nr_blocks;
  struct MB_block *mb_blocks;

  // Patch related bits: similar to multi domain, but with a little more info
  int nr_global_patches;         //Number of patches on all procs
  struct mrc_patch *patches;      // Array of global patches
  struct mrc_patch_info *patch_info; // bunch of pre-calculated patch information

  int nr_patches;                // Number of local patches
  int gpatch_off;                // Offset of this procs patches in global array
                                 // (takes the place of loc_patches list)

  // FIXME: These don't really belong in the domain.
  // The idea of carrying around values for each SW bothers me, let's see if
  // we can do without it
  /*
  int mb_N[MAX_SW];
  int mb_loc_N[MAX_SW];
  int mb_loc_start; // offset of local part into global vector
  */
};

extern struct mrc_domain_ops mrc_domain_mb_ops;

#define MB_foreach_block(mb, block)					\
  for (struct MB_block *block = mrc_domain_mb(mb)->mb_blocks;		\
       block != mrc_domain_mb(mb)->mb_blocks + mrc_domain_mb(mb)->nr_blocks; block++) { \
    struct MB_block __attribute__((unused)) __mrc_block = *block;	\
    
#define MB_foreach_block_end }

    
#define MB_foreach_patch(mb, p)					\
  {								\
  for (int p = 0; p < mrc_domain_mb(mb)->nr_patches; p++)	\
    {

#define MB_foreach_patch_end }}

#define mrc_patch_foreach(mb, patch, jx, jy, jz, l, r)			\
  struct mrc_patch_info _info;						\
  mrc_domain_get_local_patch_info(mb, patch, &_info);			\
  for (int jz = -(l); jz < _info.ldims[2]+(r); jz++) {			\
    for (int jy = -(l); jy < _info.ldims[1]+(r); jy++) {			\
      for (int jx = -(l); jx < _info.ldims[0]+(r); jx++)			\

#define mrc_patch_foreach_2d(mb, patch, jx, jy, jz, l, r)			\
  struct mrc_patch_info _info;						\
  mrc_domain_get_local_patch_info(mb, patch, &_info);			\
  for (int jz = 0; jz < _info.ldims[2]; jz++) {				\
    for (int jy = -(l); jy < _info.ldims[1]+(r); jy++) {			\
      for (int jx = -(l); jx < _info.ldims[0]+(r); jx++)			\

#define mrc_patch_foreach_end }}


#endif
