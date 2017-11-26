
#ifndef MRC_DOMAIN_H
#define MRC_DOMAIN_H

#include <mrc_common.h>
#include <mrc_obj.h>
#include <mrc_crds.h>

#include <mpi.h>

BEGIN_C_DECLS

enum {
  SW_0,
  SW_1,
  SW_2,
  MAX_SW,
};

#define MAX_BS (27)

enum {
  FACE_LEFT,
  FACE_RIGHT,
  FACE_BOTTOM,
  FACE_TOP,
  FACE_BACK,
  FACE_FRONT,
  NR_FACES,
};

// Btypes, needed visible by domain_mb. If these can't be made general,
// they need to go away.
enum {
  BTYPE_NONE,  // none specified, typically block-internal
  BTYPE_PATCH, // patch-internal boundary
  BTYPE_OUTER, // used by MB_CreateSimple(), also cylindrical / butterfly
  BTYPE_SP,    // cylindrical, butterfly
  BTYPE_SPC,   // butterfly
  BTYPE_USER = 10,
};


struct mrc_patch {
  int ldims[3];
  int off[3];
};

struct MB_pface {      // patch face info
  int pf_patch;
  int pf_face;
  int pf_map[3];
  int pf_btype;
};

struct mrc_patch_info {
  int rank;
  int patch;
  int global_patch;
  int level;
  int idx3[3];
  int ldims[3];
  int off[3];
  // stuff needed by mb
  // (could trim a bit)
  int p_bmx[3];            // size of parent block (cache)
  int p_ix[3];             // index of this patch in block
  int p_block;             // parent block
  struct MB_pface p_pface[NR_FACES];
};


struct MB_face {
  int block; // Which block do we exchange points for?
  int face;  // Which face on that block has my points?
  int map[3]; // Which dimension [x, y, z] on that face maps to each of my dimensions [0,1,2]?
  int btype; // What boundary condition should I use when exchanging data?
};

struct MB_block {
  int mx[3];
  int nr_block; // Sort of pointless to have this, but it's expected at the moment.
  struct MB_face faces[NR_FACES];  
  struct mrc_crds_gen *coord_gen[3]; // What generators are used for the 0,1,2 coordinates
  double xl[3]; // Lower bounds of this block in coordinate space
  double xh[3]; // Upper bounds of this block in coordinate space
};

MRC_CLASS_DECLARE(mrc_domain, struct mrc_domain);

void mrc_domain_get_global_dims(struct mrc_domain *domain, int *dims);
void mrc_domain_get_bc(struct mrc_domain *domain, int *bc);
void mrc_domain_get_nr_procs(struct mrc_domain *domain, int *nr_procs);
void mrc_domain_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches);
void mrc_domain_get_global_patch_info(struct mrc_domain *domain, int gpatch,
				      struct mrc_patch_info *info);
void mrc_domain_get_local_patch_info(struct mrc_domain *domain, int patch,
				     struct mrc_patch_info *info);
void mrc_domain_get_level_idx3_patch_info(struct mrc_domain *domain, int level, int idx[3],
					  struct mrc_patch_info *info);
void mrc_domain_get_nr_levels(struct mrc_domain *domain, int *p_nr_levels);
void mrc_domain_plot(struct mrc_domain *domain);
int  mrc_domain_get_neighbor_rank(struct mrc_domain *domain, int shift[3]);
void mrc_domain_get_neighbor_rank_patch(struct mrc_domain *domain, int p, int dir[3],
					int *nei_rank, int *nei_patch);
void mrc_domain_add_patch(struct mrc_domain *domain, int l, int idx3[3]);

int mrc_domain_nr_patches(struct mrc_domain *domain);
struct mrc_patch *mrc_domain_get_patches(struct mrc_domain *domain, int *nr_patches);
struct mrc_crds *mrc_domain_get_crds(struct mrc_domain *domain);
struct mrc_ddc *mrc_domain_get_ddc(struct mrc_domain *domain);

struct mrc_fld *mrc_domain_f1_create(struct mrc_domain *domain);
struct mrc_m3 *mrc_domain_m3_create(struct mrc_domain *domain);
struct mrc_fld *mrc_domain_m1_create(struct mrc_domain *domain);
struct mrc_fld *mrc_domain_fld_create(struct mrc_domain *domain, int sw, const char *comps);

struct mrc_ddc *mrc_domain_create_ddc(struct mrc_domain *domain);

END_C_DECLS

#endif
