
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
  int is_setup;
  struct mrc_crds *crds;
  struct mrc_ddc *ddc;
};

struct mrc_domain_ops {
  MRC_SUBCLASS_OPS(struct mrc_domain);
  struct mrc_patch *(*get_patches)(struct mrc_domain *domain, int *nr_patches);
  void (*get_global_dims)(struct mrc_domain *domain, int *dims);
  void (*get_nr_procs)(struct mrc_domain *domain, int *nr_procs);
  void (*get_bc)(struct mrc_domain *domain, int *bc);
  int  (*get_neighbor_rank)(struct mrc_domain *, int shift[3]);
  void (*get_nr_global_patches)(struct mrc_domain *domain, int *nr_global_patches);
  void (*get_global_patch_info)(struct mrc_domain *domain, int gp,
				struct mrc_patch_info *info);
  void (*get_local_patch_info)(struct mrc_domain *domain, int p,
			       struct mrc_patch_info *info);
  void (*get_idx3_patch_info)(struct mrc_domain *domain, int idx[3],
			      struct mrc_patch_info *info);
  void (*plot)(struct mrc_domain *domain);
  struct mrc_ddc *(*create_ddc)(struct mrc_domain *);
  int* (*get_offset)(struct mrc_domain* domain);
};

// ======================================================================
// mrc_domain_simple

struct mrc_domain_simple {
  int gdims[3];
  struct mrc_patch patch;
  int nr_procs[3];
  int bc[3];

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
};

void sfc_setup(struct mrc_sfc *sfc, int *np);
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
  int bc[3];
  int *gpatch_off_all;
  int *ldims[3]; //< patch dims for all patches by direction
  int *off[3]; //< offsets for all patches by direction
  struct mrc_sfc sfc;

  // map
  struct bintree g_patches;	//Provides a mapping gpatch -> gpatchinfo / patches
  int* gp;	//Maps [0..nr_gpatches] -> gpatch

  struct bitfield3d* p_activepatches;	//Only used as a parameter. Will be invalid after setup()
  struct bitfield3d activepatches; 	//Index of which patches are active
};

extern struct mrc_domain_ops mrc_domain_multi_ops;
extern struct mrc_domain_ops mrc_domain_dynamic_ops;

#endif
