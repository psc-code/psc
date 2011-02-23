
#ifndef MRC_DOMAIN_PRIVATE_H
#define MRC_DOMAIN_PRIVATE_H

#include <mrc_domain.h>

#include <assert.h>

// ======================================================================
// mrc_domain

struct mrc_domain {
  struct mrc_obj obj;
  int rank, size;
  int is_setup;
  struct mrc_crds *crds;
};

struct mrc_domain_ops {
  MRC_OBJ_OPS;
  struct mrc_patch *(*get_patches)(struct mrc_domain *domain, int *nr_patches);
  void (*get_local_idx)(struct mrc_domain *domain, int *idx);
  void (*get_global_dims)(struct mrc_domain *domain, int *dims);
  void (*get_nr_procs)(struct mrc_domain *domain, int *nr_procs);
  void (*get_bc)(struct mrc_domain *domain, int *bc);
  int  (*get_neighbor_rank)(struct mrc_domain *, int shift[3]);
  struct mrc_ddc * (*create_ddc)(struct mrc_domain *, struct mrc_ddc_params *ddc_par,
				 struct mrc_ddc_ops *ddc_ops);
};

void libmrc_domain_register(struct mrc_domain_ops *ops);

// ======================================================================
// mrc_domain_simple

struct mrc_domain_simple {
  int gdims[3];
  struct mrc_patch patch;
  int nr_procs[3];
  int bc[3];

  int proc[3];
};

void libmrc_domain_register_simple(void);

// ======================================================================
// mrc_domain_multi

struct mrc_domain_multi {
  int gdims[3];
  struct mrc_patch patch;
  int nr_procs[3];
  int bc[3];

  int proc[3];
};

void libmrc_domain_register_multi(void);

// ======================================================================

static inline struct mrc_domain *
to_mrc_domain(struct mrc_obj *obj)
{
  assert(obj->class == &mrc_class_mrc_domain);
  return container_of(obj, struct mrc_domain, obj);
}


#endif

