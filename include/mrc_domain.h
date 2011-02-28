
#ifndef MRC_DOMAIN_H
#define MRC_DOMAIN_H

#include <mrc_common.h>
#include <mrc_obj.h>
#include <mrc_crds.h>

#include <mpi.h>

enum {
  SW_0,
  SW_1,
  SW_2,
};

extern struct mrc_class mrc_class_mrc_domain;

MRC_OBJ_DEFINE_STANDARD_METHODS(mrc_domain, struct mrc_domain)
void mrc_domain_setup(struct mrc_domain *domain);
void mrc_domain_get_global_dims(struct mrc_domain *domain, int *dims);
void mrc_domain_get_local_idx(struct mrc_domain *domain, int *idx);
void mrc_domain_get_nr_procs(struct mrc_domain *domain, int *nr_procs);
void mrc_domain_get_bc(struct mrc_domain *domain, int *bc);
int  mrc_domain_get_neighbor_rank(struct mrc_domain *domain, int shift[3]);
bool mrc_domain_is_setup(struct mrc_domain *domain);
struct mrc_patch *mrc_domain_get_patches(struct mrc_domain *domain, int *nr_patches);
struct mrc_crds *mrc_domain_get_crds(struct mrc_domain *domain);

struct mrc_f3 *mrc_domain_f3_create(struct mrc_domain *domain, int bnd);
struct mrc_m3 *mrc_domain_m3_create(struct mrc_domain *domain, int bnd);

struct mrc_ddc_params;
struct mrc_ddc_ops;

struct mrc_ddc *
mrc_domain_create_ddc(struct mrc_domain *domain, struct mrc_ddc_params *params,
		      struct mrc_ddc_ops *ddc_ops);

struct mrc_patch {
  int off[3];
  int ldims[3];
};

#endif
