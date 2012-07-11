
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

struct mrc_patch {
  int ldims[3];
  int off[3];
};

struct mrc_patch_info {
  int rank;
  int patch;
  int global_patch;
  int level;
  int idx3[3];
  int ldims[3];
  int off[3];
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
bool mrc_domain_is_setup(struct mrc_domain *domain);
void mrc_domain_add_patch(struct mrc_domain *domain, int l, int idx3[3]);

struct mrc_patch *mrc_domain_get_patches(struct mrc_domain *domain, int *nr_patches);
struct mrc_crds *mrc_domain_get_crds(struct mrc_domain *domain);
struct mrc_ddc *mrc_domain_get_ddc(struct mrc_domain *domain);

struct mrc_f1 *mrc_domain_f1_create(struct mrc_domain *domain);
struct mrc_f3 *mrc_domain_f3_create(struct mrc_domain *domain, int bnd);
struct mrc_m3 *mrc_domain_m3_create(struct mrc_domain *domain);
struct mrc_m1 *mrc_domain_m1_create(struct mrc_domain *domain);
struct mrc_a3 *mrc_domain_a3_create(struct mrc_domain *domain);

struct mrc_ddc *mrc_domain_create_ddc(struct mrc_domain *domain);

#endif
