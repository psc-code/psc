
#ifndef MRCTEST_H
#define MRCTEST_H

#include <mrc_domain.h>
#include <mrc_mod.h>

struct mrc_fld_init_values_entry {
  int m;
  float (*ini)(float xx, float yy, float zz);
};

struct mrc_fld_init_values_info {
  struct mrc_fld_init_values_entry *ini_flds;
};

void mrc_fld_init_values(struct mrc_fld *fld, struct mrc_fld_init_values_info *iv_info);

void mrctest_init(int *argc, char ***argv);
void mrctest_finalize();

// ----------------------------------------------------------------------

void mrctest_fld_compare(struct mrc_fld *fld1, struct mrc_fld *fld2, float eps);
void mrctest_m1_compare(struct mrc_fld *fld1, struct mrc_fld *fld2, float eps);
void mrctest_m3_compare(struct mrc_m3 *m3_1, struct mrc_m3 *m3_2);
void mrctest_crds_compare(struct mrc_crds *crds1, struct mrc_crds *crds2);

// ----------------------------------------------------------------------
// mrctest_domain

struct mrctest_domain_params {
  int gdims[3];
  int nproc[3];
  bool use_diagsrv;
};

void mrctest_domain_init(struct mrctest_domain_params *par);
struct mrc_domain *mrctest_create_domain(MPI_Comm comm, struct mrctest_domain_params *par);
struct mrc_domain *mrctest_create_domain_rectilinear(MPI_Comm comm,
						     struct mrctest_domain_params *par);
void mrctest_domain_init_values_0(struct mrc_fld *f);
void mrctest_domain(void (*mod_domain)(struct mrc_mod *mod, void *arg));
void mrctest_set_crds_rectilinear_1(struct mrc_domain *domain);
struct mrc_fld *mrctest_create_field_1(struct mrc_domain *domain);
struct mrc_fld *mrctest_create_field_2(struct mrc_domain *domain);
struct mrc_fld *mrctest_create_m1_1(struct mrc_domain *domain, int dim);

// ----------------------------------------------------------------------
// mrctest_set_amr_domain

void mrctest_set_amr_domain_0(struct mrc_domain *domain);
void mrctest_set_amr_domain_1(struct mrc_domain *domain);
void mrctest_set_amr_domain_2(struct mrc_domain *domain);
void mrctest_set_amr_domain_3(struct mrc_domain *domain);
void mrctest_set_amr_domain_4(struct mrc_domain *domain);

#endif
