
#include "mrc_domain_private.h"
#include "mrc_params.h"
#include "mrc_ddc.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define TAG_SCAN_OFF (1000)

static inline struct mrc_domain_multi *
mrc_domain_multi(struct mrc_domain *domain)
{
  return domain->obj.subctx;
}

// ======================================================================

static int
idx3_to_gpatch(struct mrc_domain_multi *multi, const int p[3])
{
  int *np = multi->np;
  return (p[2] * np[1] + p[1]) * np[0] + p[0];
}

static void
gpatch_to_idx3(struct mrc_domain_multi *multi, int gpatch, int p[3])
{
  int *np = multi->np;
  p[0] = gpatch % np[0]; gpatch /= np[0];
  p[1] = gpatch % np[1]; gpatch /= np[1];
  p[2] = gpatch;
}

// ======================================================================

static void
get_nr_patches_gpatch_off(struct mrc_domain *domain, int rank,
			  int *nr_patches, int *gpatch_off)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int nr_global_patches = multi->np[0] * multi->np[1] * multi->np[2];
  int patches_per_proc = nr_global_patches / domain->size;
  int patches_per_proc_rmndr = nr_global_patches % domain->size;

  *nr_patches = patches_per_proc + (rank < patches_per_proc_rmndr);
  if (rank < patches_per_proc_rmndr) {
    *gpatch_off = rank * (patches_per_proc + 1);
  } else {
    *gpatch_off = rank * patches_per_proc + patches_per_proc_rmndr;
  }
}

static void
gpatch_to_rank_patch(struct mrc_domain *domain, int gpatch,
		     int *rank, int *patch)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int nr_global_patches = multi->np[0] * multi->np[1] * multi->np[2];
  int patches_per_proc = nr_global_patches / domain->size;
  int patches_per_proc_rmndr = nr_global_patches % domain->size;
  
  if (gpatch < (patches_per_proc + 1) * patches_per_proc_rmndr) {
    *rank = gpatch / (patches_per_proc + 1);
    *patch = gpatch % (patches_per_proc + 1);
  } else {
    int tmp = gpatch - (patches_per_proc + 1) * patches_per_proc_rmndr;
    *rank = tmp / patches_per_proc + patches_per_proc_rmndr;
    *patch = tmp % patches_per_proc;
  }
}

// ======================================================================

static void
mrc_domain_multi_view(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int proc = 0; proc < domain->size; proc++) {
    if (domain->rank == proc) {
      for (int p = 0; p < multi->nr_patches; p++) {
	struct mrc_patch *patch = &multi->patches[p];
	mprintf("  patch: ldims %dx%dx%d off %dx%dx%d\n",
		patch->ldims[0], patch->ldims[1], patch->ldims[2],
	      patch->off[0], patch->off[1], patch->off[2]);
      }
    }
    MPI_Barrier(obj->comm);
  }
}

// FIXME, get rid of this one always use idx3 one?
static void
mrc_domain_multi_get_global_patch_info(struct mrc_domain *domain, int gpatch,
				       struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  gpatch_to_rank_patch(domain, gpatch, &info->rank, &info->patch);

  int p3[3];
  gpatch_to_idx3(multi, gpatch, p3);
  for (int d = 0; d < 3; d++) {
    info->ldims[d] = multi->ldims[d][p3[d]];
    info->off[d] = multi->off[d][p3[d]];
  }
}

static void
mrc_domain_multi_setup(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);
  assert(!domain->is_setup);
  domain->is_setup = true;

  MPI_Comm_rank(domain->obj.comm, &domain->rank);
  MPI_Comm_size(domain->obj.comm, &domain->size);

  // FIXME: allow setting of desired decomposition by user?
  int *np = multi->np;
  for (int d = 0; d < 3; d++) {
    int ldims[3], rmndr[3];
    ldims[d] = multi->gdims[d] / np[d];
    rmndr[d] = multi->gdims[d] % np[d];

    multi->ldims[d] = calloc(np[d], sizeof(*multi->ldims[d]));
    multi->off[d] = calloc(np[d], sizeof(*multi->off[d]));
    for (int i = 0; i < np[d]; i++) {
      multi->ldims[d][i] = ldims[d] + (i < rmndr[d]);
      if (i > 0) {
	multi->off[d][i] = multi->off[d][i-1] + multi->ldims[d][i-1];
      }
    }
  }

  get_nr_patches_gpatch_off(domain, domain->rank,
			    &multi->nr_patches, &multi->gpatch_off);

  multi->patches = calloc(multi->nr_patches, sizeof(*multi->patches));
  for (int p = 0; p < multi->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_multi_get_global_patch_info(domain, p + multi->gpatch_off, &info);
    for (int d = 0; d < 3; d++) {
      multi->patches[p].ldims[d] = info.ldims[d];
      multi->patches[p].off[d] = info.off[d];
    }
  }
}

static void
mrc_domain_multi_destroy(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  free(multi->patches);
}

static struct mrc_patch *
mrc_domain_multi_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);
  if (nr_patches) {
    *nr_patches = multi->nr_patches;
  }
  return multi->patches;
}

static void
mrc_domain_multi_get_patch_idx3(struct mrc_domain *domain, int p, int *idx)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);
  assert(p >= 0 && p < multi->nr_patches);
  int gpatch = p + multi->gpatch_off;
  int *np = multi->np;

  idx[0] = gpatch % np[0]; gpatch /= np[0];
  idx[1] = gpatch % np[1]; gpatch /= np[1];
  idx[2] = gpatch % np[2];
}

static void
mrc_domain_multi_get_global_dims(struct mrc_domain *domain, int *dims)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    dims[d] = multi->gdims[d];
  }
}

static void
mrc_domain_multi_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    nr_procs[d] = multi->np[d];
  }
}

static void
mrc_domain_multi_get_bc(struct mrc_domain *domain, int *bc)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    bc[d] = multi->bc[d];
  }
}

static void
mrc_domain_multi_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  *nr_global_patches = multi->np[0] * multi->np[1] * multi->np[2];
}

static void
mrc_domain_multi_get_idx3_patch_info(struct mrc_domain *domain, int idx[3],
				     struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int gpatch = idx3_to_gpatch(multi, idx);
  mrc_domain_multi_get_global_patch_info(domain, gpatch, info);
}

static void
mrc_domain_multi_plot(struct mrc_domain *domain)
{
  if (domain->rank != 0) {
    return;
  }
  const char *name = mrc_domain_name(domain);
  char filename[strlen(name) + 20];
  sprintf(filename, "%s-patches.asc", name);
  FILE *file = fopen(filename, "w");
  sprintf(filename, "%s-curve.asc", name);
  FILE *file_curve = fopen(filename, "w");

  int nr_global_patches;
  mrc_domain_get_nr_global_patches(domain, &nr_global_patches);
  
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);
    fprintf(file, "%d %d %d\n", info.off[0], info.off[1], info.rank);
    fprintf(file, "%d %d %d\n", info.off[0] + info.ldims[0], info.off[1], info.rank);
    fprintf(file, "\n");
    fprintf(file, "%d %d %d\n", info.off[0], info.off[1] + info.ldims[1], info.rank);
    fprintf(file, "%d %d %d\n", info.off[0] + info.ldims[0], info.off[1] + info.ldims[1],
	    info.rank);
    fprintf(file, "\n");
    fprintf(file, "\n");

    fprintf(file_curve, "%g %g %d\n", info.off[0] + .5 * info.ldims[0],
	    info.off[1] + .5 * info.ldims[1], info.rank);
  }

  fclose(file);
  fclose(file_curve);
}

static struct mrc_ddc *
mrc_domain_multi_create_ddc(struct mrc_domain *domain)
{
  struct mrc_ddc *ddc = mrc_ddc_create(domain->obj.comm);
  mrc_ddc_set_type(ddc, "multi");
  mrc_ddc_set_domain(ddc, domain);
  return ddc;
}

static struct mrc_param_select bc_descr[] = {
  { .val = BC_NONE       , .str = "none"     },
  { .val = BC_PERIODIC   , .str = "periodic" },
  {},
};

#define VAR(x) (void *)offsetof(struct mrc_domain_multi, x)
static struct param mrc_domain_multi_params_descr[] = {
  { "m"               , VAR(gdims)           , PARAM_INT3(32, 32, 32) },
  { "np"              , VAR(np)              , PARAM_INT3(1, 1, 1)    },
  { "bcx"             , VAR(bc[0])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcy"             , VAR(bc[1])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcz"             , VAR(bc[2])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  {},
};
#undef VAR

static struct mrc_domain_ops mrc_domain_multi_ops = {
  .name                  = "multi",
  .size                  = sizeof(struct mrc_domain_multi),
  .param_descr           = mrc_domain_multi_params_descr,
  .setup                 = mrc_domain_multi_setup,
  .view                  = mrc_domain_multi_view,
  .destroy               = mrc_domain_multi_destroy,
  .get_patches           = mrc_domain_multi_get_patches,
  .get_patch_idx3        = mrc_domain_multi_get_patch_idx3,
  .get_global_dims       = mrc_domain_multi_get_global_dims,
  .get_nr_procs          = mrc_domain_multi_get_nr_procs,
  .get_bc                = mrc_domain_multi_get_bc,
  .get_nr_global_patches = mrc_domain_multi_get_nr_global_patches,
  .get_global_patch_info = mrc_domain_multi_get_global_patch_info,
  .get_idx3_patch_info   = mrc_domain_multi_get_idx3_patch_info,
  .plot                  = mrc_domain_multi_plot,
  .create_ddc            = mrc_domain_multi_create_ddc,
};

void
libmrc_domain_register_multi()
{
  libmrc_domain_register(&mrc_domain_multi_ops);
}
