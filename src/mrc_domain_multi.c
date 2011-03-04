
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

static void
mrc_domain_multi_rank2proc(struct mrc_domain *domain, int rank, int proc[3])
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);
  int *np = multi->np;

  proc[0] = rank % np[0]; rank /= np[0];
  proc[1] = rank % np[1]; rank /= np[1];
  proc[2] = rank % np[2];
}

static int
mrc_domain_multi_proc2rank(struct mrc_domain *domain, int proc[3])
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    if (proc[d] < 0 || proc[d] >= multi->np[d])
      return -1;
  }
  return ((proc[2] * multi->np[1]) + proc[1]) * multi->np[0] + proc[0];
}

static int
part2index(struct mrc_domain_multi *multi, int p[3])
{
  int *np = multi->np;
  return (p[2] * np[1] + p[1]) * np[0] + p[0];
}

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

#if 0
  assert(domain->is_setup);
  int *np = multi->nr_procs;
  int p3[3];
  for (p3[2] = 0; p3[2] < np[2]; p3[2]++) {
    for (p3[1] = 0; p3[1] < np[1]; p3[1]++) {
      for (p3[0] = 0; p3[0] < np[0]; p3[0]++) {
	struct mrc_patch *patch = &multi->all_patches[part2index(multi, p3)];
	mprintf("  patch: ldims %dx%dx%d off %dx%dx%d\n",
		patch->ldims[0], patch->ldims[1], patch->ldims[2],
		patch->off[0], patch->off[1], patch->off[2]);
      }
    }
  }
#endif
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

  // FIXME: allow setting of desired decomposition using ldims
  int ldims[3], rmndr[3], *np = multi->np;
  int nr_all_patches = 1;
  for (int d = 0; d < 3; d++) {
    nr_all_patches *= np[d];
    ldims[d] = multi->gdims[d] / np[d];
    rmndr[d] = multi->gdims[d] % np[d];
  }
  multi->all_patches = calloc(nr_all_patches, sizeof(*multi->all_patches));
  int p3[3];
  for (p3[2] = 0; p3[2] < np[2]; p3[2]++) {
    for (p3[1] = 0; p3[1] < np[1]; p3[1]++) {
      for (p3[0] = 0; p3[0] < np[0]; p3[0]++) {
	struct mrc_patch *patch = &multi->all_patches[part2index(multi, p3)];
	for (int d = 0; d < 3; d++) {
	  patch->ldims[d] = ldims[d] + (p3[d] < rmndr[d]);
	  if (p3[d] == 0) {
	    patch->off[d] = 0;
	  } else {
	    int q3[3] = { p3[0], p3[1], p3[2] };
	    q3[d]--;
	    struct mrc_patch *p_patch = &multi->all_patches[part2index(multi, q3)];
	    patch->off[d] = p_patch->off[d] + p_patch->ldims[d];
	  }
	}
      }
    }
  }

  int patches_per_proc = nr_all_patches / domain->size;
  int patches_per_proc_rmndr = nr_all_patches % domain->size;

  multi->nr_patches = patches_per_proc + (domain->rank < patches_per_proc_rmndr);
  int patches_off;
  if (domain->rank < patches_per_proc_rmndr) {
    patches_off = domain->rank * (patches_per_proc + 1);
  } else {
    patches_off = domain->rank * patches_per_proc + patches_per_proc_rmndr;
  }
  multi->patches = calloc(multi->nr_patches, sizeof(*multi->patches));
  for (int p = 0; p < multi->nr_patches; p++) {
    multi->patches[p] = multi->all_patches[p + patches_off];
  }
}

static void
mrc_domain_multi_destroy(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  free(multi->patches);
}

static int
mrc_domain_multi_get_neighbor_rank(struct mrc_domain *domain, int shift[3])
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int nei[3];
  for (int d = 0; d < 3; d++) {
    nei[d] = multi->proc[d] + shift[d];
  }

  return mrc_domain_multi_proc2rank(domain, nei);
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
mrc_domain_multi_get_local_idx(struct mrc_domain *domain, int *idx)
{
  mrc_domain_multi_rank2proc(domain, domain->rank, idx);
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
mrc_domain_multi_get_global_patch_info(struct mrc_domain *domain, int gpatch,
				       struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int nr_global_patches = multi->np[0] * multi->np[1] * multi->np[2];
  int patches_per_proc = nr_global_patches / domain->size;
  int patches_per_proc_rmndr = nr_global_patches % domain->size;
  
  if (gpatch < (patches_per_proc + 1) * patches_per_proc_rmndr) {
    info->rank = gpatch / (patches_per_proc + 1);
    info->patch = gpatch % (patches_per_proc + 1);
  } else {
    int tmp = gpatch - (patches_per_proc + 1) * patches_per_proc_rmndr;
    info->rank = tmp / patches_per_proc + patches_per_proc_rmndr;
    info->patch = tmp % patches_per_proc;
  }

  for (int d = 0; d < 3; d++) {
    info->ldims[d] = multi->all_patches[gpatch].ldims[d];
  }
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
  .get_neighbor_rank     = mrc_domain_multi_get_neighbor_rank,
  .get_patches           = mrc_domain_multi_get_patches,
  .get_local_idx         = mrc_domain_multi_get_local_idx,
  .get_global_dims       = mrc_domain_multi_get_global_dims,
  .get_nr_procs          = mrc_domain_multi_get_nr_procs,
  .get_bc                = mrc_domain_multi_get_bc,
  .get_nr_global_patches = mrc_domain_multi_get_nr_global_patches,
  .get_global_patch_info = mrc_domain_multi_get_global_patch_info,
  .create_ddc            = mrc_domain_multi_create_ddc,
};

void
libmrc_domain_register_multi()
{
  libmrc_domain_register(&mrc_domain_multi_ops);
}
