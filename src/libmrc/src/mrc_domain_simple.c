
#include "mrc_domain_private.h"
#include "mrc_params.h"
#include "mrc_ddc.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define TAG_SCAN_OFF (1000)

static inline struct mrc_domain_simple *
mrc_domain_simple(struct mrc_domain *domain)
{
  return domain->obj.subctx;
}

static void
mrc_domain_simple_rank2proc(struct mrc_domain *domain, int rank, int proc[3])
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);
  int *nr_procs = simple->nr_procs;

  proc[0] = rank % nr_procs[0]; rank /= nr_procs[0];
  proc[1] = rank % nr_procs[1]; rank /= nr_procs[1];
  proc[2] = rank % nr_procs[2];
}

static int
mrc_domain_simple_proc2rank(struct mrc_domain *domain, int proc[3])
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  for (int d = 0; d < 3; d++) {
    if (proc[d] < 0 || proc[d] >= simple->nr_procs[d])
      return -1;
  }
  return ((proc[2] * simple->nr_procs[1]) + proc[1]) * simple->nr_procs[0] + proc[0];
}

static int
mrc_domain_simple_get_neighbor_rank(struct mrc_domain *domain, int shift[3]);

static void
mrc_domain_simple_do_setup(struct mrc_domain *domain)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  MPI_Comm_rank(domain->obj.comm, &domain->rank);
  MPI_Comm_size(domain->obj.comm, &domain->size);

  if (domain->size == 1 &&
      simple->nr_procs[0] * simple->nr_procs[1] * simple->nr_procs[2] > 1) {
    // FIXME, this is just to make reading on one proc work
    for (int d = 0; d < 3; d++) {
      simple->nr_procs[d] = 1;
      simple->patch.ldims[d] = 0;
    }
  }
  int total_procs = 1;
  for (int d = 0; d < 3; d++) {
    if (simple->patch.ldims[d] == 0) {
      assert(simple->gdims[d] % simple->nr_procs[d] == 0);
      simple->patch.ldims[d] = simple->gdims[d] / simple->nr_procs[d];
    }
    total_procs *= simple->nr_procs[d];
  }
  assert(total_procs == domain->size);
  mrc_domain_simple_rank2proc(domain, domain->rank, simple->proc);

  // find local offsets into global domain
  // these are basically scans along subsets (lines).
  int *proc = simple->proc;
  for (int d = 0; d < 3; d++) {
    // get local offset from left neighbor
    if (proc[d] > 0) {
      int dir[3] = {};
      dir[d] = -1;
      int rank_nei = mrc_domain_simple_get_neighbor_rank(domain, dir);
      MPI_Recv(&simple->patch.off[d], 1, MPI_INTEGER, rank_nei, TAG_SCAN_OFF + d, domain->obj.comm,
	       MPI_STATUS_IGNORE);
    } else {
      simple->patch.off[d] = 0;
    }

    // then send next offset to the right
    if (proc[d] < simple->nr_procs[d] - 1) {
      int off_nei = simple->patch.off[d] + simple->patch.ldims[d];
      int dir[3] = {};
      dir[d] = 1;
      int rank_nei = mrc_domain_simple_get_neighbor_rank(domain, dir);
      MPI_Send(&off_nei, 1, MPI_INTEGER, rank_nei, TAG_SCAN_OFF + d, domain->obj.comm);
    }
  }
    
  int rank_last =
    mrc_domain_simple_proc2rank(domain, (int [3]) { simple->nr_procs[0] - 1,
	  simple->nr_procs[1] - 1, simple->nr_procs[2] - 1});

  if (domain->rank == rank_last) {
    for (int d = 0; d < 3; d++) {
      simple->gdims[d] = simple->patch.off[d] + simple->patch.ldims[d];
    }
  }
  MPI_Bcast(simple->gdims, 3, MPI_INTEGER, rank_last, domain->obj.comm);

  //  mprintf("off   %d %d %d\n", simple->off[0], simple->off[1], simple->off[2]);
  //  mprintf("gdims %d %d %d\n", simple->gdims[0], simple->gdims[1], simple->gdims[2]);
}

static void
mrc_domain_simple_setup(struct mrc_domain *domain)
{
  mrc_domain_simple_do_setup(domain);
  mrc_domain_setup_super(domain);
}

static void
mrc_domain_simple_read(struct mrc_domain *domain, struct mrc_io *io)
{
  mrc_domain_simple_do_setup(domain);
  mrc_domain_read_super(domain, io);
}

static int
mrc_domain_simple_get_neighbor_rank(struct mrc_domain *domain, int shift[3])
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  int nei[3];
  for (int d = 0; d < 3; d++) {
    nei[d] = simple->proc[d] + shift[d];
  }

  return mrc_domain_simple_proc2rank(domain, nei);
}

static void
mrc_domain_simple_get_neighbor_rank_patch(struct mrc_domain *domain, int p, int dir[3],
					  int *nei_rank, int *nei_patch)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);
  assert(p == 0);

  int proc_nei[3];
  for (int d = 0; d < 3; d++) {
    proc_nei[d] = simple->proc[d] + dir[d];
    if (domain->bc[d] == BC_PERIODIC) {
      if (proc_nei[d] < 0) {
	proc_nei[d] += simple->nr_procs[d];
      }
      if (proc_nei[d] >= simple->nr_procs[d]) {
	proc_nei[d] -= simple->nr_procs[d];
      }
    }
    if (proc_nei[d] < 0 || proc_nei[d] >= simple->nr_procs[d]) {
      *nei_rank = -1;
      *nei_patch = -1;
      return;
    }
  }
  *nei_rank = mrc_domain_simple_proc2rank(domain, proc_nei);
  *nei_patch = 0;
}

static struct mrc_patch *
mrc_domain_simple_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);
  if (nr_patches) {
    *nr_patches = 1;
  }
  return &simple->patch;
}

static void
mrc_domain_simple_get_global_dims(struct mrc_domain *domain, int *dims)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  for (int d = 0; d < 3; d++) {
    dims[d] = simple->gdims[d];
  }
}

static void
mrc_domain_simple_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  for (int d = 0; d < 3; d++) {
    nr_procs[d] = simple->nr_procs[d];
  }
}

static void
mrc_domain_simple_get_global_patch_info(struct mrc_domain *domain, int gpatch,
					struct mrc_patch_info *info)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  mrc_domain_simple_rank2proc(domain, gpatch, info->idx3);
  for (int d = 0; d < 3; d++) {
    assert(simple->gdims[d] % simple->nr_procs[d] == 0);
    info->ldims[d] = simple->gdims[d] / simple->nr_procs[d];
    info->off[d] = info->idx3[d] * info->ldims[d];
  }
  info->rank = gpatch;
  info->patch = 0;
  info->global_patch = gpatch;
  info->level = 0;
}

static void
mrc_domain_simple_get_local_patch_info(struct mrc_domain *domain, int patch,
				       struct mrc_patch_info *info)
{
  assert(patch == 0);
  mrc_domain_simple_get_global_patch_info(domain, domain->rank, info);
}

static void
mrc_domain_simple_get_level_idx3_patch_info(struct mrc_domain *domain, int level,
					    int idx[3], struct mrc_patch_info *info)
{
  assert(level == 0);
  int rank = mrc_domain_simple_proc2rank(domain, idx);
  mrc_domain_simple_get_global_patch_info(domain, rank, info);
}

static void
mrc_domain_simple_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  *nr_global_patches = simple->nr_procs[0] * simple->nr_procs[1] * simple->nr_procs[2];
}

static struct mrc_ddc *
mrc_domain_simple_create_ddc(struct mrc_domain *domain)
{
  struct mrc_ddc *ddc = mrc_ddc_create(domain->obj.comm);
  mrc_ddc_set_type(ddc, "multi");
  mrc_ddc_set_domain(ddc, domain);
  return ddc;
}

#define VAR(x) (void *)offsetof(struct mrc_domain_simple, x)
static struct param mrc_domain_simple_params_descr[] = {
  { "lm"              , VAR(patch.ldims)     , PARAM_INT3(0, 0, 0)    },
  { "m"               , VAR(gdims)           , PARAM_INT3(32, 32, 32),
  .help = "global number of grid points in x, y, z direction" },
  { "np"              , VAR(nr_procs)        , PARAM_INT3(1, 1, 1)    },
  {},
};
#undef VAR

struct mrc_domain_ops mrc_domain_simple_ops = {
  .name                      = "simple",
  .size                      = sizeof(struct mrc_domain_simple),
  .param_descr               = mrc_domain_simple_params_descr,
  .setup                     = mrc_domain_simple_setup,
  .read                      = mrc_domain_simple_read,
  .get_neighbor_rank         = mrc_domain_simple_get_neighbor_rank,
  .get_neighbor_rank_patch   = mrc_domain_simple_get_neighbor_rank_patch,
  .get_patches               = mrc_domain_simple_get_patches,
  .get_level_idx3_patch_info = mrc_domain_simple_get_level_idx3_patch_info,
  .get_global_dims           = mrc_domain_simple_get_global_dims,
  .get_nr_procs              = mrc_domain_simple_get_nr_procs,
  .get_local_patch_info      = mrc_domain_simple_get_local_patch_info,
  .get_global_patch_info     = mrc_domain_simple_get_global_patch_info,
  .get_nr_global_patches     = mrc_domain_simple_get_nr_global_patches,
  .create_ddc                = mrc_domain_simple_create_ddc,
};

