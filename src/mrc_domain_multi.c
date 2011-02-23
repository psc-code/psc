
#include "mrc_domain_private.h"
#include "mrc_params.h"
#include "mrc_ddc.h"

#include <stdlib.h>
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
  int *nr_procs = multi->nr_procs;

  proc[0] = rank % nr_procs[0]; rank /= nr_procs[0];
  proc[1] = rank % nr_procs[1]; rank /= nr_procs[1];
  proc[2] = rank % nr_procs[2];
}

static int
mrc_domain_multi_proc2rank(struct mrc_domain *domain, int proc[3])
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    if (proc[d] < 0 || proc[d] >= multi->nr_procs[d])
      return -1;
  }
  return ((proc[2] * multi->nr_procs[1]) + proc[1]) * multi->nr_procs[0] + proc[0];
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

  int total_procs = 1;
  for (int d = 0; d < 3; d++) {
    if (multi->patch.ldims[d] == 0) {
      assert(multi->gdims[d] % multi->nr_procs[d] == 0);
      multi->patch.ldims[d] = multi->gdims[d] / multi->nr_procs[d];
    }
    total_procs *= multi->nr_procs[d];
  }
  assert(total_procs == domain->size);
  mrc_domain_multi_rank2proc(domain, domain->rank, multi->proc);

  // find local offsets into global domain
  // these are basically scans along subsets (lines).
  int *proc = multi->proc;
  for (int d = 0; d < 3; d++) {
    // get local offset from left neighbor
    if (proc[d] > 0) {
      int dir[3] = {};
      dir[d] = -1;
      int rank_nei = mrc_domain_get_neighbor_rank(domain, dir);
      MPI_Recv(&multi->patch.off[d], 1, MPI_INTEGER, rank_nei, TAG_SCAN_OFF + d, domain->obj.comm,
	       MPI_STATUS_IGNORE);
    } else {
      multi->patch.off[d] = 0;
    }

    // then send next offset to the right
    if (proc[d] < multi->nr_procs[d] - 1) {
      int off_nei = multi->patch.off[d] + multi->patch.ldims[d];
      int dir[3] = {};
      dir[d] = 1;
      int rank_nei = mrc_domain_get_neighbor_rank(domain, dir);
      MPI_Send(&off_nei, 1, MPI_INTEGER, rank_nei, TAG_SCAN_OFF + d, domain->obj.comm);
    }
  }
    
  int rank_last =
    mrc_domain_multi_proc2rank(domain, (int [3]) { multi->nr_procs[0] - 1,
	  multi->nr_procs[1] - 1, multi->nr_procs[2] - 1});

  if (domain->rank == rank_last) {
    for (int d = 0; d < 3; d++) {
      multi->gdims[d] = multi->patch.off[d] + multi->patch.ldims[d];
    }
  }
  MPI_Bcast(multi->gdims, 3, MPI_INTEGER, rank_last, domain->obj.comm);

  //  mprintf("off   %d %d %d\n", multi->off[0], multi->off[1], multi->off[2]);
  //  mprintf("gdims %d %d %d\n", multi->gdims[0], multi->gdims[1], multi->gdims[2]);
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
    *nr_patches = 1;
  }
  return &multi->patch;
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
    nr_procs[d] = multi->nr_procs[d];
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

static struct mrc_ddc *
mrc_domain_multi_create_ddc(struct mrc_domain *domain, struct mrc_ddc_params *ddc_par,
			     struct mrc_ddc_ops *ddc_ops)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int nr_procs[3], bc[3];
  mrc_domain_get_nr_procs(domain, nr_procs);
  mrc_domain_get_bc(domain, bc);

  for (int d = 0; d < 3; d++) {
    ddc_par->ilo[d] = 0;
    ddc_par->ihi[d] = multi->patch.ldims[d];
    ddc_par->n_proc[d] = nr_procs[d];
    ddc_par->bc[d] = bc[d];
  }

  return mrc_ddc_create(domain->obj.comm, ddc_par, ddc_ops);
}

static struct mrc_param_select bc_descr[] = {
  { .val = BC_NONE       , .str = "none"     },
  { .val = BC_PERIODIC   , .str = "periodic" },
  {},
};

#define VAR(x) (void *)offsetof(struct mrc_domain_multi, x)
static struct param mrc_domain_multi_params_descr[] = {
  { "lm"              , VAR(patch.ldims)     , PARAM_INT3(0, 0, 0)    },
  { "m"               , VAR(gdims)           , PARAM_INT3(32, 32, 32) },
  { "np"              , VAR(nr_procs)        , PARAM_INT3(1, 1, 1)    },
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
  .get_neighbor_rank     = mrc_domain_multi_get_neighbor_rank,
  .get_patches           = mrc_domain_multi_get_patches,
  .get_local_idx         = mrc_domain_multi_get_local_idx,
  .get_global_dims       = mrc_domain_multi_get_global_dims,
  .get_nr_procs          = mrc_domain_multi_get_nr_procs,
  .get_bc                = mrc_domain_multi_get_bc,
  .create_ddc            = mrc_domain_multi_create_ddc,
};

void
libmrc_domain_register_multi()
{
  libmrc_domain_register(&mrc_domain_multi_ops);
}
