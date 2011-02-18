
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

void
mrc_domain_simple_set_params(struct mrc_domain *domain, struct mrc_domain_simple_params *par)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  simple->par = *par;
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

static void
mrc_domain_simple_setup(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);
  assert(!domain->is_setup);
  domain->is_setup = true;

  MPI_Comm_rank(domain->obj.comm, &domain->rank);
  MPI_Comm_size(domain->obj.comm, &domain->size);

  int total_procs = 1;
  for (int d = 0; d < 3; d++) {
    simple->ldims[d] = simple->par.ldims[d];
    simple->nr_procs[d] = simple->par.nr_procs[d];
    total_procs *= simple->par.nr_procs[d];
    simple->bc[d] = simple->par.bc[d];
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
      int rank_nei = mrc_domain_get_neighbor_rank(domain, dir);
      MPI_Recv(&simple->off[d], 1, MPI_INTEGER, rank_nei, TAG_SCAN_OFF + d, domain->obj.comm,
	       MPI_STATUS_IGNORE);
    } else {
      simple->off[d] = 0;
    }

    // then send next offset to the right
    if (proc[d] < simple->nr_procs[d] - 1) {
      int off_nei = simple->off[d] + simple->ldims[d];
      int dir[3] = {};
      dir[d] = 1;
      int rank_nei = mrc_domain_get_neighbor_rank(domain, dir);
      MPI_Send(&off_nei, 1, MPI_INTEGER, rank_nei, TAG_SCAN_OFF + d, domain->obj.comm);
    }
  }
    
  int rank_last =
    mrc_domain_simple_proc2rank(domain, (int [3]) { simple->nr_procs[0] - 1,
	  simple->nr_procs[1] - 1, simple->nr_procs[2] - 1});

  if (domain->rank == rank_last) {
    for (int d = 0; d < 3; d++) {
      simple->gdims[d] = simple->off[d] + simple->ldims[d];
    }
  }
  MPI_Bcast(simple->gdims, 3, MPI_INTEGER, rank_last, domain->obj.comm);

  //  mprintf("off   %d %d %d\n", simple->off[0], simple->off[1], simple->off[2]);
  //  mprintf("gdims %d %d %d\n", simple->gdims[0], simple->gdims[1], simple->gdims[2]);
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
mrc_domain_simple_get_local_offset_dims(struct mrc_domain *domain, int *off, int *dims)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  if (off) {
    for (int d = 0; d < 3; d++) {
      off[d] = simple->off[d];
    }
  }
  if (dims) {
    for (int d = 0; d < 3; d++) {
      dims[d] = simple->ldims[d];
    }
  }
}

static void
mrc_domain_simple_get_local_idx(struct mrc_domain *domain, int *idx)
{
  mrc_domain_simple_rank2proc(domain, domain->rank, idx);
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
mrc_domain_simple_get_bc(struct mrc_domain *domain, int *bc)
{
  struct mrc_domain_simple *simple = mrc_domain_simple(domain);

  for (int d = 0; d < 3; d++) {
    bc[d] = simple->bc[d];
  }
}

static struct mrc_ddc *
mrc_domain_simple_create_ddc(struct mrc_domain *domain, struct mrc_ddc_params *ddc_par,
			     struct mrc_ddc_ops *ddc_ops)
{
  int off[3], ldims[3], nr_procs[3], bc[3];
  mrc_domain_get_local_offset_dims(domain, off, ldims);
  mrc_domain_get_nr_procs(domain, nr_procs);
  mrc_domain_get_bc(domain, bc);

  for (int d = 0; d < 3; d++) {
    ddc_par->ilo[d] = 0;
    ddc_par->ihi[d] = ldims[d];
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

#define VAR(x) (void *)offsetof(struct mrc_domain_simple_params, x)
static struct param mrc_domain_simple_params_descr[] = {
  { "lmx"             , VAR(ldims[0])        , PARAM_INT(32)          },
  { "lmy"             , VAR(ldims[1])        , PARAM_INT(32)          },
  { "lmz"             , VAR(ldims[2])        , PARAM_INT(32)          },
  { "npx"             , VAR(nr_procs[0])     , PARAM_INT(1)           },
  { "npy"             , VAR(nr_procs[1])     , PARAM_INT(1)           },
  { "npz"             , VAR(nr_procs[2])     , PARAM_INT(1)           },
  { "bcx"             , VAR(bc[0])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcy"             , VAR(bc[1])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcz"             , VAR(bc[2])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  {},
};
#undef VAR

static struct mrc_domain_ops mrc_domain_simple_ops = {
  .name                  = "simple",
  .size                  = sizeof(struct mrc_domain_simple),
  .param_descr           = mrc_domain_simple_params_descr,
  .param_offset          = offsetof(struct mrc_domain_simple, par),
  .setup                 = mrc_domain_simple_setup,
  .get_neighbor_rank     = mrc_domain_simple_get_neighbor_rank,
  .get_local_offset_dims = mrc_domain_simple_get_local_offset_dims,
  .get_local_idx         = mrc_domain_simple_get_local_idx,
  .get_global_dims       = mrc_domain_simple_get_global_dims,
  .get_nr_procs          = mrc_domain_simple_get_nr_procs,
  .get_bc                = mrc_domain_simple_get_bc,
  .create_ddc            = mrc_domain_simple_create_ddc,
};

void
libmrc_domain_register_simple()
{
  libmrc_domain_register(&mrc_domain_simple_ops);
}
