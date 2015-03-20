
#include "mrc_decomposition_private.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// ======================================================================
// mrc_decomposition

// ----------------------------------------------------------------------
// mrc_decomposition_setup

static void
_mrc_decomposition_setup(struct mrc_decomposition *dc)
{
  MPI_Comm comm = mrc_decomposition_comm(dc);

  // So far, we only handle having the local sizes prescribed
  if (dc->n != 0 && dc->N == 0) {
    MPI_Allreduce(&dc->n, &dc->N, 1, MPI_INT, MPI_SUM, comm);
  } else {
    mprintf("FIXME n %d N %d\n", dc->n, dc->N);
    //    assert(0);
  }

  MPI_Comm_rank(comm, &dc->mpi_rank);
  MPI_Comm_size(comm, &dc->mpi_size);

  dc->offs_by_rank = calloc(dc->mpi_size + 1, sizeof(dc->offs_by_rank));
  MPI_Allgather(&dc->n, 1, MPI_INT, dc->offs_by_rank, 1, MPI_INT, comm);
  int off = 0;
  for (int r = 0; r <= dc->mpi_size; r++) {
    int cnt = dc->offs_by_rank[r];
    dc->offs_by_rank[r] = off;
    off += cnt;
  }

  dc->off = dc->offs_by_rank[dc->mpi_rank];
}

// ----------------------------------------------------------------------
// mrc_decomposition_global_to_local

int
mrc_decomposition_global_to_local(struct mrc_decomposition *dc, int gidx)
{
  assert(gidx >= dc->off && gidx < dc->off + dc->n);
  return gidx - dc->off;
}

// ----------------------------------------------------------------------
// mrc_decomposition_is_local

bool
mrc_decomposition_is_local(struct mrc_decomposition *dc, int gidx)
{
  return gidx >= dc->off && gidx < dc->off + dc->n;
}

// ----------------------------------------------------------------------
// mrc_decomposition_find_rank

int
mrc_decomposition_find_rank(struct mrc_decomposition *dc, int gidx)
{
  static int r = 0;

  // see whether we can start from last rank
  if (r >= dc->mpi_size || gidx < dc->offs_by_rank[r]) {
    // no, need to start over
    r = 0;
  }

  for (; r < dc->mpi_size; r++) {
    if (gidx < dc->offs_by_rank[r+1]) {
      return r;
    }
  }
  assert(0);
}




// ----------------------------------------------------------------------
// mrc_decomposition class definition

struct mrc_class_mrc_decomposition mrc_class_mrc_decomposition = {
  .name         = "mrc_decomposition",
  .size         = sizeof(struct mrc_decomposition),
  .setup        = _mrc_decomposition_setup,
};


