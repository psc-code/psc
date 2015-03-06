
#include "mrc_decomposition_private.h"

#include <assert.h>

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
    assert(0);
  }

  MPI_Exscan(&dc->n, &dc->off, 1, MPI_INT, MPI_SUM, comm);
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
// mrc_decomposition class definition

struct mrc_class_mrc_decomposition mrc_class_mrc_decomposition = {
  .name         = "mrc_decomposition",
  .size         = sizeof(struct mrc_decomposition),
  .setup        = _mrc_decomposition_setup,
};


