
#include "psc.h"

#include "mrc_io.h"

#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>

#if 0

// ----------------------------------------------------------------------
// psc_read_checkpoint

struct psc *
psc_read_checkpoint(MPI_Comm comm, int n)
{
  mpi_printf(MPI_COMM_WORLD, "INFO: Reading checkpoint.\n");

  char dir[30];
  sprintf(dir, "checkpoint.%08d", n);

  struct mrc_io *io = mrc_io_create(comm);
  mrc_io_set_type(io, "xdmf_serial");
  mrc_io_set_name(io, "checkpoint");
  mrc_io_set_param_string(io, "basename", "checkpoint");
  mrc_io_set_param_string(io, "outdir", dir);
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", n, 0.);
  struct psc *psc = mrc_io_read_path(io, "checkpoint", "psc", psc);
  mrc_io_close(io);
  mrc_io_destroy(io);

  ppsc = psc;
  return psc;
}

// ----------------------------------------------------------------------
// psc_write_checkpoint

void
psc_write_checkpoint(struct psc *psc)
{
  mpi_printf(psc_comm(psc), "INFO: Writing checkpoint.\n");

  char dir[30];
  sprintf(dir, "checkpoint.%08d", psc->timestep);
  int rank;
  MPI_Comm_rank(psc_comm(psc), &rank);
  if (rank == 0) {
    if (mkdir(dir, 0777)) {
      if (errno != EEXIST) {
	perror("ERROR: mkdir");
	abort();
      }
    }
  }
  MPI_Barrier(psc_comm(psc));

  struct mrc_io *io = mrc_io_create(psc_comm(psc));
  mrc_io_set_type(io, "xdmf_serial");
  mrc_io_set_name(io, "checkpoint");
  mrc_io_set_param_string(io, "basename", "checkpoint");
  mrc_io_set_param_string(io, "outdir", dir);
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", psc->timestep, psc->timestep * psc->grid().dt);
  mrc_io_write_path(io, "checkpoint", "psc", psc);
  mrc_io_close(io);
  mrc_io_destroy(io);
}

#endif
