
#include "psc.h"
#include "psc_glue.h"

#include "mrc_io.h"

// ----------------------------------------------------------------------
// psc_read_checkpoint

struct psc *
psc_read_checkpoint(MPI_Comm comm, int n)
{
  mpi_printf(MPI_COMM_WORLD, "INFO: Reading checkpoint.\n");

  struct mrc_io *io = mrc_io_create(comm);
  mrc_io_set_type(io, "xdmf_serial");
  mrc_io_set_name(io, "checkpoint");
  mrc_io_set_param_string(io, "basename", "checkpoint");
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

  struct mrc_io *io = mrc_io_create(psc_comm(psc));
  mrc_io_set_type(io, "xdmf_serial");
  mrc_io_set_name(io, "checkpoint");
  mrc_io_set_param_string(io, "basename", "checkpoint");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", psc->timestep, psc->timestep * psc->dt);
  mrc_io_write_path(io, "checkpoint", "psc", psc);
  mrc_io_close(io);
  mrc_io_destroy(io);
}

