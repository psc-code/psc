
#include "psc.h"
#include "psc_glue.h"

#include "mrc_io.h"

// ----------------------------------------------------------------------
// psc_read_checkpoint

struct psc *
psc_read_checkpoint(MPI_Comm comm)
{
  mpi_printf(MPI_COMM_WORLD, "INFO: Reading checkpoint.\n");

  struct mrc_io *io = mrc_io_create(comm);
  mrc_io_set_type(io, "xdmf_serial");
  mrc_io_set_name(io, "checkpoint");
  mrc_io_set_param_string(io, "basename", "checkpoint");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct psc *psc = psc_read(io, "psc");
  mrc_io_close(io);
  mrc_io_destroy(io);

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
  mrc_io_open(io, "w", 0, 0.);
  psc_write(psc, io);
  mrc_io_close(io);
  mrc_io_destroy(io);
}

