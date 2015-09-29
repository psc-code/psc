
#include "ggcm_mhd.h"

#include <mrc.h>
#include <mrc_io.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_from_options(mhd);
  ggcm_mhd_setup(mhd);
  ggcm_mhd_view(mhd);

  struct mrc_io *io = mrc_io_create(ggcm_mhd_comm(mhd));
  mrc_io_set_type(io, "hdf5_serial");
  mrc_io_setup(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/mhd", "mhd", mhd);
  mrc_io_close(io);

  mrc_io_destroy(io);

  ggcm_mhd_destroy(mhd);
  libmrc_finalize();
  MPI_Finalize();
}
