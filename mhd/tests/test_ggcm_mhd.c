
#include "ggcm_mhd.h"

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_destroy(mhd);

  MPI_Finalize();
}
