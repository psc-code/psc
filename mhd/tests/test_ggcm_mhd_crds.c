
#include "ggcm_mhd_crds.h"

static void
test_1()
{
  struct ggcm_mhd_crds *crds = ggcm_mhd_crds_create(MPI_COMM_WORLD);
  ggcm_mhd_crds_destroy(crds);
}

static void
test_2()
{
  struct ggcm_mhd_crds *crds = ggcm_mhd_crds_create(MPI_COMM_WORLD);
  ggcm_mhd_crds_setup(crds);
  ggcm_mhd_crds_destroy(crds);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  test_1();
  test_2();

  MPI_Finalize();
}
