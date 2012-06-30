
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_a3.h>
#include <mrc_io.h>
#include <mrctest.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int testcase = 0;
  mrc_params_get_option_int("case", &testcase);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_domain_set_type(domain, "amr");
  mrc_domain_set_param_int3(domain, "m", (int [3]) { 32, 32, 1});
  mrc_crds_set_type(crds, "multi_uniform");

  mrc_domain_set_from_options(domain);
  switch (testcase) {
  case 0:
    mrc_domain_add_patch(domain, 0, (int [3]) { 0, 0, 0 });
    break;
  case 1:
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 1, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
    break;
  case 2:
    mrc_domain_add_patch(domain, 0, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 1, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
    break;
  case 3:
    mrc_domain_add_patch(domain, 0, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 1, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
    break;
  case 4:
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
    break;
  case 5:
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
    mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 4, 6, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 4, 7, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 5, 6, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 5, 7, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 6, 6, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 6, 7, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 7, 6, 0 });
    mrc_domain_add_patch(domain, 3, (int [3]) { 7, 7, 0 });
    break;
  default:
    assert(0);
  }

  mrc_domain_setup(domain);
  mrc_domain_plot(domain);

  // create and fill a field

  struct mrc_a3 *fld = mrc_domain_a3_create(domain);
  mrc_a3_set_from_options(fld);
  mrc_a3_setup(fld);

  mrc_a3_foreach_patch(fld, p) {
    struct mrc_a3_patch *a3p = mrc_a3_patch_get(fld, p);
    double dx = .1, dy = .1;
    mrc_a3_foreach(a3p, ix,iy,iz, 0, 0) {
      double xx = ix * dx, yy = iy * dy;
      MRC_A3(a3p, 0, ix,iy,iz) = sin(xx) * cos(yy);
    } mrc_a3_foreach_end;
  }

  // write field to disk

  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_a3_write(fld, io);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_a3_destroy(fld);

  mrc_domain_destroy(domain);

  MPI_Finalize();
}
