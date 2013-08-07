
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_fld.h>
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
  mrc_domain_set_param_int3(domain, "m", (int [3]) { 16, 16, 1});
  mrc_crds_set_type(crds, "amr_uniform");
  
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

  struct mrc_fld *fld = mrc_domain_m3_create(domain);
  mrc_fld_set_name(fld, "fld");
  mrc_fld_set_from_options(fld);
  mrc_fld_setup(fld);
  mrc_fld_set_comp_name(fld, 0, "m0");

  float kx = 2. * M_PI, ky = 2. * M_PI;

  mrc_fld_foreach_patch(fld, p) {
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(fld, p);
    mrc_m3_foreach(m3p, ix,iy,iz, 0, 0) {
      float xx = MRC_MCRDX(crds, ix, p), yy = MRC_MCRDY(crds, iy, p);
      MRC_M3(m3p, 0, ix,iy,iz) = sin(kx * xx) * cos(ky * yy);
    } mrc_m3_foreach_end;
    mrc_fld_patch_put(fld);
  }

  // write field to disk

  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_type(io, "ascii");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_fld_write(fld, io);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_fld_destroy(fld);

  mrc_domain_destroy(domain);

  MPI_Finalize();
}
