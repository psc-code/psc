
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_io.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>

static void
set_m3(struct mrc_m3 *m3)
{
  struct mrc_patch *patches = mrc_domain_get_patches(m3->domain, NULL);
  struct mrc_crds *crds = mrc_domain_get_crds(m3->domain);

  mrc_m3_foreach_patch(m3, p) {
    struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
    mrc_crds_patch_get(crds, p);
    int *off = patches[p].off;
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      MRC_M3(m3p, 0, ix,iy,iz) =
	(iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]);
      MRC_M3(m3p, 1, ix,iy,iz) = MRC_MCRD(crds, 0, ix);
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(m3);
    mrc_crds_patch_put(crds);
  }
}

static void
check_m3(struct mrc_m3 *m3)
{
  struct mrc_patch *patches = mrc_domain_get_patches(m3->domain, NULL);

  mrc_m3_foreach_patch(m3, p) {
    struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
    int *off = patches[p].off;
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      assert(MRC_M3(m3p, 0, ix,iy,iz) ==
	     (iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]));
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(m3);
  }
}

static void
write_m3(struct mrc_m3 *m3)
{
  struct mrc_io *io = mrc_io_create(mrc_m3_comm(m3));

  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);

  mrc_io_open(io, "w", 0, 0.);
  mrc_m3_write(m3, io);
  mrc_io_close(io);

  mrc_io_open(io, "w", 1, 1.);
  mrc_m3_write(m3, io);
  mrc_io_close(io);

  mrc_io_destroy(io);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "multi");
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "multi_uniform");
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);

  struct mrc_m3 *m3 = mrc_domain_m3_create(domain, 0);
  mrc_m3_set_name(m3, "test_m3");
  mrc_m3_set_param_int(m3, "nr_comps", 2);
  mrc_m3_setup(m3);
  m3->name[0] = strdup("fld0");
  m3->name[1] = strdup("fld1");
  mrc_m3_view(m3);

  set_m3(m3);
  check_m3(m3);
  write_m3(m3);

  mrc_m3_destroy(m3);

  mrc_domain_destroy(domain);

  MPI_Finalize();
}
