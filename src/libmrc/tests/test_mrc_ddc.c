
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrc_ddc.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>

// ======================================================================

static void
set_m3(struct mrc_fld *m3)
{
  struct mrc_patch *patches = mrc_domain_get_patches(m3->_domain, NULL);

  mrc_fld_foreach_patch(m3, p) {
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    int *off = patches[p].off;
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      MRC_M3(m3p, 0, ix,iy,iz) =
	(iz + off[2]) * 10000 + (iy + off[1]) * 100 + (ix + off[0]);
    } mrc_m3_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

static void
check_m3(struct mrc_fld *m3)
{
  int bc[3], gdims[3];
  mrc_domain_get_bc(m3->_domain, bc);
  mrc_domain_get_global_dims(m3->_domain, gdims);

  struct mrc_patch *patches = mrc_domain_get_patches(m3->_domain, NULL);
  mrc_fld_foreach_patch(m3, p) {
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    int *off = patches[p].off;
    mrc_m3_foreach_bnd(m3p, ix,iy,iz) {
      int jx = ix + off[0];
      int jy = iy + off[1];
      int jz = iz + off[2];
      if (jx < 0 || jx >= gdims[0]) {
	if (bc[0] == BC_PERIODIC) {
	  jx = (jx + gdims[0]) % gdims[0];
	} else {
	  continue;
	}
      }
      if (jy < 0 || jy >= gdims[1]) {
	if (bc[1] == BC_PERIODIC) {
	  jy = (jy + gdims[1]) % gdims[1];
	} else {
	  continue;
	}
      }
      if (jz < 0 || jz >= gdims[2]) {
	if (bc[2] == BC_PERIODIC) {
	  jz = (jz + gdims[2]) % gdims[2];
	} else {
	  continue;
	}
      }
#if 1
      if (MRC_M3(m3p, 0, ix,iy,iz) != jz * 10000 + jy * 100 + jx) {
	printf("ixyz %d %d %d jxyz %d %d %d : %d val %g\n",
	       ix,iy,iz, jx,jy,jz, jz * 10000 + jy * 100 + jx,
	       MRC_M3(m3p, 0, ix,iy,iz));
      }
#endif
      assert(MRC_M3(m3p, 0, ix,iy,iz) == jz * 10000 + jy * 100 + jx);
    } mrc_m3_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "multi");
  mrc_domain_set_param_int(domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(domain, "bcy", BC_PERIODIC);
  mrc_domain_set_param_int(domain, "bcz", BC_PERIODIC);
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);
  mrc_domain_plot(domain);

  struct mrc_fld *m3 = mrc_domain_m3_create(domain);
  mrc_fld_set_name(m3, "test_m3");
  mrc_fld_set_param_int(m3, "nr_comps", 2);
  mrc_fld_set_param_int(m3, "nr_ghosts", 1);
  mrc_fld_set_from_options(m3);
  mrc_fld_setup(m3);
  mrc_fld_set_comp_name(m3, 0, "fld0");
  mrc_fld_set_comp_name(m3, 1, "fld1");
  mrc_fld_view(m3);

  set_m3(m3);

  struct mrc_ddc *ddc = mrc_domain_create_ddc(domain);
  mrc_ddc_setup(ddc);
  mrc_ddc_view(ddc);
  mrc_ddc_fill_ghosts_fld(ddc, 0, 2, m3);
  mrc_ddc_destroy(ddc);

  check_m3(m3);

  mrc_fld_destroy(m3);

  mrc_domain_destroy(domain);

  MPI_Finalize();
}
