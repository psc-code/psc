
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include <mrc_fld.h>
#include <mrc_ddc.h>
#include <mrc_domain.h>
#include <mrc_profile.h>

// FIXME, 0-based offsets and ghost points don't match well (not pretty anyway)

static void
test(bool periodic)
{
  const int bnd = SW_2;

  int bc[3] = { BC_NONE, BC_NONE, BC_NONE };
  if (periodic) {
    bc[0] = BC_PERIODIC;
    bc[1] = BC_PERIODIC;
    bc[2] = BC_PERIODIC;
  }

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_param_int3(domain, "lm", (int [3]) { 4, 8, 16 });
  mrc_domain_set_param_int3(domain, "np", (int [3]) { 2, 1, 1 });
  mrc_domain_set_param_int(domain, "bcx", bc[0]);
  mrc_domain_set_param_int(domain, "bcy", bc[1]);
  mrc_domain_set_param_int(domain, "bcz", bc[2]);
  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);
  mrc_domain_view(domain);

  struct mrc_ddc_params ddc_params = {
    .ibn = { bnd, bnd, bnd },
    .max_n_fields = 2,
    .size_of_type = sizeof(float),
  };

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &nr_patches);
  assert(nr_patches == 1);
  int *off = patches[0].off;

  struct mrc_ddc *ddc = mrc_domain_create_ddc(domain, &ddc_params, &mrc_ddc_ops_f3);

  struct mrc_f3 *fld = mrc_domain_f3_create(domain, bnd);
  mrc_f3_set_param_int(fld, "nr_comps", 2);
  mrc_f3_setup(fld);
  mrc_f3_view(fld);

  mrc_f3_foreach(fld, ix,iy,iz, bnd, bnd) {
    int jx = ix - bnd + off[0];
    int jy = iy - bnd + off[1];
    int jz = iz - bnd + off[2];
    MRC_F3(fld,0, ix,iy,iz) = jx * 10000 + jy * 100 + jz;
    MRC_F3(fld,1, ix,iy,iz) = - (jx * 10000 + jy * 100 + jz);
  } mrc_f3_foreach_end;

  for (int i = 0; i < 10000; i++) {
    mrc_ddc_fill_ghosts(ddc, 0, 2, fld);
  }

  int N[3];
  mrc_domain_get_global_dims(domain, N);
  mrc_f3_foreach(fld, ix,iy,iz, 0, 0) {
    int jx = ix - bnd + off[0];
    int jy = iy - bnd + off[1];
    int jz = iz - bnd + off[2];
    if (periodic) { 
      jx = (jx + N[0]) % N[0];
      jy = (jy + N[1]) % N[1];
      jz = (jz + N[2]) % N[2];
    }
    if (jx < 0 || jx >= N[0] ||
	jy < 0 || jy >= N[1] ||
	jz < 0 || jz >= N[2]) {
      continue;
    }
    if (MRC_F3(fld,0, ix,iy,iz) != jx * 10000 + jy * 100 + jz) {
      mprintf("!!!0 [%d,%d,%d] : %g -- %g\n", ix, iy, iz,
	      MRC_F3(fld,0, ix,iy,iz), (float) jx * 10000 + jy * 100 + jz);
      assert(0);
    }
    if (MRC_F3(fld,1, ix,iy,iz) != - (jx * 10000 + jy * 100 + jz)) {
      mprintf("!!!1 [%d,%d,%d] : %g -- %g\n", ix, iy, iz,
	      MRC_F3(fld,1, ix,iy,iz), -((float) jx * 10000 + jy * 100 + jz));
      assert(0);
    }
  } mrc_f3_foreach_end;

  mrc_f3_destroy(fld);
  mrc_ddc_destroy(ddc);
  mrc_domain_destroy(domain);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  test(false); // non periodic
  test(true);  // periodic

  prof_print();
  MPI_Finalize();
}
