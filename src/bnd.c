
#include "psc.h"
#include "util/profile.h"
#include "util/ddc.h"

#include <mpi.h>
#include <string.h>

// ======================================================================
// C bnd

static void
copy_to_buf(int fld_nr, int ilo[3], int ihi[3], f_real *buf)
{
  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	DDC_BUF(buf, ix,iy,iz) = FF3(fld_nr, ix,iy,iz);
      }
    }
  }
}

static void
add_from_buf(int fld_nr, int ilo[3], int ihi[3], f_real *buf)
{
  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	FF3(fld_nr, ix,iy,iz) += DDC_BUF(buf, ix,iy,iz);
      }
    }
  }
}

// ======================================================================
  
static void
c_add_ghosts(int m)
{
  if (!psc.bnd_data) {
    struct ddc_params prm = {
      .comm   = MPI_COMM_WORLD,
      .n_proc = { psc.domain.nproc[0], psc.domain.nproc[1], psc.domain.nproc[2] },
      .ilo    = { psc.ilo[0], psc.ilo[1], psc.ilo[2] },
      .ihi    = { psc.ihi[0], psc.ihi[1], psc.ihi[2] },
      .ibn    = { psc.ibn[0], psc.ibn[1], psc.ibn[2] },
      .copy_to_buf   = copy_to_buf,
      .add_from_buf  = add_from_buf,
    };
    for (int d = 0; d < 3; d++) {
      if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC) {
	assert(psc.domain.bnd_fld_hi[d] == BND_FLD_PERIODIC);
	prm.bc[d] = DDC_BC_PERIODIC;
      }
    }
    psc.bnd_data = ddc_create(&prm);
  }

  static int pr;
  if (!pr) {
    pr = prof_register("c_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);

  struct ddc_subdomain *ddc = psc.bnd_data;
  ddc_add_ghosts(ddc, m);

  prof_stop(pr);
}

struct psc_bnd_ops psc_bnd_ops_c = {
  .name       = "c",
  .add_ghosts = c_add_ghosts,
};
