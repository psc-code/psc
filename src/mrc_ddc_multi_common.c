
#include <mrc_domain.h>
#include <mrc_ddc.h>

// ======================================================================
// mrc_ddc_funcs_fld

#include <mrc_fld.h>

static void
mrc_fld_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
		   void *_buf, void *ctx)
{
  //mprintf("to %d:%d x %d:%d x %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
  struct mrc_fld *fld = ctx;
  mrc_fld_data_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf,m - mb, ix,iy,iz) = M3(fld, m, ix,iy,iz, p);
	}
      }
    }
  }
}

static void
mrc_fld_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
		     void *_buf, void *ctx)
{
  //mprintf("from %d:%d x %d:%d x %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
  struct mrc_fld *fld = ctx;
  mrc_fld_data_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  M3(fld, m, ix,iy,iz, p) = MRC_DDC_BUF3(buf,m - mb, ix,iy,iz);
	}
      }
    }
  }
}

struct mrc_ddc_funcs mrc_ddc_funcs_fld = {
  .copy_to_buf   = mrc_fld_copy_to_buf,
  .copy_from_buf = mrc_fld_copy_from_buf,
};

