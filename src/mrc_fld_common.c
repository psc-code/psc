
#include <mrc_ddc.h>

// This file is included multiple times to generate versions for various
// TYPE's (e.g., float, double, ...)

// ----------------------------------------------------------------------
// mrc_fld_TYPE_ddc_copy_to_buf

void
mrc_fld_TYPE_ddc_copy_to_buf(struct mrc_fld *fld, int mb, int me, int p,
			     int ilo[3], int ihi[3], void *_buf)
{
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

// ----------------------------------------------------------------------
// mrc_fld_TYPE_ddc_copy_from_buf

void
mrc_fld_TYPE_ddc_copy_from_buf(struct mrc_fld *fld, int mb, int me, int p,
			       int ilo[3], int ihi[3], void *_buf)
{
  //mprintf("from %d:%d x %d:%d x %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
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
