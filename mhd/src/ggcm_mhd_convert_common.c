
#include <ggcm_mhd_defs.h>

#include <mrc_fld.h>

void
copy_sc_ggcm_to_sc(struct mrc_fld *_ff, struct mrc_fld *_f)
{
  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int m = 0; m < _NR_FLDS; m++) {
    if (m == _B1X || m == _B2X) {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(f, m, ix,iy,iz) = F3(ff, m, ix-1,iy,iz);
      } mrc_fld_foreach_end;
    } else if (m == _B1Y || m == _B2Y) {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(f, m, ix,iy,iz) = F3(ff, m, ix,iy-1,iz);
      } mrc_fld_foreach_end;
    } else if (m == _B1Z || m == _B2Z) {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(f, m, ix,iy,iz) = F3(ff, m, ix,iy,iz-1);
      } mrc_fld_foreach_end;
    } else {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(f, m, ix,iy,iz) = F3(ff, m, ix,iy,iz);
      } mrc_fld_foreach_end;
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}

void
copy_sc_to_sc_ggcm(struct mrc_fld *_f, struct mrc_fld *_ff)
{
  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int m = 0; m < _NR_FLDS; m++) {
    if (m == _B1X || m == _B2X) {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(ff, m, ix,iy,iz) = F3(f, m, ix+1,iy,iz);
      } mrc_fld_foreach_end;
    } else if (m == _B1Y || m == _B2Y) {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(ff, m, ix,iy,iz) = F3(f, m, ix,iy+1,iz);
      } mrc_fld_foreach_end;
    } else if (m == _B1Z || m == _B2Z) {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(ff, m, ix,iy,iz) = F3(f, m, ix,iy,iz+1);
      } mrc_fld_foreach_end;
    } else {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	F3(ff, m, ix,iy,iz) = F3(f, m, ix,iy,iz);
      } mrc_fld_foreach_end;
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}

