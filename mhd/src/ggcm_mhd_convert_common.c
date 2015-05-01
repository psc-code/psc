
#include <ggcm_mhd_defs.h>
#include <ggcm_mhd_defs_extra.h>

#include <mrc_fld.h>
#include <mrc_domain.h>

void
copy_sc_ggcm_to_sc(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me)
{
  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    for (int m = mb; m < me; m++) {
      if (m == _B1X || m == _B2X) {
	mrc_fld_foreach(f, ix,iy,iz, 1, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix-1,iy,iz, p);
	} mrc_fld_foreach_end;
      } else if (m == _B1Y || m == _B2Y) {
	mrc_fld_foreach(f, ix,iy,iz, 1, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix,iy-1,iz, p);
      } mrc_fld_foreach_end;
      } else if (m == _B1Z || m == _B2Z) {
	mrc_fld_foreach(f, ix,iy,iz, 1, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix,iy,iz-1, p);
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix,iy,iz, p);
	} mrc_fld_foreach_end;
      }
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}

void
copy_sc_ggcm_to_fc(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me)
{
  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    for (int m = mb; m < me; m++) {
      if (m == _B1X || m == _B2X) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix-1,iy,iz, p);
	} mrc_fld_foreach_end;
      } else if (m == _B1Y || m == _B2Y) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix,iy-1,iz, p);
	} mrc_fld_foreach_end;
      } else if (m == _B1Z || m == _B2Z) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix,iy,iz-1, p);
	} mrc_fld_foreach_end;
      } else if (m == _UU1) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
	  float b2  = (sqr(.5f * (M3(ff, _B1X, ix,iy,iz, p) + M3(ff, _B1X, ix+1,iy  ,iz  , p))) +
		       sqr(.5f * (M3(ff, _B1Y, ix,iy,iz, p) + M3(ff, _B1Y, ix  ,iy+1,iz  , p))) +
		       sqr(.5f * (M3(ff, _B1Z, ix,iy,iz, p) + M3(ff, _B1Z, ix  ,iy  ,iz+1, p))));
	  M3(f, _UU1, ix,iy,iz, p) = M3(ff, _UU1, ix,iy,iz, p) + .5 * b2;
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(f, m, ix,iy,iz, p) = M3(ff, m, ix,iy,iz, p);
	} mrc_fld_foreach_end;
      }
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}

void
copy_sc_to_sc_ggcm(struct mrc_fld *_f, struct mrc_fld *_ff, int mb, int me)
{
  int gdims[3];
  mrc_domain_get_global_dims(_ff->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    for (int m = mb; m < me; m++) {
      if (m == _B1X || m == _B2X) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix+dx,iy,iz, p);
	} mrc_fld_foreach_end;
      } else if (m == _B1Y || m == _B2Y) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy+dy,iz, p);
	} mrc_fld_foreach_end;
      } else if (m == _B1Z || m == _B2Z) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy,iz+dz, p);
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy,iz, p);
	} mrc_fld_foreach_end;
      }
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}

void
copy_fc_to_sc_ggcm(struct mrc_fld *_f, struct mrc_fld *_ff, int mb, int me)
{
  int gdims[3];
  mrc_domain_get_global_dims(_ff->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    for (int m = mb; m < me; m++) {
      if (m == _B1X || m == _B2X) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix+dx,iy,iz, p);
	} mrc_fld_foreach_end;
      } else if (m == _B1Y || m == _B2Y) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy+dy,iz, p);
	} mrc_fld_foreach_end;
      } else if (m == _B1Z || m == _B2Z) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy,iz+dz, p);
	} mrc_fld_foreach_end;
      } else if (m == _UU1) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  float b2  = (sqr(.5f * (M3(f, _B1X, ix,iy,iz, p) + M3(f, _B1X, ix+dx,iy   ,iz   , p))) +
		       sqr(.5f * (M3(f, _B1Y, ix,iy,iz, p) + M3(f, _B1Y, ix   ,iy+dy,iz   , p))) +
		       sqr(.5f * (M3(f, _B1Z, ix,iy,iz, p) + M3(f, _B1Z, ix   ,iy   ,iz+dz, p))));
	  M3(ff, _UU1, ix,iy,iz, p) = M3(f, _UU1, ix,iy,iz, p) - .5 * b2;
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy,iz, p);
	} mrc_fld_foreach_end;
      }
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}

void
copy_fc_to_sc(struct mrc_fld *_f, struct mrc_fld *_ff, int mb, int me)
{
  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    for (int m = mb; m < me; m++) {
      if (m == _UU1) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  float b2  = (sqr(.5f * (M3(f, _B1X, ix,iy,iz, p) + M3(f, _B1X, ix+1,iy  ,iz  , p))) +
		       sqr(.5f * (M3(f, _B1Y, ix,iy,iz, p) + M3(f, _B1Y, ix  ,iy+1,iz  , p))) +
		       sqr(.5f * (M3(f, _B1Z, ix,iy,iz, p) + M3(f, _B1Z, ix  ,iy  ,iz+1, p))));
	M3(ff, _UU1, ix,iy,iz, p) = M3(f, _UU1, ix,iy,iz, p) - .5 * b2;
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy,iz, p);
	} mrc_fld_foreach_end;
      }
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}

void
copy_sc_to_fc(struct mrc_fld *_f, struct mrc_fld *_ff, int mb, int me)
{
  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);
  struct mrc_fld *ff = mrc_fld_get_as(_ff, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    for (int m = mb; m < me; m++) {
      if (m == _UU1) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
	  float b2  = (sqr(.5f * (M3(f, _B1X, ix,iy,iz, p) + M3(f, _B1X, ix+1,iy  ,iz  , p))) +
		       sqr(.5f * (M3(f, _B1Y, ix,iy,iz, p) + M3(f, _B1Y, ix  ,iy+1,iz  , p))) +
		       sqr(.5f * (M3(f, _B1Z, ix,iy,iz, p) + M3(f, _B1Z, ix  ,iy  ,iz+1, p))));
	  M3(ff, _UU1, ix,iy,iz, p) = M3(f, _UU1, ix,iy,iz, p) + .5 * b2;
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  M3(ff, m, ix,iy,iz, p) = M3(f, m, ix,iy,iz, p);
	} mrc_fld_foreach_end;
      }
    }
  }

  mrc_fld_put_as(f, _f);
  mrc_fld_put_as(ff, _ff);
}
