
#include "ggcm_mhd_bnd_private.h"

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>

#include <mrc_domain.h>
#include <mrc_bits.h>
#include <math.h>

enum {
  FIXED_RR,
  FIXED_VX,
  FIXED_VY,
  FIXED_VZ,
  FIXED_PP,
  FIXED_BX,
  FIXED_BY,
  FIXED_BZ,
  FIXED_NR,
};


// FIXME, consolidate with ggcm_mhd_iono

// ======================================================================
// ggcm_mhd_bnd subclass "sphere"

struct ggcm_mhd_bnd_sphere {
  // params
  double radius;
  double bnvals[FIXED_NR];  // constant values to set
  int test; // for testing, set to != 0
  int radial_velocity; // 0 : float, 1: reflect, 2: reflect if outflow

  // state
  struct ggcm_mhd_bnd_sphere_map map;
};

#define ggcm_mhd_bnd_sphere(bnd) mrc_to_subobj(bnd, struct ggcm_mhd_bnd_sphere)

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_setup

static void
ggcm_mhd_bnd_sphere_setup(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  ggcm_mhd_bnd_sphere_map_setup(map, bnd->mhd, sub->radius);
  ggcm_mhd_bnd_sphere_map_setup_flds(map);
  ggcm_mhd_bnd_setup_member_objs_sub(bnd);
  ggcm_mhd_bnd_sphere_map_setup_cc(map);
  ggcm_mhd_bnd_sphere_map_setup_ec(map);
  ggcm_mhd_bnd_sphere_map_setup_fc(map);
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts

static void
sphere_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld, int m)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;
  struct ggcm_mhd *mhd = bnd->mhd;

  double gamm = mhd->par.gamm;
  double bnvals[FIXED_NR];
  bnvals[FIXED_RR] = sub->bnvals[FIXED_RR] / mhd->rrnorm;
  bnvals[FIXED_VX] = sub->bnvals[FIXED_VX] / mhd->vvnorm;
  bnvals[FIXED_VY] = sub->bnvals[FIXED_VY] / mhd->vvnorm;
  bnvals[FIXED_VZ] = sub->bnvals[FIXED_VZ] / mhd->vvnorm;
  bnvals[FIXED_PP] = sub->bnvals[FIXED_PP] / mhd->ppnorm;
  bnvals[FIXED_BX] = sub->bnvals[FIXED_BX] / mhd->bbnorm;
  bnvals[FIXED_BY] = sub->bnvals[FIXED_BY] / mhd->bbnorm;
  bnvals[FIXED_BZ] = sub->bnvals[FIXED_BZ] / mhd->bbnorm;

  double rrbn = bnvals[FIXED_RR];
  double rvx = rrbn * bnvals[FIXED_VX];
  double rvy = rrbn * bnvals[FIXED_VY];
  double rvz = rrbn * bnvals[FIXED_VZ];

  double vvbn = sqr(bnvals[FIXED_VX]) + sqr(bnvals[FIXED_VY]) + sqr(bnvals[FIXED_VZ]);
  double uubn = .5f * (rrbn*vvbn) + bnvals[FIXED_PP] / (gamm - 1.f);
  double b2bn = sqr(bnvals[FIXED_BX]) + sqr(bnvals[FIXED_BY]) + sqr(bnvals[FIXED_BZ]);
  double eebn = uubn + .5 * b2bn;

  for (int i = 0; i < map->cc_n_map; i++) {
    int ix = MRC_I2(map->cc_imap, 0, i);
    int iy = MRC_I2(map->cc_imap, 1, i);
    int iz = MRC_I2(map->cc_imap, 2, i);
    int p  = MRC_I2(map->cc_imap, 3, i);

    M3 (fld, m + RR,  ix,iy,iz, p) = rrbn;
    M3 (fld, m + RVX, ix,iy,iz, p) = rvx;
    M3 (fld, m + RVY, ix,iy,iz, p) = rvy;
    M3 (fld, m + RVZ, ix,iy,iz, p) = rvz;

#if 0
    // FIXME, this is still kinda specific / hacky to ganymede
    // to avoid cutting off the initial perturbation from e.g., the mirror dipole,
    // let's just keep B untouched
    if (MT == MT_FULLY_CONSERVATIVE_CC) {
      M3(fld, m + BX , ix,iy,iz, p) = bnvals[FIXED_BX];
      M3(fld, m + BY , ix,iy,iz, p) = bnvals[FIXED_BY];
      M3(fld, m + BZ , ix,iy,iz, p) = bnvals[FIXED_BZ];
    } else {
      // we'd need to have a face-centered map to do this right,
      // but for now we'll do E field instead, anyway...
    }
#endif

    if (MT == MT_SEMI_CONSERVATIVE ||
        MT == MT_SEMI_CONSERVATIVE_GGCM) {
      M3(fld, m + UU , ix,iy,iz, p) = uubn;
    } else if (MT == MT_FULLY_CONSERVATIVE) {
      M3(fld, m + EE , ix,iy,iz, p) = eebn;
    } else if (MT == MT_FULLY_CONSERVATIVE_CC) {
      M3(fld, m + EE , ix,iy,iz, p) = uubn
	+ .5 * (sqr(M3(fld, m + BX, ix,iy,iz, p)) +
		sqr(M3(fld, m + BY, ix,iy,iz, p)) +
		sqr(M3(fld, m + BZ, ix,iy,iz, p)));
    } else {
      assert(0);
    }
  }
}

static void
sphere_fill_ghosts_test_1(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld, int m)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix - (1-bndp),iy,iz, p) = 0.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy - (1-bndp),iz, p) = 0.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz - (1-bndp), p) = 0.;
      }
    }
  }
  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix - (1-bndp),iy,iz, p) += 2.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy - (1-bndp),iz, p) += 2.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz - (1-bndp), p) += 2.;
      }
    }
  }
}

static void
sphere_fill_ghosts_test_2(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld, int m)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  struct mrc_fld *cnt = ggcm_mhd_get_3d_fld(bnd->mhd, 1);

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix-(1-bndp),iy,iz, p) = 0.;
	M3(cnt, 0, ix-(1-bndp),iy,iz, p) = 0.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy-(1-bndp),iz, p) = 0.;
	M3(cnt, 0, ix,iy-(1-bndp),iz, p) = 0.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz-(1-bndp), p) = 0.;
	M3(cnt, 0, ix,iy,iz-(1-bndp), p) = 0.;
      }
    }
  }
  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix-(1-bndp),iy,iz, p) += M3(fld, 0, ix-bndp,iy,iz, p);
	M3(cnt, 0, ix-(1-bndp),iy,iz, p) += 1.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy-(1-bndp),iz, p) += M3(fld, 0, ix,iy-bndp,iz, p);
	M3(cnt, 0, ix,iy-(1-bndp),iz, p) += 1.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz-(1-bndp), p) += M3(fld, 0, ix,iy,iz-bndp, p);
	M3(cnt, 0, ix,iy,iz-(1-bndp), p) += 1.;
      }
    }
  }

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int p  = MRC_I2(map->fc_imap[d], 3, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      if (d == 0) {
	M3(fld, 0, ix-(1-bndp),iy,iz, p) /= M3(cnt, 0, ix-(1-bndp),iy,iz, p);
	M3(cnt, 0, ix-(1-bndp),iy,iz, p) = 1.;
      } else if (d == 1) {
	M3(fld, 0, ix,iy-(1-bndp),iz, p) /= M3(cnt, 0, ix,iy-(1-bndp),iz, p);
	M3(cnt, 0, ix,iy-(1-bndp),iz, p) = 1.;
      } else if (d == 2) {
	M3(fld, 0, ix,iy,iz-(1-bndp), p) /= M3(cnt, 0, ix,iy,iz-(1-bndp), p);
	M3(cnt, 0, ix,iy,iz-(1-bndp), p) = 1.;
      }
    }
  }
  ggcm_mhd_put_3d_fld(bnd->mhd, cnt);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_fill_ghosts

static void
ggcm_mhd_bnd_sphere_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld_base,
			      int m, float bntim)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  if (map->cc_n_map == 0) {
    return;
  }

  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);
  assert(mhd_type == MT);
  assert(m == 0 || m == 8);

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);
  if (sub->test == 0) {
    sphere_fill_ghosts(bnd, fld, m);
  } else if (sub->test == 1) {
    sphere_fill_ghosts_test_1(bnd, fld, m);
  } else if (sub->test == 2) {
    sphere_fill_ghosts_test_2(bnd, fld, m);
  } else {
    assert(0);
  }
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts_E

static void
sphere_fill_ghosts_E(struct ggcm_mhd_bnd *bnd, struct mrc_fld *E)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->ec_n_map[d]; i++) {
      int ix = MRC_I2(map->ec_imap[d], 0, i);
      int iy = MRC_I2(map->ec_imap[d], 1, i);
      int iz = MRC_I2(map->ec_imap[d], 2, i);
      int p  = MRC_I2(map->ec_imap[d], 3, i);

      M3(E, d, ix,iy,iz, p) = 0.f;
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_fill_ghosts_E

static void _mrc_unused
ggcm_mhd_bnd_sphere_fill_ghosts_E(struct ggcm_mhd_bnd *bnd, struct mrc_fld *E_base)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  if (map->ec_n_map[0] + map->ec_n_map[1] + map->ec_n_map[2] == 0) {
    return;
  }

  struct mrc_fld *E = mrc_fld_get_as(E_base, FLD_TYPE);
  sphere_fill_ghosts_E(bnd, E);
  mrc_fld_put_as(E, E_base);
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts_reconstr

static void
sphere_fill_ghosts_reconstr(struct ggcm_mhd_bnd *bnd, struct mrc_fld *U_l[],
			    struct mrc_fld *U_r[], int p)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;
  struct mrc_crds *crds = mrc_domain_get_crds(bnd->mhd->domain);

  for (int d = 0; d < 3; d++) {
    for (int i = 0; i < map->fc_n_map[d]; i++) {
      // FIXME, we really should have maps per patch, or some other way
      // to not go through the entire list here
      if (p != MRC_I2(map->fc_imap[d], 3, i)) {
	continue;
      }

      int ix = MRC_I2(map->fc_imap[d], 0, i);
      int iy = MRC_I2(map->fc_imap[d], 1, i);
      int iz = MRC_I2(map->fc_imap[d], 2, i);
      int bndp = MRC_I2(map->fc_imap[d], 4, i);

      float crd_fc[3];
      mrc_crds_at_fc(crds, ix,iy,iz, p, d, crd_fc);

      int n_comps = mrc_fld_nr_comps(U_l[0]);
      // find true (inside) face values
      mrc_fld_data_t U_ghost[n_comps], U_true[n_comps];
      for (int m = 0; m < n_comps; m++) {
	if (bndp) {
	  U_true[m] = M3(U_l[d], m, ix,iy,iz, p);
	} else {
	  U_true[m] = M3(U_r[d], m, ix,iy,iz, p);
	}
      }

      // from BATSRUS / ganymede b.c.
      // for inflow float everything 
      for (int m = 0; m < n_comps; m++) {
	U_ghost[m] = U_true[m];
      }

      if (sub->radial_velocity != 0) {
	// calculate r^2
	mrc_fld_data_t r2 = 0.;
	for (int m = 0; m < 3; m++) {
	  r2 += sqr(crd_fc[m]);
	}
	
	// for outflow reflect radial velocity: uG = u - 2*(u.r)*r/r^2
	mrc_fld_data_t UdotR = 0.;
	for (int m = 0; m < 3; m++) {
	  UdotR += U_true[RVX + m] * crd_fc[m];
	}
	if (sub->radial_velocity == 1 || UdotR > 0.) {
	  for (int m = 0; m < 3; m++) {
	    U_ghost[RVX + m] = U_true[RVX + m] - 2. * UdotR / r2 * crd_fc[m];
	  }
	}
      }

      // store ghost values back
      for (int m = 0; m < n_comps; m++) {
	if (bndp) {
	  M3(U_r[d], m, ix, iy,iz, p) = U_ghost[m];
	} else {
	  M3(U_l[d], m, ix, iy,iz, p) = U_ghost[m];
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_fill_ghosts_reconstr

static void
ggcm_mhd_bnd_sphere_fill_ghosts_reconstr(struct ggcm_mhd_bnd *bnd, struct mrc_fld *U_l_base[],
					 struct mrc_fld *U_r_base[], int p)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd_bnd_sphere_map *map = &sub->map;

  if (map->fc_n_map[0] + map->fc_n_map[1] + map->fc_n_map[2] == 0) {
    return;
  }

  struct mrc_fld *U_l[3], *U_r[3];
  for (int d = 0; d < 3; d++) {
    U_l[d] = mrc_fld_get_as(U_l_base[d], FLD_TYPE);
    U_r[d] = mrc_fld_get_as(U_r_base[d], FLD_TYPE);
  }

  for (int p = 0; p < mrc_fld_nr_patches(U_l[0]); p++) {
    sphere_fill_ghosts_reconstr(bnd, U_l, U_r, p);
  }

  for (int d = 0; d < 3; d++) {
    mrc_fld_put_as(U_l[d], U_l_base[d]);
    mrc_fld_put_as(U_r[d], U_r_base[d]);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd "sphere" subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bnd_sphere, x)
static struct param ggcm_mhd_bnd_sphere_descr[] = {
  { "radius"          , VAR(radius)          , PARAM_DOUBLE(1.)          },
  { "rr"              , VAR(bnvals[FIXED_RR]), PARAM_DOUBLE(1.)          },
  { "pp"              , VAR(bnvals[FIXED_PP]), PARAM_DOUBLE(1.)          },
  { "vx"              , VAR(bnvals[FIXED_VX]), PARAM_DOUBLE(0.)          },
  { "vy"              , VAR(bnvals[FIXED_VY]), PARAM_DOUBLE(0.)          },
  { "vz"              , VAR(bnvals[FIXED_VZ]), PARAM_DOUBLE(0.)          },
  { "bx"              , VAR(bnvals[FIXED_BX]), PARAM_DOUBLE(0.)          },
  { "by"              , VAR(bnvals[FIXED_BY]), PARAM_DOUBLE(0.)          },
  { "bz"              , VAR(bnvals[FIXED_BZ]), PARAM_DOUBLE(0.)          },
  { "test"            , VAR(test),             PARAM_INT(0)              },
  { "radial_velocity" , VAR(radial_velocity),  PARAM_INT(0)              },

  { "dr"              , VAR(map.dr)          , PARAM_DOUBLE(.01)         },
  { "extra_dr"        , VAR(map.extra_dr)    , PARAM_DOUBLE(0.)          },
  { "min_dr"          , VAR(map.min_dr)      , MRC_VAR_DOUBLE            },
  { "radius"          , VAR(map.radius)      , MRC_VAR_DOUBLE            },
  { "r1"              , VAR(map.r1)          , MRC_VAR_DOUBLE            },
  { "cc_n_map"        , VAR(map.cc_n_map)    , MRC_VAR_INT               },
  { "cc_imap"         , VAR(map.cc_imap)     , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_n_map[0]"     , VAR(map.ec_n_map[0]) , MRC_VAR_INT               },
  { "ec_n_map[1]"     , VAR(map.ec_n_map[1]) , MRC_VAR_INT               },
  { "ec_n_map[2]"     , VAR(map.ec_n_map[2]) , MRC_VAR_INT               },
  { "ec_imap[0]"      , VAR(map.ec_imap[0])  , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_imap[1]"      , VAR(map.ec_imap[1])  , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_imap[2]"      , VAR(map.ec_imap[2])  , MRC_VAR_OBJ(mrc_fld)      },
  { "fc_n_map[0]"     , VAR(map.fc_n_map[0]) , MRC_VAR_INT               },
  { "fc_n_map[1]"     , VAR(map.fc_n_map[1]) , MRC_VAR_INT               },
  { "fc_n_map[2]"     , VAR(map.fc_n_map[2]) , MRC_VAR_INT               },
  { "fc_imap[0]"      , VAR(map.fc_imap[0])  , MRC_VAR_OBJ(mrc_fld)      },
  { "fc_imap[1]"      , VAR(map.fc_imap[1])  , MRC_VAR_OBJ(mrc_fld)      },
  { "fc_imap[2]"      , VAR(map.fc_imap[2])  , MRC_VAR_OBJ(mrc_fld)      },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_bnd subclass "sphere"

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere = {
  .name             = ggcm_mhd_bnd_sub_name,
  .size             = sizeof(struct ggcm_mhd_bnd_sphere),
  .param_descr      = ggcm_mhd_bnd_sphere_descr,
  .setup            = ggcm_mhd_bnd_sphere_setup,
  .fill_ghosts      = ggcm_mhd_bnd_sphere_fill_ghosts,
  //  .fill_ghosts_E    = ggcm_mhd_bnd_sphere_fill_ghosts_E,
  .fill_ghosts_reconstr = ggcm_mhd_bnd_sphere_fill_ghosts_reconstr,
};

