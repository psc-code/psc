
#include "ggcm_mhd_bnd_private.h"

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>

#include <mrc_domain.h>
#include <mrc_bits.h>
#include <math.h>

enum {
  FIXED_RR,
  FIXED_PP,
  FIXED_VX,
  FIXED_VY,
  FIXED_VZ,
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
  double *bnvals = sub->bnvals;

  double rvx = bnvals[FIXED_RR] * bnvals[FIXED_VX];
  double rvy = bnvals[FIXED_RR] * bnvals[FIXED_VY];
  double rvz = bnvals[FIXED_RR] * bnvals[FIXED_VZ];

  double vvbn  = sqr(bnvals[FIXED_VX]) + sqr(bnvals[FIXED_VY]) + sqr(bnvals[FIXED_VZ]);
  double uubn  = .5f * (bnvals[FIXED_RR]*vvbn) + bnvals[FIXED_PP] / (gamm - 1.f);
  double b2bn  = sqr(bnvals[FIXED_BX]) + sqr(bnvals[FIXED_BY]) + sqr(bnvals[FIXED_BZ]);
  double eebn = uubn + .5 * b2bn;

  for (int i = 0; i < map->cc_n_map; i++) {
    int ix = MRC_I2(map->cc_imap, 0, i);
    int iy = MRC_I2(map->cc_imap, 1, i);
    int iz = MRC_I2(map->cc_imap, 2, i);
    int p  = MRC_I2(map->cc_imap, 3, i);

    M3 (fld, m + RR,  ix,iy,iz, p) = bnvals[FIXED_RR];
    M3 (fld, m + RVX, ix,iy,iz, p) = rvx;
    M3 (fld, m + RVY, ix,iy,iz, p) = rvy;
    M3 (fld, m + RVZ, ix,iy,iz, p) = rvz;
    if (MT == MT_SEMI_CONSERVATIVE ||
        MT == MT_SEMI_CONSERVATIVE_GGCM) {
      M3(fld, m + UU , ix,iy,iz, p) = uubn;
    } else if (MT == MT_FULLY_CONSERVATIVE) {
      M3(fld, m + EE , ix,iy,iz, p) = eebn;
    } else {
      assert(0);
    }
#if 0
    // we'd need to have a face-centered map to do this right,
    // but for now we'll do E field instead, anyway
    M3(fld, m + BX , ix,iy,iz, p) = bnvals[FIXED_BX];
    M3(fld, m + BY , ix,iy,iz, p) = bnvals[FIXED_BY];
    M3(fld, m + BZ , ix,iy,iz, p) = bnvals[FIXED_BZ];
#endif
  }
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
  sphere_fill_ghosts(bnd, fld, m);
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

static void
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

  { "min_dr"          , VAR(map.min_dr)      , MRC_VAR_DOUBLE            },
  { "r1"              , VAR(map.r1)          , MRC_VAR_DOUBLE            },
  { "r2"              , VAR(map.r2)          , MRC_VAR_DOUBLE            },
  { "cc_n_map"        , VAR(map.cc_n_map)    , MRC_VAR_INT               },
  { "cc_imap"         , VAR(map.cc_imap)     , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_n_map[0]"     , VAR(map.ec_n_map[0]) , MRC_VAR_INT               },
  { "ec_n_map[1]"     , VAR(map.ec_n_map[1]) , MRC_VAR_INT               },
  { "ec_n_map[2]"     , VAR(map.ec_n_map[2]) , MRC_VAR_INT               },
  { "ec_imap[0]"      , VAR(map.ec_imap[0])  , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_imap[1]"      , VAR(map.ec_imap[1])  , MRC_VAR_OBJ(mrc_fld)      },
  { "ec_imap[2]"      , VAR(map.ec_imap[2])  , MRC_VAR_OBJ(mrc_fld)      },

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
  .fill_ghosts_E    = ggcm_mhd_bnd_sphere_fill_ghosts_E,
};

