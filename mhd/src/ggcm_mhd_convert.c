
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include <ggcm_mhd_gkeyll.h>

#include <mrc_domain.h>
#include <mrc_fld_as_double.h>

#include "pde/pde_defs.h"
#include "pde/pde_mhd_convert.c"

#include <string.h>
#include <stdlib.h>
#include <execinfo.h>

// ----------------------------------------------------------------------
// ggcm_mhd_convert_sc_ggcm_from_primitive
//
// converts from primitive variables to semi-conservative in-place.
// No ghost points are set.

static void
ggcm_mhd_convert_sc_ggcm_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 1) {
      mrc_fld_data_t prim[8], state[8];
      for (int m = 0; m < 8; m++) {
	prim[m] = M3(fld, m, ix,iy,iz, p);
      }
      convert_scons_from_prim(state, prim);
      for (int m = 0; m < 5; m++) {
	M3(fld, m, ix,iy,iz, p) = state[m];
      }
    } mrc_fld_foreach_end;
  }
  
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_sc_from_primitive
//
// converts from primitive variables to semi-conservative alt B in-place.
// No ghost points are set.

static void
ggcm_mhd_convert_sc_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t prim[8], state[8];
      for (int m = 0; m < 8; m++) {
	prim[m] = M3(fld, m, ix,iy,iz, p);
      }
      convert_scons_from_prim(state, prim);
      for (int m = 0; m < 5; m++) {
	M3(fld, m, ix,iy,iz, p) = state[m];
      }
    } mrc_fld_foreach_end;
  }
  
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_fc_from_primitive
//
// converts from primitive variables to fully-conservative in-place.
// No ghost points are set, the staggered B fields need to exist on all faces
// (that means one more than cell-centered dims)

static void
ggcm_mhd_convert_fc_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  // don't go into ghost cells in invariant directions
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] == 1) ? 0 : 1;
  int dy = (gdims[1] == 1) ? 0 : 1;
  int dz = (gdims[2] == 1) ? 0 : 1;

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t prim[8], state[8];
      for (int m = 0; m < 5; m++) {
	prim[m] = M3(fld, m, ix,iy,iz, p);
      }
      prim[BX] = .5f*(BX_(fld, ix,iy,iz, p) + BX_(fld, ix+dx,iy,iz, p));
      prim[BY] = .5f*(BY_(fld, ix,iy,iz, p) + BY_(fld, ix,iy+dy,iz, p));
      prim[BZ] = .5f*(BZ_(fld, ix,iy,iz, p) + BZ_(fld, ix,iy,iz+dz, p));
      convert_fcons_from_prim(state, prim);
      for (int m = 0; m < 5; m++) {
	M3(fld, m, ix,iy,iz, p) = state[m];
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_fc_cc_from_primitive
//
// converts from primitive variables to fully-conservative cell-centered in-place.
// No ghost points are set, the staggered B fields need to exist on all faces
// (that means one more than cell-centered dims)

static void
ggcm_mhd_convert_fc_cc_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  // don't go into ghost cells in invariant directions
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] == 1) ? 0 : 1;
  int dy = (gdims[1] == 1) ? 0 : 1;
  int dz = (gdims[2] == 1) ? 0 : 1;

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t prim[8], state[8];
      for (int m = 0; m < 5; m++) {
	prim[m] = M3(fld, m, ix,iy,iz, p);
      }
      prim[BX] = .5f*(BX_(fld, ix,iy,iz, p) + BX_(fld, ix+dx,iy,iz, p));
      prim[BY] = .5f*(BY_(fld, ix,iy,iz, p) + BY_(fld, ix,iy+dy,iz, p));
      prim[BZ] = .5f*(BZ_(fld, ix,iy,iz, p) + BZ_(fld, ix,iy,iz+dz, p));
      convert_fcons_from_prim(state, prim);
      for (int m = 0; m < 8; m++) {
	M3(fld, m, ix,iy,iz, p) = state[m];
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_gkeyll_from_primitive
//
// converts from primitive MHD variables to multi-fluid moment variables
// No ghost points are set.

// pointwise conversion from primitive mhd quantities to 5m fluid quantities
void
primitive_to_gkeyll_5m_fluids_point(struct mrc_fld *fld, int nr_fluids, int idx[],
    float mass[], float charge[], float pressure_ratios[],
    float gamma_m1, int ix, int iy, int iz, int p)
{
  mrc_fld_data_t rr  = M3(fld, RR, ix,iy,iz, p);
  mrc_fld_data_t rvx = rr * M3(fld, VX, ix,iy,iz, p);
  mrc_fld_data_t rvy = rr * M3(fld, VY, ix,iy,iz, p);
  mrc_fld_data_t rvz = rr * M3(fld, VZ, ix,iy,iz, p);
  mrc_fld_data_t pp  = M3(fld, PP, ix,iy,iz, p);

  float mass_ratios[nr_fluids];
  float momentum_ratios[nr_fluids];
  float mass_total = 0.;
  for (int s = 0; s < nr_fluids; s++) {
     mass_total += mass[s];
  }
  for (int s = 0; s < nr_fluids; s++) {
     mass_ratios[s] = mass[s]/mass_total;
     // assume no current and all plasma co-move
     momentum_ratios[s] = mass_ratios[s];
  }

  // each species
  for (int s = 0; s < nr_fluids; s++) {
    M3(fld, idx[s] + G5M_RRS,  ix,iy,iz, p) = rr   * mass_ratios[s];
    M3(fld, idx[s] + G5M_RVXS, ix,iy,iz, p) = rvx * momentum_ratios[s];
    M3(fld, idx[s] + G5M_RVYS, ix,iy,iz, p) = rvy * momentum_ratios[s];
    M3(fld, idx[s] + G5M_RVZS, ix,iy,iz, p) = rvz * momentum_ratios[s];
    M3(fld, idx[s] + G5M_UUS,  ix,iy,iz, p) = 
      pp * pressure_ratios[s] / gamma_m1
      + .5 * (sqr(M3(fld, idx[s] + G5M_RVXS, ix,iy,iz, p))
          +   sqr(M3(fld, idx[s] + G5M_RVYS, ix,iy,iz, p))
          +   sqr(M3(fld, idx[s] + G5M_RVZS, ix,iy,iz, p)))
             / M3(fld, idx[s] + G5M_RRS, ix,iy,iz,p);
  }
}

// pointwise conversion from primitive mhd quantities to em fields
void
primitive_to_gkeyll_em_fields_point(struct mrc_fld *fld, int idx_em,
    int dx, int dy, int dz, int ix, int iy, int iz, int p)
{
  mrc_fld_data_t vx = M3(fld, VX, ix,iy,iz, p);
  mrc_fld_data_t vy = M3(fld, VY, ix,iy,iz, p);
  mrc_fld_data_t vz = M3(fld, VZ, ix,iy,iz, p);

  // staggered to cell-center B fields
  M3(fld, idx_em + GK_BX, ix,iy,iz, p) = 
    .5 * (M3(fld, BX, ix,iy,iz, p) + M3(fld, BX, ix+dx,iy,iz, p));
  M3(fld, idx_em + GK_BY, ix,iy,iz, p) = 
    .5 * (M3(fld, BY, ix,iy,iz, p) + M3(fld, BY, ix,iy+dy,iz, p));
  M3(fld, idx_em + GK_BZ, ix,iy,iz, p) = 
    .5 * (M3(fld, BZ, ix,iy,iz, p) + M3(fld, BZ, ix,iy,iz+dz, p));

  // E=-vxB, i.e., only convection E field
  M3(fld, idx_em + GK_EX, ix,iy,iz, p) =
    - vy * M3(fld, idx_em + GK_BZ, ix,iy,iz, p)
    + vz * M3(fld, idx_em + GK_BY, ix,iy,iz, p);
  M3(fld, idx_em + GK_EY, ix,iy,iz, p) =
    - vz * M3(fld, idx_em + GK_BX, ix,iy,iz, p)
    + vx * M3(fld, idx_em + GK_BZ, ix,iy,iz, p);
  M3(fld, idx_em + GK_EZ, ix,iy,iz, p) =
    - vx * M3(fld, idx_em + GK_BY, ix,iy,iz, p)
    + vy * M3(fld, idx_em + GK_BX, ix,iy,iz, p);

  // FIXME sensible correction potentials?
  // e.g., copy ghost vals for inoutflow bnd
  M3(fld, idx_em + GK_PHI, ix,iy,iz, p) = .0;
  M3(fld, idx_em + GK_PSI, ix,iy,iz, p) = .0;
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_primitive_gkeyll_5m_point
//
// pointwise conversion from primitive mhd quantities to 5m-em quantities

void
ggcm_mhd_convert_primitive_gkeyll_5m_point(struct mrc_fld *fld, int nr_fluids,
    int idx[], float mass[], float charge[],
    float pressure_ratios[], float gamma_m1, int idx_em, int dx, int dy, int dz,
    int ix, int iy, int iz, int p)
{
  // em fields should be calculated before V used for E=-VxB are overwritten
  primitive_to_gkeyll_em_fields_point(fld, idx_em, dx, dy, dz, ix, iy, iz, p);
  primitive_to_gkeyll_5m_fluids_point(fld, nr_fluids, idx, mass, charge,
      pressure_ratios, gamma_m1, ix, iy, iz, p);
}

static void
ggcm_mhd_convert_gkeyll_from_primitive(struct ggcm_mhd *mhd, 
    struct mrc_fld *fld_base)
{
  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  int nr_moments = ggcm_mhd_gkeyll_nr_moments(mhd);
  int nr_fluids = ggcm_mhd_gkeyll_nr_fluids(mhd);

  int idx[nr_fluids];
  ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);
  int idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);

  float *mass = ggcm_mhd_gkeyll_mass(mhd);
  float *charge = ggcm_mhd_gkeyll_charge(mhd);
  float *pressure_ratios = ggcm_mhd_gkeyll_pressure_ratios(mhd);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] == 1) ? 0 : 1;
  int dy = (gdims[1] == 1) ? 0 : 1;
  int dz = (gdims[2] == 1) ? 0 : 1;

  if (nr_moments == 5) {
    float gamma_m1 = mhd->par.gamm - 1.;
    for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
      mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
        ggcm_mhd_convert_primitive_gkeyll_5m_point(fld, nr_fluids, idx,
            mass, charge, pressure_ratios, gamma_m1, idx_em,
            dx,dy,dz, ix,iy,iz, p);
      } mrc_fld_foreach_end;
    }
  } else if (nr_moments == 10) {
    // TODO
  }
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_from_primitive
//
// converts from primitive variables to the appropriate fully-conservative /
// semi-conservative state vector.

void
ggcm_mhd_convert_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);

  if (mhd_type == MT_SCONS_FC_GGCM) {
    return ggcm_mhd_convert_sc_ggcm_from_primitive(mhd, fld_base);
  } else if (mhd_type == MT_FCONS_FC) {
    return ggcm_mhd_convert_fc_from_primitive(mhd, fld_base);
  } else if (mhd_type == MT_SCONS_FC) {
    return ggcm_mhd_convert_sc_from_primitive(mhd, fld_base);
  } else if (mhd_type == MT_FCONS_CC) {
    return ggcm_mhd_convert_fc_cc_from_primitive(mhd, fld_base);
  } else if (mhd_type == MT_GKEYLL) {
    return ggcm_mhd_convert_gkeyll_from_primitive(mhd, fld_base);
  } else {
    assert(0);
  }
}

// ======================================================================

void copy_sc_ggcm_to_sc_float(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_sc_ggcm_to_fc_float(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_sc_to_sc_ggcm_float(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_fc_to_sc_ggcm_float(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_fc_to_sc_float(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_sc_to_fc_float(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);

void copy_sc_ggcm_to_sc_double(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_sc_ggcm_to_fc_double(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_sc_to_sc_ggcm_double(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_fc_to_sc_ggcm_double(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_fc_to_sc_double(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);
void copy_sc_to_fc_double(struct mrc_fld *_ff, struct mrc_fld *_f, int mb, int me);

// ----------------------------------------------------------------------
// ggcm_mhd_fld_get_as
//
// get mhd::fld with specified data type and mhd_type

struct mrc_fld *
ggcm_mhd_fld_get_as(struct mrc_fld *fld_base, const char *type,
		    int mhd_type, int mb, int me)
{
  assert(me <= mrc_fld_nr_comps(fld_base));

  int mhd_type_base;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type_base);

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, type);
  if (fld != fld_base) {
    mrc_fld_dict_add_int(fld, "mhd_type", mhd_type_base);
  }

  if (mhd_type == mhd_type_base) {
    return fld;
  }

#if 0
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      mprintf("XXXXXX ggcm_mhd_fld_get_as\n");
      void* callstack[128];
      int frames = backtrace(callstack, 128);
      char** strs = backtrace_symbols(callstack, frames);
      for (int i = 0; i < frames; i++) {
	mprintf("%s\n", strs[i]);
      }
      free(strs);
    }
#endif

  struct mrc_fld *fld2 = mrc_fld_create(mrc_fld_comm(fld_base));
  mrc_fld_set_type(fld2, type);
  mrc_fld_set_param_obj(fld2, "domain", fld->_domain);
  mrc_fld_set_param_int(fld2, "nr_spatial_dims", 3);
  mrc_fld_set_param_int(fld2, "nr_comps", me);
  mrc_fld_set_param_int(fld2, "nr_ghosts", 2);
  mrc_fld_dict_add_int(fld2, "mhd_type", mhd_type);
  mrc_fld_setup(fld2);

  if (mhd_type_base == MT_SCONS_FC_GGCM && mhd_type == MT_SCONS_FC) {
    if (strcmp(type, "float") == 0) {
      copy_sc_ggcm_to_sc_float(fld, fld2, mb, me);
    } else {
      assert(0);
    }
    mrc_fld_put_as(fld, fld_base);
    return fld2;
  } else if (mhd_type_base == MT_SCONS_FC && mhd_type == MT_SCONS_FC_GGCM) {
    if (strcmp(type, "float") == 0) {
      copy_sc_to_sc_ggcm_float(fld, fld2, mb, me);
    } else if (strcmp(type, "double") == 0) {
      copy_sc_to_sc_ggcm_double(fld, fld2, mb, me);
    } else {
      assert(0);
    }
    mrc_fld_put_as(fld, fld_base);
    return fld2;
  } else if (mhd_type_base == MT_FCONS_FC && mhd_type == MT_SCONS_FC_GGCM) {
    if (strcmp(type, "float") == 0) {
      copy_fc_to_sc_ggcm_float(fld, fld2, mb, me);
    } else if (strcmp(type, "double") == 0) {
      copy_fc_to_sc_ggcm_double(fld, fld2, mb, me);
    } else {
      assert(0);
    }
    mrc_fld_put_as(fld, fld_base);
    return fld2;
  } else if (mhd_type_base == MT_FCONS_FC && mhd_type == MT_SCONS_FC) {
    if (strcmp(type, "float") == 0) {
      copy_fc_to_sc_float(fld, fld2, mb, me);
    } else if (strcmp(type, "double") == 0) {
      copy_fc_to_sc_double(fld, fld2, mb, me);
    } else {
      assert(0);
    }
    mrc_fld_put_as(fld, fld_base);
    return fld2;
  } else {
    mprintf("ggcm_mhd_fld_get_as %d -> %d\n", mhd_type_base, mhd_type);
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_fld_put_as
//
// called when one is done with the field from _get_as()

void
ggcm_mhd_fld_put_as(struct mrc_fld *fld, struct mrc_fld *fld_base, int mb, int me)
{
  int mhd_type_base;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type_base);

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  if (mhd_type == mhd_type_base) {
    mrc_fld_put_as(fld, fld_base);
    return;
  }

  if (mhd_type_base == MT_SCONS_FC && mhd_type == MT_SCONS_FC_GGCM) {
    if (strcmp(mrc_fld_type(fld), "float") == 0) {
      copy_sc_ggcm_to_sc_float(fld, fld_base, mb, me);
    } else if (strcmp(mrc_fld_type(fld), "double") == 0) {
      copy_sc_ggcm_to_sc_double(fld, fld_base, mb, me);
    } else {
      assert(0);
    }
    mrc_fld_destroy(fld);
  } else if (mhd_type_base == MT_FCONS_FC && mhd_type == MT_SCONS_FC_GGCM) {
    if (strcmp(mrc_fld_type(fld), "float") == 0) {
      copy_sc_ggcm_to_fc_float(fld, fld_base, mb, me);
    } else if (strcmp(mrc_fld_type(fld), "double") == 0) {
      copy_sc_ggcm_to_fc_double(fld, fld_base, mb, me);
    } else {
      assert(0);
    }
    mrc_fld_destroy(fld);
  } else if (mhd_type_base == MT_FCONS_FC && mhd_type == MT_SCONS_FC) {
    if (strcmp(mrc_fld_type(fld), "float") == 0) {
      copy_sc_to_fc_float(fld, fld_base, mb, me);
    } else if (strcmp(mrc_fld_type(fld), "double") == 0) {
      copy_sc_to_fc_double(fld, fld_base, mb, me);
    } else {
      assert(0);
    }
    mrc_fld_destroy(fld);
  } else {
    mprintf("ggcm_mhd_fld_put_as %d <- %d\n", mhd_type_base, mhd_type);
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_fld_as_fortran
//
// get mhd::fld as Fortran common block field

struct mrc_fld *
ggcm_mhd_get_fld_as_fortran(struct mrc_fld *fld_base)
{
#if 0
  MHERE;
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    mprintf("XXXXXXX get_fld_as_fortran\n");
    void* callstack[128];
    int frames = backtrace(callstack, 128);
    char** strs = backtrace_symbols(callstack, frames);
    for (int i = 0; i < frames; i++) {
      mprintf("%s\n", strs[i]);
    }
    free(strs);
  }
#endif
  return ggcm_mhd_fld_get_as(fld_base, "float", MT_SCONS_FC_GGCM, 0, _NR_FLDS);
}

// ----------------------------------------------------------------------
// ggcm_mhd_put_fld_as_fortran

void
ggcm_mhd_put_fld_as_fortran(struct mrc_fld *fld, struct mrc_fld *fld_base)
{
  return ggcm_mhd_fld_put_as(fld, fld_base, 0, _NR_FLDS);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_free

// FIXME, never actually called
void
ggcm_mhd_convert_free(void)
{
  pde_free();
}
