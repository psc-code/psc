
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "psc_fields_cuda.h"
#include "bnd_fields.hxx"

#include "psc.h"
#include <mrc_profile.h>

struct BndFieldsCuda : BndFieldsBase
{
  // ----------------------------------------------------------------------
  // fill_ghosts_E

  void fill_ghosts_E(MfieldsBase& mflds_base) override
  {
    const auto& grid = mflds_base.grid();
    if (grid.bc.fld_lo[0] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[1] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[2] == BND_FLD_PERIODIC) {
      return;
    }

    auto& mflds = mflds_base.get_as<MfieldsCuda>(EX, EX + 3);
    if (grid.bc.fld_lo[0] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
	grid.bc.fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
	grid.bc.fld_lo[2] == BND_FLD_PERIODIC) {
      int d = 1;
      for (int p = 0; p < mflds.n_patches(); p++) {
	if (grid.atBoundaryLo(p, d)) {
	  cuda_conducting_wall_E_lo_y(mflds.cmflds(), p);
	}
	if (grid.atBoundaryHi(p, d)) {
	  cuda_conducting_wall_E_hi_y(mflds.cmflds(), p);
	}
      }
    } else {
      assert(0);
    }
    mflds_base.put_as(mflds, EX, EX + 3);
  }

  // ----------------------------------------------------------------------
  // fill_ghosts_H

  void fill_ghosts_H(MfieldsBase& mflds_base) override
  {
    const auto& grid = mflds_base.grid();
    if (grid.bc.fld_lo[0] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[1] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[2] == BND_FLD_PERIODIC) {
      return;
    }

    auto& mflds = mflds_base.get_as<MfieldsCuda>(HX, HX + 3);
    if (grid.bc.fld_lo[0] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
	grid.bc.fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
	grid.bc.fld_lo[2] == BND_FLD_PERIODIC) {
      int d = 1;
      for (int p = 0; p < mflds.n_patches(); p++) {
	if (grid.atBoundaryLo(p, d)) {
	  cuda_conducting_wall_H_lo_y(mflds.cmflds(), p);
	}
	if (grid.atBoundaryHi(p, d)) {
	  cuda_conducting_wall_H_hi_y(mflds.cmflds(), p);

	}
      }
    } else {
      assert(0);
    }
    mflds_base.put_as(mflds, HX, HX + 3);
  }

  // ----------------------------------------------------------------------
  // add_ghosts_J

  void add_ghosts_J(MfieldsBase& mflds_base) override
  {
    const auto& grid = mflds_base.grid();
    if (grid.bc.fld_lo[0] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[1] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[2] == BND_FLD_PERIODIC) {
      return;
    }

    auto& mflds = mflds_base.get_as<MfieldsCuda>(JXI, JXI + 3);
    if (grid.bc.fld_lo[0] == BND_FLD_PERIODIC &&
	grid.bc.fld_lo[1] == BND_FLD_CONDUCTING_WALL &&
	grid.bc.fld_hi[1] == BND_FLD_CONDUCTING_WALL &&
	grid.bc.fld_lo[2] == BND_FLD_PERIODIC) {
      int d = 1;
      for (int p = 0; p < mflds.n_patches(); p++) {
	if (grid.atBoundaryLo(p, d)) {
	  cuda_conducting_wall_J_lo_y(mflds.cmflds(), p);
	}
	if (grid.atBoundaryHi(p, d)) {
	  cuda_conducting_wall_J_hi_y(mflds.cmflds(), p);
	}
      }
    } else {
      assert(0);
    }
    mflds_base.put_as(mflds, JXI, JXI + 3);
  }
};

