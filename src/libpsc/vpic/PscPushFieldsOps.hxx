
#pragma once

// ======================================================================
// PscPushFieldOps

template<typename MfieldsState, typename LocalOps, typename RemoteOps>
struct PscPushFieldsOps
{
  using Grid = typename MfieldsState::Grid;
  using SfaParams = typename MfieldsState::SfaParams;
  using MaterialCoefficient = typename MfieldsState::MaterialCoefficient;
  using F3D = Field3D<typename MfieldsState::Patch>;

  // FIXME, uses enum-based components vs struct-based components, should
  // settle on one or the other
  
  // ----------------------------------------------------------------------
  // advance_b
  
#define CBX MfieldsState::CBX
#define CBY MfieldsState::CBY
#define CBZ MfieldsState::CBZ
#define EX MfieldsState::EX
#define EY MfieldsState::EY
#define EZ MfieldsState::EZ

#define UPDATE_CBX() F(CBX, i,j,k) -= (py*(F(EZ, i,j+1,k) - F(EZ ,i,j,k)) - pz*(F(EY, i,j,k+1) - F(EY, i,j,k)))
#define UPDATE_CBY() F(CBY, i,j,k) -= (pz*(F(EX, i,j,k+1) - F(EX ,i,j,k)) - px*(F(EZ, i+1,j,k) - F(EZ, i,j,k)))
#define UPDATE_CBZ() F(CBZ, i,j,k) -= (px*(F(EY, i+1,j,k) - F(EY, i,j,k)) - py*(F(EX, i,j+1,k) - F(EX, i,j,k)))

  static void advance_b(MfieldsState& mflds, double frac)
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    int nx = g->nx, ny = g->ny, nz = g->nz;

    // FIXME, invariant should be based on global dims
    const float px   = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;
    const float py   = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;
    const float pz   = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0;

    // bulk
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  UPDATE_CBX(); UPDATE_CBY(); UPDATE_CBZ();
	}
      }
    }

    // leftover bx
    { int i = nx + 1;
      for (int k = 1; k <= nz; k++) {
	for (int j = 1; j <= ny; j++) {
	  UPDATE_CBX();
	}
      }
    }
    
    // leftover by
    { int j = ny + 1;
      for (int k = 1; k <= nz; k++) {
	for (int i = 1; i <= nx; i++) {
	  UPDATE_CBY();
	}
      }
    }
    
    // leftover bz
    { int k = nz + 1;
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  UPDATE_CBZ();
	}
      }
    }
    
    LocalOps::local_adjust_norm_b(mflds);
  }

#undef CBX
#undef CBY
#undef CBZ
#undef EX
#undef EY
#undef EZ

  // ----------------------------------------------------------------------
  // vacuum_advance_e
  
  static void vacuum_advance_e(MfieldsState& mflds, double frac)
  {
    // Update interior fields
    // Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
    // Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
    // Note: ez all (1:nx+1,1:ny+1,1:nz ) interior (1:nx,1:ny,2:nz)

    struct AdvanceE
    {
      AdvanceE(typename MfieldsState::Patch& fa, const Grid *g, const MaterialCoefficient* m,
	       const double damp_)
	: F(fa),
	  decayx(m->decayx),
	  decayy(m->decayy),
	  decayz(m->decayz),
	  drivex(m->drivex),
	  drivey(m->drivey),
	  drivez(m->drivez),
	  px_muz(g->nx > 1 ? (1+damp_)*g->cvac*g->dt*g->rdx*m->rmuz : 0),
	  px_muy(g->nx > 1 ? (1+damp_)*g->cvac*g->dt*g->rdx*m->rmuy : 0),
	  py_mux(g->ny > 1 ? (1+damp_)*g->cvac*g->dt*g->rdy*m->rmux : 0),
	  py_muz(g->ny > 1 ? (1+damp_)*g->cvac*g->dt*g->rdy*m->rmuz : 0),
	  pz_muy(g->nz > 1 ? (1+damp_)*g->cvac*g->dt*g->rdz*m->rmuy : 0),
	  pz_mux(g->nz > 1 ? (1+damp_)*g->cvac*g->dt*g->rdz*m->rmux : 0),
	  damp(damp_),
	  cj(g->dt / g->eps0)
      {
      }

      void x(int i, int j, int k)
      {
	F(i,j,k).tcax = ((py_muz*(F(i,j,k).cbz - F(i,j-1,k).cbz) -
			  pz_muy*(F(i,j,k).cby - F(i,j,k-1).cby)) -
			 damp*F(i,j,k).tcax);
	F(i,j,k).ex   = decayx*F(i,j,k).ex + drivex*(F(i,j,k).tcax - cj*F(i,j,k).jfx);
      }
      
      void y(int i, int j, int k)
      {
	F(i,j,k).tcay = ((pz_mux*(F(i,j,k).cbx - F(i,j,k-1).cbx) -
			  px_muz*(F(i,j,k).cbz - F(i-1,j,k).cbz)) -
			 damp*F(i,j,k).tcay);
	F(i,j,k).ey   = decayy*F(i,j,k).ey + drivey*(F(i,j,k).tcay - cj*F(i,j,k).jfy);
      }
      
      void z(int i, int j, int k)
      {
	F(i,j,k).tcaz = ((px_muy*(F(i,j,k).cby - F(i-1,j,k).cby) -
			  py_mux*(F(i,j,k).cbx - F(i,j-1,k).cbx)) -
			 damp*F(i,j,k).tcaz);
	F(i,j,k).ez   = decayz*F(i,j,k).ez + drivez*(F(i,j,k).tcaz - cj*F(i,j,k).jfz);
      }

      Field3D<typename MfieldsState::Patch> F;
      const float decayx, decayy, decayz, drivex, drivey, drivez;
      const float px_muz, px_muy, py_mux, py_muz, pz_muy, pz_mux;
      const float damp, cj;
    };

    auto& fa = mflds.getPatch(0);
    assert(frac == 1.);

    SfaParams& prm = mflds.params();
    assert(prm.size() == 1);
    const MaterialCoefficient* m = prm[0];

    AdvanceE advanceE(fa, mflds.vgrid(), m, prm.damp);

    RemoteOps::begin_remote_ghost_tang_b(mflds);

    LocalOps::local_ghost_tang_b(mflds);
    foreach_ec_interior(advanceE, mflds.vgrid());

    RemoteOps::end_remote_ghost_tang_b(mflds);

    foreach_ec_boundary(advanceE, mflds.vgrid());
    LocalOps::local_adjust_tang_e(mflds);
  }

  // ----------------------------------------------------------------------
  // advance_e
  
  static void advance_e(MfieldsState& mflds, double frac)
  {
    // FIXME vacuum hardcoded
    return vacuum_advance_e(mflds, frac);
  }
};

