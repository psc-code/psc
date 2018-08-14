
#ifndef PSC_FIELD_ARRAY_H
#define PSC_FIELD_ARRAY_H

#include <mrc_bits.h>

// ======================================================================
// PscCleanDivOps

template<typename MfieldsState, typename LocalOps, typename RemoteOps>
struct PscCleanDivOps
{
  using Grid = typename MfieldsState::Grid;
  using MaterialCoefficient = typename MfieldsState::MaterialCoefficient;
  using F3D = Field3D<typename MfieldsState::Patch>;
  
  // ----------------------------------------------------------------------
  // clear_rhof

  static void clear_rhof(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    const int nv = mflds.vgrid()->nv;

    for (int v = 0; v < nv; v++) {
      fa[v].rhof = 0;
    }
  }

  // ----------------------------------------------------------------------
  // synchronize_rho
  
  // Note: synchronize_rho assumes that rhof has _not_ been adjusted at
  // the local domain boundary to account for partial cells but that
  // rhob _has_.  Specifically it is very expensive to accumulate rhof
  // and doing the adjustment for each particle is adds even more
  // expense.  Worse, if we locally corrected it after each species,
  // we cannot accumulate the next species in the same unless we use
  // (running sum of locally corrected results and thw current species
  // rhof being accumulated).  Further, rhof is always accumulated from
  // scratch so we don't have to worry about whether or not the previous
  // rhof values were in a locally corrected form.  Thus, after all
  // particles have accumulated to rhof, we correct it for partial cells
  // and remote cells for use with divergence cleaning and this is
  // the function that does the correction.
  //
  // rhob is another story though.  rhob is continuously incrementally
  // accumulated over time typically through infrequent surface area
  // scaling processes.  Like rho_f, after synchronize_rhob, rhob _must_
  // be corrected for partial and remote celle for the benefit of
  // divergence cleaning. And like rho_f, since we don't want to have
  // to keep around two versions of rhob (rhob contributions since last
  // sync and rhob as of last sync), we have no choice but to do the
  // charge accumulation per particle to rhob in a locally corrected
  // form.

  // ----------------------------------------------------------------------
  // CommRho

  template<class G, class F3D>
  struct CommRho : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;
    
    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommRho(Grid *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = 2 * (nx_[Y] + 1) * (nx_[Z] + 1);
      }
    }

    void begin_send(int X, int side, float* p, F3D& F)
    {
      int face = side ? nx_[X] + 1 : 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) { *p++ = F(x,y,z).rhof; });
      foreach_node(g_, X, face, [&](int x, int y, int z) { *p++ = F(x,y,z).rhob; });
    }
    
    void end_recv(int X, int side, float* p, F3D& F)
    {
      int face = side ? 1 : nx_[X] + 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) { F(x,y,z).rhof += *p++; });
      foreach_node(g_, X, face, [&](int x, int y, int z) { F(x,y,z).rhob = .5f*(F(x,y,z).rhob + *p++); });
    }
  };

  static void synchronize_rho(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommRho<Grid, F3D> comm(mflds.vgrid());

    LocalOps::local_adjust_rhof(mflds);
    LocalOps::local_adjust_rhob(mflds);

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, F);
      comm.end(dir, F);
    }
  }

  // ----------------------------------------------------------------------
  // compute_div_e_err
  
  static void compute_div_e_err(MfieldsState& mflds)
  {
    struct CalcDivE {
      CalcDivE(typename MfieldsState::Patch& fa, const MaterialCoefficient* m)
	: F(fa),
	  nc(m->nonconductive),
	  px(fa.grid()->nx > 1 ? fa.grid()->eps0 * fa.grid()->rdx : 0),
	  py(fa.grid()->ny > 1 ? fa.grid()->eps0 * fa.grid()->rdy : 0),
	  pz(fa.grid()->nz > 1 ? fa.grid()->eps0 * fa.grid()->rdz : 0),
	  cj(1. / fa.grid()->eps0)
      {
      }

      void operator()(int i, int j, int k)
      {
	F(i,j,k).div_e_err = nc*(px * (F(i,j,k).ex - F(i-1,j,k).ex) +
				 py * (F(i,j,k).ey - F(i,j-1,k).ey) +
				 pz * (F(i,j,k).ez - F(i,j,k-1).ez) -
				 cj * (F(i,j,k).rhof + F(i,j,k).rhob));
      }

      F3D F;
      const float nc, px, py, pz, cj;
    };

    auto& fa = mflds.getPatch(0);
    auto& prm = mflds.params();
    assert(prm.size() == 1);
    const MaterialCoefficient* m = prm[0];

    CalcDivE updater(fa, m);
    
    // Begin setting normal e ghosts
    RemoteOps::begin_remote_ghost_norm_e(mflds);

    // Overlap local computation
    LocalOps::local_ghost_norm_e(mflds);
    foreach_nc_interior(updater, mflds.vgrid());

    // Finish setting normal e ghosts
    RemoteOps::end_remote_ghost_norm_e(mflds);

    // Now do points on boundary
    foreach_nc_boundary(updater, mflds.vgrid());

    LocalOps::local_adjust_div_e(mflds);
  }

  // ----------------------------------------------------------------------
  // compute_rms_div_e_err
  //
  // OPT: doing this at the same time as compute_div_e_err might be
  // advantageous (faster)
  
  static double compute_rms_div_e_err(MfieldsState& mflds)
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    // Interior points
    // FIXME, it's inconsistent to calc the square in single prec here, but
    // double prec later
    double err = 0;
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	for (int i = 2; i <= nx; i++) {
	  err += sqr(F(i,j,k).div_e_err);
	}
      }
    }

    int x, y, z;
    // Exterior faces

    for( y=2; y<=ny; y++ ) {
      for( z=2; z<=nz; z++ ) {
	err += 0.5*sqr((double) (F(1   ,y,z).div_e_err));
	err += 0.5*sqr((double) (F(nx+1,y,z).div_e_err));
      }
    }

    for( z=2; z<=nz; z++ ) {
      for( x=2; x<=nx; x++ ) {
	err += 0.5*sqr((double) (F(x,1   ,z).div_e_err));
	err += 0.5*sqr((double) (F(x,ny+1,z).div_e_err));
      }
    }
    
    for( x=2; x<=nx; x++ ) {
      for( y=2; y<=ny; y++ ) {
	err += 0.5*sqr((double) (F(x,y,1   ).div_e_err));
	err += 0.5*sqr((double) (F(x,y,nz+1).div_e_err));
      }
    }

    // Exterior edges

    for( x=2; x<=nx; x++ ) {
      err += 0.25*sqr((double) (F(x,1   ,1   ).div_e_err));
      err += 0.25*sqr((double) (F(x,ny+1,1   ).div_e_err));
      err += 0.25*sqr((double) (F(x,1   ,nz+1).div_e_err));
      err += 0.25*sqr((double) (F(x,ny+1,nz+1).div_e_err));
    }

    for( y=2; y<=ny; y++ ) {
      err += 0.25*sqr((double) (F(1   ,y,1   ).div_e_err));
      err += 0.25*sqr((double) (F(1   ,y,nz+1).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,y,1   ).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,y,nz+1).div_e_err));
    }

    for( z=2; z<=nz; z++ ) {
      err += 0.25*sqr((double) (F(1   ,1   ,z).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,1   ,z).div_e_err));
      err += 0.25*sqr((double) (F(1   ,ny+1,z).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,ny+1,z).div_e_err));
    }

    // Exterior corners

    err += 0.125*sqr((double) (F(1   ,1   ,   1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,1   ,   1).div_e_err));
    err += 0.125*sqr((double) (F(1   ,ny+1,   1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,ny+1,   1).div_e_err));
    err += 0.125*sqr((double) (F(1   ,1   ,nz+1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,1   ,nz+1).div_e_err));
    err += 0.125*sqr((double) (F(1   ,ny+1,nz+1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,ny+1,nz+1).div_e_err));
    
    double local[2], _global[2]; // FIXME, name clash with global macro
    local[0] = err * g->dV;
    local[1] = (nx * ny * nz) * g->dV;
    MPI_Allreduce(local, _global, 2, MPI_DOUBLE, MPI_SUM, psc_comm_world);
    return g->eps0 * sqrt(_global[0]/_global[1]);
  }

  // ----------------------------------------------------------------------
  // clean_div_e
  
#define MARDER_EX(i,j,k) F(i,j,k).ex += px * (F(i+1,j,k).div_e_err - F(i,j,k).div_e_err)
#define MARDER_EY(i,j,k) F(i,j,k).ey += py * (F(i,j+1,k).div_e_err - F(i,j,k).div_e_err)
#define MARDER_EZ(i,j,k) F(i,j,k).ez += pz * (F(i,j,k+1).div_e_err - F(i,j,k).div_e_err)

  static void vacuum_clean_div_e(MfieldsState& mflds)
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    auto& prm = mflds.params();
    assert(prm.size() == 1);
    const MaterialCoefficient* m = prm[0];

    const int nx = g->nx, ny = g->ny, nz = g->nz;

    const float _rdx = (nx>1) ? g->rdx : 0;
    const float _rdy = (ny>1) ? g->rdy : 0;
    const float _rdz = (nz>1) ? g->rdz : 0;
    const float alphadt = 0.3888889/(_rdx*_rdx + _rdy*_rdy + _rdz*_rdz);
    const float px = (alphadt*_rdx) * m->drivex;
    const float py = (alphadt*_rdy) * m->drivey;
    const float pz = (alphadt*_rdz) * m->drivez;
                     
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  MARDER_EX(i,j,k);
	  MARDER_EY(i,j,k);
	  MARDER_EZ(i,j,k);
	}
      }
    }
    
    int x, y, z;
  
    // Do left over ex
    for(y=1; y<=ny+1; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_EX(x,y,nz+1);
      }
    }
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_EX(x,ny+1,z);
      }
    }

    // Do left over ey
    for(z=1; z<=nz+1; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_EY(nx+1,y,z);
      }
    }
    for(y=1; y<=ny; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_EY(x,y,nz+1);
      }
    }

    // Do left over ez
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx+1; x++) {
	MARDER_EZ(x,ny+1,z);
      }
    }
    for(z=1; z<=nz; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_EZ(nx+1,y,z);
      }
    }

    LocalOps::local_adjust_tang_e(mflds);
  }

  static void clean_div_e(MfieldsState& mflds)
  {
    vacuum_clean_div_e(mflds);
  }

  // ----------------------------------------------------------------------
  // compute_div_b_err
  
  static void compute_div_b_err(MfieldsState& mflds)
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    const int nx = g->nx, ny = g->ny, nz = g->nz;
    const float px = (nx>1) ? g->rdx : 0;  // FIXME, should be based on global dims
    const float py = (ny>1) ? g->rdy : 0;
    const float pz = (nz>1) ? g->rdz : 0;

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  F(i,j,k).div_b_err = (px * (F(i+1,j,k).cbx - F(i,j,k).cbx) +
				py * (F(i,j+1,k).cby - F(i,j,k).cby) +
				pz * (F(i,j,k+1).cbz - F(i,j,k).cbz));
	}
      }
    }
  }

  // ----------------------------------------------------------------------
  // compute_rms_div_b_err
  //
  // OPT: doing that at the same time as div_b should be faster

  static double compute_rms_div_b_err(MfieldsState& mflds)
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    const int nx = g->nx, ny = g->ny, nz = g->nz;

    double err = 0;
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  err += sqr(F(i,j,k).div_b_err);
	}
      }
    }
    
    double local[2], _global[2]; // FIXME, name clash with global macro
    local[0] = err * g->dV;
    local[1] = (nx * ny * nz) * g->dV;
    MPI_Allreduce(local, _global, 2, MPI_DOUBLE, MPI_SUM, psc_comm_world);
    return g->eps0 * sqrt(_global[0]/_global[1]);
  }

  // ----------------------------------------------------------------------
  // clean_div_b
  
#define MARDER_CBX(i,j,k) F(i,j,k).cbx += px * (F(i,j,k).div_b_err - F(i-1,j,k).div_b_err)
#define MARDER_CBY(i,j,k) F(i,j,k).cby += py * (F(i,j,k).div_b_err - F(i,j-1,k).div_b_err)
#define MARDER_CBZ(i,j,k) F(i,j,k).cbz += pz * (F(i,j,k).div_b_err - F(i,j,k-1).div_b_err)

  static void clean_div_b(MfieldsState& mflds)
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    const int nx = g->nx, ny = g->ny, nz = g->nz;
    float px = (nx>1) ? g->rdx : 0;
    float py = (ny>1) ? g->rdy : 0;
    float pz = (nz>1) ? g->rdz : 0;
    float alphadt = 0.3888889/(px*px + py*py + pz*pz);
    px *= alphadt;
    py *= alphadt;
    pz *= alphadt;

    // Have pipelines do Marder pass in interior.  The host handles
    // stragglers.

    // Begin setting ghosts
    RemoteOps::begin_remote_ghost_div_b(mflds);
    LocalOps::local_ghost_div_b(mflds);

    // Interior
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	for (int i = 2; i <= nx; i++) {
	  MARDER_CBX(i,j,k);
	  MARDER_CBY(i,j,k);
	  MARDER_CBZ(i,j,k);
	}
      }
    }
  
    int x, y, z;

    // Do left over interior bx
    for(y=1; y<=ny; y++) {
      for(x=2; x<=nx; x++) {
	MARDER_CBX(x,y,1);
      }
    }
    for(z=2; z<=nz; z++) {
      for(x=2; x<=nx; x++) {
	MARDER_CBX(x,1,z);
      }
    }

    // Left over interior by
    for(z=1; z<=nz; z++) {
      for(y=2; y<=ny; y++) {
	MARDER_CBY(1,y,z);
      }
    }
    for(y=2; y<=ny; y++) {
      for(x=2; x<=nx; x++) {
	MARDER_CBY(x,y,1);
      }
    }

    // Left over interior bz
    for(z=2; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBZ(x,1,z);
      }
    }
    for(z=2; z<=nz; z++) {
      for(y=2; y<=ny; y++) {
	MARDER_CBZ(1,y,z);
      }
    }

    // Finish setting derr ghosts
  
    RemoteOps::end_remote_ghost_div_b(mflds);

    // Do Marder pass in exterior

    // Exterior bx
    for(z=1; z<=nz; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_CBX(1,y,z);
      }
    }
    for(z=1; z<=nz; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_CBX(nx+1,y,z);
      }
    }

    // Exterior by
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBY(x,1,z);
      }
    }
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBY(x,ny+1,z);
      }
    }

    // Exterior bz
    for(y=1; y<=ny; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBZ(x,y,1);
      }
    }
    for(y=1; y<=ny; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBZ(x,y,nz+1);
      }
    }

    LocalOps::local_adjust_norm_b(mflds);
  }

  // ----------------------------------------------------------------------
  // synchronize_tang_e_norm_b

  template<class G, class F3D>
  struct CommTangENormB : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;
    
    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommTangENormB(Grid *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = (nx_[Y] * nx_[Z] +
			2 * nx_[Y] * (nx_[Z] + 1) +
			2 * nx_[Z] * (nx_[Y] + 1));
      }
      err = 0.;
    }

    void begin_send(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? nx_[X] + 1 : 1;
      foreach_face(g_, X, face, [&](int x, int y, int z)
		   { *p++ = (&F(x,y,z).cbx)[X]; });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z)
		   { *p++ = (&F(x,y,z).ex)[Y]; *p++ = (&F(x,y,z).tcax)[Y]; });
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z)
		   { *p++ = (&F(x,y,z).ex)[Z]; *p++ = (&F(x,y,z).tcax)[Z]; });
    }
    
    void end_recv(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? 1 : nx_[X] + 1;
      foreach_face(g_, X, face, [&](int x, int y, int z)
		   { double w1 = *p++, w2 = (&F(x,y,z).cbx)[X];
		     (&F(x,y,z).cbx)[X] = .5 * (w1+w2);
		     err += sqr(w1-w2); });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z)
		   { double w1 = *p++, w2 = (&F(x,y,z).ex)[Y];
		     (&F(x,y,z).ex)[Y] = .5 * (w1+w2);
		     err += sqr(w1-w2);
		     w1 = *p++, w2 = (&F(x,y,z).tcax)[Y];
		     (&F(x,y,z).tcax)[Y] = .5 * (w1+w2); });
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z)
		   { double w1 = *p++, w2 = (&F(x,y,z).ex)[Z];
		     (&F(x,y,z).ex)[Z] = .5 * (w1+w2);
		     err += sqr(w1-w2);
		     w1 = *p++, w2 = (&F(x,y,z).tcax)[Z];
		     (&F(x,y,z).tcax)[Z] = .5 * (w1+w2); });
    }

    double err;
  };

  static double synchronize_tang_e_norm_b(MfieldsState& mflds)
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommTangENormB<Grid, F3D> comm(mflds.vgrid());
    
    LocalOps::local_adjust_tang_e(mflds);
    LocalOps::local_adjust_norm_b(mflds);

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, F);
      comm.end(dir, F);
    }
    
    double gerr;
    MPI_Allreduce(&comm.err, &gerr, 1, MPI_DOUBLE, MPI_SUM, psc_comm_world);
    return gerr;
  }

};

// ======================================================================
// PscAccumulateOps

template<typename MfieldsState, typename LocalOps, typename RemoteOps>
struct PscAccumulateOps
{
  using Grid = typename MfieldsState::Grid;
  using MaterialCoefficient = typename MfieldsState::MaterialCoefficient;
  using F3D = Field3D<typename MfieldsState::Patch>;
  
  // ----------------------------------------------------------------------
  // clear_jf

  static void clear_jf(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    const int nv = mflds.vgrid()->nv;

    for (int v = 0; v < nv; v++) {
      fa[v].jfx = 0;
      fa[v].jfy = 0;
      fa[v].jfz = 0;
    }
  }

  // ----------------------------------------------------------------------
  // CommJf

  template<class G, class F3D>
  struct CommJf : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;

    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommJf(Grid *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = nx_[Y] * (nx_[Z] + 1) + nx_[Z] * (nx_[Y] + 1);
      }
    }

    void begin_send(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? nx_[X] + 1 : 1;
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).jfx)[Y]; });
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).jfx)[Z]; });
    }
    
    void end_recv(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? 1 : nx_[X] + 1;
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { (&F(x,y,z).jfx)[Y] += *p++; });
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { (&F(x,y,z).jfx)[Z] += *p++; });
    }
  };

  // ----------------------------------------------------------------------
  // synchronize_jf
  
  static void synchronize_jf(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    Field3D<typename MfieldsState::Patch> F(fa);
    CommJf<Grid, Field3D<typename MfieldsState::Patch>> comm(mflds.vgrid());

    LocalOps::local_adjust_jf(mflds);

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, F);
      comm.end(dir, F);
    }
  }

  // ----------------------------------------------------------------------
  // compute_rhob

  static void compute_rhob(MfieldsState& mflds)
  {
    struct CalcRhoB {
      CalcRhoB(typename MfieldsState::Patch& fa, const MaterialCoefficient* m)
	: F(fa),
	  nc(m->nonconductive),
	  px(fa.grid()->nx > 1 ? fa.grid()->eps0 * m->epsx * fa.grid()->rdx : 0),
	  py(fa.grid()->ny > 1 ? fa.grid()->eps0 * m->epsy * fa.grid()->rdy : 0),
	  pz(fa.grid()->nz > 1 ? fa.grid()->eps0 * m->epsz * fa.grid()->rdz : 0)
      {
      }

      void operator()(int i, int j, int k)
      {
	F(i,j,k).rhob = nc*(px * (F(i,j,k).ex - F(i-1,j,k).ex) +
			    py * (F(i,j,k).ey - F(i,j-1,k).ey) +
			    pz * (F(i,j,k).ez - F(i,j,k-1).ez) -
			    F(i,j,k).rhof);
      }

      F3D F;
      const float nc, px, py, pz;
    };

    auto& fa = mflds.getPatch(0);
    auto& prm = mflds.params();
    assert(prm.size() == 1);
    const MaterialCoefficient* m = prm[0];

    CalcRhoB updater(fa, m);

    // Begin setting normal e ghosts
    RemoteOps::begin_remote_ghost_norm_e(mflds);

    // Overlap local computation
    LocalOps::local_ghost_norm_e(mflds);
    foreach_nc_interior(updater, mflds.vgrid());
    
    // Finish setting normal e ghosts
    RemoteOps::end_remote_ghost_norm_e(mflds);

    // Now do points on boundary
    foreach_nc_boundary(updater, mflds.vgrid());

    LocalOps::local_adjust_rhob(mflds);
  }

  // ----------------------------------------------------------------------
  // compute_curl_b
  
  static void vacuum_compute_curl_b(MfieldsState& mflds)
  {
    // Update interior fields
    // Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
    // Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
    // Note: ez all (1:nx+1,1:ny+1,1:nz ) interior (1:nx,1:ny,2:nz)

    struct CurlB {
      CurlB(typename MfieldsState::Patch& fa, const Grid *g, const MaterialCoefficient* m)
	: F(fa),
	  px_muz(g->nx > 1 ? g->cvac*g->dt*g->rdx*m->rmuz : 0),
	  px_muy(g->nx > 1 ? g->cvac*g->dt*g->rdx*m->rmuy : 0),
	  py_mux(g->ny > 1 ? g->cvac*g->dt*g->rdy*m->rmux : 0),
	  py_muz(g->ny > 1 ? g->cvac*g->dt*g->rdy*m->rmuz : 0),
	  pz_muy(g->nz > 1 ? g->cvac*g->dt*g->rdz*m->rmuy : 0),
	  pz_mux(g->nz > 1 ? g->cvac*g->dt*g->rdz*m->rmux : 0)
      {
      }

      void x(int i, int j, int k)
      {
	F(i,j,k).tcax = (py_muz*(F(i,j,k).cbz - F(i,j-1,k).cbz) -
			 pz_muy*(F(i,j,k).cby - F(i,j,k-1).cby));
      }
      
      void y(int i, int j, int k)
      {
	F(i,j,k).tcay = (pz_mux*(F(i,j,k).cbx - F(i,j,k-1).cbx) -
			 px_muz*(F(i,j,k).cbz - F(i-1,j,k).cbz));
      }
      
      void z(int i, int j, int k)
      {
	F(i,j,k).tcaz = (px_muy*(F(i,j,k).cby - F(i-1,j,k).cby) -
			 py_mux*(F(i,j,k).cbx - F(i,j-1,k).cbx));
      }

      F3D F;
      const float px_muz, px_muy, py_mux, py_muz, pz_muy, pz_mux;
    };

    auto& fa = mflds.getPatch(0);
    auto& prm = mflds.params();
    assert(prm.size() == 1);
    const MaterialCoefficient* m = prm[0];

    CurlB curlB(fa, mflds.vgrid(), m);
      
    RemoteOps::begin_remote_ghost_tang_b(mflds);

    LocalOps::local_ghost_tang_b(mflds);
    foreach_ec_interior(curlB, mflds.vgrid());

    RemoteOps::end_remote_ghost_tang_b(mflds);

    foreach_ec_boundary(curlB, mflds.vgrid());
    LocalOps::local_adjust_tang_e(mflds); // FIXME, is this right here?
  }

  static void compute_curl_b(MfieldsState& mflds)
  {
    vacuum_compute_curl_b(mflds);
  }
};

// ======================================================================
// PscDiagOps

template<typename MfieldsState>
struct PscDiagOps
{
  using Grid = typename MfieldsState::Grid;
  using F3D = Field3D<typename MfieldsState::Patch>;
  
  // ----------------------------------------------------------------------
  // vacuum_energy_f

  static void vacuum_energy_f(MfieldsState& mflds, double global[6])
  {
    const Grid* g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    auto& prm = mflds.params();
    assert(prm.size() == 1);
    auto m = prm[0];

    const int nx = g->nx, ny = g->ny, nz = g->nz;

    const float qepsx = 0.25*m->epsx;
    const float qepsy = 0.25*m->epsy;
    const float qepsz = 0.25*m->epsz;
    const float hrmux = 0.5*m->rmux;
    const float hrmuy = 0.5*m->rmuy;
    const float hrmuz = 0.5*m->rmuz;
    double en[6] = {};

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  en[0] += qepsx*(sqr(F(i  ,j  ,k  ).ex) +
			  sqr(F(i  ,j+1,k  ).ex) +
			  sqr(F(i  ,j  ,k+1).ex) +
			  sqr(F(i  ,j+1,k+1).ex));
	  en[1] += qepsy*(sqr(F(i  ,j  ,k  ).ey) +
			  sqr(F(i  ,j  ,k+1).ey) +
			  sqr(F(i+1,j  ,k  ).ey) +
			  sqr(F(i+1,j  ,k+1).ey));
	  en[2] += qepsz*(sqr(F(i  ,j  ,k  ).ez) +
			  sqr(F(i+1,j  ,k  ).ez) +
			  sqr(F(i  ,j+1,k  ).ez) +
			  sqr(F(i+1,j+1,k  ).ez));
	  en[3] += hrmux*(sqr(F(i  ,j  ,k  ).cbx) +
			  sqr(F(i+1,j  ,k  ).cbx));
	  en[4] += hrmuy*(sqr(F(i  ,j  ,k  ).cby) +
			  sqr(F(i  ,j+1,k  ).cby));
	  en[5] += hrmuz*(sqr(F(i  ,j  ,k  ).cbz) +
			  sqr(F(i  ,j  ,k+1).cbz));
	}
      }
    }

    // Convert to physical units
    double v0 = 0.5 * g->eps0 * g->dV;
    for (int m = 0; m < 6; m++) {
      en[m] *= v0;
    }

    // Reduce results between nodes
    MPI_Allreduce(en, global, 6, MPI_DOUBLE, MPI_SUM, psc_comm_world);
  }

  // ----------------------------------------------------------------------
  // energy_f

  static void energy_f(MfieldsState& mflds, double en[6])
  {
    vacuum_energy_f(mflds, en);
  }
  
};

#include "PscFieldArrayRemoteOps.h" // FIXME, only because of Comm stuff

#endif
