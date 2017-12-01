
#ifndef PSC_FIELD_ARRAY_H
#define PSC_FIELD_ARRAY_H

#include <mrc_bits.h>

#include "PscFieldArrayRemoteOps.h" // FIXME, only because of Comm stuff

#include "vpic.h"

#define IN_sfa
#include "field_advance/standard/sfa_private.h"

// ======================================================================
// PscFieldArray

template<class B, class FieldArrayLocalOps, class FieldArrayRemoteOps>
struct PscFieldArray : B, FieldArrayLocalOps, FieldArrayRemoteOps
{
  typedef B Base;
  typedef PscFieldArray<B, FieldArrayLocalOps, FieldArrayRemoteOps> FieldArray;

  using Base::Base;

  using Base::g;
  using Base::params;
  
  // ----------------------------------------------------------------------
  // foreach
  // FIXME, move somewhere else

  template<class F>
  static void foreach(F f, int ib, int ie, int jb, int je, int kb, int ke)
  {
    for (int k = kb; k <= ke; k++) {
      for (int j = jb; j <= je; j++) {
	for (int i = ib; i <= ie; i++) {
	  f(i,j,k);
	}
      }
    }
  };

  template<class F>
  static void foreach_nc_interior(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    foreach(f, 2, nx, 2, ny, 2, nz);
  }

  template<class F>
  static void foreach_nc_boundary(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    // z faces, x edges, y edges and all corners
    foreach(f, 1   , nx+1, 1   , ny+1, 1   , 1   );
    foreach(f, 1   , nx+1, 1   , ny+1, nz+1, nz+1);

    // y faces, z edges
    foreach(f, 1   , nx+1, 1   , 1   ,    2, nz  );
    foreach(f, 1   , nx+1, ny+1, ny+1,    2, nz  );

    // x faces
    foreach(f, 1   , 1   , 2   , ny  ,    2, nz  );
    foreach(f, nx+1, nx+1, 2   , ny  ,    2, nz  );
  }

  template<class F>
  static void foreach_ec_interior(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	for (int i = 2; i <= nx; i++) {
	  f.x(i,j,k);
	  f.y(i,j,k);
	  f.z(i,j,k);
	}
      }
    }

    // Do leftover interior ex
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	f.x(1,j,k);
      }
    }

    // Do leftover interior ey
    for (int k = 2; k <= nz; k++) {
      for (int i = 2; i <= nx; i++) {
	f.y(i,1,k);
      }
    }

    // Do leftover interior ez
    for (int j = 2; j <= ny; j++) {
      for (int i = 2; i <= nx; i++) {
	f.z(i,j,1);
      }
    }
  }

  template<class F>
  static void foreach_ec_boundary(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    for (int j = 1; j <= ny+1; j++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,j,1);
      }
    }
    for (int j = 1; j <= ny+1; j++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,j,nz+1);
      }
    }
    for (int k = 2; k <= nz; k++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,1,k);
      }
    }
    for (int k = 2; k <= nz; k++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,ny+1,k);
      }
    }

    // Do exterior ey
    for (int k = 1; k <= nz+1; k++) {
      for (int j = 1; j <= ny; j++) {
	f.y(1,j,k);
      }
    }
    for (int k = 1; k <= nz+1; k++) {
      for (int j = 1; j <= ny; j++) {
	f.y(nx+1,j,k);
      }
    }
    for (int j = 1; j <= ny; j++) {
      for (int i = 2; i <= nx; i++) {
	f.y(i,j,1);
      }
    }
    for (int j = 1; j <= ny; j++) {
      for (int i = 2; i <= nx; i++) {
	f.y(i,j,nz+1);
      }
    }

    // Do exterior ez
    for (int k = 1; k <= nz; k++) {
      for (int i = 1; i <= nx+1; i++) {
	f.z(i,1,k);
      }
    }
    for (int k = 1; k <= nz; k++) {
      for (int i = 1; i <= nx+1; i++) {
	f.z(i,ny+1,k);
      }
    }
    for (int k = 1; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	f.z(1,j,k);
      }
    }
    for (int k = 1; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	f.z(nx+1,j,k);
      }
    }
  }
  
  // ----------------------------------------------------------------------
  // advance_b
  
#define CBX FieldArray::CBX
#define CBY FieldArray::CBY
#define CBZ FieldArray::CBZ
#define EX FieldArray::EX
#define EY FieldArray::EY
#define EZ FieldArray::EZ

#define UPDATE_CBX() F(CBX, i,j,k) -= (py*(F(EZ, i,j+1,k) - F(EZ ,i,j,k)) - pz*(F(EY, i,j,k+1) - F(EY, i,j,k)))
#define UPDATE_CBY() F(CBY, i,j,k) -= (pz*(F(EX, i,j,k+1) - F(EX ,i,j,k)) - px*(F(EZ, i+1,j,k) - F(EZ, i,j,k)))
#define UPDATE_CBZ() F(CBZ, i,j,k) -= (px*(F(EY, i+1,j,k) - F(EY, i,j,k)) - py*(F(EX, i,j+1,k) - F(EX, i,j,k)))

  void advance_b(double frac)
  {
    Field3D<FieldArray> F(*this);
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
    
    this->local_adjust_norm_b(*this);
  }

#undef CBX
#undef CBY
#undef CBZ
#undef EX
#undef EY
#undef EZ

  // ----------------------------------------------------------------------
  // advance_e
  
  void vacuum_advance_e(double frac)
  {
    sfa_params_t* prm = static_cast<sfa_params_t*>(params);
    assert(prm->n_mc == 1);
    assert(frac == 1.);

    // Update interior fields
    // Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
    // Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
    // Note: ez all (1:nx+1,1:ny+1,1:nz ) interior (1:nx,1:ny,2:nz)

    const material_coefficient_t* m = prm->mc;

    struct AdvanceE {
      AdvanceE(FieldArray& fa, const grid_t *g, const material_coefficient_t *m,
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

      Field3D<FieldArray> F;
      const float decayx, decayy, decayz, drivex, drivey, drivez;
      const float px_muz, px_muy, py_mux, py_muz, pz_muy, pz_mux;
      const float damp, cj;
    };

    AdvanceE advanceE(*this, g, m, prm->damp);

    this->begin_remote_ghost_tang_b(*this);

    this->local_ghost_tang_b(*this);
    foreach_ec_interior(advanceE, g);

    this->end_remote_ghost_tang_b(*this);

    foreach_ec_boundary(advanceE, g);
    this->local_adjust_tang_e(*this);
  }

  void advance_e(double frac)
  {
    vacuum_advance_e(frac);
  }

  // ----------------------------------------------------------------------
  // energy_f

  void vacuum_energy_f(double global[6])
  {
    sfa_params_t* prm = static_cast<sfa_params_t*>(params);
    assert(prm->n_mc == 1);
    const material_coefficient_t* m = prm->mc;

    Field3D<FieldArray> F(*this);
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
    mp_allsum_d(en, global, 6 );
  }

  void energy_f(double en[6])
  {
    vacuum_energy_f(en);
  }

  // ----------------------------------------------------------------------
  // clear_jf

  void clear_jf()
  {
    const int nv = g->nv;

    for (int v = 0; v < nv; v++) {
      (*this)[v].jfx = 0;
      (*this)[v].jfy = 0;
      (*this)[v].jfz = 0;
    }
  }

  // ----------------------------------------------------------------------
  // clear_rhof

  void clear_rhof()
  {
    const int nv = g->nv;
    for (int v = 0; v < nv; v++) {
      (*this)[v].rhof = 0;
    }
  }

  // ----------------------------------------------------------------------
  // synchronize_jf
  
  // ----------------------------------------------------------------------
  // CommJf

  template<class F3D>
  struct CommJf : Comm<F3D>
  {
    typedef Comm<F3D> Base;
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

  void synchronize_jf()
  {
    Field3D<FieldArray> F(*this);
    CommJf<Field3D<FieldArray>> comm(this->getGrid());

    this->local_adjust_jf(*this);

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, F);
      comm.end(dir, F);
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

  template<class F3D>
  struct CommRho : Comm<F3D>
  {
    typedef Comm<F3D> Base;
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

  void synchronize_rho()
  {
    Field3D<FieldArray> F(*this);
    CommRho<Field3D<FieldArray>> comm(this->getGrid());

    this->local_adjust_rhof(*this);
    this->local_adjust_rhob(*this);

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, F);
      comm.end(dir, F);
    }
  }

  // ----------------------------------------------------------------------
  // compute_rhob

  void compute_rhob()
  {
    sfa_params_t* prm = static_cast<sfa_params_t*>(params);
    assert(prm->n_mc == 1);
    const material_coefficient_t* m = prm->mc;

    struct CalcRhoB {
      CalcRhoB(FieldArray& fa, const material_coefficient_t *m)
	: F(fa),
	  nc(m->nonconductive),
	  px(fa.g->nx > 1 ? fa.g->eps0 * m->epsx * fa.g->rdx : 0),
	  py(fa.g->ny > 1 ? fa.g->eps0 * m->epsy * fa.g->rdy : 0),
	  pz(fa.g->nz > 1 ? fa.g->eps0 * m->epsz * fa.g->rdz : 0)
      {
      }

      void operator()(int i, int j, int k)
      {
	F(i,j,k).rhob = nc*(px * (F(i,j,k).ex - F(i-1,j,k).ex) +
			    py * (F(i,j,k).ey - F(i,j-1,k).ey) +
			    pz * (F(i,j,k).ez - F(i,j,k-1).ez) -
			    F(i,j,k).rhof);
      }

      Field3D<FieldArray> F;
      const float nc, px, py, pz;
    };

    CalcRhoB updater(*this, m);

    // Begin setting normal e ghosts
    this->begin_remote_ghost_norm_e(*this);

    // Overlap local computation
    this->local_ghost_norm_e(*this);
    foreach_nc_interior(updater, g);
    
    // Finish setting normal e ghosts
    this->end_remote_ghost_norm_e(*this);

    // Now do points on boundary
    foreach_nc_boundary(updater, g);

    this->local_adjust_rhob(*this);
  }

  // ----------------------------------------------------------------------
  // compute_curl_b
  
  void vacuum_compute_curl_b()
  {
    sfa_params_t* prm = static_cast<sfa_params_t*>(params);
    assert(prm->n_mc == 1);
    const material_coefficient_t* m = prm->mc;

    // Update interior fields
    // Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
    // Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
    // Note: ez all (1:nx+1,1:ny+1,1:nz ) interior (1:nx,1:ny,2:nz)

    struct CurlB {
      CurlB(FieldArray& fa, const grid_t *g, const material_coefficient_t *m)
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

      Field3D<FieldArray> F;
      const float px_muz, px_muy, py_mux, py_muz, pz_muy, pz_mux;
    };

    CurlB curlB(*this, g, m);
      
    this->begin_remote_ghost_tang_b(*this);

    this->local_ghost_tang_b(*this);
    foreach_ec_interior(curlB, g);

    this->end_remote_ghost_tang_b(*this);

    foreach_ec_boundary(curlB, g);
    this->local_adjust_tang_e(*this); // FIXME, is this right here?
  }

  void compute_curl_b()
  {
    vacuum_compute_curl_b();
  }

  // ----------------------------------------------------------------------
  // synchronize_tang_e_norm_b

  template<class F3D>
  struct CommTangENormB : Comm<F3D>
  {
    typedef Comm<F3D> Base;
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

  double synchronize_tang_e_norm_b()
  {
    Field3D<FieldArray> F(*this);
    CommTangENormB<Field3D<FieldArray>> comm(this->getGrid());
    
    this->local_adjust_tang_e(*this);
    this->local_adjust_norm_b(*this);

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, F);
      comm.end(dir, F);
    }
    
    double gerr;
    mp_allsum_d(&comm.err, &gerr, 1);
    return gerr;
  }

  // ----------------------------------------------------------------------
  // compute_div_e_err
  
  void compute_div_e_err()
  {
    sfa_params_t* prm = static_cast<sfa_params_t*>(params);
    assert(prm->n_mc == 1);
    const material_coefficient_t* m = prm->mc;

    struct CalcDivE {
      CalcDivE(FieldArray& fa, const material_coefficient_t *m)
	: F(fa),
	  nc(m->nonconductive),
	  px(fa.g->nx > 1 ? fa.g->eps0 * fa.g->rdx : 0),
	  py(fa.g->ny > 1 ? fa.g->eps0 * fa.g->rdy : 0),
	  pz(fa.g->nz > 1 ? fa.g->eps0 * fa.g->rdz : 0),
	  cj(1. / fa.g->eps0)
      {
      }

      void operator()(int i, int j, int k)
      {
	F(i,j,k).div_e_err = nc*(px * (F(i,j,k).ex - F(i-1,j,k).ex) +
				 py * (F(i,j,k).ey - F(i,j-1,k).ey) +
				 pz * (F(i,j,k).ez - F(i,j,k-1).ez) -
				 cj * (F(i,j,k).rhof + F(i,j,k).rhob));
      }

      Field3D<FieldArray> F;
      const float nc, px, py, pz, cj;
    };

    CalcDivE updater(*this, m);
    
    // Begin setting normal e ghosts
    this->begin_remote_ghost_norm_e(*this);

    // Overlap local computation
    this->local_ghost_norm_e(*this);
    foreach_nc_interior(updater, g);

    // Finish setting normal e ghosts
    this->end_remote_ghost_norm_e(*this);

    // Now do points on boundary
    foreach_nc_boundary(updater, g);

    this->local_adjust_div_e(*this);
  }

  // ----------------------------------------------------------------------
  // compute_rms_div_e_err
  //
  // OPT: doing this at the same time as compute_div_e_err might be
  // advantageous (faster)
  
  double compute_rms_div_e_err()
  {
    Field3D<FieldArray> F(*this);
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
    mp_allsum_d(local, _global, 2);
    return g->eps0 * sqrt(_global[0]/_global[1]);
  }

  // ----------------------------------------------------------------------
  // clean_div_e
  
#define MARDER_EX(i,j,k) F(i,j,k).ex += px * (F(i+1,j,k).div_e_err - F(i,j,k).div_e_err)
#define MARDER_EY(i,j,k) F(i,j,k).ey += py * (F(i,j+1,k).div_e_err - F(i,j,k).div_e_err)
#define MARDER_EZ(i,j,k) F(i,j,k).ez += pz * (F(i,j,k+1).div_e_err - F(i,j,k).div_e_err)

  void vacuum_clean_div_e()
  {
    sfa_params_t* prm = static_cast<sfa_params_t*>(params);
    assert(prm->n_mc == 1);

    Field3D<FieldArray> F(*this);
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    const float _rdx = (nx>1) ? g->rdx : 0;
    const float _rdy = (ny>1) ? g->rdy : 0;
    const float _rdz = (nz>1) ? g->rdz : 0;
    const float alphadt = 0.3888889/(_rdx*_rdx + _rdy*_rdy + _rdz*_rdz);
    const float px = (alphadt*_rdx) * prm->mc[0].drivex;
    const float py = (alphadt*_rdy) * prm->mc[0].drivey;
    const float pz = (alphadt*_rdz) * prm->mc[0].drivez;
                     
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

    this->local_adjust_tang_e(*this);
  }

  void clean_div_e()
  {
    vacuum_clean_div_e();
  }

  // ----------------------------------------------------------------------
  // compute_div_b_err
  
  void compute_div_b_err()
  {
    Field3D<FieldArray> F(*this);
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

  double compute_rms_div_b_err()
  {
    Field3D<FieldArray> F(*this);
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
    mp_allsum_d(local, _global, 2);
    return g->eps0 * sqrt(_global[0]/_global[1]);
  }

  // ----------------------------------------------------------------------
  // clean_div_b
  
#define MARDER_CBX(i,j,k) F(i,j,k).cbx += px * (F(i,j,k).div_b_err - F(i-1,j,k).div_b_err)
#define MARDER_CBY(i,j,k) F(i,j,k).cby += py * (F(i,j,k).div_b_err - F(i,j-1,k).div_b_err)
#define MARDER_CBZ(i,j,k) F(i,j,k).cbz += pz * (F(i,j,k).div_b_err - F(i,j,k-1).div_b_err)

  void clean_div_b()
  {
    Field3D<FieldArray> F(*this);

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
    this->begin_remote_ghost_div_b(*this);
    this->local_ghost_div_b(*this);

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
  
    this->end_remote_ghost_div_b(*this);

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

    this->local_adjust_norm_b(*this);
  }

};
  

#endif
