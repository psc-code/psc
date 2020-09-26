
#ifndef PSC_FIELD_ARRAY_H
#define PSC_FIELD_ARRAY_H

#include <mrc_bits.h>

// ======================================================================
// PscAccumulateOps

template <typename MfieldsState, typename LocalOps, typename RemoteOps>
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
    const int nv = mflds.vgrid().nv;

    for (int v = 0; v < nv; v++) {
      fa[v].jfx = 0;
      fa[v].jfy = 0;
      fa[v].jfz = 0;
    }
  }

  // ----------------------------------------------------------------------
  // CommJf

  template <class G, class F3D>
  struct CommJf : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;

    using Base::begin;
    using Base::end;

    using Base::buf_size_;
    using Base::g_;
    using Base::nx_;

    CommJf(Grid& g) : Base(g)
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
      foreach_edge(g_, Y, Z, face,
                   [&](int x, int y, int z) { *p++ = (&F(x, y, z).jfx)[Y]; });
      foreach_edge(g_, Z, Y, face,
                   [&](int x, int y, int z) { *p++ = (&F(x, y, z).jfx)[Z]; });
    }

    void end_recv(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? 1 : nx_[X] + 1;
      foreach_edge(g_, Y, Z, face,
                   [&](int x, int y, int z) { (&F(x, y, z).jfx)[Y] += *p++; });
      foreach_edge(g_, Z, Y, face,
                   [&](int x, int y, int z) { (&F(x, y, z).jfx)[Z] += *p++; });
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
    struct CalcRhoB
    {
      CalcRhoB(typename MfieldsState::Patch& fa, const MaterialCoefficient* m)
        : F(fa),
          nc(m->nonconductive),
          px(fa.grid()->nx > 1 ? fa.grid()->eps0 * m->epsx * fa.grid()->rdx
                               : 0),
          py(fa.grid()->ny > 1 ? fa.grid()->eps0 * m->epsy * fa.grid()->rdy
                               : 0),
          pz(fa.grid()->nz > 1 ? fa.grid()->eps0 * m->epsz * fa.grid()->rdz : 0)
      {}

      void operator()(int i, int j, int k)
      {
        F(i, j, k).rhob =
          nc * (px * (F(i, j, k).ex - F(i - 1, j, k).ex) +
                py * (F(i, j, k).ey - F(i, j - 1, k).ey) +
                pz * (F(i, j, k).ez - F(i, j, k - 1).ez) - F(i, j, k).rhof);
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

    struct CurlB
    {
      CurlB(typename MfieldsState::Patch& fa, const Grid* g,
            const MaterialCoefficient* m)
        : F(fa),
          px_muz(g->nx > 1 ? g->cvac * g->dt * g->rdx * m->rmuz : 0),
          px_muy(g->nx > 1 ? g->cvac * g->dt * g->rdx * m->rmuy : 0),
          py_mux(g->ny > 1 ? g->cvac * g->dt * g->rdy * m->rmux : 0),
          py_muz(g->ny > 1 ? g->cvac * g->dt * g->rdy * m->rmuz : 0),
          pz_muy(g->nz > 1 ? g->cvac * g->dt * g->rdz * m->rmuy : 0),
          pz_mux(g->nz > 1 ? g->cvac * g->dt * g->rdz * m->rmux : 0)
      {}

      void x(int i, int j, int k)
      {
        F(i, j, k).tcax = (py_muz * (F(i, j, k).cbz - F(i, j - 1, k).cbz) -
                           pz_muy * (F(i, j, k).cby - F(i, j, k - 1).cby));
      }

      void y(int i, int j, int k)
      {
        F(i, j, k).tcay = (pz_mux * (F(i, j, k).cbx - F(i, j, k - 1).cbx) -
                           px_muz * (F(i, j, k).cbz - F(i - 1, j, k).cbz));
      }

      void z(int i, int j, int k)
      {
        F(i, j, k).tcaz = (px_muy * (F(i, j, k).cby - F(i - 1, j, k).cby) -
                           py_mux * (F(i, j, k).cbx - F(i, j - 1, k).cbx));
      }

      F3D F;
      const float px_muz, px_muy, py_mux, py_muz, pz_muy, pz_mux;
    };

    auto& fa = mflds.getPatch(0);
    auto& prm = mflds.params();
    assert(prm.size() == 1);
    const MaterialCoefficient* m = prm[0];

    CurlB curlB(fa, &mflds.vgrid(), m);

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

template <typename MfieldsState>
struct PscDiagOps
{
  using Grid = typename MfieldsState::Grid;
  using F3D = Field3D<typename MfieldsState::Patch>;

  // ----------------------------------------------------------------------
  // vacuum_energy_f

  static void vacuum_energy_f(MfieldsState& mflds, double global[6])
  {
    const auto& g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    auto& prm = mflds.params();
    assert(prm.size() == 1);
    auto m = prm[0];

    const int nx = g.nx, ny = g.ny, nz = g.nz;

    const float qepsx = 0.25 * m->epsx;
    const float qepsy = 0.25 * m->epsy;
    const float qepsz = 0.25 * m->epsz;
    const float hrmux = 0.5 * m->rmux;
    const float hrmuy = 0.5 * m->rmuy;
    const float hrmuz = 0.5 * m->rmuz;
    double en[6] = {};

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          en[0] +=
            qepsx * (sqr(F(i, j, k).ex) + sqr(F(i, j + 1, k).ex) +
                     sqr(F(i, j, k + 1).ex) + sqr(F(i, j + 1, k + 1).ex));
          en[1] +=
            qepsy * (sqr(F(i, j, k).ey) + sqr(F(i, j, k + 1).ey) +
                     sqr(F(i + 1, j, k).ey) + sqr(F(i + 1, j, k + 1).ey));
          en[2] +=
            qepsz * (sqr(F(i, j, k).ez) + sqr(F(i + 1, j, k).ez) +
                     sqr(F(i, j + 1, k).ez) + sqr(F(i + 1, j + 1, k).ez));
          en[3] += hrmux * (sqr(F(i, j, k).cbx) + sqr(F(i + 1, j, k).cbx));
          en[4] += hrmuy * (sqr(F(i, j, k).cby) + sqr(F(i, j + 1, k).cby));
          en[5] += hrmuz * (sqr(F(i, j, k).cbz) + sqr(F(i, j, k + 1).cbz));
        }
      }
    }

    // Convert to physical units
    double v0 = 0.5 * g.eps0 * g.dV;
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
