
#pragma once

#ifdef USE_VPIC

// ======================================================================
// MarderVpicOpsWrap

template <typename _Mparticles, typename _MfieldsState>
struct MarderVpicOpsWrap
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;

  void clear_rhof(MfieldsState& mflds)
  {
    mflds.fa()->kernel->clear_rhof(mflds.fa());
  }
  void synchronize_rho(MfieldsState& mflds)
  {
    mflds.fa()->kernel->synchronize_rho(mflds.fa());
  }
  void compute_div_e_err(MfieldsState& mflds)
  {
    mflds.fa()->kernel->compute_div_e_err(mflds.fa());
  }
  double compute_rms_div_e_err(MfieldsState& mflds)
  {
    return mflds.fa()->kernel->compute_rms_div_e_err(mflds.fa());
  }
  void clean_div_e(MfieldsState& mflds)
  {
    mflds.fa()->kernel->clean_div_e(mflds.fa());
  }
  void compute_div_b_err(MfieldsState& mflds)
  {
    mflds.fa()->kernel->compute_div_b_err(mflds.fa());
  }
  double compute_rms_div_b_err(MfieldsState& mflds)
  {
    return mflds.fa()->kernel->compute_rms_div_b_err(mflds);
  }
  void clean_div_b(MfieldsState& mflds)
  {
    mflds.fa()->kernel->clean_div_b(mflds.fa());
  }
  double synchronize_tang_e_norm_b(MfieldsState& mflds)
  {
    return mflds.fa()->kernel->synchronize_tang_e_norm_b(mflds.fa());
  }

  void accumulate_rho_p(Mparticles& mprts, MfieldsState& mflds)
  {
    auto prts = mprts[0];
    for (auto sp = prts.begin(); sp != prts.end(); ++sp) {
      TIC ::accumulate_rho_p(mflds.fa(), &*sp);
      TOC(accumulate_rho_p, 1);
    }
  }
};

#endif

// ======================================================================
// MarderVpicOps

template <typename _Mparticles, typename _MfieldsState>
struct MarderVpicOps
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Grid = typename MfieldsState::Grid;
  using MaterialCoefficient = typename MfieldsState::MaterialCoefficient;
  using LocalOps = PscFieldArrayLocalOps<MfieldsState>;
  using RemoteOps = PscFieldArrayRemoteOps<MfieldsState>;
  using F3D = Field3D<typename MfieldsState::Patch>;

  // ----------------------------------------------------------------------
  // clear_rhof

  void clear_rhof(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    const int nv = mflds.vgrid().nv;

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

  template <class G, class F3D>
  struct CommRho : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;

    using Base::begin;
    using Base::end;

    using Base::buf_size_;
    using Base::g_;
    using Base::nx_;

    CommRho(Grid& g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
        int Y = (X + 1) % 3, Z = (X + 2) % 3;
        buf_size_[X] = 2 * (nx_[Y] + 1) * (nx_[Z] + 1);
      }
    }

    void begin_send(int X, int side, float* p, F3D& F)
    {
      int face = side ? nx_[X] + 1 : 1;
      foreach_node(g_, X, face,
                   [&](int x, int y, int z) { *p++ = F(x, y, z).rhof; });
      foreach_node(g_, X, face,
                   [&](int x, int y, int z) { *p++ = F(x, y, z).rhob; });
    }

    void end_recv(int X, int side, float* p, F3D& F)
    {
      int face = side ? 1 : nx_[X] + 1;
      foreach_node(g_, X, face,
                   [&](int x, int y, int z) { F(x, y, z).rhof += *p++; });
      foreach_node(g_, X, face, [&](int x, int y, int z) {
        F(x, y, z).rhob = .5f * (F(x, y, z).rhob + *p++);
      });
    }
  };

  void synchronize_rho(MfieldsState& mflds)
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

  void compute_div_e_err(MfieldsState& mflds)
  {
    struct CalcDivE
    {
      CalcDivE(typename MfieldsState::Patch& fa, const MaterialCoefficient* m)
        : F(fa),
          nc(m->nonconductive),
          px(fa.grid()->nx > 1 ? fa.grid()->eps0 * fa.grid()->rdx : 0),
          py(fa.grid()->ny > 1 ? fa.grid()->eps0 * fa.grid()->rdy : 0),
          pz(fa.grid()->nz > 1 ? fa.grid()->eps0 * fa.grid()->rdz : 0),
          cj(1. / fa.grid()->eps0)
      {}

      void operator()(int i, int j, int k)
      {
        F(i, j, k).div_e_err = nc * (px * (F(i, j, k).ex - F(i - 1, j, k).ex) +
                                     py * (F(i, j, k).ey - F(i, j - 1, k).ey) +
                                     pz * (F(i, j, k).ez - F(i, j, k - 1).ez) -
                                     cj * (F(i, j, k).rhof + F(i, j, k).rhob));
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

  double compute_rms_div_e_err(MfieldsState& mflds)
  {
    const auto& g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    const int nx = g.nx, ny = g.ny, nz = g.nz;

    // Interior points
    // FIXME, it's inconsistent to calc the square in single prec here, but
    // double prec later
    double err = 0;
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
        for (int i = 2; i <= nx; i++) {
          err += sqr(F(i, j, k).div_e_err);
        }
      }
    }

    int x, y, z;
    // Exterior faces

    for (y = 2; y <= ny; y++) {
      for (z = 2; z <= nz; z++) {
        err += 0.5 * sqr((double)(F(1, y, z).div_e_err));
        err += 0.5 * sqr((double)(F(nx + 1, y, z).div_e_err));
      }
    }

    for (z = 2; z <= nz; z++) {
      for (x = 2; x <= nx; x++) {
        err += 0.5 * sqr((double)(F(x, 1, z).div_e_err));
        err += 0.5 * sqr((double)(F(x, ny + 1, z).div_e_err));
      }
    }

    for (x = 2; x <= nx; x++) {
      for (y = 2; y <= ny; y++) {
        err += 0.5 * sqr((double)(F(x, y, 1).div_e_err));
        err += 0.5 * sqr((double)(F(x, y, nz + 1).div_e_err));
      }
    }

    // Exterior edges

    for (x = 2; x <= nx; x++) {
      err += 0.25 * sqr((double)(F(x, 1, 1).div_e_err));
      err += 0.25 * sqr((double)(F(x, ny + 1, 1).div_e_err));
      err += 0.25 * sqr((double)(F(x, 1, nz + 1).div_e_err));
      err += 0.25 * sqr((double)(F(x, ny + 1, nz + 1).div_e_err));
    }

    for (y = 2; y <= ny; y++) {
      err += 0.25 * sqr((double)(F(1, y, 1).div_e_err));
      err += 0.25 * sqr((double)(F(1, y, nz + 1).div_e_err));
      err += 0.25 * sqr((double)(F(nx + 1, y, 1).div_e_err));
      err += 0.25 * sqr((double)(F(nx + 1, y, nz + 1).div_e_err));
    }

    for (z = 2; z <= nz; z++) {
      err += 0.25 * sqr((double)(F(1, 1, z).div_e_err));
      err += 0.25 * sqr((double)(F(nx + 1, 1, z).div_e_err));
      err += 0.25 * sqr((double)(F(1, ny + 1, z).div_e_err));
      err += 0.25 * sqr((double)(F(nx + 1, ny + 1, z).div_e_err));
    }

    // Exterior corners

    err += 0.125 * sqr((double)(F(1, 1, 1).div_e_err));
    err += 0.125 * sqr((double)(F(nx + 1, 1, 1).div_e_err));
    err += 0.125 * sqr((double)(F(1, ny + 1, 1).div_e_err));
    err += 0.125 * sqr((double)(F(nx + 1, ny + 1, 1).div_e_err));
    err += 0.125 * sqr((double)(F(1, 1, nz + 1).div_e_err));
    err += 0.125 * sqr((double)(F(nx + 1, 1, nz + 1).div_e_err));
    err += 0.125 * sqr((double)(F(1, ny + 1, nz + 1).div_e_err));
    err += 0.125 * sqr((double)(F(nx + 1, ny + 1, nz + 1).div_e_err));

    double local[2], _global[2]; // FIXME, name clash with global macro
    local[0] = err * g.dV;
    local[1] = (nx * ny * nz) * g.dV;
    MPI_Allreduce(local, _global, 2, MPI_DOUBLE, MPI_SUM, psc_comm_world);
    return g.eps0 * sqrt(_global[0] / _global[1]);
  }

  // ----------------------------------------------------------------------
  // clean_div_e

#define MARDER_EX(i, j, k)                                                     \
  F(i, j, k).ex += px * (F(i + 1, j, k).div_e_err - F(i, j, k).div_e_err)
#define MARDER_EY(i, j, k)                                                     \
  F(i, j, k).ey += py * (F(i, j + 1, k).div_e_err - F(i, j, k).div_e_err)
#define MARDER_EZ(i, j, k)                                                     \
  F(i, j, k).ez += pz * (F(i, j, k + 1).div_e_err - F(i, j, k).div_e_err)

  void vacuum_clean_div_e(MfieldsState& mflds)
  {
    const auto& g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    auto& prm = mflds.params();
    assert(prm.size() == 1);
    const MaterialCoefficient* m = prm[0];

    const int nx = g.nx, ny = g.ny, nz = g.nz;

    const float _rdx = (nx > 1) ? g.rdx : 0;
    const float _rdy = (ny > 1) ? g.rdy : 0;
    const float _rdz = (nz > 1) ? g.rdz : 0;
    const float alphadt = 0.3888889 / (_rdx * _rdx + _rdy * _rdy + _rdz * _rdz);
    const float px = (alphadt * _rdx) * m->drivex;
    const float py = (alphadt * _rdy) * m->drivey;
    const float pz = (alphadt * _rdz) * m->drivez;

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          MARDER_EX(i, j, k);
          MARDER_EY(i, j, k);
          MARDER_EZ(i, j, k);
        }
      }
    }

    int x, y, z;

    // Do left over ex
    for (y = 1; y <= ny + 1; y++) {
      for (x = 1; x <= nx; x++) {
        MARDER_EX(x, y, nz + 1);
      }
    }
    for (z = 1; z <= nz; z++) {
      for (x = 1; x <= nx; x++) {
        MARDER_EX(x, ny + 1, z);
      }
    }

    // Do left over ey
    for (z = 1; z <= nz + 1; z++) {
      for (y = 1; y <= ny; y++) {
        MARDER_EY(nx + 1, y, z);
      }
    }
    for (y = 1; y <= ny; y++) {
      for (x = 1; x <= nx; x++) {
        MARDER_EY(x, y, nz + 1);
      }
    }

    // Do left over ez
    for (z = 1; z <= nz; z++) {
      for (x = 1; x <= nx + 1; x++) {
        MARDER_EZ(x, ny + 1, z);
      }
    }
    for (z = 1; z <= nz; z++) {
      for (y = 1; y <= ny; y++) {
        MARDER_EZ(nx + 1, y, z);
      }
    }

    LocalOps::local_adjust_tang_e(mflds);
  }

  void clean_div_e(MfieldsState& mflds)
  {
    vacuum_clean_div_e(mflds); // FIXME eventually?
  }

  // ----------------------------------------------------------------------
  // compute_div_b_err

  void compute_div_b_err(MfieldsState& mflds)
  {
    const auto& g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    const int nx = g.nx, ny = g.ny, nz = g.nz;
    const float px =
      (nx > 1) ? g.rdx : 0; // FIXME, should be based on global dims
    const float py = (ny > 1) ? g.rdy : 0;
    const float pz = (nz > 1) ? g.rdz : 0;

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          F(i, j, k).div_b_err = (px * (F(i + 1, j, k).cbx - F(i, j, k).cbx) +
                                  py * (F(i, j + 1, k).cby - F(i, j, k).cby) +
                                  pz * (F(i, j, k + 1).cbz - F(i, j, k).cbz));
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // compute_rms_div_b_err
  //
  // OPT: doing that at the same time as div_b should be faster

  double compute_rms_div_b_err(MfieldsState& mflds)
  {
    const auto& g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    const int nx = g.nx, ny = g.ny, nz = g.nz;

    double err = 0;
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          err += sqr(F(i, j, k).div_b_err);
        }
      }
    }

    double local[2], _global[2]; // FIXME, name clash with global macro
    local[0] = err * g.dV;
    local[1] = (nx * ny * nz) * g.dV;
    MPI_Allreduce(local, _global, 2, MPI_DOUBLE, MPI_SUM, psc_comm_world);
    return g.eps0 * sqrt(_global[0] / _global[1]);
  }

  // ----------------------------------------------------------------------
  // clean_div_b

#define MARDER_CBX(i, j, k)                                                    \
  F(i, j, k).cbx += px * (F(i, j, k).div_b_err - F(i - 1, j, k).div_b_err)
#define MARDER_CBY(i, j, k)                                                    \
  F(i, j, k).cby += py * (F(i, j, k).div_b_err - F(i, j - 1, k).div_b_err)
#define MARDER_CBZ(i, j, k)                                                    \
  F(i, j, k).cbz += pz * (F(i, j, k).div_b_err - F(i, j, k - 1).div_b_err)

  void clean_div_b(MfieldsState& mflds)
  {
    const auto& g = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    F3D F(fa);

    const int nx = g.nx, ny = g.ny, nz = g.nz;
    float px = (nx > 1) ? g.rdx : 0;
    float py = (ny > 1) ? g.rdy : 0;
    float pz = (nz > 1) ? g.rdz : 0;
    float alphadt = 0.3888889 / (px * px + py * py + pz * pz);
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
          MARDER_CBX(i, j, k);
          MARDER_CBY(i, j, k);
          MARDER_CBZ(i, j, k);
        }
      }
    }

    int x, y, z;

    // Do left over interior bx
    for (y = 1; y <= ny; y++) {
      for (x = 2; x <= nx; x++) {
        MARDER_CBX(x, y, 1);
      }
    }
    for (z = 2; z <= nz; z++) {
      for (x = 2; x <= nx; x++) {
        MARDER_CBX(x, 1, z);
      }
    }

    // Left over interior by
    for (z = 1; z <= nz; z++) {
      for (y = 2; y <= ny; y++) {
        MARDER_CBY(1, y, z);
      }
    }
    for (y = 2; y <= ny; y++) {
      for (x = 2; x <= nx; x++) {
        MARDER_CBY(x, y, 1);
      }
    }

    // Left over interior bz
    for (z = 2; z <= nz; z++) {
      for (x = 1; x <= nx; x++) {
        MARDER_CBZ(x, 1, z);
      }
    }
    for (z = 2; z <= nz; z++) {
      for (y = 2; y <= ny; y++) {
        MARDER_CBZ(1, y, z);
      }
    }

    // Finish setting derr ghosts

    RemoteOps::end_remote_ghost_div_b(mflds);

    // Do Marder pass in exterior

    // Exterior bx
    for (z = 1; z <= nz; z++) {
      for (y = 1; y <= ny; y++) {
        MARDER_CBX(1, y, z);
      }
    }
    for (z = 1; z <= nz; z++) {
      for (y = 1; y <= ny; y++) {
        MARDER_CBX(nx + 1, y, z);
      }
    }

    // Exterior by
    for (z = 1; z <= nz; z++) {
      for (x = 1; x <= nx; x++) {
        MARDER_CBY(x, 1, z);
      }
    }
    for (z = 1; z <= nz; z++) {
      for (x = 1; x <= nx; x++) {
        MARDER_CBY(x, ny + 1, z);
      }
    }

    // Exterior bz
    for (y = 1; y <= ny; y++) {
      for (x = 1; x <= nx; x++) {
        MARDER_CBZ(x, y, 1);
      }
    }
    for (y = 1; y <= ny; y++) {
      for (x = 1; x <= nx; x++) {
        MARDER_CBZ(x, y, nz + 1);
      }
    }

    LocalOps::local_adjust_norm_b(mflds);
  }

  // ----------------------------------------------------------------------
  // synchronize_tang_e_norm_b

  template <class G, class F3D>
  struct CommTangENormB : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;

    using Base::begin;
    using Base::end;

    using Base::buf_size_;
    using Base::g_;
    using Base::nx_;

    CommTangENormB(Grid& g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
        int Y = (X + 1) % 3, Z = (X + 2) % 3;
        buf_size_[X] = (nx_[Y] * nx_[Z] + 2 * nx_[Y] * (nx_[Z] + 1) +
                        2 * nx_[Z] * (nx_[Y] + 1));
      }
      err = 0.;
    }

    void begin_send(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? nx_[X] + 1 : 1;
      foreach_face(g_, X, face,
                   [&](int x, int y, int z) { *p++ = (&F(x, y, z).cbx)[X]; });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) {
        *p++ = (&F(x, y, z).ex)[Y];
        *p++ = (&F(x, y, z).tcax)[Y];
      });
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) {
        *p++ = (&F(x, y, z).ex)[Z];
        *p++ = (&F(x, y, z).tcax)[Z];
      });
    }

    void end_recv(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? 1 : nx_[X] + 1;
      foreach_face(g_, X, face, [&](int x, int y, int z) {
        double w1 = *p++, w2 = (&F(x, y, z).cbx)[X];
        (&F(x, y, z).cbx)[X] = .5 * (w1 + w2);
        err += sqr(w1 - w2);
      });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) {
        double w1 = *p++, w2 = (&F(x, y, z).ex)[Y];
        (&F(x, y, z).ex)[Y] = .5 * (w1 + w2);
        err += sqr(w1 - w2);
        w1 = *p++, w2 = (&F(x, y, z).tcax)[Y];
        (&F(x, y, z).tcax)[Y] = .5 * (w1 + w2);
      });
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) {
        double w1 = *p++, w2 = (&F(x, y, z).ex)[Z];
        (&F(x, y, z).ex)[Z] = .5 * (w1 + w2);
        err += sqr(w1 - w2);
        w1 = *p++, w2 = (&F(x, y, z).tcax)[Z];
        (&F(x, y, z).tcax)[Z] = .5 * (w1 + w2);
      });
    }

    double err;
  };

  double synchronize_tang_e_norm_b(MfieldsState& mflds)
  {
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

  // ----------------------------------------------------------------------
  // accumulate_rho_p

  void accumulate_rho_p(MfieldsState& mflds,
                        const typename Mparticles::Species& sp)
  {
    const auto* RESTRICT ALIGNED(128) p = sp.p;

    const auto& g = sp.vgrid();
    const float q_8V = sp.q * g.r8V;
    const int np = sp.np;
    const int sy = g.sy;
    const int sz = g.sz;

    float w0, w1, w2, w3, w4, w5, w6, w7, dz;

    int n, v;

    // Load the grid data

    for (n = 0; n < np; n++) {

      // Load the particle data

      w0 = p[n].dx;
      w1 = p[n].dy;
      dz = p[n].dz;
      v = p[n].i;
      w7 = p[n].w * q_8V;

      // Compute the trilinear weights
      // Though the PPE should have hardware fma/fmaf support, it was
      // measured to be more efficient _not_ to use it here.  (Maybe the
      // compiler isn't actually generating the assembly for it.

#define FMA(x, y, z) ((z) + (x) * (y))
#define FNMS(x, y, z) ((z) - (x) * (y))
      w6 = FNMS(w0, w7, w7); // q(1-dx)
      w7 = FMA(w0, w7, w7);  // q(1+dx)
      w4 = FNMS(w1, w6, w6);
      w5 = FNMS(w1, w7, w7); // q(1-dx)(1-dy), q(1+dx)(1-dy)
      w6 = FMA(w1, w6, w6);
      w7 = FMA(w1, w7, w7); // q(1-dx)(1+dy), q(1+dx)(1+dy)
      w0 = FNMS(dz, w4, w4);
      w1 = FNMS(dz, w5, w5);
      w2 = FNMS(dz, w6, w6);
      w3 = FNMS(dz, w7, w7);
      w4 = FMA(dz, w4, w4);
      w5 = FMA(dz, w5, w5);
      w6 = FMA(dz, w6, w6);
      w7 = FMA(dz, w7, w7);
#undef FNMS
#undef FMA

      // Reduce the particle charge to rhof

      auto& fa = mflds.getPatch(0);
      fa[v].rhof += w0;
      fa[v + 1].rhof += w1;
      fa[v + sy].rhof += w2;
      fa[v + sy + 1].rhof += w3;
      fa[v + sz].rhof += w4;
      fa[v + sz + 1].rhof += w5;
      fa[v + sz + sy].rhof += w6;
      fa[v + sz + sy + 1].rhof += w7;
    }
  }

  void accumulate_rho_p(Mparticles& mprts, MfieldsState& mflds)
  {
    for (auto& sp : mprts[0]) {
      TIC accumulate_rho_p(mflds, sp);
      TOC(accumulate_rho_p, 1);
    }
  }
};

// ======================================================================
// MarderVpic_

template <typename Ops>
struct MarderVpic_ : Ops
{
  using Mparticles = typename Ops::Mparticles;
  using MfieldsState = typename Ops::MfieldsState;
  using real_t = typename MfieldsState::real_t;

  MarderVpic_(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : comm_{grid.comm()}
  {}

  // ----------------------------------------------------------------------
  // psc_marder_vpic_clean_div_e

  void psc_marder_vpic_clean_div_e(MfieldsState& mflds, Mparticles& mprts)
  {
    mpi_printf(comm_, "Divergence cleaning electric field\n");

    this->clear_rhof(mflds);
    this->accumulate_rho_p(mprts, mflds);
    this->synchronize_rho(mflds);

    for (int round = 0; round < num_div_e_round_; round++) {
      this->compute_div_e_err(mflds);
      if (round == 0 || round == num_div_e_round_ - 1) {
        double err = this->compute_rms_div_e_err(mflds);
        mpi_printf(comm_, "%s rms error = %e (charge/volume)\n",
                   round == 0 ? "Initial" : "Cleaned", err);
      }
      this->clean_div_e(mflds);
    }
  }

  // ----------------------------------------------------------------------
  // psc_marder_vpic_clean_div_b

  void psc_marder_vpic_clean_div_b(MfieldsState& mflds)
  {
    mpi_printf(comm_, "Divergence cleaning magnetic field\n");

    for (int round = 0; round < num_div_b_round_; round++) {
      this->compute_div_b_err(mflds);
      if (round == 0 || round == num_div_b_round_ - 1) {
        double err = this->compute_rms_div_b_err(mflds);
        mpi_printf(comm_, "%s rms error = %e (charge/volume)\n",
                   round == 0 ? "Initial" : "Cleaned", err);
      }
      this->clean_div_b(mflds);
    }
  }

  // ----------------------------------------------------------------------
  // run

  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    const auto& grid = mflds.grid();
    bool clean_div_e = (clean_div_e_interval_ > 0 &&
                        grid.timestep() % clean_div_e_interval_ == 0);
    bool clean_div_b = (clean_div_b_interval_ > 0 &&
                        grid.timestep() % clean_div_b_interval_ == 0);
    bool sync_shared = (sync_shared_interval_ > 0 &&
                        grid.timestep() % sync_shared_interval_ == 0);

    if (!(clean_div_e || clean_div_b || sync_shared)) {
      return;
    }

    // Divergence clean e
    if (clean_div_e) {
      // needs E, rhof, rhob, material
      psc_marder_vpic_clean_div_e(mflds, mprts);
      // upates E, rhof, div_e_err
    }

    // Divergence clean b
    if (clean_div_b) {
      // needs B
      psc_marder_vpic_clean_div_b(mflds);
      // updates B, div_b_err
    }

    // Synchronize the shared faces
    if (sync_shared) {
      // needs E, B, TCA
      mpi_printf(comm_, "Synchronizing shared tang e, norm b\n");
      double err = synchronize_tang_e_norm_b(mflds);
      mpi_printf(comm_, "Domain desynchronization error = %e (arb units)\n",
                 err);
      // updates E, B, TCA
    }
  }

  // ----------------------------------------------------------------------
  // the following are profiling wrappers around the methods inherited through
  // Ops

  void clear_rhof(MfieldsState& mflds)
  {
    TIC Ops::clear_rhof(mflds);
    TOC(clear_rhof, 1);
  }

  void synchronize_rho(MfieldsState& mflds)
  {
    TIC Ops::synchronize_rho(mflds);
    TOC(synchronize_rho, 1);
  }

  void compute_div_e_err(MfieldsState& mflds)
  {
    TIC Ops::compute_div_e_err(mflds);
    TOC(compute_div_e_err, 1);
  }

  double compute_rms_div_e_err(MfieldsState& mflds)
  {
    double err;
    TIC err = Ops::compute_rms_div_e_err(mflds);
    TOC(compute_rms_div_e_err, 1);
    return err;
  }

  void clean_div_e(MfieldsState& mflds)
  {
    TIC Ops::clean_div_e(mflds);
    TOC(clean_div_e, 1);
  }

  void compute_div_b_err(MfieldsState& mflds)
  {
    TIC Ops::compute_div_b_err(mflds);
    TOC(compute_div_b_err, 1);
  }

  double compute_rms_div_b_err(MfieldsState& mflds)
  {
    double err;
    TIC err = Ops::compute_rms_div_b_err(mflds);
    TOC(compute_rms_div_e_err, 1);
    return err;
  }

  void clean_div_b(MfieldsState& mflds)
  {
    TIC Ops::clean_div_b(mflds);
    TOC(clean_div_e, 1);
  }

  double synchronize_tang_e_norm_b(MfieldsState& mflds)
  {
    double err;
    TIC err = Ops::synchronize_tang_e_norm_b(mflds);
    TOC(synchronize_tang_e_norm_b, 1);
    return err;
  }

  void accumulate_rho_p(Mparticles& mprts, MfieldsState& mflds)
  {
    Ops::accumulate_rho_p(mprts, mflds);
  }

private:
  MPI_Comm comm_;
  int clean_div_e_interval_ = 10; // FIXME, hardcoded...
  int clean_div_b_interval_ = 10;
  int sync_shared_interval_ = 10;
  int num_div_e_round_ = 3;
  int num_div_b_round_ = 3;
};

#ifdef USE_VPIC
template <typename Mparticles, typename MfieldsState>
using MarderVpicWrap = MarderVpic_<MarderVpicOpsWrap<Mparticles, MfieldsState>>;
#endif

template <typename Mparticles, typename MfieldsState>
using MarderVpic = MarderVpic_<MarderVpicOps<Mparticles, MfieldsState>>;
