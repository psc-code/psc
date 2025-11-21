
#ifndef GRID_HXX
#define GRID_HXX

#include <balance.hxx>
#include "kg/Vec3.h"
#include "psc_bits.h"
#include "mrc_domain.hxx"
#include "grid/BC.h"
#include "grid/domain.hxx"
#include <mrc_ddc.h>

#include <vector>
#include <string>
#include <cmath>

// ======================================================================
// Grid_

template <class T>
struct Grid_
{
  using real_t = T;
  using Real3 = Vec3<real_t>;

  struct NormalizationParams;
  struct Normalization;

  using Domain = psc::grid::Domain<real_t>;

  struct Kind;
  using Kinds = std::vector<Kind>;

  struct Patch
  {
    Patch() {}

    Patch(const Int3& _off, const Real3& _xb, const Real3& _xe, const Real3& dx)
      : off(_off), xb(_xb), xe(_xe), dx_(dx)
    {}

    real_t x_nc(int i) const { return xb[0] + i * dx_[0]; }
    real_t y_nc(int j) const { return xb[1] + j * dx_[1]; }
    real_t z_nc(int k) const { return xb[2] + k * dx_[2]; }
    real_t get_nc(int i, int axis) const { return xb[axis] + i * dx_[axis]; }

    real_t x_cc(int i) const { return xb[0] + (i + .5f) * dx_[0]; }
    real_t y_cc(int j) const { return xb[1] + (j + .5f) * dx_[1]; }
    real_t z_cc(int k) const { return xb[2] + (k + .5f) * dx_[2]; }
    real_t get_cc(int i, int axis) const
    {
      return xb[axis] + (i + .5f) * dx_[axis];
    }

    Int3 off;
    Real3 xb;
    Real3 xe;

  private:
    Real3 dx_;
  };

  // ----------------------------------------------------------------------
  // ctor

  Grid_() {}

  Grid_(const Domain& domain, const psc::grid::BC& bc, const Kinds& kinds,
        const Normalization& norm, double dt, int n_patches = -1, Int3 ibn = {})
    : domain{domain},
      ldims{domain.ldims},
      bc{bc},
      kinds{kinds},
      norm{norm},
      dt{dt},
      ibn{ibn},
      mrc_domain_{domain, bc, n_patches}
  {
#if 0
    mpi_printf(MPI_COMM_WORLD, "::: dt      = %g\n", dt);
    mpi_printf(MPI_COMM_WORLD, "::: dx      = %g %g %g\n", domain.dx[0], domain.dx[1], domain.dx[2]);
#endif

    for (auto off : mrc_domain_.offs()) {
      patches.push_back(Patch(
        off, Vec3<double>(off) * domain.dx + domain.corner,
        Vec3<double>(off + ldims) * domain.dx + domain.corner, domain.dx));
    }

    for (int d = 0; d < 3; d++) {
      if (isInvar(d)) {
        // if invariant in this direction: set bnd to periodic (FIXME?)
        this->bc.fld_lo[d] = BND_FLD_PERIODIC;
        this->bc.fld_hi[d] = BND_FLD_PERIODIC;
        this->bc.prt_lo[d] = BND_PRT_PERIODIC;
        this->bc.prt_hi[d] = BND_PRT_PERIODIC;
      }
    }

    reset_ddc();
  }

  void view() const
  {
    mprintf("Grid_::view()\n");
    domain.view();
    mprintf("Grid_: ldims %d x %d x %d\n", ldims[0], ldims[1], ldims[2]);
  }

  int n_patches() const { return patches.size(); }

  bool isInvar(int d) const { return domain.isInvar(d); }

  bool atBoundaryLo(int p, int d) const { return patches[p].off[d] == 0; }
  bool atBoundaryHi(int p, int d) const
  {
    return patches[p].off[d] + ldims[d] == domain.gdims[d];
  }

  int timestep() const { return timestep_; }

  template <typename FUNC>
  void Foreach_3d(int l, int r, FUNC F) const
  {
    int __ilo[3] = {isInvar(0) ? 0 : -l, isInvar(1) ? 0 : -l,
                    isInvar(2) ? 0 : -l};
    int __ihi[3] = {ldims[0] + (isInvar(0) ? 0 : r),
                    ldims[1] + (isInvar(1) ? 0 : r),
                    ldims[2] + (isInvar(2) ? 0 : r)};
    for (int k = __ilo[2]; k < __ihi[2]; k++) {
      for (int j = __ilo[1]; j < __ihi[1]; j++) {
        for (int i = __ilo[0]; i < __ihi[0]; i++) {
          F(i, j, k);
        }
      }
    }
  }

  Int3 ldims;
  Domain domain;
  psc::grid::BC bc;
  Normalization norm;
  // FIXME? this default might be a bit dangerous, as they're useful
  // for testing but might hide if one forgets to init dt correctly
  // for a real simulation
  real_t dt = {1.};
  std::vector<Patch> patches;
  std::vector<Kind> kinds;
  Int3 ibn; // FIXME
  int timestep_ = 0;
  MrcDomain mrc_domain_;
  mutable mrc_ddc* ddc_ = {};
  mutable int balance_generation_cnt_;

  const MrcDomain& mrc_domain() const { return mrc_domain_; }

  mrc_ddc* ddc() const
  {
    if (balance_generation_cnt_ != psc_balance_generation_cnt) {
      reset_ddc();
    }
    assert(ddc_);
    return ddc_;
  }

  MPI_Comm comm() const { return MPI_COMM_WORLD; }
  int nGlobalPatches() const { return mrc_domain().nGlobalPatches(); }
  void neighborRankPatch(int p, int dir[3], int* nei_rank, int* nei_patch) const
  {
    mrc_domain().neighborRankPatch(p, dir, nei_rank, nei_patch);
  }
  mrc_patch_info globalPatchInfo(int p) const
  {
    return mrc_domain().globalPatchInfo(p);
  }
  mrc_patch_info localPatchInfo(int p) const
  {
    return mrc_domain().localPatchInfo(p);
  }
  mrc_ddc* create_ddc() const { return mrc_domain().create_ddc(); }

  void reset_ddc() const
  {
    if (ddc_) {
      mrc_ddc_destroy(ddc_);
    }
    if (ibn[0] > 0 || ibn[1] > 0 || ibn[2] > 0) {
      ddc_ = create_ddc();
      mrc_ddc_set_param_int3(ddc_, "ibn", ibn);
      mrc_ddc_set_param_int(ddc_, "max_n_fields", 24);
      mrc_ddc_setup(ddc_);
    }
    balance_generation_cnt_ = psc_balance_generation_cnt;
  }
};

// ======================================================================
// Grid::NormalizationParams

template <class T>
struct Grid_<T>::NormalizationParams
{
  static NormalizationParams dimensionless()
  {
    auto prm = NormalizationParams{};

    prm.cc = 1.;
    prm.qq = 1.;
    prm.mm = 1.;
    prm.tt = 1.;
    prm.eps0 = 1.;
    prm.e0 = 1.;
    prm.lw = 2. * M_PI;
    prm.i0 = 0.;
    prm.n0 = 1.;

    return prm;
  }

  NormalizationParams()
  {
    qq = 1.6021e-19;
    mm = 9.1091e-31;
    tt = 1.6021e-16;
    cc = 3.0e8;

    eps0 = 8.8542e-12;
    lw = 3.2e-6;
    i0 = 1e21;
    n0 = 1e26;
    e0 = 0.;

    nicell = 0;
  }

public:
  double cc;   // speed of light
  double qq;   // elemental charge
  double mm;   // mass
  double tt;   // some measurement for energy ? (default is 1keV in fortran)
  double eps0; // vacuum permittivity
  double e0;   // field intensity

  double lw; // normalization coefficient for laser wavelength (omega)
  double i0; // laser intensity
  double n0; // electron density

  int nicell; // number of particles per gridpoint to represent a normalized
              // density of 1
};

// ======================================================================
// Grid::Normalization

template <class T>
struct Grid_<T>::Normalization
{
  // ----------------------------------------------------------------------
  // ctor

  Normalization() = default; // FIXME

  Normalization(const NormalizationParams& prm)
  {
    assert(prm.nicell > 0);
    cc = prm.cc;
    eps0 = prm.eps0;

    double wl = 2. * M_PI * cc / prm.lw;
    double ld = cc / wl;
    assert(ld == 1.); // FIXME, not sure why? (calculation of fnqs?)
    e0 = prm.e0;
    if (e0 == 0.) {
      e0 = sqrt(2.0 * prm.i0 / prm.eps0 / cc) / prm.lw / 1.0e6;
    }
    b0 = e0 / cc;
    rho0 = prm.eps0 * wl * b0;
    phi0 = ld * prm.e0;
    a0 = e0 / wl;

    double vos = prm.qq * prm.e0 / (prm.mm * wl);
    double vt = sqrt(prm.tt / prm.mm);
    double wp = sqrt(sqr(prm.qq) * prm.n0 / prm.eps0 / prm.mm);

    prts_per_unit_density = prm.nicell;
    cori = 1. / prm.nicell;
    double alpha_ = wp / wl;
    beta = vt / cc;
    eta = vos / cc;
    fnqs = sqr(alpha_) * cori / eta;
  }

  real_t cc;
  real_t eps0;
  real_t fnqs = {1.};
  real_t eta = {1.};
  real_t beta = {1.};
  real_t cori = {1.};
  real_t prts_per_unit_density;

  real_t e0;
  real_t b0;
  real_t rho0;
  real_t phi0;
  real_t a0;
};

// ======================================================================
// Grid::Kind

template <class T>
struct Grid_<T>::Kind
{
  Kind() = default;

  Kind(real_t q, real_t m, const std::string& name) : q(q), m(m), name(name) {}

  real_t q = {};
  real_t m = {};
  std::string name;
};

using Grid_t = Grid_<double>;

#endif
