
#ifndef GRID_HXX
#define GRID_HXX

#include "vec3.hxx"
#include "psc_bits.h"
#include "mrc_domain.hxx"

#include <vector>
#include <cstring>
#include <cmath>

///Possible boundary conditions for fields
enum {
  BND_FLD_OPEN,
  BND_FLD_PERIODIC,
  BND_FLD_CONDUCTING_WALL,
  BND_FLD_ABSORBING,
};

///Possible boundary conditions for particles
enum {
  BND_PRT_REFLECTING,
  BND_PRT_PERIODIC,
  BND_PRT_ABSORBING,
  BND_PRT_OPEN,
};

///Describes the spatial domain to operate on.
///
///This struct describes the spatial dimension of the simulation-box
///@note Here, you can also set the dimensionality by eliminating a dimension. Example: To simulate in xy only, set
///\verbatim psc_domain.gdims[2]=1 \endverbatim
///Also, set the boundary conditions for the eliminated dimensions to BND_FLD_PERIODIC or you'll get invalid \a dt and \a dx

struct GridBc
{
  GridBc() = default; // FIXME
  
  GridBc(Int3 fld_lo, Int3 fld_hi, Int3 prt_lo, Int3 prt_hi)
    : fld_lo(fld_lo), fld_hi(fld_hi),
      prt_lo(prt_lo), prt_hi(prt_hi)
  {}
  
  Int3 fld_lo;	///<Boundary conditions of the fields. Can be any value of BND_FLD.
  Int3 fld_hi;	///<Boundary conditions of the fields. Can be any value of BND_FLD.
  Int3 prt_lo;	///<Boundary conditions of the particles. Can be any value of BND_PART.
  Int3 prt_hi;  ///<Boundary conditions of the particles. Can be any value of BND_PART.
};

// ======================================================================
// Grid_

template<class T>
struct Grid_
{
  using real_t = T;
  using Real3 = Vec3<real_t>;

  struct NormalizationParams;
  struct Normalization;

  struct Domain;
  
  struct Kind;
  using Kinds = std::vector<Kind>;
  
  struct Patch
  {
    Patch(const Int3& _off, const Real3& _xb, const Real3& _xe, const Real3& dx)
      : off(_off), xb(_xb), xe(_xe), dx_(dx)
    {}

    real_t x_nc(int i) const { return xb[0] + i * dx_[0]; }
    real_t y_nc(int j) const { return xb[1] + j * dx_[1]; }
    real_t z_nc(int k) const { return xb[2] + k * dx_[2]; }
    
    real_t x_cc(int i) const { return xb[0] + (i + .5f) * dx_[0]; }
    real_t y_cc(int j) const { return xb[1] + (j + .5f) * dx_[1]; }
    real_t z_cc(int k) const { return xb[2] + (k + .5f) * dx_[2]; }
    
    Int3 off;
    Real3 xb;
    Real3 xe;
  private:
    Real3 dx_;
  };

  // Construct a single≈≈ patch covering the whole domain
  // Args: global dimensions, and length of the domain in all 3 dims
  // -- mostly useful ≈for testing
  //
  // Maybe the named ctor idiom would be good here (but right now
  // can't be done since the copy ctor is deleted.
  Grid_(const Domain& domain)
    : domain{domain},
      ldims{domain.ldims}
  {
    patches.emplace_back(Patch({ 0, 0, 0 }, { 0., 0., 0.}, domain.length, domain.dx));
  }

  Grid_(const Domain& domain, const std::vector<Int3>& offs)
    : domain{domain},
      ldims{domain.ldims}
  {
    for (auto off : offs) {
      patches.push_back(Patch(off,
			      Vec3<double>(off        ) * domain.dx + domain.corner,
			      Vec3<double>(off + ldims) * domain.dx + domain.corner,
			      domain.dx));
    }
  }

  int n_patches() const { return patches.size(); }

  bool isInvar(int d) const { return domain.isInvar(d); }

  template<typename FUNC>
  void Foreach_3d(int l, int r, FUNC F) const
  {
    int __ilo[3] = { isInvar(0) ? 0 : -l,
		     isInvar(1) ? 0 : -l,
		     isInvar(2) ? 0 : -l };
    int __ihi[3] = { ldims[0] + (isInvar(0) ? 0 : r),
		     ldims[1] + (isInvar(1) ? 0 : r),
		     ldims[2] + (isInvar(2) ? 0 : r) };
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
  GridBc bc;
  Normalization norm;
  // FIXME? this default might be a bit dangerous, as they're useful
  // for testing but might hide if one forgets to init dt correctly
  // for a real simulation
  real_t dt = { 1. };
  std::vector<Patch> patches;
  std::vector<Kind> kinds;

  MrcDomain mrc_domain() const { assert(mrc_domain_); return {mrc_domain_}; }

  MPI_Comm comm() const { return mrc_domain().comm(); }
  int nGlobalPatches() const { return mrc_domain().nGlobalPatches(); }
  void neighborRankPatch(int p, int dir[3], int* nei_rank, int* nei_patch) const { mrc_domain().neighborRankPatch(p, dir, nei_rank, nei_patch); }
  mrc_patch_info globalPatchInfo(int p) const { return mrc_domain().globalPatchInfo(p); }
  mrc_patch_info localPatchInfo(int p) const { return mrc_domain().localPatchInfo(p); }
  mrc_ddc* create_ddc() const { return mrc_domain().create_ddc(); }

  struct mrc_domain* mrc_domain_ = {};
};

// ======================================================================
// Grid::Domain

template<class T>
struct Grid_<T>::Domain
{
  Domain(Int3 gdims, Real3 length, Real3 corner = {0., 0., 0.}, Int3 np = {1, 1, 1})
    : gdims(gdims), length(length), corner(corner), np(np)
  {
    for (int d = 0; d < 3; d++) {
      assert(gdims[d] % np[d] == 0);
      ldims[d] = gdims[d] / np[d];
    }
    dx = length / Real3(gdims);
  }

  bool isInvar(int d) const { return gdims[d] == 1; }
  
  Int3 gdims;		///<Number of grid-points in each dimension
  Real3 length;	///<The physical size of the simulation-box 
  Real3 corner;
  Int3 np;		///<Number of patches in each dimension
  
  Int3 ldims;
  Real3 dx;
};

// ======================================================================
// Grid::NormalizationParams

template<class T>
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
    prm.lw = 2.*M_PI;
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
    lw   = 3.2e-6;
    i0   = 1e21;
    n0   = 1e26;
    e0   = 0.;
    
    nicell = 0;
  }

public:
  double cc;   // speed of light
  double qq;   // elemental charge 
  double mm;   // mass
  double tt;   // some measurement for energy ? (default is 1keV in fortran) 
  double eps0; // vacuum permittivity
  double e0;   // field intensity

  double lw;   // normalization coefficient for laser wavelength (omega)
  double i0;   // laser intensity
  double n0;   // electron density

  int nicell;  // number of particles per gridpoint to represent a normalized density of 1 
};

// ======================================================================
// Grid::Normalization

template<class T>
struct Grid_<T>::Normalization
{
  // ----------------------------------------------------------------------
  // ctor

  Normalization() // FIXME
  {}
  
  Normalization(NormalizationParams& prm)
  {
    assert(prm.nicell > 0);
    cc = prm.cc;

    double wl = 2. * M_PI * cc / prm.lw;
    double ld = cc / wl;
    assert(ld == 1.); // FIXME, not sure why? (calculation of fnqs?)
    if (prm.e0 == 0.) {
      prm.e0 = sqrt(2.0 * prm.i0 / prm.eps0 / cc) /
	prm.lw / 1.0e6;
    }
    b0 = prm.e0 / cc;
    rho0 = prm.eps0 * wl * b0;
    phi0 = ld * prm.e0;
    a0 = prm.e0 / wl;
    
    double vos = prm.qq * prm.e0 / (prm.mm * wl);
    double vt = sqrt(prm.tt / prm.mm);
    double wp = sqrt(sqr(prm.qq) * prm.n0 / prm.eps0 / prm.mm);
    
    cori = 1. / prm.nicell;
    double alpha_ = wp / wl;
    beta = vt / cc;
    eta = vos / cc;
    fnqs = sqr(alpha_) * cori / eta;
  }

  real_t cc;
  real_t fnqs = { 1. };
  real_t eta = { 1. };
  real_t beta = { 1. };
  real_t cori = { 1. };

  real_t b0;
  real_t rho0;
  real_t phi0;
  real_t a0;
};

// ======================================================================
// Grid::Kind

template<class T>
struct Grid_<T>::Kind
{
  Kind() // FIXME, do we want to keep this ctor?
  {}
  
  Kind(real_t q_, real_t m_, const char *name_)
    : q(q_), m(m_), name(strdup(name_))
  {};
  
  Kind(const Kind& k)
    : q(k.q), m(k.m), name(strdup(k.name))
  {}
  
  Kind(Kind&& k)
    : q(k.q), m(k.m), name(k.name)
  {
    k.name = nullptr;
  }

  Kind& operator=(const Kind& other)
  {
    free((void*) name);
    q = other.q;
    m = other.m;
    name = strdup(other.name);
    return *this;
  }
  
  ~Kind()
  {
    free((void*) name);
  }
  
  real_t q;
  real_t m;
  const char *name;
};

using Grid_t = Grid_<double>;

#endif

