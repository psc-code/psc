
#ifndef GRID_HXX
#define GRID_HXX

#include "vec3.hxx"
#include <vector>
#include <cstring>

///Possible boundary conditions for fields
enum {
  BND_FLD_OPEN,
  BND_FLD_PERIODIC,
  BND_FLD_UPML,
  BND_FLD_TIME,
  BND_FLD_CONDUCTING_WALL,
  BND_FLD_ABSORBING,
};

///Possible boundary conditions for particles
enum {
  BND_PART_REFLECTING,
  BND_PART_PERIODIC,
  BND_PART_ABSORBING,
  BND_PART_OPEN,
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

struct GridParams
{
  using Double3 = Vec3<double>;
  
  Double3 length;	///<The physical size of the simulation-box 
  Double3 corner;
  Int3 gdims;		///<Number of grid-points in each dimension
  Int3 np;		///<Number of patches in each dimension
  Int3 bs;
  Int3 bc_fld_lo;	///<Boundary conditions of the fields. Can be any value of BND_FLD.
  Int3 bc_fld_hi;	///<Boundary conditions of the fields. Can be any value of BND_FLD.
  Int3 bc_prt_lo; 	///<Boundary conditions of the particles. Can be any value of BND_PART.
  Int3 bc_prt_hi;       ///<Boundary conditions of the particles. Can be any value of BND_PART.
};

// ======================================================================
// Grid_

template<class T>
struct Grid_
{
  using real_t = T;
  using Real3 = Vec3<real_t>;

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

  // Construct a single patch covering the whole domain
  // Args: global dimensions, and length of the domain in all 3 dims
  // -- mostly useful for testing
  //
  // Maybe the named ctor idiom would be good here (but right now
  // can't be done since the copy ctor is deleted.
  Grid_(const Int3& _gdims, const Real3& length)
    : gdims(_gdims)
  {
    ldims = gdims;
    dx = length / Real3(gdims);

    patches.emplace_back(Patch({ 0, 0, 0 }, { 0., 0., 0.}, length, dx));

    for (int d = 0; d < 3; d++) {
      assert(ldims[d] % bs[d] == 0); // FIXME, % operator for Vec3
    }
  }

  Grid_(const Int3& _gdims, const Int3& _ldims, const Real3& length,
	const Real3& corner, const std::vector<Int3>& offs)
    : gdims(_gdims),
      ldims(_ldims),
      dx(length / Real3(gdims))
  {
    for (auto off : offs) {
      patches.push_back(Patch(off,
			      Vec3<double>(off        ) * dx + corner,
			      Vec3<double>(off + ldims) * dx + corner, dx));
    }
  }
  
  int n_patches() const { return patches.size(); }

  Int3 gdims;
  Int3 ldims;
  Real3 dx;
  GridBc bc;
  // FIXME? these defaults, in particular for dt, might be a bit
  // dangerous, as they're useful for testing but might hide if one
  // forgets to init them correctly for an real simulation
  real_t fnqs = { 1. };
  real_t eta = { 1. };
  real_t dt = { 1. };
  std::vector<Patch> patches;
  std::vector<Kind> kinds;
  Int3 bs = { 1, 1, 1 };
};

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

