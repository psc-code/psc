
#ifndef GRID_HXX
#define GRID_HXX

#include "vec3.hxx"
#include <vector>
#include <cstring>

// ======================================================================
// Grid_

template<class T>
struct Grid_
{
  using real_t = T;
  using Real3 = Vec3<real_t>;
  
  struct Patch
  {
    // FIXME, default ctor should go away
    Patch() = default;
    
    Patch(const Real3& _xb, const Real3& _xe)
      : xb(_xb), xe(_xe)
    {}
    
    Real3 xb;
    Real3 xe;
  };

  struct Kind
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

    ~Kind()
    {
      free((void*) name);
    }
    
    real_t q;
    real_t m;
    const char *name;
  };

  // FIXME, default constructor should maybe go away, since it doesn't
  // guarantee a consistent state
  Grid_() = default;

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

    patches.emplace_back(Patch({ 0., 0., 0.}, length));

    for (int d = 0; d < 3; d++) {
      assert(ldims[d] % bs[d] == 0); // FIXME, % operator for Vec3
    }
  }
  
  Grid_(const Grid_& grid) = delete;

  Int3 gdims;
  Int3 ldims;
  Real3 dx;
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

using Grid_t = Grid_<double>;

#endif

