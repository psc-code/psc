
#include <array>

#include "VariableByPatch.h"

// ======================================================================
// Descr<Grid_t::Domain>

// FIXME, this should be templated by Grid_<T>::Domain, but can't do that...

template <>
class kg::io::Descr<Grid_t::Domain>
{
public:
  using value_type = typename Grid_t::Domain;

  static void put(kg::io::Engine& writer, const value_type& domain,
                  const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    writer.put("gdims", domain.gdims, launch);
    writer.put("length", domain.length, launch);
    writer.put("corner", domain.corner, launch);
    writer.put("np", domain.np, launch);
    writer.put("ldims", domain.ldims, launch);
    writer.put("dx", domain.dx, launch);
  }

  static void get(Engine& reader, value_type& domain,
                  const Mode launch = Mode::NonBlocking)
  {
    reader.get("gdims", domain.gdims, launch);
    reader.get("length", domain.length, launch);
    reader.get("corner", domain.corner, launch);
    reader.get("np", domain.np, launch);
    reader.get("ldims", domain.ldims, launch);
    reader.get("dx", domain.dx, launch);
  }
};

// ======================================================================
// Descr<psc::grid::BC>

template <>
class kg::io::Descr<psc::grid::BC>
{
public:
  using value_type = psc::grid::BC;

  static void put(kg::io::Engine& writer, const value_type& bc,
                  const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    writer.put("fld_lo", bc.fld_lo, launch);
    writer.put("fld_hi", bc.fld_hi, launch);
    writer.put("prt_lo", bc.prt_lo, launch);
    writer.put("prt_hi", bc.prt_hi, launch);
  }

  static void get(Engine& reader, value_type& bc,
                  const Mode launch = Mode::NonBlocking)
  {
    reader.get("fld_lo", bc.fld_lo, launch);
    reader.get("fld_hi", bc.fld_hi, launch);
    reader.get("prt_lo", bc.prt_lo, launch);
    reader.get("prt_hi", bc.prt_hi, launch);
  }
};

// ======================================================================
// Descr<Normalization>

template <>
class kg::io::Descr<Grid_t::Normalization>
{
public:
  using value_type = Grid_t::Normalization;

  static void put(kg::io::Engine& writer, const Grid_t::Normalization& norm,
                  const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    writer.put("cc", norm.cc, launch);
    writer.put("fnqs", norm.fnqs, launch);
    writer.put("eta", norm.eta, launch);
    writer.put("beta", norm.beta, launch);
    writer.put("cori", norm.cori, launch);
    writer.put("b0", norm.b0, launch);
    writer.put("rho0", norm.rho0, launch);
    writer.put("phi0", norm.phi0, launch);
    writer.put("a0", norm.a0, launch);
  }

  static void get(Engine& reader, Grid_t::Normalization& norm,
                  const Mode launch = Mode::NonBlocking)
  {
    reader.get("cc", norm.cc, launch);
    reader.get("fnqs", norm.fnqs, launch);
    reader.get("eta", norm.eta, launch);
    reader.get("beta", norm.beta, launch);
    reader.get("cori", norm.cori, launch);
    reader.get("b0", norm.b0, launch);
    reader.get("rho0", norm.rho0, launch);
    reader.get("phi0", norm.phi0, launch);
    reader.get("a0", norm.a0, launch);
  }
};

// ======================================================================
// Descr<Grid::Kinds>

template <>
class kg::io::Descr<Grid_t::Kinds>
{
  using real_t = Grid_t::real_t;

public:
  using value_type = Grid_t::Kinds;

  static void put(kg::io::Engine& writer, const Grid_t::Kinds& kinds,
                  const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    auto n_kinds = kinds.size();
    auto names = std::vector<std::string>(n_kinds);
    auto q = std::vector<real_t>(n_kinds);
    auto m = std::vector<real_t>(n_kinds);
    for (int kind = 0; kind < n_kinds; kind++) {
      q[kind] = kinds[kind].q;
      m[kind] = kinds[kind].m;
      names[kind] = kinds[kind].name;
    }

    writer.put("names", names, kg::io::Mode::Blocking);
    writer.put("q", q, kg::io::Mode::Blocking);
    writer.put("m", m, kg::io::Mode::Blocking);
  }

  static void get(Engine& reader, Grid_t::Kinds& kinds,
                  const Mode launch = Mode::NonBlocking)
  {
    auto q = std::vector<real_t>{};
    auto m = std::vector<real_t>{};
    auto names = std::vector<std::string>{};
    reader.get("names", names, kg::io::Mode::Blocking);
    reader.get("q", q, kg::io::Mode::Blocking);
    reader.get("m", m, kg::io::Mode::Blocking);
    reader.get("names", names);

    kinds.resize(q.size());
    for (int kind = 0; kind < q.size(); kind++) {
      kinds[kind].q = q[kind];
      kinds[kind].m = m[kind];
      kinds[kind].name = strdup(names[kind].c_str());
    }
  }
};

// ======================================================================
// DescrGrid_<T>>

template <typename T>
class kg::io::Descr<Grid_<T>>
{
  using Grid = Grid_<T>;
  using real_t = typename Grid::real_t;
  using Real3 = typename Grid::Real3;

public:
  using value_type = Grid;

  void put(kg::io::Engine& writer, const Grid& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    writer.put("ldims", grid.ldims);
    writer.put("domain", grid.domain, launch);
    writer.put("bc", grid.bc, launch);
    writer.put("norm", grid.norm, launch);
    writer.put("dt", grid.dt);

    size_t patches_n_local = grid.patches.size();
    writer.putLocal("n_local", patches_n_local);

    auto patches_off = std::vector<Int3>(patches_n_local);
    auto patches_xb = std::vector<Real3>(patches_n_local);
    auto patches_xe = std::vector<Real3>(patches_n_local);
    for (int p = 0; p < patches_n_local; p++) {
      auto& patch = grid.patches[p];
      patches_off[p] = patch.off;
      patches_xb[p] = patch.xb;
      patches_xe[p] = patch.xe;
    }

    writer.put<VariableByPatch>("off", patches_off, grid, launch);
    writer.put<VariableByPatch>("xb", patches_xb, grid, launch);
    writer.put<VariableByPatch>("xe", patches_xe, grid, launch);

    writer.put("kinds", grid.kinds, launch);
    writer.put("ibn", grid.ibn);
    writer.put("timestep", grid.timestep_);

    writer.performPuts(); // because we're writing temp local vars (the
                          // patches_*)
  }

  void get(kg::io::Engine& reader, Grid& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    reader.get("ldims", grid.ldims);
    reader.get("domain", grid.domain, launch);
    reader.get("bc", grid.bc, launch);
    reader.get("norm", grid.norm, launch);
    reader.get("dt", grid.dt);

    size_t patches_n_local;
    reader.getLocal("n_local", patches_n_local, launch);

    reader.performGets(); // need patches_n_local, domain, bc to be read
    grid.mrc_domain_ =
      MrcDomain(grid.domain, grid.bc, patches_n_local);

    grid.patches.resize(patches_n_local);
    auto patches_off = std::vector<Int3>(patches_n_local);
    auto patches_xb = std::vector<Real3>(patches_n_local);
    auto patches_xe = std::vector<Real3>(patches_n_local);
    reader.get<VariableByPatch>("off", patches_off, grid, launch);
    reader.get<VariableByPatch>("xb", patches_xb, grid, launch);
    reader.get<VariableByPatch>("xe", patches_xe, grid, launch);

    reader.performGets(); // need to actually read the temp local vars
    for (int p = 0; p < patches_n_local; p++) {
      auto& patch = grid.patches[p];
      patch.off = patches_off[p];
      patch.xb = patches_xb[p];
      patch.xe = patches_xe[p];
    }

    reader.get("kinds", grid.kinds, launch);
    reader.get("ibn", grid.ibn);
    reader.get("timestep", grid.timestep_);

    grid.reset_ddc();
  }
};
