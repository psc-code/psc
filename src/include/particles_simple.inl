
#include "VariableByPatch.h"

// ======================================================================
// VariableByParticle

template <typename T>
struct VariableByParticle;

template <typename T>
struct VariableByParticle<std::vector<T>>
{
  using value_type = std::vector<T>;

  void put(kg::io::Engine& writer, const value_type& vec, const Grid_t& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    unsigned long n = vec.size(), N, off = 0;
    MPI_Allreduce(&n, &N, 1, MPI_UNSIGNED_LONG, MPI_SUM, grid.comm());
    MPI_Exscan(&n, &off, 1, MPI_UNSIGNED_LONG, MPI_SUM, grid.comm());
    kg::io::Dims shape = {size_t(N)};
    kg::io::Dims start = {size_t(off)};
    kg::io::Dims count = {size_t(n)};
    writer.putVariable(vec.data(), launch, shape, {start, count});
  }

#if 0
  void get(kg::io::Engine& reader, value_type& vec, const Grid_t& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    kg::io::Dims shape = {static_cast<size_t>(grid.nGlobalPatches())};
    kg::io::Dims start = {
      static_cast<size_t>(grid.localPatchInfo(0).global_patch)};
    kg::io::Dims count = {static_cast<size_t>(grid.n_patches())};
    assert(reader.variableShape<T>() == shape);
    vec.resize(count[0]);
    reader.getVariable(vec.data(), launch, {start, count});
  }
#endif
};

// ======================================================================
// Variable<MparticlesSimple>

template <typename Mparticles>
class PutComponent
{
public:
  PutComponent(kg::io::Engine& writer, const Mparticles& mprts)
    : writer_{writer}, mprts_{mprts}
  {}

  template <typename FUNC>
  void operator()(const std::string& name, FUNC&& func)
  {
    using Ret = decltype(func(mprts_[0][0]));
    std::vector<Ret> vec(mprts_.size());
    auto it = vec.begin();
    for (int p = 0; p < mprts_.n_patches(); p++) {
      auto prts = mprts_[p];
      for (int n = 0; n < prts.size(); n++) {
        *it++ = func(prts[n]);
      }
    }

    writer_.put<VariableByParticle>(name, vec, mprts_.grid(),
                                    kg::io::Mode::Blocking);
  }

private:
  kg::io::Engine& writer_;
  const Mparticles& mprts_;
};


template <typename R>
class kg::io::Descr<MparticlesSimple<R>>
{
public:
  using Mparticles = MparticlesSimple<R>;
  using Particle = typename Mparticles::Particle;

  void put(kg::io::Engine& writer, const Mparticles& mprts,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    auto& grid = mprts.grid();
    auto sizeByPatch = mprts.sizeByPatch();
    writer.put<VariableByPatch>("sizeByPatch", sizeByPatch, grid,
                                Mode::NonBlocking);

    PutComponent<Mparticles> put_component{writer, mprts};
    DoComponents<Particle>::run(put_component);
    
    writer.performPuts();
  }

  void get(kg::io::Engine& reader, Mparticles& mprts,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {}

private:
  template <typename F>
  void putComponent(kg::io::Engine& writer, const Mparticles& mprts,
                    const Grid_t& grid, const std::string& name, F&& func)
  {
    using Ret = decltype(func(mprts[0][0]));
    std::vector<Ret> vec(mprts.size());
    auto it = vec.begin();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto prts = mprts[p];
      for (int n = 0; n < prts.size(); n++) {
        *it++ = func(prts[n]);
      }
    }
    writer.put<VariableByParticle>(name, vec, grid, Mode::Blocking);
  }
};
