
#include "particles_simple.hxx"

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

  void get(kg::io::Engine& reader, value_type& vec, const Grid_t& grid,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    unsigned long n = vec.size(), N, off = 0;
    MPI_Allreduce(&n, &N, 1, MPI_UNSIGNED_LONG, MPI_SUM, grid.comm());
    MPI_Exscan(&n, &off, 1, MPI_UNSIGNED_LONG, MPI_SUM, grid.comm());
    kg::io::Dims shape = {size_t(N)};
    kg::io::Dims start = {size_t(off)};
    kg::io::Dims count = {size_t(n)};
    reader.getVariable(vec.data(), launch, {start, count});
  }
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
    using Ret = typename std::remove_const_t<
      typename std::remove_pointer_t<decltype(func(mprts_.at(0, 0)))>>;
    std::vector<Ret> vec(mprts_.size());
    auto it = vec.begin();
    for (int p = 0; p < mprts_.n_patches(); p++) {
      for (int n = 0; n < mprts_.size(p); n++) {
        *it++ = *func(mprts_.at(p, n));
      }
    }

    writer_.put<VariableByParticle>(name, vec, mprts_.grid(),
                                    kg::io::Mode::Blocking);
  }

private:
  kg::io::Engine& writer_;
  const Mparticles& mprts_;
};

template <typename Mparticles>
class GetComponent
{
public:
  GetComponent(kg::io::Engine& reader, Mparticles& mprts)
    : reader_{reader}, mprts_{mprts}
  {}

  template <typename FUNC>
  void operator()(const std::string& name, FUNC&& func)
  {
    using Ret =
      typename std::remove_pointer<decltype(func(mprts_.at(0, 0)))>::type;
    std::vector<Ret> vec(mprts_.size());
    reader_.get<VariableByParticle>(name, vec, mprts_.grid(),
                                    kg::io::Mode::Blocking);
    auto it = vec.begin();
    for (int p = 0; p < mprts_.n_patches(); p++) {
      for (int n = 0; n < mprts_.size(p); n++) {
        *func(mprts_.at(p, n)) = *it++;
      }
    }
  }

private:
  kg::io::Engine& reader_;
  Mparticles& mprts_;
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

    auto size_by_patch = mprts.sizeByPatch();
    writer.put<VariableByPatch>("size_by_patch", size_by_patch, grid,
                                Mode::NonBlocking);

    PutComponent<Mparticles> put_component{writer, mprts};
    ForComponents<Particle>::run(put_component);

    writer.performPuts();
  }

  void get(kg::io::Engine& reader, Mparticles& mprts,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    auto& grid = mprts.grid();

    auto size_by_patch = std::vector<uint>(mprts.n_patches());
    reader.get<VariableByPatch>("size_by_patch", size_by_patch, grid,
                                Mode::Blocking);
    mprts.reserve_all(size_by_patch);
    mprts.resize_all(size_by_patch);

    GetComponent<Mparticles> get_component{reader, mprts};
    ForComponents<Particle>::run(get_component);
  }
};
