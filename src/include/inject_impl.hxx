

#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include <bnd.hxx>
#include <fields.hxx>
#include <fields_item.hxx>
#include <inject.hxx>

#include <stdlib.h>
#include <string>

// ======================================================================
// Inject_

template <typename Target_t, typename ItemMoment, typename Mparticles>
class Functor
{
public:
  Functor(const Grid_t& grid, Target_t& target, int HE_population,
          int base_population, double HE_ratio, ItemMoment& moment_n,
          Mparticles& mprts, int interval, double tau)
    : target(target),
      HE_population(HE_population),
      base_population(base_population),
      HE_ratio(HE_ratio),
      mf_n(grid, grid.kinds.size(), grid.ibn)
  {
    moment_n.update(mprts);
    mf_n = evalMfields(moment_n);

    fac = (interval * grid.dt / tau) / (1. + interval * grid.dt / tau);
  }

  void operator()(int kind, Double3 pos, int p, Int3 idx, psc_particle_npt& npt)
  {
    if (target.is_inside(pos)) {

      if (kind == HE_population) {
        target.init_npt(base_population, pos, npt);
        npt.n -= mf_n[p](base_population, idx[0], idx[1], idx[2]);
        npt.n *= HE_ratio;
      } else {
        // should electrons inject by moments and then scale to (1-HE_ratio)?
        target.init_npt(kind, pos, npt);
        npt.n -= mf_n[p](kind, idx[0], idx[1], idx[2]);
        if (kind == base_population)
          npt.n *= (1 - HE_ratio);
      }
      if (npt.n < 0) {
        npt.n = 0;
      }
      npt.n *= fac;
    }
  }

private:
  Target_t& target;
  int HE_population;
  int base_population;
  double HE_ratio;
  MfieldsC mf_n;
  double fac;
};

template <typename _Mparticles, typename _Mfields, typename Target_t,
          typename _ItemMoment>
struct Inject_ : InjectBase
{
  using Mfields = _Mfields;
  using Mparticles = _Mparticles;
  using real_t = typename Mparticles::real_t;
  using ItemMoment_t = _ItemMoment;
  using SetupParticles = ::SetupParticles<Mparticles>;

  // ----------------------------------------------------------------------
  // ctor

  Inject_(const Grid_t& grid, int interval, int tau,
          SetupParticles& setup_particles)
    : InjectBase{interval, tau},
      setup_particles_{setup_particles}
  {}

  // ----------------------------------------------------------------------
  // operator()

  template <typename F>
  void operator()(Mparticles& mprts, F& func)
  {
    static int pr, pr_1;
    if (!pr) {
      pr_1 = prof_register("inject_setup_prts", 1., 0, 0);
    }

    prof_start(pr);
    const auto& grid = mprts.grid();

    prof_barrier("inject_barrier");

    prof_start(pr_1);
    setup_particles_.setupParticles(mprts, func);
    prof_stop(pr_1);

    prof_stop(pr);
  }

private:
  SetupParticles setup_particles_;
};

// ======================================================================
// InjectSelector
//
// FIXME, should go away eventually

template <typename Mparticles, typename InjectShape, typename Dim,
          typename Enable = void>
struct InjectSelector
{
  using Inject =
    Inject_<Mparticles, MfieldsC, InjectShape,
            Moment_n_1st<Mparticles, MfieldsC>>; // FIXME, shouldn't
                                                 // always use MfieldsC
};

#ifdef USE_CUDA

#include "../libpsc/cuda/fields_item_moments_1st_cuda.hxx"
#include <psc_fields_single.h>

// This not particularly pretty template arg specializes InjectSelector for all
// CUDA Mparticles types
template <typename Mparticles, typename InjectShape, typename Dim>
struct InjectSelector<Mparticles, InjectShape, Dim,
                      typename std::enable_if<Mparticles::is_cuda::value>::type>
{
  using Inject = Inject_<Mparticles, MfieldsSingle, InjectShape,
                         Moment_n_1st_cuda<Mparticles, Dim>>;
};

#endif
