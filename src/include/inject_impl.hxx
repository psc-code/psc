

#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include <bnd.hxx>
#include <fields.hxx>
#include <fields_item.hxx>
#include <inject.hxx>

#include <stdlib.h>
#include <string>

// ======================================================================
// Inject_

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

  Inject_(const Grid_t& grid, int interval, int tau, int kind_n,
          Target_t target, SetupParticles& setup_particles)
    : InjectBase{interval, tau, kind_n},
      target_{target},
      moment_n_{grid},
      setup_particles_{setup_particles}
  {}

  // ----------------------------------------------------------------------
  // operator()

  void operator()(Mparticles& mprts)
  {
    static int pr, pr_1, pr_2, pr_3, pr_4;
    if (!pr) {
      pr = prof_register("inject_impl", 1., 0, 0);
      pr_1 = prof_register("inject_moment_n", 1., 0, 0);
      pr_2 = prof_register("inject_eval", 1., 0, 0);
      pr_3 = prof_register("inject_get_as", 1., 0, 0);
      pr_4 = prof_register("inject_setup_prts", 1., 0, 0);
    }

    prof_start(pr);
    const auto& grid = mprts.grid();

    prof_barrier("inject_barrier");

    prof_start(pr_1);
    moment_n_.update(mprts);
    prof_stop(pr_1);
    
    prof_start(pr_2);
    auto mres = evalMfields(moment_n_);
    prof_stop(pr_2);
    
    prof_start(pr_3);
    auto& mf_n = mres.template get_as<Mfields>(kind_n, kind_n + 1);
    prof_stop(pr_3);
    
    real_t fac = (interval * grid.dt / tau) / (1. + interval * grid.dt / tau);

    auto lf_init_npt = [&](int kind, Double3 pos, int p, Int3 idx,
                           psc_particle_npt& npt) {
      if (target_.is_inside(pos)) {
        target_.init_npt(kind, pos, npt);
        npt.n -= mf_n[p](kind_n, idx[0], idx[1], idx[2]);
        if (npt.n < 0) {
          npt.n = 0;
        }
        npt.n *= fac;
      }
    };

    prof_start(pr_4);
    setup_particles_.setupParticles(mprts, lf_init_npt);
    prof_stop(pr_4);
    
    mres.put_as(mf_n, 0, 0);
    prof_stop(pr);
  }

private:
  Target_t target_;
  ItemMoment_t moment_n_;
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
