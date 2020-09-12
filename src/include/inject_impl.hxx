

#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
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
#ifdef USE_CUDA
  using MFields = HMFields;  
#else
  using MFields = MFieldsC;
#endif

  // ----------------------------------------------------------------------
  // ctor
  template <typename FUNC>
  Inject_( Grid_t& grid, int interval,
          SetupParticles& setup_particles,
          FUNC init_npt)
          //std::function<void (int, int, psc_particle_npt&, MFields)> init_npt)
    : InjectBase{interval},
      moment_n_{grid},
      setup_particles_{setup_particles},
      init_npt_{init_npt}
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
    auto mf_n = evalMfields(moment_n_);
    prof_stop(pr_2);
    
    auto lf_init_npt = [&](int kind, Double3 pos, int p, Int3 idx,
                           psc_particle_npt& npt) {

        init_npt_(kind, pos, p, idx, npt, mf_n);
    };
    prof_start(pr_4);
    setup_particles_.setupParticles(mprts, lf_init_npt);
    prof_stop(pr_4);
    
    prof_stop(pr);
  }

private:
  ItemMoment_t moment_n_;
  SetupParticles setup_particles_;
  std::function<void (int, Double3, int, Int3,  psc_particle_npt&, MFields)> init_npt_;
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
