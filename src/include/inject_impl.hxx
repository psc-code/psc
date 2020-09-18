

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
  using ItemMoment_t = _ItemMoment;
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
