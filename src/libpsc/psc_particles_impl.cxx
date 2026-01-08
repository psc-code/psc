
#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "particle_with_id.h"

template <typename MP_FROM, typename MP_TO>
struct Convert
{
  using MparticlesFrom = MP_FROM;
  using MparticlesTo = MP_TO;
  using particle_from_t = typename MP_FROM::Particle;
  using particle_to_t = typename MP_TO::Particle;

  particle_to_t operator()(const particle_from_t& prt_from, const Grid_t& grid)
  {
    using real_t = typename MparticlesTo::real_t;
    using Real3 = typename MparticlesTo::Real3;

    return {Real3(prt_from.x), Real3(prt_from.u), real_t(prt_from.qni_wni),
            prt_from.kind,     prt_from.id(),     prt_from.tag()};
  }
};

template <typename MP_FROM, typename MP_TO>
void psc_mparticles_copy(MP_FROM& mp_from, MP_TO& mp_to)
{
  Convert<MP_FROM, MP_TO> convert;
  auto n_prts_by_patch = mp_from.sizeByPatch();
  mp_to.reserve_all(n_prts_by_patch);
  mp_to.clear();

  for (int p = 0; p < mp_to.n_patches(); p++) {
    auto&& prts_from = mp_from[p];
    auto&& prts_to = mp_to[p];
    int n_prts = mp_from.size(p);
    for (int n = 0; n < n_prts; n++) {
      mp_to.push_back(p, convert(prts_from[n], mp_from.grid()));
    }
  }
}

template <typename MP_FROM, typename MP_TO>
void psc_mparticles_copy_to(MparticlesBase& mp_from, MparticlesBase& mp_to)
{
  psc_mparticles_copy<MP_FROM, MP_TO>(dynamic_cast<MP_FROM&>(mp_from),
                                      dynamic_cast<MP_TO&>(mp_to));
}

template <typename MP_FROM, typename MP_TO>
void psc_mparticles_copy_from(MparticlesBase& mp_from, MparticlesBase& mp_to)
{
  psc_mparticles_copy<MP_TO, MP_FROM>(dynamic_cast<MP_TO&>(mp_to),
                                      dynamic_cast<MP_FROM&>(mp_from));
}

// ======================================================================
// psc_mparticles: subclass "single"

// ----------------------------------------------------------------------
// conversion to/from "double"

template <>
const MparticlesBase::Convert MparticlesSingle::convert_to_ = {
  {std::type_index(typeid(MparticlesDouble)),
   psc_mparticles_copy_to<MparticlesSingle, MparticlesDouble>},
};

template <>
const MparticlesBase::Convert MparticlesSingle::convert_from_ = {
  {std::type_index(typeid(MparticlesDouble)),
   psc_mparticles_copy_from<MparticlesSingle, MparticlesDouble>},
};

// ======================================================================
// psc_mparticles: subclass "double"

template <>
const MparticlesBase::Convert MparticlesDouble::convert_to_ = {};
template <>
const MparticlesBase::Convert MparticlesDouble::convert_from_ = {};

// ======================================================================
// MparticlesSimple<ParticleWithId<float>>

template <>
const MparticlesBase::Convert
  MparticlesSimple<ParticleWithId<float>>::convert_to_ = {};
template <>
const MparticlesBase::Convert
  MparticlesSimple<ParticleWithId<float>>::convert_from_ = {};

// ======================================================================
// MparticlesSimple<ParticleWithId<double>>

template <>
const MparticlesBase::Convert
  MparticlesSimple<ParticleWithId<double>>::convert_to_ = {};
template <>
const MparticlesBase::Convert
  MparticlesSimple<ParticleWithId<double>>::convert_from_ = {};
