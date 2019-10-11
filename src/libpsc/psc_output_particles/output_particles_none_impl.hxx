
#include "output_particles.hxx"

struct OutputParticlesNone : OutputParticlesParams, OutputParticlesBase
{
  OutputParticlesNone(const Grid_t& grid, const OutputParticlesParams& params) {}

  template <typename Mparticles>
  void operator()(Mparticles& mprts_base) {}
};
