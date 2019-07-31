
#include "output_particles.hxx"

struct OutputParticlesNone : OutputParticlesParams, OutputParticlesBase
{
  OutputParticlesNone(const Grid_t& grid, const OutputParticlesParams& params) {}

  void operator()(MparticlesBase& mprts_base) {}
};
