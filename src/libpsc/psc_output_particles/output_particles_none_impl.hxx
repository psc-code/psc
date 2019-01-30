
#include "output_particles.hxx"

struct OutputParticlesNone : OutputParticlesParams, OutputParticlesBase
{
  OutputParticlesNone(const Grid_t& grid, const OutputParticlesParams& params) {}

  void run(MparticlesBase& mprts_base) override {}
};
