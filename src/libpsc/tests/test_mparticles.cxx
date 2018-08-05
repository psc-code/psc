
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc_particles_single.h"

TEST(Mparticles, Constructor)
{
  using Mparticles_t = Mparticles<particle_single_t>;
  
  Grid_t grid = make_grid();
  grid.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));

  Mparticles_t mprts(grid);
}

