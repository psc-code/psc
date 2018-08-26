
#include "psc.h"

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================

#include "particles.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "../libpsc/cuda/psc_particles_cuda.h"
#ifdef HAVE_VPIC
#include "psc_particles_vpic.h"
#endif

void MparticlesBase::convert(MparticlesBase& mp_from, MparticlesBase& mp_to)
{
  // FIXME, implementing == wouldn't hurt
  assert(&mp_from.grid() == &mp_to.grid());
  
  auto convert_to = mp_from.convert_to().find(std::type_index(typeid(mp_to)));
  if (convert_to != mp_from.convert_to().cend()) {
    convert_to->second(mp_from, mp_to);
    return;
  }
  
  auto convert_from = mp_to.convert_from().find(std::type_index(typeid(mp_from)));
  if (convert_from != mp_to.convert_from().cend()) {
    convert_from->second(mp_to, mp_from);
    return;
  }

  fprintf(stderr, "ERROR: no conversion known from %s to %s!\n",
	  typeid(mp_from).name(), typeid(mp_to).name());
  assert(0);
}

// ======================================================================
// psc_mparticles base class

void
psc_mparticles_check(MparticlesBase& mprts_base)
{
  int fail_cnt = 0;

  auto& mprts = mprts_base.get_as<MparticlesDouble>();
  const auto& grid = mprts.grid();

  for (int p = 0; p < grid.n_patches(); p++) {
    auto& patch = grid.patches[p];
    auto& prts = mprts[p];

    double xb[3], xe[3];
    
    // New-style boundary requirements.
    // These will need revisiting when it comes to non-periodic domains.
    
    for (int d = 0; d < 3; d++) {
      xb[d] = patch.xb[d];
      xe[d] = patch.xb[d] + grid.ldims[d] * grid.domain.dx[d];
    }

    for (auto prt : prts) {
      if (prt.x[0] < 0.f || prt.x[0] >= xe[0] - xb[0] || // FIXME xz only!
	  prt.x[1] < 0.f || prt.x[2] >= xe[2] - xb[2]) {
	if (fail_cnt++ < 10) {
	  mprintf("FAIL: xi %g [%g:%g]\n", prt.x[0], 0., xe[0] - xb[0]);
	  mprintf("      zi %g [%g:%g]\n", prt.x[2], 0., xe[2] - xb[2]);
	}
      }
    }
  }
  assert(fail_cnt == 0);

  mprts_base.put_as(mprts, MP_DONT_COPY);
}

