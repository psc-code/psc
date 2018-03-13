
#include "psc.h"

#include "psc_particles_as_double.h" // only for checking...

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================

#include "particles.hxx"
#include "particles_impl.hxx"

#include "psc_particles_single.h"
#include "psc_particles_cuda.h"
#include "psc_particles_vpic.h"

template PscMparticlesSingle PscMparticles<MparticlesBase>::get_as(uint flags);
template PscMparticlesCuda PscMparticles<MparticlesBase>::get_as(uint flags);
template PscMparticlesVpic PscMparticles<MparticlesBase>::get_as(uint flags);
template void PscMparticles<MparticlesSingle>::put_as(PscMparticlesBase mprts_base, uint flags);
template void PscMparticles<MparticlesCuda>::put_as(PscMparticlesBase mprts_base, uint flags);
template void PscMparticles<MparticlesVpic>::put_as(PscMparticlesBase mprts_base, uint flags);

// ======================================================================
// psc_mparticles base class

static void
_psc_mparticles_view(struct psc_mparticles *_mprts)
{
  MPI_Comm comm = psc_mparticles_comm(_mprts);
  PscMparticlesBase mprts(_mprts);
  mpi_printf(comm, "  n_patches    = %d\n", mprts->n_patches());
  mpi_printf(comm, "  n_prts_total = %d\n", mprts->get_n_prts());

  uint n_prts_by_patch[mprts->n_patches()];
  mprts->get_size_all(n_prts_by_patch);

  for (int p = 0; p < mprts->n_patches(); p++) {
    mpi_printf(comm, "  p %d: n_prts = %d\n", p, n_prts_by_patch[p]);
  }  
}

void
psc_mparticles_check(struct psc_mparticles *_mprts_base)
{
  auto mprts_base = PscMparticlesBase{_mprts_base};
  int fail_cnt = 0;

  mparticles_t mprts = mprts_base.get_as<PscMparticlesDouble>();
  const Grid_t& grid = ppsc->grid();
  
  psc_foreach_patch(ppsc, p) {
    auto& patch = ppsc->grid().patches[p];
    mparticles_t::patch_t& prts = mparticles_t(mprts)[p];

    f_real xb[3], xe[3];
    
    // New-style boundary requirements.
    // These will need revisiting when it comes to non-periodic domains.
    
    for (int d = 0; d < 3; d++) {
      xb[d] = patch.xb[d];
      xe[d] = patch.xb[d] + grid.ldims[d] * grid.dx[d];
    }

    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      particle_t *part = &*prt_iter;
      if (part->xi < 0.f || part->xi >= xe[0] - xb[0] || // FIXME xz only!
	  part->zi < 0.f || part->zi >= xe[2] - xb[2]) {
	if (fail_cnt++ < 10) {
	  mprintf("FAIL: xi %g [%g:%g]\n", part->xi, 0., xe[0] - xb[0]);
	  mprintf("      zi %g [%g:%g]\n", part->zi, 0., xe[2] - xb[2]);
	}
      }
    }
  }
  assert(fail_cnt == 0);

  mprts.put_as(mprts_base, 0);
}

// ======================================================================

extern struct psc_mparticles_ops psc_mparticles_single_ops;
extern struct psc_mparticles_ops psc_mparticles_double_ops;
extern struct psc_mparticles_ops psc_mparticles_sse2_ops;
extern struct psc_mparticles_ops psc_mparticles_cbe_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda2_ops;
extern struct psc_mparticles_ops psc_mparticles_acc_ops;
extern struct psc_mparticles_ops psc_mparticles_vpic_ops;
extern struct psc_mparticles_ops psc_mparticles_single_by_kind_ops;

static void
psc_mparticles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_single_by_kind_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_cuda_ops);
#endif
#ifdef USE_CUDA2
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_cuda2_ops);
#endif
#ifdef USE_ACC
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_acc_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_vpic_ops);
#endif
}

struct mrc_class_psc_mparticles_ : mrc_class_psc_mparticles {
  mrc_class_psc_mparticles_() {
    name             = "psc_mparticles";
    size             = sizeof(struct psc_mparticles);
    init             = psc_mparticles_init;
    view             = _psc_mparticles_view;
  }
} mrc_class_psc_mparticles;



