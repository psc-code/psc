
#include "psc.h"
#include "psc_method.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_marder.h"
#include "psc_bnd_particles.h"
#include "psc_collision.h"
#include "psc_sort.h"
#include "psc_event_generator.h"
#include "psc_balance.h"
#include "psc_checks.h"
#include "balance.hxx"
#include "particles.hxx"
#include "push_particles.hxx"
#include "push_fields.hxx"
#include "sort.hxx"
#include "collision.hxx"
#include "bnd_particles.hxx"
#include "checks.hxx"
#include "marder.hxx"

#include <mrc_common.h>
#include <mrc_profile.h>

int st_time_output;
int st_time_comm;
int st_time_particle;
int st_time_field;

#define psc_ops(psc) ((struct psc_ops *)((psc)->obj.ops))

// ----------------------------------------------------------------------
// psc_print_profiling

void
psc_print_profiling(struct psc *psc)
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (1||(size > 1 && !psc->prm.detailed_profiling)) {
    prof_print_mpi(MPI_COMM_WORLD);
  } else {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (int i = 0; i < size; i++) {
      if (i == rank) {
	mprintf("profile\n");
	prof_print();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

// ----------------------------------------------------------------------
// psc_output

void
psc_output(struct psc *psc)
{
  psc_method_output(psc->method, psc);
}

/////////////////////////////////////////////////////////////////////////
/// psc_step
///

// This measures the time spent pushing particles and fields, exclusive of
// communication.
// Only works correctly for push_fields "variant 1"!
int pr_time_step_no_comm; // FIXME, don't like globals

void
psc_step(struct psc *psc)
{
  if (!pr_time_step_no_comm) {
    pr_time_step_no_comm = prof_register("time step w/o comm", 1., 0, 0);
  }

  // default psc_step() implementation

#if 0
  mpi_printf(psc_comm(psc), "**** Step %d / %d, Time %g\n", psc->timestep + 1,
	     psc->prm.nmax, psc->timestep * psc->dt);
#endif

  // x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}

  PscMparticlesBase mprts(psc->particles);
  PscMfieldsBase mflds(psc->flds);
  PscPushParticlesBase pushp(psc->push_particles);
  PscPushFieldsBase pushf(psc->push_fields);
  PscSortBase sort(psc->sort);
  PscCollisionBase collision(psc->collision);
  PscBndParticlesBase bndp(psc->bnd_particles);

  auto balance = PscBalanceBase{psc->balance};
  balance(psc, mprts);

  prof_start(pr_time_step_no_comm);
  prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

  sort(mprts);
  collision(mprts);
  
  //psc_bnd_particles_open_calc_moments(psc->bnd_particles, psc->particles);

  PscChecksBase{psc->checks}.continuity_before_particle_push(psc);

  // particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
  pushp(mprts, mflds);
  // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}
    
  // field propagation B^{n+1/2} -> B^{n+1}
  pushf.advance_H(mflds, .5);
  // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

  bndp(*mprts.sub());
  
  psc_event_generator_run(psc->event_generator, psc->particles, psc->flds);
  
  // field propagation E^{n+1/2} -> E^{n+3/2}
  pushf.advance_b2(mflds);
  // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}

  // field propagation B^{n+1} -> B^{n+3/2}
  pushf.advance_a(mflds);
  // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}

  PscChecksBase{psc->checks}.continuity_after_particle_push(psc);

  // E at t^{n+3/2}, particles at t^{n+3/2}
  // B at t^{n+3/2} (Note: that is not it's natural time,
  // but div B should be == 0 at any time...)
  PscMarderBase{psc->marder}(mflds, mprts);
    
  PscChecksBase{psc->checks}.gauss(psc);

  psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
}

