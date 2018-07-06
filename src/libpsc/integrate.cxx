
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
psc_output(struct psc *psc, PscMparticlesBase mprts)
{
  psc_method_output(psc->method, psc, mprts);
}

// This measures the time spent pushing particles and fields, exclusive of
// communication.
// Only works correctly for push_fields "variant 1"!
int pr_time_step_no_comm; // FIXME, don't like globals
