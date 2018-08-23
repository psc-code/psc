
#include "psc.h"
#include "balance.hxx"
#include "particles.hxx"
#include "push_particles.hxx"
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

// This measures the time spent pushing particles and fields, exclusive of
// communication.
// Only works correctly for push_fields "variant 1"!
int pr_time_step_no_comm; // FIXME, don't like globals
