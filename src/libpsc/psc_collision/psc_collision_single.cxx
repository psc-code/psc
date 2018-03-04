
#include <psc_fields_as_single.h>
#include <psc_particles_as_single.h>

#define psc_output_fields_item_coll_stats_ops psc_output_fields_item_coll_stats_single_ops
#define psc_output_fields_item_coll_rei_ops psc_output_fields_item_coll_rei_single_ops

#include "psc_collision_common.cxx"

// ======================================================================
// psc_collision: subclass "single"

psc_collision_ops_<Collision_<mparticles_t>> psc_collision_single_ops;


