
#include <psc_fields_as_c.h>
#include <psc_particles_as_double.h>

#define psc_output_fields_item_coll_stats_ops psc_output_fields_item_coll_stats_double_ops
#define psc_output_fields_item_coll_rei_ops psc_output_fields_item_coll_rei_double_ops

#include "psc_collision_common.cxx"

// ======================================================================
// psc_collision: subclass "double"

psc_collision_ops_<Collision_<mparticles_t>> psc_collision_double_ops;

