
#include "psc_collision_impl.hxx"
#include "psc_fields_single.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_particles_double.h"
#include "fields_item.hxx"

#include <string>

// ======================================================================
// CollisionNone

struct CollisionNone : CollisionBase
{
  constexpr static const char* name = "none";

  CollisionNone(MPI_Comm comm, int interval, double nu) {}

  void operator()(MparticlesBase& mprts_base) override {}
};

void* global_collision; // FIXME

// ======================================================================
// psc_collision: subclass "single"/"double"

template<typename Collision>
struct Item_coll_stats
{
  using MfieldsState = typename Collision::MfieldsState;
  using Mfields = typename Collision::Mfields;

  constexpr static const char* name = "coll_stats";
  constexpr static int n_comps = Collision::NR_STATS;
  constexpr static fld_names_t fld_names() { return { "coll_nudt_min", "coll_nudt_med", "coll_nudt_max",
	"coll_nudt_large", "coll_ncoll" }; }
  constexpr static int flags = 0;

  static void run(MfieldsState& mflds, Mfields& mres)
  {
    assert(global_collision);
    Collision& coll = *reinterpret_cast<Collision*>(global_collision);
    
    for (int m = 0; m < coll.NR_STATS; m++) {
      // FIXME, copy could be avoided (?)
      mres.copy_comp(m, coll.mflds_stats_, m);
    }
  }
};

template<typename Collision>
struct Item_coll_rei
{
  using MfieldsState = typename Collision::MfieldsState;
  using Mfields = typename Collision::Mfields;

  constexpr static const char* name = "coll_rei";
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "coll_rei_x", "coll_rei_y", "coll_rei_z" }; }
  constexpr static int flags = 0;

  static void run(MfieldsState& mflds, Mfields& mres)
  {
    assert(global_collision);
    Collision& coll = *reinterpret_cast<Collision*>(global_collision);
    
    for (int m = 0; m < 3; m++) {
      // FIXME, copy could be avoided (?)
      mres.copy_comp(m, coll.mflds_rei_, m);
    }
  }
};
  
// ======================================================================
// psc_output_fields_item: subclass "coll_stats" / "coll_rei"

using CollisionSingle = Collision_<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>;
FieldsItemOps<FieldsItemFields<Item_coll_stats<CollisionSingle>>> psc_output_fields_item_coll_stats_single_ops;
FieldsItemOps<FieldsItemFields<Item_coll_rei<CollisionSingle>>> psc_output_fields_item_coll_rei_single_ops;

using CollisionDouble = Collision_<MparticlesDouble, MfieldsStateDouble, MfieldsC>;
FieldsItemOps<FieldsItemFields<Item_coll_stats<CollisionDouble>>> psc_output_fields_item_coll_stats_double_ops;
FieldsItemOps<FieldsItemFields<Item_coll_rei<CollisionDouble>>> psc_output_fields_item_coll_rei_double_ops;

