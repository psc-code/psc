
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

  template <typename Mparticles>
  void operator()(Mparticles& mprts_base) {}
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
  static std::vector<std::string> fld_names() {
    return { "coll_nudt_min", "coll_nudt_med", "coll_nudt_max",
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
  static std::vector<std::string> fld_names() {
    return { "coll_rei_x", "coll_rei_y", "coll_rei_z" };
  }
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
  
