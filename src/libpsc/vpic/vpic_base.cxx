
#include "vpic_config.h"
#include "vpic_iface.h"

#include <mrc_common.h>
#include <cassert>

vpic_simulation* simulation; // FIXME, never initialized, just a hack to make
                             // things link at least

// ----------------------------------------------------------------------
// vpic_base_init

void vpic_base_init(int* pargc, char*** pargv)
{
  static bool vpic_base_inited = false;

  if (vpic_base_inited) {
    return;
  }
  vpic_base_inited = true;

  //  boot_services( &argc, &argv );
  {
    // Start up the checkpointing service.  This should be first.
#ifdef USE_VPIC
    boot_checkpt(pargc, pargv);

    serial.boot(pargc, pargv);
    thread.boot(pargc, pargv);

    // Boot up the communications layer
    // See note above about thread-core-affinity

    boot_mp(pargc, pargv);
#else
    // MPI_Init(pargc, pargv);
    int provided;
    MPI_Init_thread(pargc, pargv, MPI_THREAD_MULTIPLE, &provided);
#endif

    MPI_Comm_dup(MPI_COMM_WORLD, &psc_comm_world);
    MPI_Comm_rank(psc_comm_world, &psc_world_rank);
    MPI_Comm_size(psc_comm_world, &psc_world_size);

    MPI_Barrier(psc_comm_world);
#ifdef USE_VPIC
    _boot_timestamp = 0;
    _boot_timestamp = uptime();
#endif
  }
  LOG_INFO("vpic_base_init() done\n");
}

// ======================================================================

void Simulation_set_region_resistive_harris(vpic_harris_params* prm,
                                            globals_physics* phys, double dx[3],
                                            double thickness,
                                            struct material* resistive)
{
  // Define resistive layer surrounding boundary --> set thickness=0
  // to eliminate this feature
#define resistive_layer                                                        \
  ((prm->open_bc_x && x < dx[0] * thickness) ||                                \
   (prm->open_bc_x && x > phys->Lx - dx[0] * thickness) ||                     \
   z < -phys->Lz / 2 + dx[2] * thickness ||                                    \
   z > phys->Lz / 2 - dx[2] * thickness)

  if (thickness > 0) {
    mprintf("Setting resistive layer of thickness %g", thickness);
    // FIXME!!!
    assert(0);
#if 0
#define field vpic->field
    grid_t *grid = sim->grid_->g_;
    set_region_material(resistive_layer, resistive, resistive);
#undef field
#endif
  }
}
