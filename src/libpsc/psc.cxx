
#include "psc.h"
#include "psc_fields_as_c.h"
#include "fields.hxx"
#include "setup_fields.hxx"
#ifdef USE_CUDA
#include "cuda_base.hxx"
#endif

#include <mrc_common.h>
#include <mrc_params.h>
#include <mrc_io.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <array>

#if 0
// ----------------------------------------------------------------------
// _psc_write

static void
_psc_write(struct psc *psc, struct mrc_io *io)
{
  mrc_io_write_int(io, psc, "timestep", psc->timestep);
  mrc_io_write_int(io, psc, "nr_kinds", psc->nr_kinds_);

  for (int k = 0; k < psc->nr_kinds_; k++) {
    char s[20];
    sprintf(s, "kind_q%d", k);
    mrc_io_write_double(io, psc, s, psc->kinds_[k].q);
    sprintf(s, "kind_m%d", k);
    mrc_io_write_double(io, psc, s, psc->kinds_[k].m);
    mrc_io_write_string(io, psc, s, psc->kinds_[k].name);
  }
  mrc_io_write_ref(io, psc, "mrc_domain", psc->mrc_domain_);
  mrc_io_write_ref(io, psc, "mparticles", psc->particles_);
  mrc_io_write_ref(io, psc, "mfields", psc->flds);
}

// ----------------------------------------------------------------------
// _psc_read

static void
_psc_read(struct psc *psc, struct mrc_io *io)
{
  //psc_setup_coeff(psc->norm_params);

  mrc_io_read_int(io, psc, "timestep", &psc->timestep);
  mrc_io_read_int(io, psc, "nr_kinds", &psc->nr_kinds_);
  psc->kinds_ = new psc_kind[psc->nr_kinds_]();
  for (int k = 0; k < psc->nr_kinds_; k++) {
    char s[20];
    sprintf(s, "kind_q%d", k);
    mrc_io_read_double(io, psc, s, &psc->kinds_[k].q);
    sprintf(s, "kind_m%d", k);
    mrc_io_read_double(io, psc, s, &psc->kinds_[k].m);
    mrc_io_read_string(io, psc, s, &psc->kinds_[k].name);
  }

  psc->mrc_domain_ = mrc_io_read_ref(io, psc, "mrc_domain", mrc_domain);
  //psc_setup_domain(psc, psc->domain_, psc->bc_, psc->kinds_);
#ifdef USE_FORTRAN
  psc_setup_fortran(psc);
#endif

  //psc->particles_ = mrc_io_read_ref(io, psc, "mparticles", psc_mparticles);
  //psc->flds = mrc_io_read_ref(io, psc, "mfields", psc_mfields);

  //psc_read_member_objs(psc, io);
}

#endif

// FIXME
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

void psc_init(int& argc, char**& argv)
{
#if 1
  vpic_base_init(&argc, &argv);
#else
  MPI_Init(&argc, &argv);
#endif
  libmrc_params_init(argc, argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);
#ifdef USE_CUDA
  cuda_base_init();
#endif

  // FIXME, we should use RngPool consistently throughout
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srandom(rank);
}

void psc_finalize()
{
  libmrc_params_finalize();
  MPI_Finalize();
}
