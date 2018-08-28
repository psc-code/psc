
#include "psc.h"
#include "psc_diag.h"
#include "psc_output_particles.h"
#include "psc_fields_as_c.h"
#include "fields.hxx"
#include "setup_fields.hxx"

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

using Mfields_t = MfieldsC;
using Fields = Fields3d<Mfields_t::fields_t>;

Grid_t* ggrid;

// ----------------------------------------------------------------------
// _psc_write

static void
_psc_write(struct psc *psc, struct mrc_io *io)
{
#if 0
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
#endif
}

// ----------------------------------------------------------------------
// _psc_read

static void
_psc_read(struct psc *psc, struct mrc_io *io)
{
  //psc_setup_coeff(psc->norm_params);

#if 0
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
#endif
  //psc_setup_domain(psc, psc->domain_, psc->bc_, psc->kinds_);
#ifdef USE_FORTRAN
  psc_setup_fortran(psc);
#endif

  //psc->particles_ = mrc_io_read_ref(io, psc, "mparticles", psc_mparticles);
  //psc->flds = mrc_io_read_ref(io, psc, "mfields", psc_mfields);

  //psc_read_member_objs(psc, io);
}

#ifdef USE_VPIC
// FIXME
void vpic_base_init(int *pargc, char ***pargv);
#endif

void psc_init(int& argc, char**& argv)
{
#ifdef USE_VPIC
  vpic_base_init(&argc, &argv);
#else
  MPI_Init(&argc, &argv);
#endif
  libmrc_params_init(argc, argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);
}
