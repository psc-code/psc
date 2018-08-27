
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
// psc_setup_domain

Grid_t* psc_setup_domain(const Grid_t::Domain& domain, GridBc& bc, const Grid_t::Kinds& kinds,
			 const Grid_t::Normalization& norm, double dt, Int3 ibn)
{
#if 0
  mpi_printf(MPI_COMM_WORLD, "::: dt      = %g\n", dt);
  mpi_printf(MPI_COMM_WORLD, "::: dx      = %g %g %g\n", domain.dx[0], domain.dx[1], domain.dx[2]);
#endif

  assert(domain.dx[0] > 0.);
  assert(domain.dx[1] > 0.);
  assert(domain.dx[2] > 0.);
  
  for (int d = 0; d < 3; d++) {
    if (ibn[d] != 0) {
      continue;
    }
    // FIXME, old-style particle pushers need 3 ghost points still
    if (domain.gdims[d] == 1) {
      // no ghost points
      ibn[d] = 0;
    } else {
      ibn[d] = 2;
    }
  }

  ggrid = new Grid_t{domain, bc, kinds, norm, dt};
  ggrid->ibn = ibn;

  return ggrid;
}

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

