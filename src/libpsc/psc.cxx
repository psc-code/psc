
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

struct psc *ppsc;

#define VAR(x) (void *)(offsetof(struct psc, x))
#define VAR_(x, n) (void *)(offsetof(struct psc, x) + n*sizeof(int))

static mrc_param_select bnd_fld_descr[3];
static mrc_param_select bnd_prt_descr[4];

static struct select_init {
  select_init() {
    bnd_fld_descr[0].str = "open";
    bnd_fld_descr[0].val = BND_FLD_OPEN;
    bnd_fld_descr[1].str = "periodic";
    bnd_fld_descr[1].val = BND_FLD_PERIODIC;
    bnd_fld_descr[2].str = "conducting_wall";
    bnd_fld_descr[2].val = BND_FLD_CONDUCTING_WALL;

    bnd_prt_descr[0].str = "reflecting";
    bnd_prt_descr[0].val = BND_PRT_REFLECTING;
    bnd_prt_descr[1].str = "periodic";
    bnd_prt_descr[1].val = BND_PRT_PERIODIC;
    bnd_prt_descr[2].str = "absorbing";
    bnd_prt_descr[2].val = BND_PRT_ABSORBING;
    bnd_prt_descr[3].str = "open";
    bnd_prt_descr[3].val = BND_PRT_OPEN;
  }
} select_initializer;

static struct param psc_descr[] = {
  { "diag"                    , VAR(diag)                    , MRC_VAR_OBJ(psc_diag) },
  { "output_particles"        , VAR(output_particles)        , MRC_VAR_OBJ(psc_output_particles) },

  {},
};

#undef VAR

// ----------------------------------------------------------------------
// psc_create

static void
_psc_create(struct psc *psc)
{
  assert(!ppsc);
  ppsc = psc;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srandom(rank);
}

// ======================================================================
// psc_setup

// ----------------------------------------------------------------------
// psc_setup_mrc_domain

struct mrc_domain *
psc_setup_mrc_domain(const Grid_t::Domain& grid_domain, const GridBc& grid_bc, int nr_patches)
{
  // FIXME, should be split to create, set_from_options, setup time?
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  // create a very simple domain decomposition
  int bc[3] = {};
  for (int d = 0; d < 3; d++) {
    if (grid_bc.fld_lo[d] == BND_FLD_PERIODIC && grid_domain.gdims[d] > 1) {
      bc[d] = BC_PERIODIC;
    }
  }

  mrc_domain_set_type(domain, "multi");
  mrc_domain_set_param_int3(domain, "m", grid_domain.gdims);
  mrc_domain_set_param_int(domain, "bcx", bc[0]);
  mrc_domain_set_param_int(domain, "bcy", bc[1]);
  mrc_domain_set_param_int(domain, "bcz", bc[2]);
  mrc_domain_set_param_int(domain, "nr_patches", nr_patches);
  mrc_domain_set_param_int3(domain, "np", grid_domain.np);

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", 2);
  mrc_crds_set_param_double3(crds, "l", grid_domain.corner);
  mrc_crds_set_param_double3(crds, "h", grid_domain.corner + grid_domain.length);

  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);

  return domain;
}

// ----------------------------------------------------------------------
// psc_make_grid

Grid_t* psc::make_grid(const MrcDomain& mrc_domain, const Grid_t::Domain& domain, const GridBc& bc,
		       const Grid_t::Kinds& kinds, Grid_t::Normalization coeff, double dt)
{
  int n_patches;
  auto patches = mrc_domain.getPatches(&n_patches);
  assert(n_patches > 0);
  Int3 ldims = patches[0].ldims;
  std::vector<Int3> offs;
  for (int p = 0; p < n_patches; p++) {
    assert(ldims == Int3(patches[p].ldims));
    offs.push_back(patches[p].off);
  }

  Grid_t *grid = new Grid_t(domain, offs);

  grid->kinds = kinds;

  grid->bc = bc;
  for (int d = 0; d < 3; d++) {
    if (grid->isInvar(d)) {
      // if invariant in this direction: set bnd to periodic (FIXME?)
      grid->bc.fld_lo[d] = BND_FLD_PERIODIC;
      grid->bc.fld_hi[d] = BND_FLD_PERIODIC;
      grid->bc.prt_lo[d] = BND_PRT_PERIODIC;
      grid->bc.prt_hi[d] = BND_PRT_PERIODIC;
    }
  }
  
  grid->norm = coeff;
  grid->dt = dt;

  grid->mrc_domain_ = &mrc_domain;

  return grid;
}

// ----------------------------------------------------------------------
// psc_setup_domain

Grid_t* psc_setup_domain(struct psc *psc, const Grid_t::Domain& domain, GridBc& bc, const Grid_t::Kinds& kinds,
			 const Grid_t::Normalization& norm, double dt)
{
#if 0
  mpi_printf(MPI_COMM_WORLD, "::: dt      = %g\n", dt);
  mpi_printf(MPI_COMM_WORLD, "::: dx      = %g %g %g\n", domain.dx[0], domain.dx[1], domain.dx[2]);
#endif

  assert(domain.dx[0] > 0.);
  assert(domain.dx[1] > 0.);
  assert(domain.dx[2] > 0.);
  
  for (int d = 0; d < 3; d++) {
    if (psc->ibn[d] != 0) {
      continue;
    }
    // FIXME, old-style particle pushers need 3 ghost points still
    if (domain.gdims[d] == 1) {
      // no ghost points
      psc->ibn[d] = 0;
    } else {
      psc->ibn[d] = 2;
    }
  }

  psc->mrc_domain_.reset(psc_setup_mrc_domain(domain, bc, -1));
  psc->grid_ = psc->make_grid(psc->mrc_domain_, domain, bc, kinds, norm, dt);

  // make sure that np isn't overridden on the command line
  Int3 np;
  psc->mrc_domain_.get_param_int3("np", np);
  assert(np == domain.np);

  return psc->grid_;
}

// ----------------------------------------------------------------------
// psc_destroy

static void
_psc_destroy(struct psc *psc)
{
  ppsc = NULL;
}

// ----------------------------------------------------------------------
// _psc_write

static void
_psc_write(struct psc *psc, struct mrc_io *io)
{
  mrc_io_write_int(io, psc, "timestep", psc->timestep);
#if 0
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
  assert(!ppsc);
  ppsc = psc;

  //psc_setup_coeff(psc->norm_params);

  mrc_io_read_int(io, psc, "timestep", &psc->timestep);
#if 0
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

  psc_read_member_objs(psc, io);
}

// ----------------------------------------------------------------------
// _psc_view

static void
_psc_view(struct psc *psc)
{
  psc->mrc_domain_.view();

  MPI_Comm comm = psc_comm(psc);
  const auto& kinds = psc->grid().kinds;
  mpi_printf(comm, "%20s|\n", "particle kinds");
  for (int k = 0; k < kinds.size(); k++) {
    mpi_printf(comm, "%19s | q = %g m = %g\n",
	       kinds[k].name, kinds[k].q, kinds[k].m);
  }
}

// ======================================================================
// psc class

struct mrc_class_psc_ : mrc_class_psc {
  mrc_class_psc_() {
    name             = "psc";
    size             = sizeof(struct psc);
    param_descr      = psc_descr;
    create           = _psc_create;
    view             = _psc_view;
    destroy          = _psc_destroy;
    write            = _psc_write;
    read             = _psc_read;
  }
} mrc_class_psc;

// ======================================================================
// helpers

