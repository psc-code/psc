
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_collision.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>


// ======================================================================
// psc_bnd "amr"

#include "psc_bnd_private.h"

#include <psc_fields_as_c.h>
#include <mrc_ddc.h>

struct psc_bnd_amr {
};

// ----------------------------------------------------------------------
// psc_bnd_amr_create_ddc

static void
psc_bnd_amr_create_ddc(struct psc_bnd *bnd)
{
  // FIXME, this ddc is fake, not really usable (need E, H separately)
  struct mrc_ddc *ddc = mrc_domain_create_ddc(bnd->psc->mrc_domain);
  mrc_ddc_set_param_int(ddc, "sw", 3);
  mrc_ddc_set_param_int(ddc, "size_of_type", sizeof(fields_real_t));
  mrc_ddc_setup(ddc);
  bnd->ddc = ddc;
}

// ----------------------------------------------------------------------
// psc_bnd_amr_setup

static void
psc_bnd_amr_setup(struct psc_bnd *bnd)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->create_ddc);
  ops->create_ddc(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_amr_unsetup

static void
psc_bnd_amr_unsetup(struct psc_bnd *bnd)
{
}

// ----------------------------------------------------------------------
// psc_bnd_amr_fill_ghosts

static void
psc_bnd_amr_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  mprintf("fill mb %d me %d\n", mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_amr_add_ghosts

static void
psc_bnd_amr_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  mprintf("add mb %d me %d\n", mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_amr_exchange_particles

static void
psc_bnd_amr_exchange_particles(struct psc_bnd *bnd, mparticles_base_t *particles_base)
{
  MHERE;
}

// ======================================================================
// psc_bnd: subclass "amr"

struct psc_bnd_ops psc_bnd_amr_ops = {
  .name                    = "amr",
  .size                    = sizeof(struct psc_bnd_amr),
  .setup                   = psc_bnd_amr_setup,
  .unsetup                 = psc_bnd_amr_unsetup,
  .create_ddc              = psc_bnd_amr_create_ddc,
  .fill_ghosts             = psc_bnd_amr_fill_ghosts,
  .add_ghosts              = psc_bnd_amr_add_ghosts,
  .exchange_particles      = psc_bnd_amr_exchange_particles,
};

// ======================================================================
// psc_test_fdtd

// ----------------------------------------------------------------------
// psc_test_fdtd_create

static void
psc_test_fdtd_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 100;
  psc->prm.cfl = 1.;

  psc->domain.length[0] = 1.;
  psc->domain.length[1] = 1.;
  psc->domain.length[2] = 1.;

  psc->domain.gdims[0] = 8;
  psc->domain.gdims[1] = 8;
  psc->domain.gdims[2] = 1;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;

  psc_bnd_set_type(psc->bnd, "amr");
}

// ----------------------------------------------------------------------
// psc_test_fdtd_init_field

static double
psc_test_fdtd_init_field(struct psc *psc, double x[3], int m)
{
  double kx = 2. * M_PI, ky = 2. * M_PI;

  switch (m) {
  case EX: return   1./sqrtf(2.) * sin(kx * x[0] + ky * x[1]);
  case EY: return - 1./sqrtf(2.) * sin(kx * x[0] + ky * x[1]);
  case HZ: return sin(kx * x[0] + ky * x[1]);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_fdtd_setup_mrc_domain

static struct mrc_domain *
psc_test_fdtd_setup_mrc_domain(struct psc *psc, int nr_patches)
{
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  // create a very simple domain decomposition
  int bc[3] = {};
  for (int d = 0; d < 3; d++) {
    if (psc->domain.bnd_fld_lo[d] == BND_FLD_PERIODIC &&
	psc->domain.gdims[d] > 1) {
      bc[d] = BC_PERIODIC;
    }
  }

  mrc_domain_set_type(domain, "amr");
  mrc_domain_set_param_int3(domain, "m", psc->domain.gdims);
  // FIXME, these bc's aren't yet honored by libmrc
  mrc_domain_set_param_int(domain, "bcx", bc[0]);
  mrc_domain_set_param_int(domain, "bcy", bc[1]);
  mrc_domain_set_param_int(domain, "bcz", bc[2]);

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_crds_set_type(crds, "amr_uniform");
  mrc_crds_set_param_int(crds, "sw", 3);
  mrc_crds_set_param_float3(crds, "l",  (float[3]) { psc->domain.corner[0],
	psc->domain.corner[1], psc->domain.corner[2] });
  mrc_crds_set_param_float3(crds, "h",  (float[3]) {
      psc->domain.corner[0] + psc->domain.length[0],
      psc->domain.corner[1] + psc->domain.length[1],
      psc->domain.corner[2] + psc->domain.length[2] });

  mrc_domain_add_patch(domain, 0, (int [3]) { 0, 0, 0 });

  mrc_domain_set_from_options(domain);
  mrc_domain_setup(domain);

  return domain;
}

// ======================================================================
// psc_test_fdtd_ops

struct psc_ops psc_test_fdtd_ops = {
  .name             = "test_fdtd",
  .create           = psc_test_fdtd_create,
  .init_field       = psc_test_fdtd_init_field,
  .setup_mrc_domain = psc_test_fdtd_setup_mrc_domain,
};

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_amr_ops);

  return psc_main(&argc, &argv, &psc_test_fdtd_ops);
}
