
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_collision.h>
#include <psc_bnd.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

#if 0

// ======================================================================
// psc_bnd "amr"

#include "psc_bnd_private.h"

#include <psc_fields_as_c.h>
#include <mrc_ddc.h>

struct psc_bnd_amr {
  struct mrc_ddc *ddc_E;
  struct mrc_ddc *ddc_H;
};

#define psc_bnd_amr(bnd) mrc_to_subobj(bnd, struct psc_bnd_amr)

// ----------------------------------------------------------------------

#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))

static struct mrc_ddc_amr_stencil_entry stencil_coarse_EX[2] = {
  // FIXME, 3D
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 0, 1, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_EY[2] = {
  // FIXME, 3D
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 1, 0, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_EZ[4] = {
  { .dx = { 0, 0, 0 }, .val = .25f },
  { .dx = { 1, 0, 0 }, .val = .25f },
  { .dx = { 0, 1, 0 }, .val = .25f },
  { .dx = { 1, 1, 0 }, .val = .25f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_HX[2] = {
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 1, 0, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_HY[2] = {
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 0, 1, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_HZ[2] = {
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 0, 0, 1 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil stencils_coarse[] = {
  [EX] = { stencil_coarse_EX, ARRAY_SIZE(stencil_coarse_EX) },
  [EY] = { stencil_coarse_EY, ARRAY_SIZE(stencil_coarse_EY) },
  [EZ] = { stencil_coarse_EZ, ARRAY_SIZE(stencil_coarse_EZ) },
  [HX] = { stencil_coarse_HX, ARRAY_SIZE(stencil_coarse_HX) },
  [HY] = { stencil_coarse_HY, ARRAY_SIZE(stencil_coarse_HY) },
  [HZ] = { stencil_coarse_HZ, ARRAY_SIZE(stencil_coarse_HZ) },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_EX[6] = {
  // FIXME, 3D
  { .dx = {  0, -1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = { +1, -1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = {  0,  0,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = { +1,  0,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = {  0, +1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = { +1, +1,  0 }, .val = (1.f/8.f) * 1.f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_EY[6] = {
  // FIXME, 3D
  { .dx = { -1,  0,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = {  0,  0,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = { +1,  0,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = { -1, +1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = {  0, +1,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = { +1, +1,  0 }, .val = (1.f/8.f) * 1.f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_EZ[9] = {
  // FIXME, 3D
  { .dx = { -1, -1,  0 }, .val = (2.f/8.f) * .25f },
  { .dx = {  0, -1,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = { +1, -1,  0 }, .val = (2.f/8.f) * .25f },
  { .dx = { -1,  0,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f  },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = { -1, +1,  0 }, .val = (2.f/8.f) * .25f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * .25f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_HX[6] = {
  // FIXME, 3D
  { .dx = { -1,  0,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = { -1, +1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * .5f },
};
	  
static struct mrc_ddc_amr_stencil_entry stencil_fine_HY[6] = {
  // FIXME, 3D
  { .dx = {  0, -1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = { +1, -1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * .5f },
};
	  
static struct mrc_ddc_amr_stencil_entry stencil_fine_HZ[4] = {
  // FIXME, 3D
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * 1.f },
};

static struct mrc_ddc_amr_stencil stencils_fine[] = {
  [EX] = { stencil_fine_EX, ARRAY_SIZE(stencil_fine_EX) },
  [EY] = { stencil_fine_EY, ARRAY_SIZE(stencil_fine_EY) },
  [EZ] = { stencil_fine_EZ, ARRAY_SIZE(stencil_fine_EZ) },
  [HX] = { stencil_fine_HX, ARRAY_SIZE(stencil_fine_HX) },
  [HY] = { stencil_fine_HY, ARRAY_SIZE(stencil_fine_HY) },
  [HZ] = { stencil_fine_HZ, ARRAY_SIZE(stencil_fine_HZ) },
};

// ----------------------------------------------------------------------
// psc_bnd_amr_create_ddc

static void
psc_bnd_amr_create_ddc(struct psc_bnd *bnd)
{
  struct psc_bnd_amr *bnd_amr = psc_bnd_amr(bnd);

  bnd_amr->ddc_E = mrc_domain_create_ddc(bnd->psc->mrc_domain);
  mrc_ddc_set_param_int(bnd_amr->ddc_E, "sw", 3);
  mrc_ddc_set_param_int(bnd_amr->ddc_E, "size_of_type", sizeof(fields_t::real_t));
  mrc_ddc_setup(bnd_amr->ddc_E);
  mrc_ddc_amr_set_by_stencil(bnd_amr->ddc_E, EX, 2, (int[]) { 0, 1, 1 }, &stencils_coarse[EX], &stencils_fine[EX]);
  mrc_ddc_amr_set_by_stencil(bnd_amr->ddc_E, EY, 2, (int[]) { 1, 0, 1 }, &stencils_coarse[EY], &stencils_fine[EY]);
  mrc_ddc_amr_set_by_stencil(bnd_amr->ddc_E, EZ, 2, (int[]) { 1, 1, 0 }, &stencils_coarse[EZ], &stencils_fine[EZ]);
  mrc_ddc_amr_assemble(bnd_amr->ddc_E);

  bnd_amr->ddc_H = mrc_domain_create_ddc(bnd->psc->mrc_domain);
  mrc_ddc_set_param_int(bnd_amr->ddc_H, "sw", 3);
  mrc_ddc_set_param_int(bnd_amr->ddc_H, "size_of_type", sizeof(fields_t::real_t));
  mrc_ddc_setup(bnd_amr->ddc_H);
  mrc_ddc_amr_set_by_stencil(bnd_amr->ddc_H, HX, 2, (int[]) { 1, 0, 0 }, &stencils_coarse[HX], &stencils_fine[HX]);
  mrc_ddc_amr_set_by_stencil(bnd_amr->ddc_H, HY, 2, (int[]) { 0, 1, 0 }, &stencils_coarse[HY], &stencils_fine[HY]);
  mrc_ddc_amr_set_by_stencil(bnd_amr->ddc_H, HZ, 2, (int[]) { 0, 0, 1 }, &stencils_coarse[HZ], &stencils_fine[HZ]);
  mrc_ddc_amr_assemble(bnd_amr->ddc_H);

  // FIXME, this psc_bnd::ddc is fake, not really usable (need E, H separately)
  struct mrc_ddc *ddc = mrc_domain_create_ddc(bnd->psc->mrc_domain);
  mrc_ddc_set_param_int(ddc, "sw", 3);
  mrc_ddc_set_param_int(ddc, "size_of_type", sizeof(fields_t::real_t));
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

static void _mrc_unused
psc_bnd_amr_unsetup(struct psc_bnd *bnd)
{
  // FIXME!!! this gets called after bnd_amr is already freed!
  /* struct psc_bnd_amr *bnd_amr = psc_bnd_amr(bnd); */

  /* mrc_ddc_destroy(bnd_amr->ddc_E); */
  /* mrc_ddc_destroy(bnd_amr->ddc_H); */
}

// ----------------------------------------------------------------------
// psc_bnd_amr_destroy

static void
psc_bnd_amr_destroy(struct psc_bnd *bnd)
{
}

// ----------------------------------------------------------------------
// psc_bnd_amr_fill_ghosts

static void
psc_bnd_amr_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds, int mb, int me)
{
  mfields_t mf(mflds);
  struct psc_bnd_amr *bnd_amr = psc_bnd_amr(bnd);

  fields_t::real_t **fldp = malloc(mflds->nr_patches * sizeof(*fldp));
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = mf[p];
    fldp[p] = flds.data;
  }

  if (mb == EX && me == EX + 3) {
    mrc_ddc_fill_ghosts(bnd_amr->ddc_E, -1, -1, fldp);
  } else if (mb == HX && me == HX + 3) {
    mrc_ddc_fill_ghosts(bnd_amr->ddc_H, -1, -1, fldp);
  } else {
    mprintf("fill mb %d me %d\n", mb, me);
  }
  free(fldp);
}

// ----------------------------------------------------------------------
// psc_bnd_amr_add_ghosts

static void
psc_bnd_amr_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base, int mb, int me)
{
  mprintf("add mb %d me %d\n", mb, me);
}

// ======================================================================
// psc_bnd: subclass "amr"

struct psc_bnd_ops psc_bnd_amr_ops = {
  .name                    = "amr",
  .size                    = sizeof(struct psc_bnd_amr),
  .setup                   = psc_bnd_amr_setup,
  .destroy                 = psc_bnd_amr_destroy,
  .create_ddc              = psc_bnd_amr_create_ddc,
  .fill_ghosts             = psc_bnd_amr_fill_ghosts,
  .add_ghosts              = psc_bnd_amr_add_ghosts,
};

#endif

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
  mrc_crds_set_param_double3(crds, "l",  (double[3]) { psc->domain.corner[0],
	psc->domain.corner[1], psc->domain.corner[2] });
  mrc_crds_set_param_double3(crds, "h",  (double[3]) {
      psc->domain.corner[0] + psc->domain.length[0],
      psc->domain.corner[1] + psc->domain.length[1],
      psc->domain.corner[2] + psc->domain.length[2] });

#define AMR_DOMAIN 3

#if AMR_DOMAIN == 0
  mrc_domain_add_patch(domain, 0, (int [3]) { 0, 0, 0 });
#elif AMR_DOMAIN == 1
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 1, 0 });
#elif AMR_DOMAIN == 2
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
#elif AMR_DOMAIN == 3
  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 1, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 3, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 5, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 5, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
#endif

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
#if 0
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_amr_ops);
#endif

  return psc_main(&argc, &argv, &psc_test_fdtd_ops);
}
