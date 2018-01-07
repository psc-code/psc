
#include <psc_bnd_particles_private.h>

#include <mrc_io.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>
#include <mrc_domain_private.h> // FIXME

// ----------------------------------------------------------------------
// psc_bnd_particles_destroy

static void
_psc_bnd_particles_destroy(struct psc_bnd_particles *bnd)
{
  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(bnd);
  assert(ops->unsetup);
  ops->unsetup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_write

static void
_psc_bnd_particles_write(struct psc_bnd_particles *bnd, struct mrc_io *io)
{
  mrc_io_write_ref(io, bnd, "psc", bnd->psc);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_read

static void
_psc_bnd_particles_read(struct psc_bnd_particles *bnd, struct mrc_io *io)
{
  bnd->psc = mrc_io_read_ref(io, bnd, "psc", psc);

  psc_bnd_particles_setup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_set_psc

void
psc_bnd_particles_set_psc(struct psc_bnd_particles *bnd, struct psc *psc)
{
  bnd->psc = psc;
}

// ----------------------------------------------------------------------
// check_domain
//
// check if the underlying mrc_domain changed since setup(),
// which might happen, e.g., through rebalancing.
// In this case, do setup() over.

void
psc_bnd_particles_check_domain(struct psc_bnd_particles *bnd)
{
  if (!bnd->ddcp) {
    return;
  }

  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(bnd);

  ops->unsetup(bnd);
  ops->setup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_exchange

void
psc_bnd_particles_exchange(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  static int pr;
  if (!pr) {
    pr = prof_register("xchg_prts", 1., 0, 0);
  }

  //  psc_bnd_particles_check_domain(bnd);

  prof_start(pr);
  psc_stats_start(st_time_comm);
  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(bnd);
  assert(ops->exchange_particles);
  ops->exchange_particles(bnd, mprts);
  psc_stats_stop(st_time_comm);
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_open_calc_moments

void
psc_bnd_particles_open_calc_moments(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  // no need to calculate moments if we're not having any open boundary
  if (!(ppsc->domain.bnd_part_lo[0] == BND_PART_OPEN || ppsc->domain.bnd_part_hi[0] == BND_PART_OPEN ||
	ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN || ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN ||
	ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN || ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN)) {
    return;
  }

  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(bnd);
  assert(ops->open_calc_moments);
  ops->open_calc_moments(bnd, mprts);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_init

static void
psc_bnd_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_auto_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_single2_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_double_omp_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_fortran_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_cuda_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_vpic_ops);
#endif
}

// ======================================================================

#define VAR(x) (void *)offsetof(struct psc_bnd_particles, x)
static struct param psc_bnd_particles_descr[] = {
  { "time_relax" , VAR(time_relax), PARAM_DOUBLE(.1), },

  {},
};
#undef VAR

// ======================================================================
// psc_bnd_particles class

struct mrc_class_psc_bnd_particles mrc_class_psc_bnd_particles = {
  .name             = "psc_bnd_particles",
  .size             = sizeof(struct psc_bnd_particles),
  .param_descr      = psc_bnd_particles_descr,
  .init             = psc_bnd_particles_init,
  .destroy          = _psc_bnd_particles_destroy,
  .write            = _psc_bnd_particles_write,
  .read             = _psc_bnd_particles_read,
};

