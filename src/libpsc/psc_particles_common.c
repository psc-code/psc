
#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#if PSC_PARTICLES_AS_DOUBLE

#define PFX(x) psc_particles_double_ ## x
#define MPFX(x) psc_mparticles_double_ ## x
#define psc_particles_sub psc_particles_double
#define psc_mparticles_sub psc_mparticles_double

#elif PSC_PARTICLES_AS_SINGLE

#define PFX(x) psc_particles_single_ ## x
#define MPFX(x) psc_mparticles_single_ ## x
#define psc_particles_sub psc_particles_single
#define psc_mparticles_sub psc_mparticles_single

#elif PSC_PARTICLES_AS_SINGLE_BY_BLOCK

#define PFX(x) psc_particles_single_by_block_ ## x
#define MPFX(x) psc_mparticles_single_by_block_ ## x
#define psc_particles_sub psc_particles_single_by_block
#define psc_mparticles_sub psc_mparticles_single_by_block

#elif PSC_PARTICLES_AS_C

#define PFX(x) psc_particles_c_ ## x
#define MPFX(x) psc_mparticles_c_ ## x
#define psc_particles_sub psc_particles_c
#define psc_mparticles_sub psc_mparticles_c

#elif PSC_PARTICLES_AS_FORTRAN

#define PFX(x) psc_particles_fortran_ ## x
#define MPFX(x) psc_mparticles_fortran_ ## x
#define psc_particles_sub psc_particles_fortran
#define psc_mparticles_sub psc_mparticles_fortran

#endif

// ======================================================================
// psc_particles "single" / "double" / "c"

// ----------------------------------------------------------------------
// psc_particles_sub_setup

static void
PFX(setup)(struct psc_particles *prts)
{
  struct psc_particles_sub *sub = psc_particles_sub(prts);

  int n_alloced = psc_particles_size(prts) * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  sub->particles = calloc(n_alloced, sizeof(*sub->particles));

#if PSC_PARTICLES_AS_SINGLE
  sub->particles_alt = calloc(n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  sub->particles_alt = calloc(n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
  sub->b_off = calloc(sub->nr_blocks + 2, sizeof(*sub->b_off));
#endif
}

// ----------------------------------------------------------------------
// psc_particles_sub_destroy

static void
PFX(destroy)(struct psc_particles *prts)
{
  struct psc_particles_sub *sub = psc_particles_sub(prts);

  free(sub->particles);

#if PSC_PARTICLES_AS_SINGLE
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
  free(sub->b_off);
#endif
}

// ----------------------------------------------------------------------
// psc_particles_sub_realloc

void
PFX(realloc)(struct psc_particles *prts, int new_n_part)
{
  struct psc_particles_sub *sub = psc_particles_sub(prts);

  if (new_n_part <= psc_mparticles_n_alloced(prts->mprts, prts->p))
    return;

  int n_alloced = new_n_part * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  sub->particles = realloc(sub->particles, n_alloced * sizeof(*sub->particles));

#if PSC_PARTICLES_AS_SINGLE
  sub->b_idx = realloc(sub->b_idx, n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(n_alloced * sizeof(*sub->particles_alt));
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  sub->b_idx = realloc(sub->b_idx, n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(n_alloced * sizeof(*sub->particles_alt));
#endif
}

// ----------------------------------------------------------------------
// psc_particles: subclass ops

struct psc_particles_ops PFX(ops) = {
  .name                    = PARTICLE_TYPE,
  .size                    = sizeof(struct psc_particles_sub),
  .setup                   = PFX(setup),
  .destroy                 = PFX(destroy),
};

// ======================================================================
// psc_mparticles

#if (PSC_PARTICLES_AS_DOUBLE || PSC_PARTICLES_AS_SINGLE) && HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_mparticles_sub_write

static void
MPFX(write)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  assert(sizeof(particle_t) / sizeof(particle_real_t) == 8);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      // in a rather ugly way, we write the int "kind" member as a float / double
      hsize_t hdims[2] = { n_prts, 8 };
#if PSC_PARTICLES_AS_DOUBLE
      ierr = H5LTmake_dataset_double(pgroup, "data", 2, hdims,
				    (double *) particle_iter_deref(prts.begin)); CE;
#elif PSC_PARTICLES_AS_SINGLE
      ierr = H5LTmake_dataset_float(pgroup, "data", 2, hdims,
				    (float *) particle_iter_deref(prts.begin)); CE;
#else
      assert(0);
#endif
    }
    ierr = H5Gclose(pgroup); CE;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_read

static void
MPFX(read)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  // FIXME those should be superclass bits
  mprts->domain = mrc_io_read_ref(io, mprts, "domain", mrc_domain);
  mrc_domain_get_patches(mprts->domain, &mprts->nr_patches);
  mrc_io_read_int(io, mprts, "flags", (int *) &mprts->flags);

  mprts->prts = calloc(mprts->nr_patches, sizeof(*mprts->prts));
  mprts->mpatch = calloc(mprts->nr_patches, sizeof(*mprts->mpatch));

  for (int p = 0; p < mprts->nr_patches; p++) {
    char name[20]; sprintf(name, "prts%d", p);
    mprts->prts[p] = psc_particles_create(MPI_COMM_NULL);
    psc_particles_set_type(mprts->prts[p], PARTICLE_TYPE);
    psc_particles_set_name(mprts->prts[p], name);
    mprts->prts[p]->p = p;

    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    psc_particles_set_n_prts(mprts->prts[p], n_prts);
    psc_particles_setup(mprts->prts[p]);
    
    if (n_prts > 0) {
#if PSC_PARTICLES_AS_SINGLE
      ierr = H5LTread_dataset_float(pgroup, "data",
				    (float *) particle_iter_deref(prts.begin)); CE;
#elif PSC_PARTICLES_AS_DOUBLE
      ierr = H5LTread_dataset_double(pgroup, "data",
				    (double *) particle_iter_deref(prts.begin)); CE;
#else
      assert(0);
#endif
    }
    ierr = H5Gclose(pgroup); CE;
  }
  ierr = H5Gclose(group); CE;
}

#else

static void
MPFX(write)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  assert(0);
}

static void
MPFX(read)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  assert(0);
}

#endif

struct psc_mparticles_ops MPFX(ops) = {
  .name                    = PARTICLE_TYPE,
  .methods                 = MPFX(methods),
  .write                   = MPFX(write),
  .read                    = MPFX(read),
};

