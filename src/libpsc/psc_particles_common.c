
#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// ======================================================================
// psc_mparticles "single" / "double" / "c"

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup_patch

static void
PFX(setup_patch)(struct psc_mparticles *mprts, int p, int n_alloced)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  struct PFX(patch) *patch = &sub->patch[p];

  patch->n_prts = 0;
  patch->n_alloced = n_alloced;
  patch->prt_array = calloc(n_alloced, sizeof(*patch->prt_array));

#if PSC_PARTICLES_AS_SINGLE
  patch->prt_array_alt = calloc(n_alloced, sizeof(*patch->prt_array_alt));
  patch->b_idx = calloc(n_alloced, sizeof(*patch->b_idx));
  patch->b_ids = calloc(n_alloced, sizeof(*patch->b_ids));

  int *b_mx = patch->b_mx;
  for (int d = 0; d < 3; d++) {
    b_mx[d] = ppsc->patch[p].ldims[d];
    patch->b_dxi[d] = 1.f / ppsc->patch[p].dx[d];
  }
  patch->nr_blocks = b_mx[0] * b_mx[1] * b_mx[2];
  patch->b_cnt = calloc(patch->nr_blocks + 1, sizeof(*patch->b_cnt));
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  patch->prt_array_alt = calloc(n_alloced, sizeof(*patch->prt_array_alt));
  patch->b_idx = calloc(n_alloced, sizeof(*patch->b_idx));
  patch->b_ids = calloc(n_alloced, sizeof(*patch->b_ids));

  for (int d = 0; d < 3; d++) {
    patch->b_mx[d] = ppsc->patch[p].ldims[d];
    patch->b_dxi[d] = 1.f / ppsc->patch[p].dx[d];
  }
  patch->nr_blocks = patch->b_mx[0] * patch->b_mx[1] * patch->b_mx[2];
  patch->b_cnt = calloc(patch->nr_blocks + 1, sizeof(*patch->b_cnt));
  patch->b_off = calloc(patch->nr_blocks + 2, sizeof(*patch->b_off));
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_destroy_patch

static void
PFX(destroy_patch)(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  struct PFX(patch) *patch = &sub->patch[p];

  free(patch->prt_array);

#if PSC_PARTICLES_AS_SINGLE
  free(patch->prt_array_alt);
  free(patch->b_idx);
  free(patch->b_ids);
  free(patch->b_cnt);
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  free(patch->prt_array_alt);
  free(patch->b_idx);
  free(patch->b_ids);
  free(patch->b_cnt);
  free(patch->b_off);
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup

static void
PFX(setup)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  psc_mparticles_setup_super(mprts);
  sub->patch = calloc(mprts->nr_patches, sizeof(*sub->patch));
}

// ----------------------------------------------------------------------
// psc_mparticls_sub_write/read
  
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
PFX(write)(struct psc_mparticles *mprts, struct mrc_io *io)
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
    int n_prts = particle_range_size(prts);
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
PFX(read)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  psc_mparticles_read_super(mprts, io);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);

  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    PFX(setup_patch)(mprts, p, n_prts);
    
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
PFX(write)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  assert(0);
}

static void
PFX(read)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  assert(0);
}

#endif

static void
PFX(destroy)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    PFX(destroy_patch)(mprts, p);
  }
  free(sub->patch);
}

static void
PFX(reserve)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    PFX(setup_patch)(mprts, p, n_prts_by_patch[p]);
  }
}

static void
PFX(get_n_prts_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    n_prts_by_patch[p] = sub->patch[p].n_prts;
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_ops

struct psc_mparticles_ops PFX(ops) = {
  .name                    = PARTICLE_TYPE,
  .size                    = sizeof(struct psc_mparticles_sub),
  .methods                 = PFX(methods),
  .setup                   = PFX(setup),
  .destroy                 = PFX(destroy),
  .write                   = PFX(write),
  .read                    = PFX(read),
  .reserve                 = PFX(reserve),
  .patch_resize            = PFX(patch_resize),
  .get_n_prts_all          = PFX(get_n_prts_all),
};

