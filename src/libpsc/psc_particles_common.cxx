
#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

using real_t = mparticles_t::real_t;

// ======================================================================
// psc_mparticles "single" / "double" / "c"

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
  assert(sizeof(particle_t) / sizeof(real_t) == 8);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mprts->nr_patches; p++) {
    mparticles_t::patch_t& prts = mparticles_t(mprts)[p];
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = prts.size();
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      // in a rather ugly way, we write the int "kind" member as a float / double
      hsize_t hdims[2];
      hdims[0] = n_prts;
      hdims[1] = 8;
#if PSC_PARTICLES_AS_DOUBLE
      ierr = H5LTmake_dataset_double(pgroup, "data", 2, hdims,
				    (double *) &*prts.begin()); CE;
#elif PSC_PARTICLES_AS_SINGLE
      ierr = H5LTmake_dataset_float(pgroup, "data", 2, hdims,
				    (float *) &*prts.begin()); CE;
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

  PFX(setup)(mprts);
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);

  for (int p = 0; p < mprts->nr_patches; p++) {
    mparticles_t::patch_t& prts = mparticles_t(mprts)[p];
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    mparticles_t(mprts)[p].reserve(n_prts);
    
    if (n_prts > 0) {
#if PSC_PARTICLES_AS_SINGLE
      ierr = H5LTread_dataset_float(pgroup, "data",
				    (float *) &*prts.begin()); CE;
#elif PSC_PARTICLES_AS_DOUBLE
      ierr = H5LTread_dataset_double(pgroup, "data",
				    (double *) &*prts.begin()); CE;
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

// ----------------------------------------------------------------------
// psc_mparticles_ops

