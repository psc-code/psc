
#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

using real_t = mparticles_t::real_t;

// ======================================================================
// psc_mparticles "single" / "double" / "c"

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup

static void
PFX(setup)(struct psc_mparticles *_mprts)
{
  mparticles_t mprts(_mprts);
  new(mprts.sub_) mparticles_t::sub_t(ppsc->grid);

  psc_mparticles_setup_super(_mprts);
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

static void
PFX(reserve_all)(struct psc_mparticles *_mprts, uint *n_prts_by_patch)
{
  mparticles_t mprts(_mprts);

  mprts->reserve_all(n_prts_by_patch);
}

static void
PFX(resize_all)(struct psc_mparticles *_mprts, uint *n_prts_by_patch)
{
  mparticles_t mprts(_mprts);

  mprts->resize_all(n_prts_by_patch);
}

static void
PFX(get_size_all)(struct psc_mparticles *_mprts, uint *n_prts_by_patch)
{
  mparticles_t mprts(_mprts);

  mprts->get_size_all(n_prts_by_patch);
}

static unsigned int
PFX(get_nr_particles)(struct psc_mparticles *_mprts)
{
  mparticles_t mprts(_mprts);

  return mprts->get_n_prts();
}

#if PSC_PARTICLES_AS_SINGLE

static void
PFX(inject)(struct psc_mparticles *mprts, int p,
	    const struct psc_particle_inject *new_prt)
{
  int kind = new_prt->kind;

  const Grid_t& grid = ppsc->grid;
  const Grid_t::Patch& patch = grid.patches[p];
  for (int d = 0; d < 3; d++) {
    assert(new_prt->x[d] >= patch.xb[d]);
    assert(new_prt->x[d] <= patch.xe[d]);
  }
  
  float dVi = 1.f / (grid.dx[0] * grid.dx[1] * grid.dx[2]);

  particle_t prt;
  prt.xi      = new_prt->x[0] - patch.xb[0];
  prt.yi      = new_prt->x[1] - patch.xb[1];
  prt.zi      = new_prt->x[2] - patch.xb[2];
  prt.pxi     = new_prt->u[0];
  prt.pyi     = new_prt->u[1];
  prt.pzi     = new_prt->u[2];
  prt.qni_wni_ = new_prt->w * ppsc->kinds[kind].q * dVi;
  prt.kind_   = kind;

  mparticles_t(mprts)[p].push_back(prt);
}

#endif

// ----------------------------------------------------------------------
// psc_mparticles_ops

struct PFX(OPS) : psc_mparticles_ops {
  PFX(OPS)() {
    name                    = PARTICLE_TYPE;
    size                    = sizeof(psc_mparticles_sub);
    methods                 = PFX(methods);
    setup                   = PFX(setup);
    write                   = PFX(write);
    read                    = PFX(read);
    reserve_all             = PFX(reserve_all);
    resize_all              = PFX(resize_all);
    get_size_all            = PFX(get_size_all);
    get_nr_particles        = PFX(get_nr_particles);
#if PSC_PARTICLES_AS_SINGLE
    inject                  = PFX(inject);
#endif
  }
} PFX(ops);

