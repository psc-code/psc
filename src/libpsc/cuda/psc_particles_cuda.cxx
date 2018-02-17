
#include "psc.h"
#include "cuda_iface.h"
#include "psc_particles_cuda.h"
#include "psc_push_particles.h"

#include <mrc_io.h>

// ======================================================================
// psc_mparticles "cuda"

// ----------------------------------------------------------------------
// psc_mparticles_cuda_methods

static struct mrc_obj_method psc_mparticles_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_mparticles_cuda::copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mparticles_cuda::copy_from_single),
  MRC_OBJ_METHOD("copy_to_double"  , psc_mparticles_cuda::copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_mparticles_cuda::copy_from_double),
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_cuda_setup

static void
psc_mparticles_cuda_setup(struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);

  cuda_base_init();

  new(mprts.sub()) psc_mparticles_cuda(ppsc->grid);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_destroy

static void
psc_mparticles_cuda_destroy(struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);
  
  mprts.sub()->~psc_mparticles_cuda();
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_reserve_all

static void
psc_mparticles_cuda_reserve_all(struct psc_mparticles *_mprts, uint *n_prts_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->reserve_all(n_prts_by_patch);
}

#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_mparticles_cuda_write

static void
psc_mparticles_cuda_write(struct psc_mparticles *_mprts, struct mrc_io *io)
{
  mparticles_cuda_t mprts(_mprts);
  int ierr;
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  uint n_prts_by_patch[mprts.n_patches()];
  mprts->get_size_all(n_prts_by_patch);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, _mprts), H5P_DEFAULT); H5_CHK(group);
  uint off = 0;
  // FIXME, reorder first if necessary
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      float4 *xi4  = (float4 *) calloc(n_prts, sizeof(*xi4));
      float4 *pxi4 = (float4 *) calloc(n_prts, sizeof(*pxi4));
      
      mprts->from_device(xi4, pxi4, n_prts, off);
      
      hsize_t hdims[2];
      hdims[0] = n_prts; hdims[1] = 4;
      ierr = H5LTmake_dataset_float(pgroup, "xi4", 2, hdims, (float *) xi4); CE;
      ierr = H5LTmake_dataset_float(pgroup, "pxi4", 2, hdims, (float *) pxi4); CE;
      
      free(xi4);
      free(pxi4);
    }
    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_read

static void
psc_mparticles_cuda_read(struct psc_mparticles *_mprts, struct mrc_io *io)
{
  mparticles_cuda_t mprts(_mprts);

  psc_mparticles_read_super(_mprts, io);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, _mprts), H5P_DEFAULT); H5_CHK(group);

  int n_prts_by_patch[mprts.n_patches()];

  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    n_prts_by_patch[p] = n_prts;
    ierr = H5Gclose(pgroup); CE;
  }

  psc_mparticles_setup(_mprts);
  psc_mparticles_reserve_all(_mprts, n_prts_by_patch);

  uint off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    if (n_prts > 0) {
      float4 *xi4  = (float4*) calloc(n_prts, sizeof(float4));
      float4 *pxi4 = (float4*) calloc(n_prts, sizeof(float4));
      
      ierr = H5LTread_dataset_float(pgroup, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(pgroup, "pxi4", (float *) pxi4); CE;
      
      mparticles_cuda_t mprts = mparticles_cuda_t(_mprts);
      mprts->to_device(xi4, pxi4, n_prts, off);
      
      free(xi4);
      free(pxi4);
    }

    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }

  ierr = H5Gclose(group); CE;
  mprts->setup_internals();
}

#endif

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_nr_particles

static uint
psc_mparticles_cuda_get_nr_particles(struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);

  return mprts->get_n_prts();
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_size_all

static void
psc_mparticles_cuda_get_size_all(struct psc_mparticles *_mprts, uint *n_prts_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->get_size_all(n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_resize_all

static void
psc_mparticles_cuda_resize_all(struct psc_mparticles *_mprts, uint *n_prts_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->resize_all(n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_inject

void
psc_mparticles_cuda_inject(struct psc_mparticles *_mprts, struct cuda_mparticles_prt *buf,
			   uint *buf_n_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->inject(buf, buf_n_by_patch);
}

// ----------------------------------------------------------------------
// mparticles_cuda_t::patch_t::get_b_mx

const int* mparticles_cuda_t::patch_t::get_b_mx() const
{
  return mp_->patch_get_b_mx(p_);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_ops_cuda : psc_mparticles_ops {
  psc_mparticles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(struct psc_mparticles_cuda);
    methods                 = psc_mparticles_cuda_methods;
    setup                   = psc_mparticles_cuda_setup;
    destroy                 = psc_mparticles_cuda_destroy;
#ifdef HAVE_LIBHDF5_HL
    read                    = psc_mparticles_cuda_read;
    write                   = psc_mparticles_cuda_write;
#endif
    reserve_all             = psc_mparticles_cuda_reserve_all;
    get_nr_particles        = psc_mparticles_cuda_get_nr_particles;
    resize_all              = psc_mparticles_cuda_resize_all;
    get_size_all            = psc_mparticles_cuda_get_size_all;
  }
} psc_mparticles_cuda_ops;

