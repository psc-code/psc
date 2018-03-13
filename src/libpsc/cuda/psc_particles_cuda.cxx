
#include "psc.h"
#include "cuda_iface.h"
#include "psc_particles_cuda.h"
#include "psc_push_particles.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"

#include <mrc_io.h>

// ======================================================================
// psc_mparticles "cuda"

// ----------------------------------------------------------------------
// psc_mparticles_cuda_methods

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
  PscMparticlesCuda mprts(_mprts);
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
  PscMparticlesCuda mprts(_mprts);

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
  mprts->reserve_all(n_prts_by_patch);

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
      
      PscMparticlesCuda mprts = PscMparticlesCuda(_mprts);
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
// psc_mparticles_cuda_inject

void
psc_mparticles_cuda_inject(struct psc_mparticles *_mprts, struct cuda_mparticles_prt *buf,
			   uint *buf_n_by_patch)
{
  PscMparticlesCuda mprts(_mprts);

  mprts->inject_buf(buf, buf_n_by_patch);
}

// ----------------------------------------------------------------------
// PscMparticlesCuda::patch_t::get_b_mx

const int* PscMparticlesCuda::patch_t::get_b_mx() const
{
  return mp_.patch_get_b_mx(p_);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
psc_mparticles_ops_<MparticlesCuda> psc_mparticles_cuda_ops;
