
#include "psc.h"
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "psc_fields_cuda.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"
#include "fields.hxx"

#include <mrc_params.h>

using FieldsH = Fields3d<fields_single_t>; // host
using FieldsS = Fields3d<fields_single_t> ;// host
using FieldsC = Fields3d<fields_c_t>; // host

// OPT, CUDA fields have too many ghostpoints, and 7 points in the invar direction!

// ======================================================================
// convert from/to "c"

static void
psc_mfields_cuda_copy_from_c(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_c,
			    int mb, int me)
{
  mfields_c_t mf_c(mflds_c);
  mfields_cuda_t mf_cuda(mflds_cuda);
  fields_single_t flds = mf_cuda->get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda->n_patches(); p++) {
    FieldsC F_c(mf_c[p]);
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    F(m, jx,jy,jz) = F_c( m, jx,jy,jz);
	  }
	}
      }
    }

    mf_cuda->copy_to_device(p, flds, mb, me);
  }
  
  flds.dtor();
}

static void
psc_mfields_cuda_copy_to_c(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_c,
			  int mb, int me)
{
  mfields_c_t mf_c(mflds_c);
  mfields_cuda_t mf_cuda(mflds_cuda);
  fields_single_t flds = mf_cuda->get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda->n_patches(); p++) {
    FieldsC F_c(mf_c[p]);
    mf_cuda->copy_from_device(p, flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    F_c(m, jx,jy,jz) = F(m, jx,jy,jz);
	  }
	}
      }
    }
  }

  flds.dtor();
}

// ======================================================================
// convert from/to "single"

static void
psc_mfields_cuda_copy_from_single(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_single,
				  int mb, int me)
{
  mfields_single_t mf_single(mflds_single);
  mfields_cuda_t mf_cuda(mflds_cuda);
  fields_single_t flds = mf_cuda->get_host_fields();
  FieldsH F(flds);
  
  for (int p = 0; p < mf_cuda->n_patches(); p++) {
    FieldsS F_s(mf_single[p]);

    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    F(m, jx,jy,jz) = F_s(m, jx,jy,jz);
	  }
	}
      }
    }

    mf_cuda->copy_to_device(p, flds, mb, me);
  }
  
  flds.dtor();
}

static void
psc_mfields_cuda_copy_to_single(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_single,
				int mb, int me)
{
  mfields_single_t mf_single(mflds_single);
  mfields_cuda_t mf_cuda(mflds_cuda);
  fields_single_t flds = mf_cuda->get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda->n_patches(); p++) {
    FieldsS F_s(mf_single[p]);
    mf_cuda->copy_from_device(p, flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    F_s(m, jx,jy,jz) = F(m, jx,jy,jz);
	  }
	}
      }
    }
  }

  flds.dtor();
}

// ======================================================================

// ----------------------------------------------------------------------
// psc_mfields_cuda_setup

static void
psc_mfields_cuda_setup(struct psc_mfields *_mflds)
{
  mfields_cuda_t mflds(_mflds);

  psc_mfields_setup_super(_mflds);

  cuda_base_init();

  assert(_mflds->grid);
  new(mflds.sub()) psc_mfields_cuda(*_mflds->grid, mflds.n_fields(), _mflds->ibn);
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_destroy

static void
psc_mfields_cuda_destroy(struct psc_mfields *_mflds)
{
  mfields_cuda_t mflds(_mflds);

  mflds->~psc_mfields_cuda();
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_zero_comp

static void
psc_mfields_cuda_zero_comp(struct psc_mfields *_mflds, int m)
{
  mfields_cuda_t mflds(_mflds);

  mflds->zero_comp(m);
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_axpy_comp

static void
psc_mfields_cuda_axpy_comp(struct psc_mfields *_mflds_y, int my, double alpha,
			   struct psc_mfields *_mflds_x, int mx)
{
  mfields_cuda_t mflds_y(_mflds_y);
  mfields_cuda_t mflds_x(_mflds_x);

  assert(ppsc->domain.gdims[0] == 1);
  mflds_y->axpy_comp_yz(my, alpha, mflds_x, mx);
}

#ifdef HAVE_LIBHDF5_HL

#include <mrc_io.h>

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_mfields_write

static void
psc_mfields_cuda_write(struct psc_mfields *_mflds, struct mrc_io *io)
{
  mfields_cuda_t mflds(_mflds);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, _mflds), H5P_DEFAULT); H5_CHK(group0);

  fields_single_t flds = mflds->get_host_fields();

  for (int p = 0; p < mflds.n_patches(); p++) {
    mflds->copy_from_device(p, flds, 0, flds.nr_comp);
    char name[20]; sprintf(name, "flds%d", p);
    hid_t group = H5Gcreate(group0, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    
    ierr = H5LTset_attribute_int(group, ".", "ib", flds.ib, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "im", flds.im, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "nr_comp", &flds.nr_comp, 1); CE;
    // write components separately instead?
    hsize_t hdims[4];
    hdims[0] = flds.nr_comp;
    hdims[1] = flds.im[2]; hdims[2] = flds.im[1]; hdims[3] = flds.im[0];
    ierr = H5LTmake_dataset_float(group, "fields_cuda", 4, hdims, flds.data); CE;
    ierr = H5Gclose(group); CE;
  }
  flds.dtor();

  ierr = H5Gclose(group0); CE;
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_read

static void
psc_mfields_cuda_read(struct psc_mfields *_mflds, struct mrc_io *io)
{
  psc_mfields_read_super(_mflds, io);
  
  psc_mfields_cuda_setup(_mflds);

  mfields_cuda_t mflds(_mflds);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, _mflds), H5P_DEFAULT); H5_CHK(group0);

  fields_single_t flds = mflds->get_host_fields();
  for (int p = 0; p < mflds.n_patches(); p++) {
    char name[20]; sprintf(name, "flds%d", p);
    hid_t group = H5Gopen(group0, name, H5P_DEFAULT); H5_CHK(group);

    int ib[3], im[3], nr_comp;
    ierr = H5LTget_attribute_int(group, ".", "ib", ib); CE;
    ierr = H5LTget_attribute_int(group, ".", "im", im); CE;
    ierr = H5LTget_attribute_int(group, ".", "nr_comp", &nr_comp); CE;
    for (int d = 0; d < 3; d++) {
      assert(ib[d] == flds.ib[d]);
      assert(im[d] == flds.im[d]);
    }
    assert(nr_comp == flds.nr_comp);

    ierr = H5LTread_dataset_float(group, "fields_cuda", flds.data); CE;
    mflds->copy_to_device(p, flds, 0, flds.nr_comp);
    ierr = H5Gclose(group); CE;
  }
  flds.dtor();
  ierr = H5Gclose(group0); CE;
}

#endif

// ======================================================================
// psc_mfields: subclass "cuda"
  
static struct mrc_obj_method psc_mfields_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_mfields_cuda_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_mfields_cuda_copy_from_c),
  MRC_OBJ_METHOD("copy_to_single"  , psc_mfields_cuda_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mfields_cuda_copy_from_single),
  {}
};

struct psc_mfields_ops_cuda : psc_mfields_ops {
  psc_mfields_ops_cuda() {
    name                  = "cuda";
    size                  = sizeof(struct psc_mfields_cuda);
    methods               = psc_mfields_cuda_methods;
    setup                 = psc_mfields_cuda_setup;
    destroy               = psc_mfields_cuda_destroy;
#ifdef HAVE_LIBHDF5_HL
    write                 = psc_mfields_cuda_write;
    read                  = psc_mfields_cuda_read;
#endif
  }
} psc_mfields_cuda_ops;

