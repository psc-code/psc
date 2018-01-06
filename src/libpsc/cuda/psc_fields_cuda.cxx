
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
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds_cuda)->cmflds;
  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  FieldsH F(flds);

  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
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

    cuda_mfields_copy_to_device(cmflds, p, flds, mb, me);
  }
  
  flds.dtor();
}

static void
psc_mfields_cuda_copy_to_c(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_c,
			  int mb, int me)
{
  mfields_c_t mf_c(mflds_c);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds_cuda)->cmflds;
  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  FieldsH F(flds);

  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
    FieldsC F_c(mf_c[p]);
    cuda_mfields_copy_from_device(cmflds, p, flds, mb, me);
  
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
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds_cuda)->cmflds;
  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  FieldsH F(flds);
  
  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
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

    cuda_mfields_copy_to_device(cmflds, p, flds, mb, me);
  }
  
  flds.dtor();
}

static void
psc_mfields_cuda_copy_to_single(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_single,
				int mb, int me)
{
  mfields_single_t mf_single(mflds_single);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds_cuda)->cmflds;
  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  FieldsH F(flds);

  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
    FieldsS F_s(mf_single[p]);
    cuda_mfields_copy_from_device(cmflds, p, flds, mb, me);
  
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
psc_mfields_cuda_setup(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda *sub = psc_mfields_cuda(mflds);

  psc_mfields_setup_super(mflds);

  struct mrc_patch *patches = mrc_domain_get_patches(mflds->domain,
						     &mflds->nr_patches);

  int im[3], ib[3];
  assert(mflds->nr_patches > 0);
  for (int p = 0; p < mflds->nr_patches; p++) {
    for (int d = 0; d < 3; d++) {
      if (p == 0) {
	ib[d] = -mflds->ibn[d];
	im[d] = patches[0].ldims[d] + 2 * mflds->ibn[d];
      } else {
	assert(patches[p].ldims[d] == patches[0].ldims[d]);
      }
    }
  }
  
  cuda_base_init();

  mrc_json_t json = mrc_json_object_new(0);

  mrc_json_t info = mrc_json_object_new(0);
  mrc_json_object_push(json, "info", info);
  mrc_json_object_push_integer(info, "n_patches", mflds->nr_patches);
  mrc_json_object_push_integer(info, "n_fields", mflds->nr_fields);
  mrc_json_object_push_integer_array(info, "ib", 3, ib);
  mrc_json_object_push_integer_array(info, "im", 3, im);
  mrc_json_object_push_integer_array(info, "ldims", 3, ppsc->patch[0].ldims);
  mrc_json_object_push_double_array(info, "dx", 3, ppsc->patch[0].dx);

  mrc_json_print(json, 0);

  sub->cmflds = cuda_mfields_create();
  cuda_mfields_ctor(sub->cmflds, json);

  // FIXME json_builder_free(obj);
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_destroy

static void
psc_mfields_cuda_destroy(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda *sub = psc_mfields_cuda(mflds);

  cuda_mfields_dtor(sub->cmflds);
  cuda_mfields_destroy(sub->cmflds);
  sub->cmflds = NULL;
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_zero_comp

static void
psc_mfields_cuda_zero_comp(struct psc_mfields *mflds, int m)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;

  assert(ppsc->domain.gdims[0] == 1);
  cuda_mfields_zero_comp_yz(cmflds, m);
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_axpy_comp

static void
psc_mfields_cuda_axpy_comp(struct psc_mfields *mflds_y, int my, double alpha,
			   struct psc_mfields *mflds_x, int mx)
{
  struct cuda_mfields *cmflds_y = psc_mfields_cuda(mflds_y)->cmflds;
  struct cuda_mfields *cmflds_x = psc_mfields_cuda(mflds_x)->cmflds;

  assert(ppsc->domain.gdims[0] == 1);
  cuda_mfields_axpy_comp_yz(cmflds_y, my, alpha, cmflds_x, mx);
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
psc_mfields_cuda_write(struct psc_mfields *mflds, struct mrc_io *io)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);

  for (int p = 0; p < mflds->nr_patches; p++) {
    cuda_mfields_copy_from_device(cmflds, p, flds, 0, flds.nr_comp);
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
psc_mfields_cuda_read(struct psc_mfields *mflds, struct mrc_io *io)
{
  psc_mfields_read_super(mflds, io);
  
  psc_mfields_cuda_setup(mflds);

  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  for (int p = 0; p < mflds->nr_patches; p++) {
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
    cuda_mfields_copy_to_device(cmflds, p, flds, 0, flds.nr_comp);
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

struct psc_mfields_ops psc_mfields_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_mfields_cuda),
  .methods               = psc_mfields_cuda_methods,
  .setup                 = psc_mfields_cuda_setup,
  .destroy               = psc_mfields_cuda_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_mfields_cuda_write,
  .read                  = psc_mfields_cuda_read,
#endif
  .zero_comp             = psc_mfields_cuda_zero_comp,
  .axpy_comp             = psc_mfields_cuda_axpy_comp,
};

