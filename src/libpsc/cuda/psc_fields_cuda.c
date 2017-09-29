
#include "psc.h"
#include "cuda_iface.h"
#include "cuda_mfields.h" // FIXME
#include "psc_fields_cuda.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"
#include "psc_cuda.h"

#include <mrc_params.h>

// OPT, CUDA fields have too many ghostpoints, and 7 points in the invar direction!

// ----------------------------------------------------------------------
// this really should just use fields_single_t

fields_cuda_t
_fields_cuda_t_mflds(struct psc_mfields *mflds, int p, fields_cuda_real_t *h_flds)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;

  fields_cuda_t flds;

  flds.data = h_flds;
  for (int d = 0; d < 3; d++) {
    flds.ib[d] = cmflds->ib[d];
    flds.im[d] = cmflds->im[d];
  }
  flds.nr_comp = mflds->nr_fields;
  flds.first_comp = mflds->first_comp;

  return flds;
}

// ----------------------------------------------------------------------
// psc_mfields_get_patch_size

static unsigned int
psc_mfields_get_patch_size(struct psc_mfields *mflds)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;

  return cmflds->im[0] * cmflds->im[1] * cmflds->im[2] * mflds->nr_fields;
}

// ======================================================================
// convert from/to "c"

static void
psc_mfields_cuda_copy_from_c(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_c,
			    int mb, int me)
{
  float *h_flds = malloc(psc_mfields_get_patch_size(mflds_cuda) * sizeof(*h_flds));

  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
    fields_cuda_t flds_cuda = _fields_cuda_t_mflds(mflds_cuda, p, h_flds);
    fields_c_t flds_c = fields_c_t_mflds(mflds_c, p);
    for (int m = mb; m < me; m++) {
      for (int jz = flds_cuda.ib[2]; jz < flds_cuda.ib[2] + flds_cuda.im[2]; jz++) {
	for (int jy = flds_cuda.ib[1]; jy < flds_cuda.ib[1] + flds_cuda.im[1]; jy++) {
	  for (int jx = flds_cuda.ib[0]; jx < flds_cuda.ib[0] + flds_cuda.im[0]; jx++) {
	    _F3_CUDA(flds_cuda, m, jx,jy,jz) = _F3_C(flds_c, m, jx,jy,jz);
	  }
	}
      }
    }

    __fields_cuda_to_device(mflds_cuda, p, h_flds, mb, me);
  }
  
  free(h_flds);
}

static void
psc_mfields_cuda_copy_to_c(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_c,
			  int mb, int me)
{
  float *h_flds = malloc(psc_mfields_get_patch_size(mflds_cuda) * sizeof(*h_flds));

  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
    fields_cuda_t flds_cuda = _fields_cuda_t_mflds(mflds_cuda, p, h_flds);
    fields_c_t flds_c = fields_c_t_mflds(mflds_c, p);
    __fields_cuda_from_device(mflds_cuda, p, h_flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds_cuda.ib[2]; jz < flds_cuda.ib[2] + flds_cuda.im[2]; jz++) {
	for (int jy = flds_cuda.ib[1]; jy < flds_cuda.ib[1] + flds_cuda.im[1]; jy++) {
	  for (int jx = flds_cuda.ib[0]; jx < flds_cuda.ib[0] + flds_cuda.im[0]; jx++) {
	    _F3_C(flds_c, m, jx,jy,jz) = _F3_CUDA(flds_cuda, m, jx,jy,jz);
	  }
	}
      }
    }
  }

  free(h_flds);
}

// ======================================================================
// convert from/to "single"

static void
psc_mfields_cuda_copy_from_single(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_single,
				  int mb, int me)
{
  float *h_flds = malloc(psc_mfields_get_patch_size(mflds_cuda) * sizeof(*h_flds));

  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
    fields_cuda_t flds_cuda = _fields_cuda_t_mflds(mflds_cuda, p, h_flds);
    fields_single_t flds_single = fields_single_t_mflds(mflds_single, p);

    for (int m = mb; m < me; m++) {
      for (int jz = flds_cuda.ib[2]; jz < flds_cuda.ib[2] + flds_cuda.im[2]; jz++) {
	for (int jy = flds_cuda.ib[1]; jy < flds_cuda.ib[1] + flds_cuda.im[1]; jy++) {
	  for (int jx = flds_cuda.ib[0]; jx < flds_cuda.ib[0] + flds_cuda.im[0]; jx++) {
	    _F3_CUDA(flds_cuda, m, jx,jy,jz) = _F3_S(flds_single, m, jx,jy,jz);
	  }
	}
      }
    }

    __fields_cuda_to_device(mflds_cuda, p, h_flds, mb, me);
  }
  
  free(h_flds);
}

static void
psc_mfields_cuda_copy_to_single(struct psc_mfields *mflds_cuda, struct psc_mfields *mflds_single,
				int mb, int me)
{
  float *h_flds = malloc(psc_mfields_get_patch_size(mflds_cuda) * sizeof(*h_flds));

  for (int p = 0; p < mflds_cuda->nr_patches; p++) {
    fields_cuda_t flds_cuda = _fields_cuda_t_mflds(mflds_cuda, p, h_flds);
    fields_single_t flds_single = fields_single_t_mflds(mflds_single, p);
    __fields_cuda_from_device(mflds_cuda, p, h_flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds_cuda.ib[2]; jz < flds_cuda.ib[2] + flds_cuda.im[2]; jz++) {
	for (int jy = flds_cuda.ib[1]; jy < flds_cuda.ib[1] + flds_cuda.im[1]; jy++) {
	  for (int jx = flds_cuda.ib[0]; jx < flds_cuda.ib[0] + flds_cuda.im[0]; jx++) {
	    _F3_S(flds_single, m, jx,jy,jz) = _F3_CUDA(flds_cuda, m, jx,jy,jz);
	  }
	}
      }
    }
  }

  free(h_flds);
}

// ======================================================================

// ----------------------------------------------------------------------
// psc_mfields_cuda_setup

static void
psc_mfields_cuda_setup(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda *sub = psc_mfields_cuda(mflds);

  struct mrc_patch *patches = mrc_domain_get_patches(mflds->domain,
						     &mflds->nr_patches);

  psc_mfields_setup_super(mflds);

  cuda_base_init();

  struct cuda_mfields *cmflds = cuda_mfields_create();
  sub->cmflds = cmflds;

  cmflds->bnd_by_patch = calloc(mflds->nr_patches, sizeof(*cmflds->bnd_by_patch));

  __psc_mfields_cuda_setup(mflds, patches);
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_destroy

static void
psc_mfields_cuda_destroy(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct cuda_mfields *cmflds = mflds_cuda->cmflds;

  __psc_mfields_cuda_destroy(mflds);

  free(cmflds->bnd_by_patch);
  cuda_mfields_destroy(cmflds);
  mflds_cuda->cmflds = NULL;
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_zero_comp

static void
psc_mfields_cuda_zero_comp(struct psc_mfields *mflds, int m)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    assert(ppsc->domain.gdims[0] == 1);
    cuda_zero_comp_yz(mflds, m, p);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_axpy_comp

static void
psc_mfields_cuda_axpy_comp(struct psc_mfields *y, int my, double alpha,
			   struct psc_mfields *x, int mx)
{
  for (int p = 0; p < y->nr_patches; p++) {
    assert(ppsc->domain.gdims[0] == 1);
    cuda_axpy_comp_yz(y, my, alpha, x, mx, p);
  }
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
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  for (int p = 0; p < mflds->nr_patches; p++) {
    float *h_flds = malloc(psc_mfields_get_patch_size(mflds) * sizeof(*h_flds));
    fields_cuda_t flds = _fields_cuda_t_mflds(mflds, p, h_flds);
    __fields_cuda_from_device(mflds, p, flds.data, 0, flds.nr_comp);
    char name[20]; sprintf(name, "flds%d", p);
    hid_t group = H5Gcreate(group0, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    
    ierr = H5LTset_attribute_int(group, ".", "ib", flds.ib, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "im", flds.im, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "nr_comp", &flds.nr_comp, 1); CE;
    // write components separately instead?
    hsize_t hdims[4] = { flds.nr_comp, flds.im[2], flds.im[1], flds.im[0] };
    ierr = H5LTmake_dataset_float(group, "fields_cuda", 4, hdims, flds.data); CE;
    free(h_flds);
    ierr = H5Gclose(group); CE;
  }

  ierr = H5Gclose(group0); CE;
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_read

static void
psc_mfields_cuda_read(struct psc_mfields *mflds, struct mrc_io *io)
{
  psc_mfields_read_super(mflds, io);
  
  psc_mfields_cuda_setup(mflds);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  for (int p = 0; p < mflds->nr_patches; p++) {
    float *h_flds = malloc(psc_mfields_get_patch_size(mflds) * sizeof(*h_flds));
    fields_cuda_t flds = _fields_cuda_t_mflds(mflds, p, h_flds);
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
    __fields_cuda_to_device(mflds, p, flds.data, 0, flds.nr_comp);
    free(h_flds);
    ierr = H5Gclose(group); CE;
  }
  ierr = H5Gclose(group0); CE;
}

#endif

// ----------------------------------------------------------------------
// psc_mfields_cuda_get_field_t
//
// FIXME, should not exist...

fields_cuda_t
psc_mfields_cuda_get_field_t(struct psc_mfields *mflds, int p)
{
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;
  fields_cuda_t flds;

  flds.data = NULL;
  for (int d = 0; d < 3; d++) {
    flds.ib[d] = cmflds->ib[d];
    flds.im[d] = cmflds->im[d];
  }
  flds.nr_comp = mflds->nr_fields;
  flds.first_comp = mflds->first_comp;

  return flds;
}

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

