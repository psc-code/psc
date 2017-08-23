
#include "psc.h"
#include "psc_fields_as_single.h"
#include "psc_fields_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "psc_fields_inc.h"

#include "psc_fields_common.c"

// FIXME, very duplicated from psc_fields_c.c

static void
psc_fields_single_scale_comp(struct psc_fields *pf, int m, double _val)
{
  fields_single_real_t val = _val;
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_S(pf, m, jx, jy, jz) *= val;
      }
    }
  }
}

static void
psc_fields_single_copy_comp(struct psc_fields *pto, int m_to, struct psc_fields *pfrom, int m_from)
{
  for (int jz = pto->ib[2]; jz < pto->ib[2] + pto->im[2]; jz++) {
    for (int jy = pto->ib[1]; jy < pto->ib[1] + pto->im[1]; jy++) {
      for (int jx = pto->ib[0]; jx < pto->ib[0] + pto->im[0]; jx++) {
	F3_S(pto, m_to, jx, jy, jz) = F3_S(pfrom, m_from, jx, jy, jz);
      }
    }
  }
}

static void
psc_fields_single_axpy_comp(struct psc_fields *y, int ym, double _a, struct psc_fields *x, int xm)
{
  fields_single_real_t a = _a;

  for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
    for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
      for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	F3_S(y, ym, jx, jy, jz) += a * F3_S(x, xm, jx, jy, jz);
      }
    }
  }
}

// ======================================================================
// convert to c

static void
psc_fields_single_copy_from_c(struct psc_fields *flds_single, struct psc_fields *flds_c,
			      int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_single->ib[2]; jz < flds_single->ib[2] + flds_single->im[2]; jz++) {
      for (int jy = flds_single->ib[1]; jy < flds_single->ib[1] + flds_single->im[1]; jy++) {
	for (int jx = flds_single->ib[0]; jx < flds_single->ib[0] + flds_single->im[0]; jx++) {
	  F3_S(flds_single, m, jx,jy,jz) = F3_C(flds_c, m, jx,jy,jz);
	}
      }
    }
  }
}

void
psc_fields_single_copy_to_c(struct psc_fields *flds_single, struct psc_fields *flds_c,
			    int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_single->ib[2]; jz < flds_single->ib[2] + flds_single->im[2]; jz++) {
      for (int jy = flds_single->ib[1]; jy < flds_single->ib[1] + flds_single->im[1]; jy++) {
	for (int jx = flds_single->ib[0]; jx < flds_single->ib[0] + flds_single->im[0]; jx++) {
	  F3_C(flds_c, m, jx,jy,jz) = F3_S(flds_single, m, jx,jy,jz);
	}
      }
    }
  }
}

// ======================================================================

#ifdef HAVE_LIBHDF5_HL

#include <mrc_io.h>

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_fields_single_write

static void
psc_fields_single_write(struct psc_fields *flds, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, flds), H5P_DEFAULT); H5_CHK(group);
  ierr = H5LTset_attribute_int(group, ".", "p", &flds->p, 1); CE;
  ierr = H5LTset_attribute_int(group, ".", "ib", flds->ib, 3); CE;
  ierr = H5LTset_attribute_int(group, ".", "im", flds->im, 3); CE;
  ierr = H5LTset_attribute_int(group, ".", "nr_comp", &flds->nr_comp, 1); CE;
  // write components separately instead?
  hsize_t hdims[4] = { flds->nr_comp, flds->im[2], flds->im[1], flds->im[0] };
  ierr = H5LTmake_dataset_float(group, "fields_single", 4, hdims, flds->data); CE;
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_fields_single_read

static void
psc_fields_single_read(struct psc_fields *flds, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, flds), H5P_DEFAULT); H5_CHK(group);
  int ib[3], im[3], nr_comp;
  ierr = H5LTget_attribute_int(group, ".", "p", &flds->p); CE;
  ierr = H5LTget_attribute_int(group, ".", "ib", ib); CE;
  ierr = H5LTget_attribute_int(group, ".", "im", im); CE;
  ierr = H5LTget_attribute_int(group, ".", "nr_comp", &nr_comp); CE;
  for (int d = 0; d < 3; d++) {
    assert(ib[d] == flds->ib[d]);
    assert(im[d] == flds->im[d]);
  }
  assert(nr_comp == flds->nr_comp);
  psc_fields_setup(flds);
  ierr = H5LTread_dataset_float(group, "fields_single", flds->data); CE;
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================
// psc_mfields: subclass "single"
  
struct psc_mfields_ops psc_mfields_single_ops = {
  .name                  = "single",
};

// ======================================================================
// psc_fields: subclass "single"
  
static struct mrc_obj_method psc_fields_single_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_fields_single_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_fields_single_copy_from_c),
  {}
};

struct psc_fields_ops psc_fields_single_ops = {
  .name                  = "single",
  .methods               = psc_fields_single_methods,
  .setup                 = psc_fields_single_setup,
  .destroy               = psc_fields_single_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_fields_single_write,
  .read                  = psc_fields_single_read,
#endif
  .zero_comp             = psc_fields_single_zero_comp,
  .set_comp              = psc_fields_single_set_comp,
  .scale_comp            = psc_fields_single_scale_comp,
  .copy_comp             = psc_fields_single_copy_comp,
  .axpy_comp             = psc_fields_single_axpy_comp,
};

