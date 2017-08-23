
#include "psc.h"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "psc_fields_inc.h"

#include "psc_fields_common.c"

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
// psc_fields_c_write

static void
psc_fields_c_write(struct psc_fields *flds, struct mrc_io *io)
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
  ierr = H5LTmake_dataset_double(group, "fields_c", 4, hdims, flds->data); CE;
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_fields_c_read

static void
psc_fields_c_read(struct psc_fields *flds, struct mrc_io *io)
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
  ierr = H5LTread_dataset_double(group, "fields_c", flds->data); CE;
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================
// psc_mfields: subclass "c"
  
struct psc_mfields_ops psc_mfields_c_ops = {
  .name                  = "c",
};

// ======================================================================
// psc_fields: subclass "c"
  
struct psc_fields_ops psc_fields_c_ops = {
  .name                  = "c",
  .setup                 = psc_fields_c_setup,
  .destroy               = psc_fields_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                  = psc_fields_c_read,
  .write                 = psc_fields_c_write,
#endif
  .zero_comp             = psc_fields_c_zero_comp,
  .set_comp              = psc_fields_c_set_comp,
  .scale_comp            = psc_fields_c_scale_comp,
  .copy_comp             = psc_fields_c_copy_comp,
  .axpy_comp             = psc_fields_c_axpy_comp,
};

