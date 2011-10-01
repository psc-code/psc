
#include "psc.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

#define VAR(x) (void *)offsetof(struct psc_mfields_c, x)

static struct param psc_mfields_descr[] = {
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
  {},
};

#undef VAR

#define MAKE_MFIELDS_METHODS(type)					\
									\
LIST_HEAD(mfields_##type##_list);					\
									\
void									\
psc_mfields_##type##_set_domain(mfields_##type##_t *flds,		\
				struct mrc_domain *domain)		\
{									\
  flds->domain = domain;						\
}									\
									\
static void								\
_psc_mfields_##type##_setup(mfields_##type##_t *flds)			\
{									\
  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,	\
						     &flds->nr_patches); \
  flds->f = calloc(flds->nr_patches, sizeof(*flds->f));			\
  for (int p = 0; p < flds->nr_patches; p++) {				\
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };	\
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],			\
		   patches[p].ldims[1] + flds->ibn[1],			\
		   patches[p].ldims[2] + flds->ibn[2] };		\
    fields_##type##_alloc(&flds->f[p], ilg, ihg, flds->nr_fields);	\
  }									\
  list_add_tail(&flds->entry, &mfields_##type##_list);			\
}									\
									\
static void								\
_psc_mfields_##type##_destroy(mfields_##type##_t *flds)		        \
{									\
  for (int p = 0; p < flds->nr_patches; p++) {				\
    fields_##type##_free(&flds->f[p]);					\
  }									\
  free(flds->f);							\
  list_del(&flds->entry);						\
}									\
									\
struct mrc_class_psc_mfields_##type mrc_class_psc_mfields_##type = {	\
  .name             = "psc_mfields_" #type,				\
  .size             = sizeof(struct psc_mfields_##type),		\
  .param_descr      = psc_mfields_descr,				\
  .setup            = _psc_mfields_##type##_setup,			\
  .destroy          = _psc_mfields_##type##_destroy,			\
};

MAKE_MFIELDS_METHODS(fortran)
//MAKE_MFIELDS_METHODS(sse2)

// ======================================================================
// psc_fields_c

#include <mrc_io.h>

LIST_HEAD(mfields_c_list);

void
psc_mfields_c_set_domain(mfields_c_t *flds, struct mrc_domain *domain)
{
  flds->domain = domain;
}

static void
_psc_mfields_c_setup(mfields_c_t *flds)
{
  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->f = calloc(flds->nr_patches, sizeof(*flds->f));
  for (int p = 0; p < flds->nr_patches; p++) {
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],
		   patches[p].ldims[1] + flds->ibn[1],
		   patches[p].ldims[2] + flds->ibn[2] };
    fields_c_alloc(&flds->f[p], ilg, ihg, flds->nr_fields);
  }
  list_add_tail(&flds->entry, &mfields_c_list);
}

static void
_psc_mfields_c_destroy(mfields_c_t *flds)
{
  psc_mfields_c_free(flds);
  list_del(&flds->entry);
}

#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

static void
_psc_mfields_c_write(mfields_c_t *mfields, struct mrc_io *io)
{
  int ierr;
  const char *path = psc_mfields_c_name(mfields);
  mrc_io_write_obj_ref(io, path, "domain", (struct mrc_obj *) mfields->domain);
  mrc_io_write_attr_int(io, path, "nr_patches", mfields->nr_patches);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mfields->nr_patches; p++) {
    fields_c_t *fields = &mfields->f[p];
    char name[10]; sprintf(name, "p%d", p);

    hid_t groupp = H5Gcreate(group, name, H5P_DEFAULT, H5P_DEFAULT,
			     H5P_DEFAULT); H5_CHK(groupp);
    ierr = H5LTset_attribute_int(groupp, ".", "ib", fields->ib, 3); CE;
    ierr = H5LTset_attribute_int(groupp, ".", "im", fields->im, 3); CE;
    ierr = H5LTset_attribute_int(groupp, ".", "nr_comp", &fields->nr_comp, 1); CE;
    int with_array = fields->with_array;
    ierr = H5LTset_attribute_int(groupp, ".", "with_array", &with_array, 1); CE;
    for (int m = 0; m < fields->nr_comp; m++) {
      char namec[10]; sprintf(namec, "m%d", m);
      char *s = fields->name[m];
      if (!s) {
	s = "(null)";
      }
      ierr = H5LTset_attribute_string(groupp, ".", namec, s); CE;
    }
    // write components separately instead?
    hsize_t hdims[4] = { fields->nr_comp, fields->im[2], fields->im[1], fields->im[0] };
    ierr = H5LTmake_dataset_double(groupp, "fields_c", 4, hdims, fields->flds); CE;
    ierr = H5Gclose(groupp); CE;
  }

  ierr = H5Gclose(group); CE;
}

static void
_psc_mfields_c_read(mfields_c_t *mfields, struct mrc_io *io)
{
  int ierr;
  const char *path = psc_mfields_c_name(mfields);
  mfields->domain = (struct mrc_domain *)
    mrc_io_read_obj_ref(io, path, "domain", &mrc_class_mrc_domain);
  mrc_io_read_attr_int(io, path, "nr_patches", &mfields->nr_patches);
  psc_mfields_c_setup(mfields);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mfields->nr_patches; p++) {
    fields_c_t *fields = &mfields->f[p];
    char name[10]; sprintf(name, "p%d", p);

    hid_t groupp = H5Gopen(group, name, H5P_DEFAULT); H5_CHK(groupp);
    int ib[3], im[3], nr_comp, with_array;
    ierr = H5LTget_attribute_int(groupp, ".", "ib", ib); CE;
    ierr = H5LTget_attribute_int(groupp, ".", "im", im); CE;
    ierr = H5LTget_attribute_int(groupp, ".", "nr_comp", &nr_comp); CE;
    ierr = H5LTget_attribute_int(groupp, ".", "with_array", &with_array); CE;
    for (int d = 0; d < 3; d++) {
      assert(ib[d] == fields->ib[d]);
      assert(im[d] == fields->im[d]);
    }
    assert(nr_comp == fields->nr_comp);
    assert(!with_array);
    for (int m = 0; m < fields->nr_comp; m++) {
      char namec[10]; sprintf(namec, "m%d", m);

      hsize_t dims;
      H5T_class_t class;
      size_t sz;
      ierr = H5LTget_attribute_info(groupp, ".", namec, &dims, &class, &sz); CE;
      char *s = malloc(sz);
      ierr = H5LTget_attribute_string(groupp, ".", namec, s); CE;
      if (strcmp(s, "(null)") != 0) {
	fields->name[m] = s;
      }
    }

    ierr = H5LTread_dataset_double(groupp, "fields_c", fields->flds); CE;
    ierr = H5Gclose(groupp); CE;
  }

  ierr = H5Gclose(group); CE;
}

#endif

struct mrc_class_psc_mfields_c mrc_class_psc_mfields_c = {
  .name             = "psc_mfields_c",
  .size             = sizeof(struct psc_mfields_c),
  .param_descr      = psc_mfields_descr,
  .setup            = _psc_mfields_c_setup,
  .destroy          = _psc_mfields_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write            = _psc_mfields_c_write,
  .read             = _psc_mfields_c_read,
#endif
};

