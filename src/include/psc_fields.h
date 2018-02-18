
#ifndef PSC_FIELDS_H
#define PSC_FIELDS_H

#include <mrc_obj.h>

#include "fields_traits.hxx"

// ----------------------------------------------------------------------
// psc_mfields class

struct psc_mfields {
  struct mrc_obj obj;

  // state
  int nr_patches;
  char **comp_name; //> name for each field component

  // parameters
  struct mrc_domain *domain;
  int nr_fields; //> number of field components
  int ibn[3]; //> number of ghost points
  int first_comp; //> The first component in this field (normally 0)

  template<typename MF>
  MF get_as(int mb, int me);
};

MRC_CLASS_DECLARE(psc_mfields, struct psc_mfields);

struct psc_mparticles;

struct psc_mfields_ops {
  MRC_SUBCLASS_OPS(struct psc_mfields);
};

typedef void (*psc_mfields_copy_func_t)(struct psc_mfields *, struct psc_mfields *,
					int, int);

void psc_mfields_set_comp_name(struct psc_mfields *flds, int m, const char *s);
const char *psc_mfields_comp_name(struct psc_mfields *flds, int m);
struct psc_mfields *psc_mfields_get_as(struct psc_mfields *mflds_base,
				       const char *type, int mb, int me);
void psc_mfields_put_as(struct psc_mfields *mflds,
			struct psc_mfields *mflds_base, int mb, int me);

void psc_mfields_write_as_mrc_fld(struct psc_mfields *mflds, struct mrc_io *io);

template<typename MF>
inline MF psc_mfields::get_as(int mb, int me)
{
  const char *type = fields_traits<typename MF::fields_t>::name;
  struct psc_mfields *mflds = psc_mfields_get_as(this, type, mb, me);
  return MF(mflds);
}

struct psc_mfields_list_entry {
  struct psc_mfields **flds_p;
  list_t entry;
};

void psc_mfields_list_add(list_t *head, struct psc_mfields **flds_p);
void psc_mfields_list_del(list_t *head, struct psc_mfields **flds_p);

#define psc_mfields_ops(flds) (struct psc_mfields_ops *) ((flds)->obj.ops)

extern list_t psc_mfields_base_list;

#endif


