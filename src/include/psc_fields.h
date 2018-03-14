
#ifndef PSC_FIELDS_H
#define PSC_FIELDS_H

#include <mrc_obj.h>

#include "grid.hxx"
#include "fields_traits.hxx"

// ----------------------------------------------------------------------
// psc_mfields class

struct psc_mfields {
  struct mrc_obj obj;

  // state
  char **comp_name; //> name for each field component

  // parameters
  int nr_fields; //> number of field components
  int ibn[3]; //> number of ghost points

  const Grid_t* grid;
};

MRC_CLASS_DECLARE(psc_mfields, struct psc_mfields);

struct psc_mparticles;

struct psc_mfields_ops {
  MRC_SUBCLASS_OPS(struct psc_mfields);
};

void psc_mfields_set_comp_name(struct psc_mfields *flds, int m, const char *s);
const char *psc_mfields_comp_name(struct psc_mfields *flds, int m);

void psc_mfields_write_as_mrc_fld(struct psc_mfields *mflds, struct mrc_io *io);

struct psc_mfields_list_entry {
  struct psc_mfields **flds_p;
  list_t entry;
};

void psc_mfields_list_add(list_t *head, struct psc_mfields **flds_p);
void psc_mfields_list_del(list_t *head, struct psc_mfields **flds_p);

#define psc_mfields_ops(flds) (struct psc_mfields_ops *) ((flds)->obj.ops)

extern list_t psc_mfields_base_list;

#endif


