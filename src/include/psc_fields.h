
#ifndef PSC_FIELDS_H
#define PSC_FIELDS_H

#include <mrc_obj.h>

// ----------------------------------------------------------------------
// psc_fields class

MRC_CLASS_DECLARE(psc_fields, struct psc_fields);

// ----------------------------------------------------------------------
// psc_mfields class

struct psc_mfields {
  struct mrc_obj obj;
  struct psc_fields **flds;
  int nr_patches;
  struct mrc_domain *domain;
  int nr_fields; //> number of field components
  char **comp_name; //> name for each field component
  int ibn[3];
  int first_comp; //> The first component in this field (normally 0)
  int ib[3]; //> lower left corner for each patch (incl. ghostpoints)
  int im[3]; //> extent for each patch (incl. ghostpoints)
  void **data;
};

MRC_CLASS_DECLARE(psc_mfields, struct psc_mfields);

struct psc_mfields_ops {
  MRC_SUBCLASS_OPS(struct psc_mfields);
  void (*zero_comp)(struct psc_mfields *mflds, int m);
  void (*set_comp)(struct psc_mfields *mflds, int m, double alpha);
  void (*scale_comp)(struct psc_mfields *mflds, int m, double alpha);
  void (*copy_comp)(struct psc_mfields *to, int mto, struct psc_mfields *from, int mfrom);
  void (*axpy_comp)(struct psc_mfields *y, int my, double alpha,
		    struct psc_mfields *x, int mx);
};

typedef void (*psc_mfields_copy_func_t)(struct psc_mfields *, struct psc_mfields *,
					int, int);

void psc_mfields_set_domain(struct psc_mfields *flds,
			    struct mrc_domain *domain);
void psc_mfields_zero_comp(struct psc_mfields *flds, int m);
void psc_mfields_zero_range(struct psc_mfields *flds, int mb, int me);
void psc_mfields_set_comp(struct psc_mfields *flds, int m, double alpha);
void psc_mfields_scale(struct psc_mfields *flds, double alpha);
void psc_mfields_copy_comp(struct psc_mfields *to, int mto,
			   struct psc_mfields *from, int mfrom);
void psc_mfields_axpy(struct psc_mfields *yf, double alpha,
		      struct psc_mfields *xf);
void psc_mfields_axpy_comp(struct psc_mfields *yf, int ym, double alpha,
			   struct psc_mfields *xf, int xm);
void psc_mfields_set_comp_name(struct psc_mfields *flds, int m, const char *s);
const char *psc_mfields_comp_name(struct psc_mfields *flds, int m);
struct psc_mfields *psc_mfields_get_as(struct psc_mfields *mflds_base,
				       const char *type, int mb, int me);
void psc_mfields_put_as(struct psc_mfields *mflds,
			struct psc_mfields *mflds_base, int mb, int me);

void psc_mfields_write_as_mrc_fld(struct psc_mfields *mflds, struct mrc_io *io);

static inline struct psc_fields *
psc_mfields_get_patch(struct psc_mfields *flds, int p)
{
  return flds->flds[p];
}

extern struct psc_mfields_ops psc_mfields_c_ops;
extern struct psc_mfields_ops psc_mfields_fortran_ops;
extern struct psc_mfields_ops psc_mfields_single_ops;
extern struct psc_mfields_ops psc_mfields_cuda_ops;

struct psc_mfields_list_entry {
  struct psc_mfields **flds_p;
  list_t entry;
};

void psc_mfields_list_add(list_t *head, struct psc_mfields **flds_p);
void psc_mfields_list_del(list_t *head, struct psc_mfields **flds_p);

#define psc_mfields_ops(flds) (struct psc_mfields_ops *) ((flds)->obj.ops)

extern list_t psc_mfields_base_list;

#endif


