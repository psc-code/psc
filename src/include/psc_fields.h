
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
  void (*zero_comp)(struct psc_mfields *mflds, int m);
  void (*set_comp)(struct psc_mfields *mflds, int m, double alpha);
  void (*scale_comp)(struct psc_mfields *mflds, int m, double alpha);
  void (*copy_comp)(struct psc_mfields *to, int mto, struct psc_mfields *from, int mfrom);
  void (*axpy_comp)(struct psc_mfields *y, int my, double alpha,
		    struct psc_mfields *x, int mx);
  double (*max_comp)(struct psc_mfields *mflds, int m);

  double (*synchronize_tang_e_norm_b)(struct psc_mfields *mflds);
  void (*compute_div_b_err)(struct psc_mfields *mflds);
  double (*compute_rms_div_b_err)(struct psc_mfields *mflds);
  void (*clean_div_b)(struct psc_mfields *mflds);
  void (*compute_div_e_err)(struct psc_mfields *mflds);
  double (*compute_rms_div_e_err)(struct psc_mfields *mflds);
  void (*clean_div_e)(struct psc_mfields *mflds);
  void (*clear_rhof)(struct psc_mfields *mflds);
  void (*accumulate_rho_p)(struct psc_mfields *mflds, struct psc_mparticles *mprts);
  void (*synchronize_rho)(struct psc_mfields *mflds);
  void (*compute_rhob)(struct psc_mfields *mflds);
  void (*compute_curl_b)(struct psc_mfields *mflds);
};

typedef void (*psc_mfields_copy_func_t)(struct psc_mfields *, struct psc_mfields *,
					int, int);

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

double psc_mfields_max_comp(struct psc_mfields *mflds, int m);

void psc_mfields_set_comp_name(struct psc_mfields *flds, int m, const char *s);
const char *psc_mfields_comp_name(struct psc_mfields *flds, int m);
struct psc_mfields *psc_mfields_get_as(struct psc_mfields *mflds_base,
				       const char *type, int mb, int me);
void psc_mfields_put_as(struct psc_mfields *mflds,
			struct psc_mfields *mflds_base, int mb, int me);

void psc_mfields_write_as_mrc_fld(struct psc_mfields *mflds, struct mrc_io *io);

double psc_mfields_synchronize_tang_e_norm_b(struct psc_mfields *mflds);
void psc_mfields_compute_div_b_err(struct psc_mfields *mflds);
double psc_mfields_compute_rms_div_b_err(struct psc_mfields *mflds);
void psc_mfields_clean_div_b(struct psc_mfields *mflds);
void psc_mfields_compute_div_e_err(struct psc_mfields *mflds);
double psc_mfields_compute_rms_div_e_err(struct psc_mfields *mflds);
void psc_mfields_clean_div_e(struct psc_mfields *mflds);
void psc_mfields_clear_rhof(struct psc_mfields *mflds);
void psc_mfields_accumulate_rho_p(struct psc_mfields *mflds, struct psc_mparticles *mprts);
void psc_mfields_synchronize_rho(struct psc_mfields *mflds);
void psc_mfields_compute_rhob(struct psc_mfields *mflds);
void psc_mfields_compute_curl_b(struct psc_mfields *mflds);

template<typename MF>
inline MF psc_mfields::get_as(int mb, int me)
{
  const char *type = fields_traits<typename MF::fields_t>::name;
  struct psc_mfields *mflds = psc_mfields_get_as(this, type, mb, me);
  return MF(mflds);
}

extern struct psc_mfields_ops psc_mfields_c_ops;
extern struct psc_mfields_ops psc_mfields_fortran_ops;
extern struct psc_mfields_ops psc_mfields_single_ops;
extern struct psc_mfields_ops psc_mfields_vpic_ops;
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


