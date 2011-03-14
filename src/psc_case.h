
#ifndef PSC_CASE_H
#define PSC_CASE_H

#include "psc.h"

struct psc_particle_npt {
  double q; ///< charge
  double m; ///< mass
  double n; ///< density
  double p[3]; ///< momentum
  double T[3]; ///< temperature
};

/////////////////////////////////////////////////////////////////////////
/// Physics cases (initial conditions and other parameters).
///
/// A "case" defines what problem the code is supposed to run.
/// Currently, we have the following cases predefined:
/// - "harris": Double Harris sheet in 2D (xz)
/// - "langmuir": A 1D Langmuir wave (z)
///
/// If you want to run a physics problem that is not included in the predefined
/// cases, you will have to create a case for it yourself. This requires the
/// following steps:
/// - Create a "case_<casename>.c" file in the "src" directory. Instead of
///   writing this from scratch, take one of the existing cases which seems close
///   and copy it, then search and replace the old case's name with your new
///   <casename>.
///   In particular, you may have to modify:
///   @param init_param Set up all kinds of parameters, in particular, you
///     can override defaults (e.g. instead of the physical units used by default,
///     you can set the electron charge to 1).
///   @param init_field Set initial condition for the fields (E, B, j).
///   @param init_npt Set up particles. For the cell located at the passed
///     arguments, set charge (q), density (n), mass (m), velocity (v) and
///     temperature (T).
/// - Add "case_<casename>.c" to "src/Makefile.am" next to the other cases.
/// - Add the psc_case_ops_<casename> variable to psc.h following how it's done
///   for the other cases.
/// - Add your psc_case_ops_<casename> to the list of available cases in 
///   "src/init_parameters.c".

struct psc_case;

struct psc_case_ops {
  const char *name; ///< Name of case.
  size_t ctx_size; ///< Size of private context (e.g., for parameters)
  struct param *ctx_descr; ///< Description of user-settable parameters
  void (*create)(struct psc_case *); ///< Function to set up needed environment.
  void (*destroy)(struct psc_case *); ///< Funtion to cleanup environment.
  void (*init_param)(struct psc_case *); ///< Initialize simulation parameters based on case.
  void (*init_field)(struct psc_case *, mfields_base_t *flds); ///< Initialize fields relevant to case.
  void (*init_npt)(struct psc_case *, int kind, double x[3],
		   struct psc_particle_npt *npt);
};

struct psc_case {
  struct psc_case_ops *ops;
  void *ctx;
};

struct psc_case *psc_case_create(const char *case_name);
void psc_case_destroy(struct psc_case *Case);

static inline void
psc_case_init_param(struct psc_case *Case)
{
  if (Case->ops->init_param) {
    Case->ops->init_param(Case);
  }
}

static inline void
psc_case_init_field(struct psc_case *Case, mfields_base_t *flds)
{
  if (Case->ops->init_field) {
    Case->ops->init_field(Case, flds);
  }
}

static inline void
psc_case_init_npt(struct psc_case *Case, int kind, double xx[3],
		  struct psc_particle_npt *npt)
{
  Case->ops->init_npt(Case, kind, xx, npt);
}

MRC_CLASS_DECLARE(_psc_case, struct _psc_case);

#endif

