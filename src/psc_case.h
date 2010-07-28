
#ifndef PSC_CASE_H
#define PSC_CASE_H

#include "psc.h"

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
///   @param init_nvt Set up particles. For the cell located at the passed
///     arguments, set charge (q), density (n), mass (m), velocity (v) and
///     temperature (T).
/// - Add "case_<casename>.c" to "src/Makefile.am" next to the other cases.
/// - Add the psc_case_ops_<casename> variable to psc.h following how it's done
///   for the other cases.
/// - Add your psc_case_ops_<casename> to the list of available cases in 
///   "src/init_parameters.c".

struct psc_case_ops {
  const char *name; ///< Name of case.
  void (*create)(void); ///< Function to set up needed environment.
  void (*destroy)(void); ///< Funtion to cleanup environment.
  void (*init_param)(void); ///< Initialize simulation parameters based on case.
  void (*init_field)(void); ///< Initialize fields relevant to case.
  void (*init_nvt)(int kind, double x[3], double *q, double *m, double *n,
		   double v[3], double T[3]);
};

struct psc_case {
  struct psc_case_ops *ops;
  void *ctx;
};

static inline void
psc_case_init_param(struct psc_case *Case)
{
  if (Case->ops->init_param) {
    Case->ops->init_param();
  }
}

static inline void
psc_case_init_field(struct psc_case *Case)
{
  if (Case->ops->init_field) {
    Case->ops->init_field();
  }
}

static inline void
psc_case_init_nvt(struct psc_case *Case, int kind, double xx[3],
		  double *q, double *m, double *n, double v[3], double T[3])
{
  Case->ops->init_nvt(kind, xx, q, m, n, v, T);
}

#endif

