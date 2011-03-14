
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

// FIXME, update
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

MRC_CLASS_DECLARE(psc_case, struct psc_case);

void psc_case_init_field(struct psc_case *_case, mfields_base_t *flds);
void psc_case_init_npt(struct psc_case *_case, int kind, double x[3],
			struct psc_particle_npt *npt);

#endif

