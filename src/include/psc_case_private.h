
#ifndef PSC_CASE_PRIVATE_H
#define PSC_CASE_PRIVATE_H

#include <psc_case.h>

///Physics cases (initial conditions and other parameters)
///
///A "case" defines what problem the code is supposed to run.
///
///If you want to run a physics problem that is not included in the predefined cases, you will have to create a case for it yourself. This requires the following steps:
///
///Create a "case_<casename>.c" file in the "src" directory. Instead of writing this from scratch, take one of the existing cases which seems close and copy it, then search and replace the old case's name with your new <casename>. In particular, you may have to modify:
///@param init_param 	Set up all kinds of parameters, in particular, you can override defaults (e.g. instead of the physical units used by default, you can set the electron charge to 1).
///@param init_field 	Set initial condition for the fields (E, B, j).
///@param init_npt 	Set up particles. For the cell located at the passed arguments, set charge (q), density (n), mass (m), velocity (v) and temperature (T).
///\sa psc_case_ops
///
///Add \p case_<casename>.c to "src/Makefile.am" next to the other cases.
///
///Add \code extern struct psc_case_ops psc_case_<casename>_ops; \endcode 
///to \p psc_case_private.h, following how it's done for the other cases.
///
///Add your psc_case_<casename>_ops to the registration by adding 
///\code mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_<casename>_ops); \endcode
///near the end of "src/psc_case.c".

struct psc_case {
  struct mrc_obj obj;
  bool seed_by_time;
  int nr_kinds;

  struct psc *psc;
};


///Interface for custom cases
///
///Take a look at the existing test-cases to see how to overload them.
///@note 
///When you have added an implementation of psc_case_ops, add its declaration
///as extern to the end of psc_case_private.h and psc_case.c
///@sa \ref custom_case
struct psc_case_ops {
  MRC_SUBCLASS_OPS(struct psc_case);
  void (*init_npt)(struct psc_case *_case, int kind, double x[3],
		   struct psc_particle_npt *npt);	///< Initizalizes particles
  void (*init_field)(struct psc_case *_case, mfields_base_t *flds);	///< Initializes fields
  void (*init_photon_np)(struct psc_case *_case, double x[3],	///< Initializes photons
			 struct psc_photon_np *np);
  void (*integrate)(struct psc_case *_case);
};

extern struct psc_case_ops psc_case_harris_ops;
extern struct psc_case_ops psc_case_test_xy_ops;
extern struct psc_case_ops psc_case_test_xz_ops;
extern struct psc_case_ops psc_case_test_yz_ops;
extern struct psc_case_ops psc_case_harris_xy_ops;
extern struct psc_case_ops psc_case_langmuir_ops;
extern struct psc_case_ops psc_case_wakefield_ops;
extern struct psc_case_ops psc_case_thinfoil_ops;
extern struct psc_case_ops psc_case_foils_ops;
extern struct psc_case_ops psc_case_curvedfoil_ops;
extern struct psc_case_ops psc_case_singlepart_ops;
extern struct psc_case_ops psc_case_collisions_ops;
extern struct psc_case_ops psc_case_cone_ops;
extern struct psc_case_ops psc_case_microsphere_ops;
extern struct psc_case_ops psc_case_photon_test_ops;
extern struct psc_case_ops psc_case_bubble_ops;

extern struct psc_case_ops psc_case_dynamic_ops;
// ======================================================================

#define psc_case_ops(c) ((struct psc_case_ops *)((c)->obj.ops))

#endif
