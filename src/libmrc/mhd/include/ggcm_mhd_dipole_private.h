
#ifndef GGCM_MHD_DIPOLE_PRIVATE_H
#define GGCM_MHD_DIPOLE_PRIVATE_H

#include "ggcm_mhd_dipole.h"

struct ggcm_mhd_dipole {
  struct mrc_obj obj;

  // params
  struct ggcm_mhd *mhd;
  double r1lim;
  // dipole strength is given as equatorial B field <dipolestrength>
  // at distance <dipolestrength_r>, e.g., 30574 (nT) at 1 (RE)
  double dipolestrength;
  double dipolestrength_r;

  // state
  struct mrc_fld *bdip;
};

// ggcm_mhd_dipole subclasses exist to provide specific implementations
// for the following purposes:
// - provide the dipole field for purpose of dipole rotation
//   (set_b_field, update_b_field)
// - provide the dipole field for purpose of setting up the initial
//   condition
//   (add_dipole)
//
// So add_dipole() is the more powerful function, which can add a
// dipole not only at the origin, but also at a given location, and it can
// zero it out sunward of xmir.
// A subclass can provide only add_dipole() and that function will be
// used for both purposes (FIXME: actual rotation does not work yet, because a
// generic update_b_field is still missing).
//
// When a subclass provides set_b_field(), this will be used instead
// of add_dipole() (if even present) for purposes of the dipole rotation.
//
// The existence of set_b_field(), update_b_field() is for retaining an interface
// to the legacy fortran functions (which aren't very good, since their dipole isn't
// even divergence-free). A modern implementation should only have to
// provide add_dipole().

struct ggcm_mhd_dipole_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_dipole);
  void (*add_dipole)(struct ggcm_mhd_dipole *mhd_dipole, struct mrc_fld *b,
		     float x0[3], float moment[3], float xmir, float keep);
  void (*set_b_field)(struct ggcm_mhd_dipole *mhd_dipole, struct mrc_fld *b,
		      float moment[3], double diptime);
  void (*update_b_field)(struct ggcm_mhd_dipole *mhd_dipole, struct mrc_fld *bdip,
			 struct mrc_fld *x, double dacttime);
  double (*vector_potential)(struct ggcm_mhd_dipole *mhd_dipole, int m, double x[3],
			     float x0[3], float moment[3], float xmir);
};

extern struct ggcm_mhd_dipole_ops ggcm_mhd_dipole_float_ops;
extern struct ggcm_mhd_dipole_ops ggcm_mhd_dipole_double_ops;
extern struct ggcm_mhd_dipole_ops ggcm_mhd_dipole_none_ops;

#endif
