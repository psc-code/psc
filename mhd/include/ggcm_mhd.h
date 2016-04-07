
#ifndef GGCM_MHD_H
#define GGCM_MHD_H

#include <mrc_obj.h>

#include <mrc_fld.h>
#include <mrc_ts.h>

// ======================================================================
// ggcm_mhd
//
// This object runs an MHD simulation

MRC_CLASS_DECLARE(ggcm_mhd, struct ggcm_mhd);

void ggcm_mhd_fill_ghosts(struct ggcm_mhd *mhd, struct mrc_fld *fld,
			  int m, float bntim);
void ggcm_mhd_fill_ghosts_E(struct ggcm_mhd *mhd, struct mrc_fld *E);
void ggcm_mhd_fill_ghosts_reconstr(struct ggcm_mhd *mhd, struct mrc_fld *U_l[],
				   struct mrc_fld *U_r[]);
void ggcm_mhd_calc_divb(struct ggcm_mhd *mhd, struct mrc_fld *fld,
			struct mrc_fld *divb);
void ggcm_mhd_calc_currcc(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m,
			struct mrc_fld *currcc);
void ggcm_mhd_get_state(struct ggcm_mhd *mhd);
void ggcm_mhd_set_state(struct ggcm_mhd *mhd);

int ggcm_mhd_ntot(struct ggcm_mhd *mhd);

void ggcm_mhd_default_box(struct ggcm_mhd *mhd);

void ggcm_mhd_convert_from_primitive(struct ggcm_mhd *mhd,
				     struct mrc_fld *fld_base);

struct mrc_fld *ggcm_mhd_fld_get_as(struct mrc_fld *fld_base, const char *type,
				    int mhd_type, int mb, int me);
void ggcm_mhd_fld_put_as(struct mrc_fld *fld, struct mrc_fld *fld_base,
			 int mb, int me);

struct mrc_fld *ggcm_mhd_get_fld_as_fortran(struct mrc_fld *fld_base);
void ggcm_mhd_put_fld_as_fortran(struct mrc_fld *fld, struct mrc_fld *fld_base);

// primitive fluid variables, face-centered B
#define MT_PRIMITIVE (0)
// primitive fluid variables, cell-centered B
#define MT_PRIMITIVE_CC (1)

// the following has B staggered the openggcm way: [-1..mx[
#define MT_SEMI_CONSERVATIVE_GGCM (2)

// the following have B staggered the "normal" way: [0..mx]
#define MT_SEMI_CONSERVATIVE (3)
#define MT_FULLY_CONSERVATIVE (4)

// cell-centered fully conservative MHD
#define MT_FULLY_CONSERVATIVE_CC (5)

// the multi-moment schemes are cell-centered for all quantities
#define MT_GKEYLL (6)

// ----------------------------------------------------------------------
// wrappers / helpers

void ggcm_mhd_wrongful_death(struct ggcm_mhd *mhd, struct mrc_fld *x, int errcode);

void ts_ggcm_mhd_step_calc_rhs(void *ctx, struct mrc_obj *_rhs, float time,
			       struct mrc_obj *_x);
void ts_ggcm_mhd_step_run(void *ctx, struct mrc_ts *ts, struct mrc_obj *_x);

int ggcm_mhd_main(int *argc, char ***argv);

// ----------------------------------------------------------------------

void ggcm_mhd_register();

#endif
