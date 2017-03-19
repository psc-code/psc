
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

void ggcm_mhd_fill_ghosts(struct ggcm_mhd *mhd, struct mrc_fld *fld, float bntim);
void ggcm_mhd_fill_ghosts_E(struct ggcm_mhd *mhd, struct mrc_fld *E);
void ggcm_mhd_fill_ghosts_reconstr(struct ggcm_mhd *mhd, struct mrc_fld *U_l[],
				   struct mrc_fld *U_r[], int p);
void ggcm_mhd_correct_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *fluxes[3]);
void ggcm_mhd_correct_E(struct ggcm_mhd *mhd, struct mrc_fld *E);
void ggcm_mhd_calc_divb(struct ggcm_mhd *mhd, struct mrc_fld *fld,
			struct mrc_fld *divb);
void ggcm_mhd_calc_currcc(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m,
			struct mrc_fld *currcc);
void ggcm_mhd_calc_rr(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld);
void ggcm_mhd_calc_v(struct ggcm_mhd *mhd, struct mrc_fld *v, struct mrc_fld *fld);
void ggcm_mhd_calc_pp(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld);
void ggcm_mhd_set_state(struct ggcm_mhd *mhd);
void ggcm_mhd_pre_step(struct ggcm_mhd *mhd, struct mrc_ts *ts, struct mrc_fld *fld);
void ggcm_mhd_post_step(struct ggcm_mhd *mhd, struct mrc_ts *ts, struct mrc_fld *fld);

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

void ggcm_mhd_setup_ts(struct ggcm_mhd *mhd, struct mrc_ts *ts);

#define MT_FORMULATION_MASK 3
#define MT_FORMULATION_PRIMITIVE 0
#define MT_FORMULATION_SCONS 1
#define MT_FORMULATION_FCONS 2
#define MT_FORMULATION_GKEYLL 3

#define MT_BGRID_MASK 12
// B staggered the "normal" way: [0..mx]
#define MT_BGRID_FC 0
// B is not staggered (cell-centered)
#define MT_BGRID_CC 4
// B staggered the openggcm way: [-1..mx[
#define MT_BGRID_FC_GGCM 8

#define MT_SCONS_FC_GGCM (MT_FORMULATION_SCONS | MT_BGRID_FC_GGCM)
#define MT_SCONS_FC (MT_FORMULATION_SCONS | MT_BGRID_FC)
#define MT_FCONS_FC (MT_FORMULATION_FCONS | MT_BGRID_FC)
#define MT_FCONS_CC (MT_FORMULATION_FCONS | MT_BGRID_CC)

// the multi-moment schemes are cell-centered for all quantities
#define MT_GKEYLL (MT_FORMULATION_GKEYLL | MT_BGRID_CC)

#define MT_FORMULATION(mhd_type) ((mhd_type) & MT_FORMULATION_MASK)
#define MT_BGRID(mhd_type) ((mhd_type) & MT_BGRID_MASK)


// ----------------------------------------------------------------------
// wrappers / helpers

void ggcm_mhd_wrongful_death(struct ggcm_mhd *mhd, struct mrc_fld *x, int errcode);

int ggcm_mhd_main(int *argc, char ***argv);

// ----------------------------------------------------------------------

void ggcm_mhd_register();

#endif
