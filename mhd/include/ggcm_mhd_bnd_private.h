
#ifndef GGCM_MHD_BND_PRIVATE_H
#define GGCM_MHD_BND_PRIVATE_H

#include "ggcm_mhd_bnd.h"

struct ggcm_mhd_bnd {
  struct mrc_obj obj;
  struct ggcm_mhd *mhd;
};

struct ggcm_mhd_bnd_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_bnd);
  void (*fill_ghosts)(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
		      float bntim);
  void (*fill_ghosts_E)(struct ggcm_mhd_bnd *bnd, struct mrc_fld *E);
  void (*fill_ghosts_reconstr)(struct ggcm_mhd_bnd *bnd, struct mrc_fld *U_l[],
			       struct mrc_fld *U_r[], int p);
};

extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_none;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_ggcm_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_sc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_fc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_fc_cc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow_gkeyll;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_sc_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_fc_float;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_sc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_sc_ggcm_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_fc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_fc_cc_double;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere_gkeyll;

#define ggcm_mhd_bnd_ops(bnd) ((struct ggcm_mhd_bnd_ops *)((bnd)->obj.ops))

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map
//
// infrastructure for internal spherical boundaries

struct ggcm_mhd_bnd_sphere_map {
  struct ggcm_mhd *mhd;
  double radius;
  double min_dr;
  double dr; // step size to use when determining r1
  double extra_dr; // extra amount (* min_dr) to subtract from final r1

  // the spherical shell where we set ghost points is r1 <= r <= radius
  double r1;

  // maps
  // for managing cell-centered ghost points
  int cc_n_map;
  struct mrc_fld *cc_imap;  // ghost cell # -> (ix,iy,iz,p)

  // for managing edge-centered ghost points
  int ec_n_map[3];
  struct mrc_fld *ec_imap[3];

  // for managing face-centered boundary
  int fc_n_map[3];
  struct mrc_fld *fc_imap[3];
};

void ggcm_mhd_bnd_sphere_map_setup(struct ggcm_mhd_bnd_sphere_map *map,
				   struct ggcm_mhd *mhd, double radius,
				   double dr, double extra_dr);
void ggcm_mhd_bnd_sphere_map_setup_flds(struct ggcm_mhd_bnd_sphere_map *map);
void ggcm_mhd_bnd_sphere_map_setup_cc(struct ggcm_mhd_bnd_sphere_map *map);
void ggcm_mhd_bnd_sphere_map_setup_ec(struct ggcm_mhd_bnd_sphere_map *map);
void ggcm_mhd_bnd_sphere_map_setup_fc(struct ggcm_mhd_bnd_sphere_map *map);


#endif
