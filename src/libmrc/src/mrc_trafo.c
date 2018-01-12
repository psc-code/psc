
#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_trafo.h>
#include <mrc_crds.h>

#include <string.h>

// ======================================================================
#define mrc_trafo_ops(trafo) ((struct mrc_trafo_ops *)trafo->obj.ops)
#define XI0(crds, m, ix) MRC_DMCRDX(crds, m, ix)
#define XI1(crds, m, ix) MRC_DMCRDY(crds, m, ix)
#define XI2(crds, m, ix) MRC_DMCRDZ(crds, m, ix)
#define CE assert(ierr == 0);


// ======================================================================
// Trafo element setup functions
// These all look over the simulation domain and feed the logical
// coordinates to the subclass functions to get the data objects
// which map the logical coords into physical space.

// If subclass functions don't exist, it does it all using a finite difference
// approximation, taken from the old mrc_trafo_curvilinear.

// ----------------------------------------------------------------------
// SetupTrafoCoord - Setup actual coordinates

static void
SetupTrafoCoord(struct mrc_trafo *trafo)
{
  // every trafo must have, at the least, a calcCRD
  assert(mrc_trafo_ops(trafo)->calcCRD);
  struct mrc_domain *mb = trafo->_domain;
  struct mrc_fld **coord = trafo->_cc;
  struct mrc_crds *crds = mrc_domain_get_crds(mb);
  int sw = coord[0]->_nr_ghosts;

  mrc_fld_foreach_patch(coord[0], patch) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mb, patch, &info);
    mrc_fld_foreach(coord[0], jx, jy, jz, sw, sw) {
      double xi[3] = { XI0(crds, jx, patch), XI1(crds, jy, patch), XI2(crds, jz, patch) };
      double xx[3];
      mrc_trafo_ops(trafo)->calcCRD(trafo, info.p_block, xi, xx);
      for (int m = 0; m < 3; m++) {
	MRC_D5(coord[m], 0, jx,jy,jz, patch) = xx[m];
      }

    } mrc_fld_foreach_end;
  };
  
}


static int
LocalToLocal(struct mrc_fld *x)
{

  //FIXME: This used to handle boundary conditions through
  // M3d_FillGhosts (doesn't anymore)
  int bs = mrc_fld_nr_comps(x);
  struct mrc_domain *mb = x->_domain;
  mrc_ddc_fill_ghosts_fld(mrc_domain_get_ddc(mb), 0, bs, x);
  return(0);
}

// FIXME: These macros require 'crds' and 'patch' to be in scope,
// but they're purely internal so I don't care.
#define iDx(jx) (2./(XI0(crds, jx+1, patch)-XI0(crds, jx-1, patch)))
#define iDy(jy) (2./(XI1(crds, jy+1, patch)-XI1(crds, jy-1, patch)))
#define i2Dx(jx) (1./(XI0(crds, jx+1, patch)-XI0(crds, jx-1, patch)))
#define i2Dy(jy) (1./(XI1(crds, jy+1, patch)-XI1(crds, jy-1, patch)))

#define D0(x,m, jx,jy,jz,patch)					\
  ((MRC_D5(x,m, jx+1,jy,jz,patch) - MRC_D5(x,m, jx-1,jy,jz,patch))*i2Dx(jx))

#define D1(x, m, jx,jy,jz,patch)					\
  ((MRC_D5(x,m, jx,jy+1,jz,patch) - MRC_D5(x,m, jx,jy-1,jz,patch))*i2Dy(jy))

// ----------------------------------------------------------------------
// SetupTrafoJAC - Setup Jacobian of the mapping

static void
SetupTrafoJAC(struct mrc_trafo *trafo)
{
  struct mrc_domain *mb = trafo->_domain;
  struct mrc_crds *crds = mrc_domain_get_crds(mb);
  struct mrc_fld **coord = trafo->_cc;
  // check if we're forwarding to subclass, or doing it
  // the hard way.
  bool have_sub_func = false;
  // FIXME: if we have to calc the jac here we can't run over the
  // full ghost region, as we need to do finite diff.
  int sw_off = 1;
  if (mrc_trafo_ops(trafo)->calcJAC) { have_sub_func = true; sw_off = 0;}
  
  int sw = trafo->_jac->_nr_ghosts - sw_off;

  mrc_fld_foreach_patch(coord[0], patch) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mb, patch, &info);
    int block = info.p_block;
    mrc_fld_foreach(coord[0], jx, jy, jz, sw, sw) {
      if (have_sub_func) {
	double xi[3] = { XI0(crds, jx, patch), XI1(crds, jy, patch), XI2(crds, jz, patch) };
	mrc_trafo_ops(trafo)->calcJAC(trafo, block, xi, &TRAFO_JAC(trafo, jx,jy,jz, patch));
      } else {
	TRAFO_JAC(trafo, jx, jy, jz, patch) =
	  D0(coord[0], 0, jx,jy,jz,patch) * D1(coord[1], 0, jx,jy,jz,patch) -
	  D0(coord[1], 0, jx,jy,jz,patch) * D1(coord[0], 0, jx,jy,jz,patch);
      }
    } mrc_fld_foreach_end;
  }

  // the metric will need the ghost points filled in
  if (!have_sub_func) {
    assert(mrc_trafo_ops(trafo)->bc_jac);
#warning Trafo boundary conditions (bc_jac) not being handled
    int ierr = LocalToLocal(trafo->_jac); CE;
  }
}

// ----------------------------------------------------------------------
// SetupTrafoEL - Setup lower index version of .. something

static void
SetupTrafoEL(struct mrc_trafo *trafo)
{
  struct mrc_domain *mb = trafo->_domain;
  struct mrc_crds *crds = mrc_domain_get_crds(mb);
  struct mrc_fld **coord = trafo->_cc;

  // check if we're forwarding to subclass, or doing it
  // the hard way.
  bool have_sub_func = false;
  // FIXME: if we have to calc here we can't run over the
  // full ghost region, as we need to do finite diff.
  int sw_off = 1;
  if (mrc_trafo_ops(trafo)->calcEL) { have_sub_func = true; sw_off = 0;}
  
  int sw = trafo->_el->_nr_ghosts - sw_off;

  mrc_fld_foreach_patch(coord[0], patch) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mb, patch, &info);
    int block = info.p_block;
    mrc_fld_foreach(coord[0], jx, jy, jz, sw, sw) {

      double J = TRAFO_JAC(trafo, jx,jy,jz, patch);
      double el[3];
      double xi[3] = { XI0(crds, jx, patch), XI1(crds, jy, patch), XI2(crds, jz, patch) };
      for (int d = 0; d < 3; d++) {
	if (have_sub_func) {
	  mrc_trafo_ops(trafo)->calcEL(trafo, block, xi, d, el);
	} else {
	  el[0] = 1./J * D0(coord[d], 0, jx,jy,jz, patch);
	  el[1] = 1./J * D1(coord[d], 0, jx,jy,jz, patch);
	  el[2] = 1./J * (d == 2 ? 1. : 0.);
	}
	for (int m = 0; m < 3; m++) {
	  TRAFO_EL(trafo, m,d, jx,jy,jz, patch) = el[m];
	}
      }
      
    } mrc_fld_foreach_end;
  }

  // the metric will need ghost points
  if (!have_sub_func) {
    assert(mrc_trafo_ops(trafo)->bc_el);
#warning Trafo boundary conditions (bc_el) not being handled
    int ierr = LocalToLocal(trafo->_el); CE;
  }
}
  

// ----------------------------------------------------------------------
// SetupTrafoEU - Setup upper index version of .. something

static void
SetupTrafoEU(struct mrc_trafo *trafo)
{
  struct mrc_domain *mb = trafo->_domain;
  struct mrc_crds *crds = mrc_domain_get_crds(mb);
  struct mrc_fld **coord = trafo->_cc;

  // check if we're forwarding to subclass, or doing it
  // the hard way.
  bool have_sub_func = false;
  // FIXME: if we have to calc here we can't run over the
  // full ghost region, as we need to do finite diff.
  int sw_off = 1;
  if (mrc_trafo_ops(trafo)->calcEU) { have_sub_func = true; sw_off = 0;}
  
  int sw = trafo->_eu->_nr_ghosts - sw_off;
    
  mrc_fld_foreach_patch(coord[0], patch) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mb, patch, &info);
    int block = info.p_block;
    mrc_fld_foreach(coord[0], jx, jy, jz, sw, sw) {

      double eu[3];
      double xi[3] = { XI0(crds, jx, patch), XI1(crds, jy, patch), XI2(crds, jz, patch) };
      if (have_sub_func) {
	for (int d = 0; d < 3; d++) {
	  mrc_trafo_ops(trafo)->calcEU(trafo, block, xi, d, eu);
	  for (int m = 0; m < 3; m++) {
	    TRAFO_EU(trafo, m,d, jx,jy,jz, patch) = eu[m];
	  }
	}
      } else {
	double J = TRAFO_JAC(trafo, jx,jy,jz, patch);
	TRAFO_EU(trafo, 0,0, jx,jy,jz, patch) =   1./J * D1(coord[1], 0, jx,jy,jz, patch);
	TRAFO_EU(trafo, 0,1, jx,jy,jz, patch) = - 1./J * D1(coord[0], 0, jx,jy,jz, patch);
	TRAFO_EU(trafo, 0,2, jx,jy,jz, patch) =   0.;
	
	TRAFO_EU(trafo, 1,0, jx,jy,jz, patch) = - 1./J * D0(coord[1], 0, jx,jy,jz, patch);
	TRAFO_EU(trafo, 1,1, jx,jy,jz, patch) =   1./J * D0(coord[0], 0,  jx,jy,jz, patch);
	TRAFO_EU(trafo, 1,2, jx,jy,jz, patch) =   0.;
	
	TRAFO_EU(trafo, 2,0, jx,jy,jz, patch) =   0.;
	TRAFO_EU(trafo, 2,1, jx,jy,jz, patch) =   0.;
	TRAFO_EU(trafo, 2,2, jx,jy,jz, patch) =   1.;
      }

    } mrc_fld_foreach_end;
  };

  // the metric will need ghost points
  if (!have_sub_func) {
    assert(mrc_trafo_ops(trafo)->bc_eu);
#warning Trafo boundary conditions (bc_eu) not being handled
    int ierr = LocalToLocal(trafo->_eu); CE;
  }
}


// patch is LOCAL patch number
static inline double
__D_kl(int patch, int k, int l, struct mrc_fld *coord, int m,
       int jx, int jy, int jz)
{
  struct mrc_domain *domain = coord->_domain;
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  if (k == l) {
    if (k == 0) {
      return iDx(jx) *
	((MRC_D5(coord,m, jx+1,jy,jz,patch)-MRC_D5(coord, m, jx ,jy,jz,patch))/(XI0(crds,jx+1,patch)-XI0(crds,jx,patch)) -
	 (MRC_D5(coord,m, jx  ,jy,jz,patch)-MRC_D5(coord,m, jx-1,jy,jz,patch))/(XI0(crds,jx,patch)-XI0(crds,jx-1,patch)));
    } else if (k == 1) {
      return iDy(jy) *
	((MRC_D5(coord,m, jx,jy+1,jz,patch)-MRC_D5(coord,m, jx,jy ,jz,patch))/(XI1(crds,jy+1,patch)-XI1(crds,jy,patch)) -
	 (MRC_D5(coord,m, jx,jy  ,jz,patch)-MRC_D5(coord,m, jx,jy-1,jz,patch))/(XI1(crds,jy,patch)-XI1(crds,jy-1,patch)));
    } else if (k == 2) {
      return 0.;
    }
  }

  if ((k == 0 && l == 1) || (k == 1 && l == 0)) {
    return ((MRC_D5(coord,m, jx+1,jy+1,jz,patch) - MRC_D5(coord,m, jx+1,jy-1,jz,patch)) -
	    (MRC_D5(coord,m, jx-1,jy+1,jz,patch) - MRC_D5(coord,m, jx-1,jy-1,jz,patch)))
      *i2Dx(jx)*i2Dy(jy);
  } else {
    return 0.;
  }
}

#define D_kl(k,l, coord,m, jx,jy,jz) \
  (__D_kl(patch, k,l, coord,m, jx,jy,jz))

// ----------------------------------------------------------------------
// SetupTrafoGAM - Setup Christoffel symbols

static void
SetupTrafoGAM(struct mrc_trafo *trafo)
{
  struct mrc_domain *mb = trafo->_domain;
  struct mrc_crds *crds = mrc_domain_get_crds(mb);
  struct mrc_fld **coord = trafo->_cc;

  // check if we're forwarding to subclass, or doing it
  // the hard way.
  bool have_sub_func = false;
  // FIXME: if we have to calc here we can't run over the
  // full ghost region, as we need to do finite diff.
  int sw_off = 1;
  if (mrc_trafo_ops(trafo)->calcGAM) { have_sub_func = true; sw_off = 0;}
  
  int sw = trafo->_gam->_nr_ghosts - sw_off;

  mrc_fld_set(trafo->_gam, 0.);

  mrc_fld_foreach_patch(coord[0], patch) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mb, patch, &info);
    int block = info.p_block;
    mrc_fld_foreach(coord[0], jx, jy, jz, sw, sw) {
      // Christoffel symbols
      double xi[3] = { XI0(crds,jx,patch), XI1(crds,jy,patch), XI2(crds,jz,patch) };
      for (int i = 0; i < 3; i++) {
	for (int k = 0; k < 3; k++) {
	  for (int l = 0; l < 3; l++) {
	    if (have_sub_func) {
	      double gam;
	      mrc_trafo_ops(trafo)->calcGAM(trafo, block, xi, i,k,l, &gam);
	      TRAFO_GAM(trafo,i,k,l, jx,jy,jz,patch) = gam;
	    } else {
	      for (int d = 0; d < 3; d++) {
		TRAFO_GAM(trafo,i,k,l, jx,jy,jz,patch) += 
		  D_kl(k,l, coord[d], 0, jx,jy,jz) * TRAFO_EU(trafo,i,d, jx,jy,jz,patch);
	      }
	    }
	  }
	}
      }

    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// SetupTrafoMetric - setup elements of upper and lower metric tensors
// (apparently these are always calculated from the pre assigned values,
// rather than being assigned analytically. Wonder why?)

static void
SetupTrafoMetric(struct mrc_trafo *trafo)
{
  int sw = trafo->_jac->_nr_ghosts;
  
  mrc_fld_set(trafo->_guu, 0.);
  mrc_fld_set(trafo->_gll, 0.);
  
  mrc_fld_foreach_patch(trafo->_guu, patch) {
    mrc_fld_foreach(trafo->_guu, jx, jy, jz, sw, sw) {
      // metric tensors
      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  double J = TRAFO_JAC(trafo,jx,jy,jz,patch);
	  for (int d = 0; d < 3; d++) {
	    TRAFO_GUU(trafo,i,j, jx,jy,jz,patch) += 
	      J * TRAFO_EU(trafo,i,d, jx,jy,jz,patch) * TRAFO_EU(trafo,j,d, jx,jy,jz,patch);
	    TRAFO_GLL(trafo,i,j, jx,jy,jz,patch) += 
	      J * TRAFO_EL(trafo,i,d, jx,jy,jz,patch) * TRAFO_EL(trafo,j,d, jx,jy,jz,patch);
	  }
	}
      }
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// MB_CorrectChristoffel

#define delta_il(i, l) (i == l ? 1 : 0)
#define D_guu(k, m, i)							\
  (k == 0 ? (TRAFO_GUU(trafo,m,i, jx+1,jy,jz,patch) - TRAFO_GUU(trafo,m,i, jx-1,jy,jz,patch))*i2Dx(jx) : \
   k == 1 ? (TRAFO_GUU(trafo,m,i, jx,jy+1,jz,patch) - TRAFO_GUU(trafo,m,i, jx,jy-1,jz,patch))*i2Dy(jy) : 0 )

static int
MB_CorrectChristoffel(struct mrc_trafo *trafo)
{
  struct mrc_domain *mb = trafo->_domain;
  struct mrc_crds *crds = mrc_domain_get_crds(mb);

  // If GAM was setup from an analytic function, don't bother correcting it.
  if (mrc_trafo_ops(trafo)->calcGAM) return(0);

  int sw = trafo->_jac->_nr_ghosts;
  mrc_fld_foreach_patch(trafo->_jac, patch) {
    mrc_fld_foreach(trafo->_jac, jx, jy, jz, sw-1, sw-1) {
      double gam[3][3][3] = {};
      for (int i = 0; i < 3; i++) {
	for (int k = 0; k < 3; k++) {
	  for (int l = 0; l < 3; l++) {
	    for (int m = 0; m < 3; m++) {
	      gam[i][k][l] +=
		-TRAFO_GLL(trafo,l,m, jx,jy,jz,patch)*D_guu(k, m, i) +
		delta_il(i, l) * TRAFO_GAM(trafo,m,k,m, jx,jy,jz,patch);
	      for (int j = 0; j < 3; j++) {
		gam[i][k][l] += 
		  -TRAFO_GLL(trafo,l,m, jx,jy,jz,patch)
		  *TRAFO_GUU(trafo,i,j, jx,jy,jz,patch)
		  * TRAFO_GAM(trafo,m,k,j, jx,jy,jz,patch);
	      }
	    }
	  }
	}
      }
      for (int i = 0; i < 3; i++) {
	for (int k = 0; k < 3; k++) {
	  for (int l = 0; l < 3; l++) {
	    TRAFO_GAM(trafo,i,k,l, jx,jy,jz,patch) = gam[i][k][l];
	  }
	}
      }
    } mrc_fld_foreach_end;
  }
  return(0);
}

// ======================================================================
// LEGACY, commented out in original source, no idea what this does.

#define delta_il(i, l) (i == l ? 1 : 0)

#if 0
#define D_Jguu(k, m, i)							\
  (k == 0 ? (.5*(CL_JAC(jx+1,jy,jz)+CL_JAC(jx,jy,jz))*			\
	     .5*(CL_GUU(m,i, jx+1,jy,jz)+CL_GUU(m,i,jx,jy,jz)) -		\
	     .5*(CL_JAC(jx-1,jy,jz)+CL_JAC(jx,jy,jz))*			\
	     .5*(CL_GUU(m,i, jx-1,jy,jz)+CL_GUU(m,i, jx,jy,jz)))/(d0) :	\
   k == 1 ? (.5*(CL_JAC(jx,jy+1,jz)+CL_JAC(jx,jy,jz))*			\
	     .5*(CL_GUU(m,i, jx,jy+1,jz)+CL_GUU(m,i, jx,jy,jz)) -		\
	     .5*(CL_JAC(jx,jy-1,jz)+CL_JAC(jx,jy,jz))*			\
	     .5*(CL_GUU(m,i, jx,jy-1,jz)+CL_GUU(m,i, jx,jy,jz)))/(d1) : 0 )

#else

#define D_Jguu(k, m, i)							\
  (k == 0 ? (CL_JAC(jx+1,jy,jz) * CL_GUU(m,i, jx+1,jy,jz) -			\
	     CL_JAC(jx-1,jy,jz) * CL_GUU(m,i, jx-1,jy,jz))*i2Dx(jx) :		\
   k == 1 ? (CL_JAC(jx,jy+1,jz) * CL_GUU(m,i, jx,jy+1,jz) -			\
	     CL_JAC(jx,jy-1,jz) * CL_GUU(m,i, jx,jy-1,jz))*i2Dy(jy) : 0 )


#endif

#if 0
// not really tested, and actually using an averaged stencil (not for 
// momentum equation pressure)

static int
MB_CorrectChristoffelAlt(struct mrc_domain *mb)
{
  pfb;
  int sw = mb->mb_trafo->mbt_jac.info->mi_sw;
  MB_foreach_patch(mb, patch) {
    mrc_patch_foreach(patch, jx,jy,jz, sw-1,sw-1) {
      double gam[3][3][3] = {};
      for (int i = 0; i < 3; i++) {
	for (int k = 0; k < 3; k++) {
	  for (int l = 0; l < 3; l++) {
	    if (i == 1) {
	      for (int m = 0; m < 3; m++) {
		gam[i][k][l] +=
		  -CL_GLL(l,m, jx,jy,jz)/CL_JAC(jx,jy,jz) * D_Jguu(k, m, i) +
		  2.*delta_il(i, l) * CL_GAM(m,k,m, jx,jy,jz);
		for (int j = 0; j < 3; j++) {
		  gam[i][k][l] += 
		    -CL_GLL(l,m, jx,jy,jz)*CL_GUU(i,j, jx,jy,jz) * CL_GAM(m,j,k, jx,jy,jz);
		}
	      }
	    } else {
	      for (int m = 0; m < 3; m++) {
		gam[i][k][l] +=
		  -CL_GLL(l,m, jx,jy,jz)*D_guu(k, m, i) +
                delta_il(i, l) * CL_GAM(m,k,m, jx,jy,jz);
		for (int j = 0; j < 3; j++) {
		  gam[i][k][l] += 
		    -CL_GLL(l,m, jx,jy,jz)*CL_GUU(i,j, jx,jy,jz) * CL_GAM(m,k,j, jx,jy,jz);
		}
	      }
	    }
	  }
	}
      }
      for (int i = 0; i < 3; i++) {
	for (int k = 0; k < 3; k++) {
	  for (int l = 0; l < 3; l++) {
	    CL_GAM(i,k,l, jx,jy,jz) = gam[i][k][l];//.5*(gam[i][k][l]+gam[i][l][k]);
	  }
	}
      }
    } mrc_patch_foreach_end;
  } MB_foreach_patch_end;
  
  pfr;
}
#endif

static void
MB_GetFld(struct mrc_domain *domain, int bs, int sw, struct mrc_fld *fld, char *name)
{
  // FIXME: There's no reason these have to be double_aos, as far as I can
  // tell, but until I'm sure I'm going to leave them this way.
  mrc_fld_set_type(fld, "double");
  mrc_fld_set_param_bool(fld, "aos", true);
  mrc_fld_set_name(fld, name);
  mrc_fld_set_param_obj(fld, "domain", domain);
  mrc_fld_set_param_int(fld, "nr_spatial_dims", 3);
  mrc_fld_set_param_int(fld, "nr_comps", bs);
  mrc_fld_set_param_int(fld, "nr_ghosts", sw);
  // FIXME; strdup is screwing the pooch is comp_name is NULL, so 
  // I'm tossing in some null strings here
  for (int d = 0; d < bs; d++) {
    char s[5];
    sprintf(s, "%d", d);
    mrc_fld_set_comp_name(fld, d, s);
  }
  mrc_fld_setup(fld);
}

static void
_mrc_trafo_setup(struct mrc_trafo *trafo)
{
  int ierr;

  for (int d = 0; d < 3; d++) {
    char s[30];
    sprintf(s, "trafo_cc[%d]", d);
    MB_GetFld(trafo->_domain,  1, SW_2, trafo->_cc[d], s);
  }
  MB_GetFld(trafo->_domain,  1, SW_1, trafo->_jac, "trafo_jac");
  MB_GetFld(trafo->_domain,  9, SW_1, trafo->_el, "trafo_el");
  MB_GetFld(trafo->_domain,  9, SW_1, trafo->_eu, "trafo_eu");
  MB_GetFld(trafo->_domain,  9, SW_1, trafo->_gll, "trafo_gll");
  MB_GetFld(trafo->_domain,  9, SW_1, trafo->_guu, "trafo_guu");
  MB_GetFld(trafo->_domain, 27, SW_1, trafo->_gam, "trafo_gam");

  // actually fill in the data for the mapping fields
  SetupTrafoCoord(trafo);
  SetupTrafoJAC(trafo);
  SetupTrafoEL(trafo);
  SetupTrafoEU(trafo);
  SetupTrafoGAM(trafo);
  SetupTrafoMetric(trafo);

  ierr = MB_CorrectChristoffel(trafo); CE;
}


// ----------------------------------------------------------------------
// mrc_trafo_read
// FIXME: This doesn't so much read as it does regenerate. That's... not good.
static void
_mrc_trafo_read(struct mrc_trafo *trafo, struct mrc_io *io)
{
  // These don't get created during standard read, and we need them
  // to be. Since trafo elements have BC bits that get filled during
  // setup, we won't read them back.(FIXME: should trafos be able to exist without domains?)
  trafo->_cc[0] = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_cc[1] = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_cc[2] = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_jac = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_eu = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_el = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_guu = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_gll = mrc_fld_create(mrc_trafo_comm(trafo));
  trafo->_gam = mrc_fld_create(mrc_trafo_comm(trafo));

  mrc_trafo_setup(trafo);
}

// ======================================================================
// mrc_trafo_block_factory_type
// Get the subclass recommendation on which block factory to use

const char *
mrc_trafo_block_factory_type(struct mrc_trafo *trafo)
{
  return mrc_trafo_ops(trafo)->_block_factory;
}

// ======================================================================
// mrc_domain_init

extern struct mrc_trafo_ops mrc_trafo_cylindrical_ops;
extern struct mrc_trafo_ops mrc_trafo_cartesian_ops;

static void
_mrc_trafo_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_trafo, &mrc_trafo_cartesian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_trafo, &mrc_trafo_cylindrical_ops);
}

// ======================================================================
// mrc_trafo class

#define VAR(x) (void *)offsetof(struct mrc_trafo, x)
static struct param mrc_trafo_descr[] = {
  { "domain"          , VAR(_domain)      , PARAM_OBJ(mrc_domain) },
  { "cc0"             , VAR(_cc[0])       , MRC_VAR_OBJ(mrc_fld)   },
  { "cc1"             , VAR(_cc[1])       , MRC_VAR_OBJ(mrc_fld)   },
  { "cc2"             , VAR(_cc[2])       , MRC_VAR_OBJ(mrc_fld)   },
  //  { "nc"              , VAR(_nc)          , MRC_VAR_OBJ(mrc_fld)   },
  { "jac"             , VAR(_jac)         , MRC_VAR_OBJ(mrc_fld)   },
  { "eu"              , VAR(_eu)          , MRC_VAR_OBJ(mrc_fld)   },
  { "el"              , VAR(_el)          , MRC_VAR_OBJ(mrc_fld)   },
  { "guu"             , VAR(_guu)         , MRC_VAR_OBJ(mrc_fld)   },
  { "gll"             , VAR(_gll)         , MRC_VAR_OBJ(mrc_fld)   },
  { "gam"             , VAR(_gam)         , MRC_VAR_OBJ(mrc_fld)   },
  {},
};
#undef VAR


struct mrc_class_mrc_trafo mrc_class_mrc_trafo = {
  .name             = "mrc_trafo",
  .size             = sizeof(struct mrc_trafo),
  .init             = _mrc_trafo_init,
  .setup            = _mrc_trafo_setup,
  .read             = _mrc_trafo_read,
  .param_descr      = mrc_trafo_descr,
};


/* int */
/* __mrc_trafo_WriteH5(struct mrc_trafo *trafo, hid_t loc, const char *path) */
/* { */
/*   int ierr; */

/*   pfb; */
/*   hid_t group = MRC_H5Gcreate(loc, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); */
/*   ierr = MRC_H5Awrite_string(group, "type", mrc_trafo_type(trafo)); CE;  ierr = MRC_H5Gclose(group); CE; */
/*   pfr; */
/* } */


/* // ====================================================================== */
/* // mrc_trafo_ReadH5 */

/* int */
/* mrc_trafo_ReadH5(struct mrc_domain *mb, hid_t loc, const char *path, */
/* 		struct mrc_trafo **ptrafo) */
/* { */
/*   int ierr; */
  
/*   pfb; */
/*   hid_t group = MRC_H5Gopen(loc, path, H5P_DEFAULT); */
/*   char *type; */
/*   ierr = MRC_H5Aread_string(group, "type", &type); CE; */

/*   struct mrc_trafo *trafo = mrc_trafo_create(mrc_domain_comm(mb)); */
/*   mrc_trafo_set_type(trafo, type); */
/*   mrc_trafo_set_param_obj(trafo, "domain", mb); */

/*   free(type); */

/*   char *names[3] = { "mbt_cc0", "mbt_cc1", "mbt_cc2" }; */
/*   for (int d = 0; d < 3; d++) { */
/*     ierr = M3d_ReadH5(trafo->_cc[d] , group, names[d] , mb); CE; */
/*   } */
/*   ierr = M3d_ReadH5(trafo->_jac, group, "mbt_jac", mb); CE; */
/*   ierr = M3d_ReadH5(trafo->_eu , group, "mbt_eu" , mb); CE; */
/*   ierr = M3d_ReadH5(trafo->_el , group, "mbt_el" , mb); CE; */
/*   ierr = M3d_ReadH5(trafo->_guu, group, "mbt_guu", mb); CE; */
/*   ierr = M3d_ReadH5(trafo->_gll, group, "mbt_gll", mb); CE; */
/*   ierr = M3d_ReadH5(trafo->_gam, group, "mbt_gam", mb); CE; */

/*   ierr = MRC_H5Gclose(group); CE; */

/*   *ptrafo = trafo; */
/*   pfr; */
/* } */

/* // ====================================================================== */
/* // mrc_trafo_WriteH5 */

/* int */
/* mrc_trafo_WriteH5(struct mrc_trafo *trafo, hid_t loc, const char *name) */
/* { */
/*   // FIXME, put group here */
/*   int ierr = __mrc_trafo_WriteH5(trafo, loc, name); CE; */
/*   hid_t group = MRC_H5Gopen(loc, name, H5P_DEFAULT); */

/*   char *names[3] = { "mbt_cc0", "mbt_cc1", "mbt_cc2" }; */
/*   for (int d = 0; d < 3; d++) { */
/*     ierr = M3d_WriteH5(trafo->_cc[d] , group, names[d] ); CE; */
/*   } */
/*   ierr = M3d_WriteH5(trafo->_jac, group, "mbt_jac"); CE; */
/*   ierr = M3d_WriteH5(trafo->_eu , group, "mbt_eu" ); CE; */
/*   ierr = M3d_WriteH5(trafo->_el , group, "mbt_el" ); CE; */
/*   ierr = M3d_WriteH5(trafo->_guu, group, "mbt_guu"); CE; */
/*   ierr = M3d_WriteH5(trafo->_gll, group, "mbt_gll"); CE; */
/*   ierr = M3d_WriteH5(trafo->_gam, group, "mbt_gam"); CE; */
/*   ierr = MRC_H5Gclose(group); CE; */

/*   return ierr; */
/* } */

