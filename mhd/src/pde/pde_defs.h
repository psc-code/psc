
#ifndef PDE_DEFS_H
#define PDE_DEFS_H

// ======================================================================
// compile time options

#define OPT_FLD1D_MRC_FLD   (1)
#define OPT_FLD1D_C_ARRAY   (2)
#define OPT_FLD1D_PTR_ARRAY (3)

#define OPT_STAGGER_REG     (0)
#define OPT_STAGGER_GGCM    (1)

#define OPT_TMP_NORMAL      (0)
#define OPT_TMP_COMPAT      (1)

#define OPT_GGCM_CRDS_NORMAL (0)
#define OPT_GGCM_CRDS_LEGACY (1)

// ======================================================================
// MHD options

// ----------------------------------------------------------------------
// these option can potentially be set at run time

// --- OPT_EQN ----------------------------------------------------------

#define OPT_EQN_MHD_FCONS     (1) // fully-conservative MHD
#define OPT_EQN_MHD_SCONS     (2) // semi-conservative MHD (w/o pressure in flux)
#define OPT_EQN_HD            (3) // hydrodynamics

// --- OPT_LIMITER ------------------------------------------------------

#define OPT_LIMITER_FLAT     (1)
#define OPT_LIMITER_MINMOD   (2)
#define OPT_LIMITER_MC       (3)
#define OPT_LIMITER_GMINMOD  (4)

// --- OPT_RESISTIVITY --------------------------------------------------

#define OPT_RESISTIVITY_NONE  (0)
#define OPT_RESISTIVITY_CONST (1)

// --- OPT_HALL ---------------------------------------------------------

#define OPT_HALL_NONE  (0)
#define OPT_HALL_CONST (1)
#define OPT_HALL_YES   (2)

// --- OPT_RIEMANN ------------------------------------------------------

#define OPT_RIEMANN_RUSANOV   (1)
#define OPT_RIEMANN_HLL       (2)
#define OPT_RIEMANN_HLLC      (3)
#define OPT_RIEMANN_HLLD      (4)

// --- OPT_DIVB ---------------------------------------------------------

#define OPT_DIVB_NONE         (0)
#define OPT_DIVB_GLM          (1)

// -- OPT_TIME_INTEGRATOR -----------------------------------------------

#define OPT_TIME_INTEGRATOR_EULER                   (1)
#define OPT_TIME_INTEGRATOR_PREDCORR                (2)
#define OPT_TIME_INTEGRATOR_TVD_RK2                 (3)

// -- OPT_BACKGROUND ----------------------------------------------------
// just true/false

// -- OPT_BC_RECONSTRUCT ------------------------------------------------
// just true/false

// -- OPT_GET_DT --------------------------------------------------------

#define OPT_GET_DT_MHD_GGCM       (1)
#define OPT_GET_DT_MHD            (2)
#define OPT_GET_DT_MHD_CT         (3)
#define OPT_GET_DT_HD             (4)

// -- OPT_MHD_PRIMVAR ---------------------------------------------------
// -- and many others
//
// in general, "MHD_C" is a straight translation into C, while "MHD_C2", 
// if it exists, maybe equivalent but rewritten, and may have small
// differences

#define OPT_MHD_C       (1)
#define OPT_MHD_FORTRAN (2)
#define OPT_MHD_C_V2    (3)

#endif
