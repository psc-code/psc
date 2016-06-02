
#ifndef PDE_DEFS_H
#define PDE_DEFS_H

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
#define OPT_LIMITER_GMINMOD  (2)

// --- OPT_RIEMANN ------------------------------------------------------

#define OPT_RIEMANN_RUSANOV   (1)
#define OPT_RIEMANN_HLL       (2)
#define OPT_RIEMANN_HLLC      (3)
#define OPT_RIEMANN_HLLD      (4)

// -- OPT_TIME_INTEGRATOR -----------------------------------------------

#define OPT_TIME_INTEGRATOR_EULER                   (1)
#define OPT_TIME_INTEGRATOR_PREDCORR                (2)

#endif
