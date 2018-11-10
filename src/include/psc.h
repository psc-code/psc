/*! @file */

#ifndef PSC_H
#define PSC_H

#include "psc_config.h"

#include "psc_bits.h"
#include "psc_particles.h"
#include "mrc_domain.hxx"
#include "grid.hxx"
#include "particles.hxx"

#include "psc_stats.h"

#include <mrc_domain.h>

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <vector>

// ----------------------------------------------------------------------

enum {
  JXI, JYI, JZI,
  EX , EY , EZ ,
  HX , HY , HZ ,
  NR_FIELDS,
};

// ----------------------------------------------------------------------
// general info / parameters for the code

///Default kinds (electrons + ions)
enum {
  KIND_ELECTRON,
  KIND_ION,
  NR_KINDS,
};

struct psc_particle_npt {
  int kind; ///< particle kind
  double q; ///< charge
  double m; ///< mass
  double n; ///< density
  double p[3]; ///< momentum
  double T[3]; ///< temperature
};

// ----------------------------------------------------------------------
// we keep this info global for now.

extern int pr_time_step_no_comm;

void psc_init(int& argc, char**& argv);

#endif
