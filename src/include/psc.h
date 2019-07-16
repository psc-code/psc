/*! @file */

#ifndef PSC_H
#define PSC_H

#include "PscConfig.h"

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

// ----------------------------------------------------------------------
// we keep this info global for now.

extern int pr_time_step_no_comm;

void psc_init(int& argc, char**& argv);

#endif
