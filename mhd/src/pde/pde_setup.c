
#include <mrc_domain.h>
#include <mrc_bits.h>

#include <stdlib.h>
#include <math.h>

// ======================================================================
// PDE/mesh parameters that we keep around statically

static int s_n_ghosts; // number of ghost points
static int s_n_comps;  // number of components in state vector
static int s_n_dims;   // number of (not invariant) dimensions

// mesh info
static int s_size_1d;  // largest local dimension (including ghosts)
static int s_ldims[3]; // local dimensions (interior only) 
static int s_sw[3];    // number of ghost points per dim
static int s_dijk[3];  // set to 1 if actual direction, 0 if invariant

// need to have the static parameters above before we can include pde_fld1d.c

#include "pde_fld1d.c"

// ----------------------------------------------------------------------
// pde_setup

static void
pde_setup(struct mrc_fld *fld)
{
  s_n_ghosts = fld->_nr_ghosts;
  s_n_comps = mrc_fld_nr_comps(fld);

  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int n_dims = 3;
  if (gdims[2] == 1) {
    n_dims--;
    if (gdims[1] == 1) {
      n_dims--;
    }
  }
  s_n_dims = n_dims;

  s_size_1d = 0;
  for (int d = 0; d < 3; d++) {
    s_ldims[d] = mrc_fld_spatial_dims(fld)[d];
    s_sw[d] = mrc_fld_spatial_sw(fld)[d];
    s_size_1d = MAX(s_size_1d, s_ldims[d] + 2 * s_sw[d]);
    s_dijk[d] = (gdims[d] > 1);
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// pde_for_each_dir

#define pde_for_each_dir(dir)			\
      for (int dir = 0; dir < 3; dir++)		\
	if (s_sw[dir] > 0)

// ----------------------------------------------------------------------
// pde_for_each_line

#define pde_for_each_line(dir, j, k)					\
  int j, k, *_i1, *_i2, _i1b, _i1e, _i2b, _i2e;				\
  if (dir == 0) {							\
    _i1 = &j; _i2 = &k; _i1b = 0; _i2b = 0; _i1e = s_ldims[1]; _i2e = s_ldims[2]; \
  } else if (dir == 1) {						\
    _i1 = &k; _i2 = &j; _i1b = 0; _i2b = 0; _i1e = s_ldims[0]; _i2e = s_ldims[2]; \
  } else if (dir == 2) {						\
    _i1 = &j; _i2 = &k; _i1b = 0; _i2b = 0; _i1e = s_ldims[0]; _i2e = s_ldims[1]; \
  } else {								\
    assert(0);								\
  }									\
  for (*_i2 = _i2b; *_i2 < _i2e; (*_i2)++)				\
    for (*_i1 = _i1b; *_i1 < _i1e; (*_i1)++)


