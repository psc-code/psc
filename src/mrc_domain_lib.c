
#include "mrc_domain_private.h"

// ======================================================================
// mrc_sfc


// ----------------------------------------------------------------------
// bydim space filling curve

static void
sfc_bydim_setup(struct mrc_sfc *sfc, int *np)
{
  for (int d = 0; d < 3; d++) {
    sfc->np[d] = np[d];
  }
}

static int
sfc_bydim_idx3_to_idx(struct mrc_sfc *sfc, const int p[3])
{
  int *np = sfc->np;
  return (p[2] * np[1] + p[1]) * np[0] + p[0];
}

static void
sfc_bydim_idx_to_idx3(struct mrc_sfc *sfc, int idx, int p[3])
{
  int *np = sfc->np;
  p[0] = idx % np[0]; idx /= np[0];
  p[1] = idx % np[1]; idx /= np[1];
  p[2] = idx;
}

// ----------------------------------------------------------------------
// morton space filling curve

static void
sfc_morton_setup(struct mrc_sfc *sfc, int *np)
{
  int *nbits = sfc->nbits;
  for (int d = 0; d < 3; d++) {
    int n = np[d];
    nbits[d] = 0;
    while (n > 1) {
      n >>= 1;
      nbits[d]++;
    }
    // each dim must be power of 2
    assert(np[d] == 1 << nbits[d]);
  }

  sfc->nbits_max = nbits[0];
  if (nbits[1] > sfc->nbits_max) sfc->nbits_max = nbits[1];
  if (nbits[2] > sfc->nbits_max) sfc->nbits_max = nbits[2];
}

static int
sfc_morton_idx3_to_idx(struct mrc_sfc *sfc, const int p[3])
{
  int *nbits = sfc->nbits;
  int nbits_max = sfc->nbits_max;

  int pos = 0;
  int idx = 0;
  for (int b = 0; b < nbits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b >= nbits[d])
	continue;

      if (p[d] & (1 << b)) {
	idx |= (1 << pos);
      }
      pos++;
    }
  }

  return idx;
}

static void
sfc_morton_idx_to_idx3(struct mrc_sfc *sfc, int idx, int p[3])
{
  int *nbits = sfc->nbits;
  int nbits_max = sfc->nbits_max;

  for (int d = 0; d < 3; d++) {
    p[d] = 0;
  }

  int pos = 0;
  for (int b = 0; b < nbits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b >= nbits[d])
	continue;

      if (idx & (1 << pos)) {
	p[d] |= (1 << b);
      }
      pos++;
    }
  }
}

// ----------------------------------------------------------------------
// hilbert space filling curve

#include "hilbert.h"

static void
sfc_hilbert_setup(struct mrc_sfc *sfc, int *np)
{
  int *nbits = sfc->nbits;
  for (int d = 0; d < 3; d++) {
    int n = np[d];
    nbits[d] = 0;
    while (n > 1) {
      n >>= 1;
      nbits[d]++;
    }
    // each dim must be power of 2
    assert(np[d] == 1 << nbits[d]);
  }

  sfc->nbits_max = nbits[0];
  if (nbits[1] > sfc->nbits_max) sfc->nbits_max = nbits[1];
  if (nbits[2] > sfc->nbits_max) sfc->nbits_max = nbits[2];
  sfc->hilbert_nr_dims = 0;
  for (int d = 0; d < 3; d++) {
    if (nbits[d] == 0)
      continue;

    sfc->hilbert_dim[sfc->hilbert_nr_dims] = d;
    sfc->hilbert_nr_dims++;
  }
  // all not invariant dimensions must be equal
  int d0 = sfc->hilbert_dim[0];
  for (int i = 0; i < sfc->hilbert_nr_dims; i++) {
    int d = sfc->hilbert_dim[i];
    assert(nbits[d] == nbits[d0]);
  }
}

static int
sfc_hilbert_idx3_to_idx(struct mrc_sfc *sfc, const int p[3])
{
  if (sfc->hilbert_nr_dims == 0)
    return 0;

  int nbits_max = sfc->nbits_max;

  bitmask_t p_bm[3];
  for (int i = 0; i < sfc->hilbert_nr_dims; i++) {
    int d = sfc->hilbert_dim[i];
    p_bm[i] = p[d];
  }
  return hilbert_c2i(sfc->hilbert_nr_dims, nbits_max, p_bm);
}

static void
sfc_hilbert_idx_to_idx3(struct mrc_sfc *sfc, int idx, int p[3])
{
  int nbits_max = sfc->nbits_max;

  bitmask_t p_bm[3];
  hilbert_i2c(sfc->hilbert_nr_dims, nbits_max, idx, p_bm);
  for (int d = 0; d < 3; d++) {
    p[d] = 0;
  }
  for (int i = 0; i < sfc->hilbert_nr_dims; i++) {
    int d = sfc->hilbert_dim[i];
    p[d] = p_bm[i];
  }
}

// ----------------------------------------------------------------------
// space filling curve

void
sfc_setup(struct mrc_sfc *sfc, int *np)
{
  switch (sfc->curve_type) {
  case CURVE_BYDIM: sfc_bydim_setup(sfc, np); break;
  case CURVE_MORTON: sfc_morton_setup(sfc, np); break;
  case CURVE_HILBERT: sfc_hilbert_setup(sfc, np); break;
  default: assert(0);
  }
}

int
sfc_idx3_to_idx(struct mrc_sfc *sfc, const int p[3])
{
  switch (sfc->curve_type) {
  case CURVE_BYDIM: return sfc_bydim_idx3_to_idx(sfc, p);
  case CURVE_MORTON: return sfc_morton_idx3_to_idx(sfc, p);
  case CURVE_HILBERT: return sfc_hilbert_idx3_to_idx(sfc, p);
  default: assert(0);
  }
}

void
sfc_idx_to_idx3(struct mrc_sfc *sfc, int idx, int p[3])
{
  switch (sfc->curve_type) {
  case CURVE_BYDIM: sfc_bydim_idx_to_idx3(sfc, idx, p); break;
  case CURVE_MORTON: sfc_morton_idx_to_idx3(sfc, idx, p); break;
  case CURVE_HILBERT: sfc_hilbert_idx_to_idx3(sfc, idx, p); break;
  default: assert(0);
  }
}




