
#ifndef PSC_PHOTON_H
#define PSC_PHOTON_H

#include "psc.h"

typedef double photon_real_t;

#define MPI_PHOTONS_REAL MPI_DOUBLE

/// a single photon

typedef struct {
  photon_real_t x[3]; //< position
  photon_real_t p[3]; //< momentum
  photon_real_t wni;  //< weight
} photon_t;

/// a list of photons

typedef struct {
  photon_t *photons; //< array holding all the photons
  int nr; //< number of photons in the array 
  int nr_alloced; //< number of photons the array has space for
} photons_t;

/// photon lists per patch

typedef struct {
  photons_t *p;
} mphotons_t;

void photons_alloc(photons_t *pp, int n_part);
void photons_realloc(photons_t *pp, int new_n_part);
void photons_free(photons_t *pp);

void mphotons_alloc(mphotons_t *mphotons);
void mphotons_destroy(mphotons_t *mphotons);

static inline photon_t *
photons_get_one(photons_t *pp, int n)
{
  return &pp->photons[n];
}

#endif
