
#include "psc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

// simple division into equal sized chunks / direction

static int
get_n_in_cell(real n)
{
  if (psc.prm.const_num_particles_per_cell) {
    return psc.prm.nicell;
  }
  return n / psc.coeff.cori + .5;
}

void
init_partition(int *n_part)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int p[3];
  int r = rank;
  p[0] = r % psc.domain.nproc[0]; r /= psc.domain.nproc[0];
  p[1] = r % psc.domain.nproc[1]; r /= psc.domain.nproc[1];
  p[2] = r;

  int m[3], n[3];
  for (int d = 0; d < 3; d++) {
    m[d] = psc.domain.ihi[d] - psc.domain.ilo[d];
    assert(m[d] % psc.domain.nproc[d] == 0);
    n[d] = m[d] / psc.domain.nproc[d];
    psc.ilo[d] = psc.domain.ilo[d] + p[d] * n[d];
    psc.ihi[d] = psc.domain.ilo[d] + (p[d] + 1) * n[d];
    int min_size = 1;
    if (p[d] == 0 && // left-most proc in this dir
	(psc.domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
      min_size += psc.pml.size;
    }
    if (p[d] == psc.domain.nproc[d] - 1 && // right-most proc in this dir
	(psc.domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
      min_size += psc.pml.size;
    }
    assert(psc.ihi[d] - psc.ilo[d] >= min_size);
  }

  MPI_Barrier(comm);
  for (int n = 0; n < size; n++) {
    if (rank == n) {
      printf("[%d] [%d:%d] x [%d:%d] x [%d:%d]\n", rank,
	     psc.ilo[0], psc.ihi[0],
	     psc.ilo[1], psc.ihi[1],
	     psc.ilo[2], psc.ihi[2]);
    }
    MPI_Barrier(comm);
  }

  int np = 0;
  for (int kind = 0; kind < 2; kind++) {
    for (int jz = psc.ilo[2]; jz < psc.ihi[2]; jz++) {
      for (int jy = psc.ilo[1]; jy < psc.ihi[1]; jy++) {
	for (int jx = psc.ilo[0]; jx < psc.ihi[0]; jx++) {
	  double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2];
	  double xx[3] = { jx * dx, jy * dy, jz * dz };
	  double q, m, n, v[3], T[3];

	  psc.Case->ops->init_nvt(kind, xx, &q, &m, &n, v, T);
	  int n_in_cell = get_n_in_cell(n);
	  np += n_in_cell;
	}
      }
    }
  }

  *n_part = np;
}

void
init_particles()
{
  double beta = psc.coeff.beta;
  int i = 0;

  for (int kind = 0; kind < 2; kind++) {
    for (int jz = psc.ilo[2]; jz < psc.ihi[2]; jz++) {
      for (int jy = psc.ilo[1]; jy < psc.ihi[1]; jy++) {
	for (int jx = psc.ilo[0]; jx < psc.ihi[0]; jx++) {
	  double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2];
	  double xx[3] = { jx * dx, jy * dy, jz * dz };
	  double q, m, n, v[3], T[3];
	  psc.Case->ops->init_nvt(kind, xx, &q, &m, &n, v, T);

	  int n_in_cell = get_n_in_cell(n);
	  for (int cnt = 0; cnt < n_in_cell; cnt++) {
	    struct f_particle *p = &psc.f_part[i++];

	    // FIXME? this gives same random numbers on all procs
	    float ran1 = random() / ((float) RAND_MAX + 1);
	    float ran2 = random() / ((float) RAND_MAX + 1);
	    float ran3 = random() / ((float) RAND_MAX + 1);
	    float ran4 = random() / ((float) RAND_MAX + 1);
	    float ran5 = random() / ((float) RAND_MAX + 1);
	    float ran6 = random() / ((float) RAND_MAX + 1);
	    
	    float px =
	      sqrtf(-2.f*T[0]/m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2)
	      + v[0];
	    float py =
	      sqrtf(-2.f*T[1]/m*sqr(beta)*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4)
	      + v[1];
	    float pz =
  	      sqrtf(-2.f*T[2]/m*sqr(beta)*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6)
	      + v[2];
	    
	    p->xi = xx[0];
	    p->yi = xx[1];
	    p->zi = xx[2];
	    p->pxi = px;
	    p->pyi = py;
	    p->pzi = pz;
	    p->qni = q;
	    p->mni = m;
	    p->lni = -1; // FIXME?
	    if (psc.prm.fortran_particle_weight_hack) {
	      p->wni = n;
	    } else {
	      p->wni = n / (n_in_cell * psc.coeff.cori);
	    }
	  }
	}
      }
    }
  }
  psc.n_part = i;
}

// ======================================================================
// Fortran glue

#define C_init_partition_F77 F77_FUNC_(c_init_partition, C_INIT_PARTITION)
#define C_init_particles_F77 F77_FUNC_(c_init_particles, C_INIT_PARTICLES)
#define SET_subdomain_F77 F77_FUNC_(set_subdomain, SET_SUBDOMAIN)
#define INIT_partition_F77 F77_FUNC_(init_partition, INIT_PARTITION)
#define INIT_idistr_F77 F77_FUNC_(init_idistr, INIT_IDISTR)
#define GET_niloc_F77 F77_FUNC_(get_niloc, GET_NILOC)

void SET_subdomain_F77(f_int *i1mn, f_int *i1mx, f_int *i2mn, f_int *i2mx,
		       f_int *i3mn, f_int *i3mx);
void INIT_partition_F77(f_int *part_label_off, f_int *rd1n, f_int *rd1x,
			f_int *rd2n, f_int *rd2x, f_int *rd3n, f_int *rd3x,
			f_int *niloc_new);
void INIT_idistr_F77(f_int *part_label_off, f_int *rd1n, f_int *rd1x,
		     f_int *rd2n, f_int *rd2x, f_int *rd3n, f_int *rd3x);
void GET_niloc_F77(f_int *niloc);

static int rd1n, rd1x, rd2n, rd2x, rd3n, rd3x;
static int part_label_offset;

void
C_init_partition_F77(f_int *part_label_off, f_int *n_part)
{
  if (!psc.Case) {
    INIT_partition_F77(&part_label_offset, &rd1n, &rd1x, &rd2n, &rd2x, &rd3n, &rd3x,
		       n_part);
    return;
  }

  init_partition(n_part);
  *part_label_off = -1; // not supported
  f_int i1mn = psc.ilo[0];
  f_int i2mn = psc.ilo[1];
  f_int i3mn = psc.ilo[2];
  f_int i1mx = psc.ihi[0] - 1;
  f_int i2mx = psc.ihi[1] - 1;
  f_int i3mx = psc.ihi[2] - 1;
  SET_subdomain_F77(&i1mn, &i1mx, &i2mn, &i2mx, &i3mn, &i3mx);
}

void
C_init_particles_F77()
{
  if (!psc.Case) {
    INIT_idistr_F77(&part_label_offset, &rd1n, &rd1x, &rd2n, &rd2x, &rd3n, &rd3x);
    GET_niloc_F77(&psc.n_part);
    return;
  }

  init_particles();
}
