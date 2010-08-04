
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
psc_init_partition(int *n_part, int *particle_label_offset)
{
  if (!psc.Case) {
    INIT_partition(n_part);
    return;
  }

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  assert(psc.domain.nproc[0] * psc.domain.nproc[1] * psc.domain.nproc[2] == size);

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
    psc.ibn[d] = psc.domain.nghost[d];
    psc.ilg[d] = psc.ilo[d] - psc.ibn[d];
    psc.ihg[d] = psc.ihi[d] + psc.ibn[d];
    psc.img[d] = psc.ihg[d] - psc.ilg[d];
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
  psc.fld_size = psc.img[0] * psc.img[1] * psc.img[2];

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
  if (psc.Case->ops->init_npt) {
    for (int kind = 0; kind < 2; kind++) {
      for (int jz = psc.ilo[2]; jz < psc.ihi[2]; jz++) {
	for (int jy = psc.ilo[1]; jy < psc.ihi[1]; jy++) {
	  for (int jx = psc.ilo[0]; jx < psc.ihi[0]; jx++) {
	    double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2];
	    double xx[3] = { jx * dx, jy * dy, jz * dz };
	    struct psc_particle_npt npt = { // init to all zero
	    };
	    psc_case_init_npt(psc.Case, kind, xx, &npt);
	    
	    int n_in_cell = get_n_in_cell(npt.n);
	    np += n_in_cell;
	  }
	}
      }
    }
  }

  *n_part = np;
  // calculate global particle label offset for unique numbering
  MPI_Exscan(n_part, particle_label_offset, 1, MPI_INT, MPI_SUM,
	   MPI_COMM_WORLD);
}

void INIT_idistr();

void
psc_init_particles(int particle_label_offset)
{
  if (!psc.Case) {
    INIT_idistr();
    GET_niloc(&psc.pp.n_part);
    return;
  }

  double beta = psc.coeff.beta;
  int i = 0;

  if (psc.Case->ops->init_npt) {
    for (int kind = 0; kind < 2; kind++) {
      for (int jz = psc.ilo[2]; jz < psc.ihi[2]; jz++) {
	for (int jy = psc.ilo[1]; jy < psc.ihi[1]; jy++) {
	  for (int jx = psc.ilo[0]; jx < psc.ihi[0]; jx++) {
	    double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2];
	    double xx[3] = { jx * dx, jy * dy, jz * dz };
	    struct psc_particle_npt npt = { // init to all zero
	    };
	    psc_case_init_npt(psc.Case, kind, xx, &npt);
	    
	    int n_in_cell = get_n_in_cell(npt.n);
	    for (int cnt = 0; cnt < n_in_cell; cnt++) {
	      particle_base_t *p = psc_particles_base_get_one(&psc.pp, i++);
	      
	      // FIXME? this gives same random numbers on all procs
	      float ran1, ran2, ran3, ran4, ran5, ran6;
	      do {
		ran1 = random() / ((float) RAND_MAX + 1);
		ran2 = random() / ((float) RAND_MAX + 1);
		ran3 = random() / ((float) RAND_MAX + 1);
		ran4 = random() / ((float) RAND_MAX + 1);
		ran5 = random() / ((float) RAND_MAX + 1);
		ran6 = random() / ((float) RAND_MAX + 1);
	      } while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f ||
		       ran4 >= 1.f || ran5 >= 1.f || ran6 >= 1.f);
	      
	      float px =
		sqrtf(-2.f*npt.T[0]/npt.m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2)
		+ npt.p[0];
	      float py =
		sqrtf(-2.f*npt.T[1]/npt.m*sqr(beta)*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4)
		+ npt.p[1];
	      float pz =
		sqrtf(-2.f*npt.T[2]/npt.m*sqr(beta)*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6)
		+ npt.p[2];
	      
	      p->xi = xx[0];
	      p->yi = xx[1];
	      p->zi = xx[2];
	      p->pxi = px;
	      p->pyi = py;
	      p->pzi = pz;
	      p->qni = npt.q;
	      p->mni = npt.m;
#if PARTICLES_BASE == PARTICLES_FORTRAN
	      p->lni = particle_label_offset + 1;
#endif
	      if (psc.prm.fortran_particle_weight_hack) {
		p->wni = npt.n;
	      } else {
		p->wni = npt.n / (n_in_cell * psc.coeff.cori);
	      }
	    }
	  }
	}
      }
    }
  }
  psc_set_n_particles(i);
}
