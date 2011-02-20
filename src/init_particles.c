#include "psc.h"
#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <time.h>

static struct param seed_by_time_descr[] = {
  { "seed_by_time"          , 0       , PARAM_BOOL(false)   },
  {},
};

static int
get_n_in_cell(real n)
{
  if (psc.prm.const_num_particles_per_cell) {
    return psc.prm.nicell;
  }
  return n / psc.coeff.cori + .5;
}

// particles must not be placed in the pml regions
// check if pml bnds are set and restrict avoid particle placement inside pml regions

static void
pml_find_bounds(int ilo[3], int ihi[3])
{
  struct psc_patch *patch = &psc.patch[0];
  for (int d = 0; d < 3; d++) {
    ilo[d] = patch->off[d];
    if (ilo[d] == 0 && // left-most proc in this dir
	(psc.domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
      ilo[d] += psc.pml.size+1;
    }
    ihi[d] = patch->off[d] + patch->ldims[d];
    if (ihi[d] == psc.domain.gdims[d] && // right-most proc in this dir
	(psc.domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
      ihi[d] -= psc.pml.size+1;
    }
  }
}

void
psc_init_partition(int *n_part, int *particle_label_offset)
{
  assert(psc.Case);

  MPI_Comm comm = MPI_COMM_WORLD;

  // create a very simple domain decomposition
  struct mrc_domain_simple_params domain_par = {};
  for (int d = 0; d < 3; d++) {
    int gdim = psc.domain.gdims[d];
    assert(gdim % psc.domain.nproc[d] == 0);
    domain_par.ldims[d] = gdim / psc.domain.nproc[d];
    domain_par.nr_procs[d] = psc.domain.nproc[d];

    if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC &&
	psc.domain.gdims[d] > 1) {
      domain_par.bc[d] = BC_PERIODIC;
    }
  }

  psc.mrc_domain = mrc_domain_create(comm);
  mrc_domain_set_type(psc.mrc_domain, "simple");
  mrc_domain_simple_set_params(psc.mrc_domain, &domain_par);

  struct mrc_crds *crds = mrc_domain_get_crds(psc.mrc_domain);
  mrc_crds_set_param_int(crds, "sw", 2);
  mrc_crds_set_param_float(crds, "xh",  psc.domain.length[0]);
  mrc_crds_set_param_float(crds, "yh",  psc.domain.length[1]);
  mrc_crds_set_param_float(crds, "zh",  psc.domain.length[2]);

  mrc_domain_setup(psc.mrc_domain);
  mrc_domain_view(psc.mrc_domain);

  // set up index bounds,
  // sanity checks for the decomposed domain
  int off[3], ldims[3], lidx[3];
  mrc_domain_get_local_offset_dims(psc.mrc_domain, off, ldims);
  mrc_domain_get_local_idx(psc.mrc_domain, lidx);
  psc.nr_patches = 1;
  for (int d = 0; d < 3; d++) {
    psc.patch[0].ldims[d] = ldims[d];
    psc.patch[0].off[d] = off[d];
    psc.patch[0].xb[d]  = off[d] * psc.dx[d];

    int min_size = 1;
    if (lidx[d] == 0 && // left-most proc in this dir
	(psc.domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
      min_size += psc.pml.size;
    }
    if (lidx[d] == psc.domain.nproc[d] - 1 && // right-most proc in this dir
	(psc.domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
      min_size += psc.pml.size;
    }
    assert(psc.patch[0].ldims[d] >= min_size);
  }

#if 0
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
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
#endif

  int ilo[3], ihi[3];
  pml_find_bounds(ilo, ihi);

  int np = 0;
  if (psc.Case->ops->init_npt) {
    for (int kind = 0; kind < 2; kind++) {
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
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
  *particle_label_offset = 0; // necessary on proc 0
  MPI_Exscan(n_part, particle_label_offset, 1, MPI_INT, MPI_SUM,
	     MPI_COMM_WORLD);
}

void
psc_init_particles(int particle_label_offset)
{
  double beta = psc.coeff.beta;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // if seed_by_time is set, different random numbers are generated in each run. 
  bool seed_by_time;
  mrc_params_parse(&seed_by_time, seed_by_time_descr, "seed by time", MPI_COMM_WORLD);
  mrc_params_print(&seed_by_time, seed_by_time_descr, "seed by time", MPI_COMM_WORLD);
  
  if (seed_by_time) {
    srandom(10*rank + time(NULL));
  } else {
    srandom(rank);
  }

  int ilo[3], ihi[3];
  pml_find_bounds(ilo, ihi);

  int i = 0;
  if (psc.Case->ops->init_npt) {
    for (int kind = 0; kind < 2; kind++) {
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {

	    double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2];
	    double xx[3] = { jx * dx, jy * dy, jz * dz };
	    struct psc_particle_npt npt = { // init to all zero
	    };
	    psc_case_init_npt(psc.Case, kind, xx, &npt);
	    
	    int n_in_cell = get_n_in_cell(npt.n);
	    for (int cnt = 0; cnt < n_in_cell; cnt++) {
	      particle_base_t *p = particles_base_get_one(&psc.pp, i++);
	      
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
