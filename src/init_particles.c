#include "psc.h"
#include "psc_case_private.h"
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
pml_find_bounds(int p, int ilo[3], int ihi[3])
{
  struct psc_patch *patch = &psc.patch[p];
  for (int d = 0; d < 3; d++) {
    ilo[d] = 0;
    if (patch->off[d] == 0 && // left-most proc in this dir
	(psc.domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
      ilo[d] += psc.pml.size+1;
    }
    ihi[d] = patch->ldims[d];
    if (ihi[d] + patch->off[d] == psc.domain.gdims[d] && // right-most proc in this dir
	(psc.domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	 psc.domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
      ihi[d] -= psc.pml.size+1;
    }
  }
}

void
psc_init_partition(int *particle_label_offset)
{
  assert(_psc_case->Case);

  MPI_Comm comm = MPI_COMM_WORLD;

  psc.mrc_domain = mrc_domain_create(comm);
  // create a very simple domain decomposition
  int bc[3] = {};
  for (int d = 0; d < 3; d++) {
    if (psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC &&
	psc.domain.gdims[d] > 1) {
      bc[d] = BC_PERIODIC;
    }
  }

  mrc_domain_set_type(psc.mrc_domain, "multi");
  mrc_domain_set_param_int3(psc.mrc_domain, "m", psc.domain.gdims);
  mrc_domain_set_param_int(psc.mrc_domain, "bcx", bc[0]);
  mrc_domain_set_param_int(psc.mrc_domain, "bcy", bc[1]);
  mrc_domain_set_param_int(psc.mrc_domain, "bcz", bc[2]);

  struct mrc_crds *crds = mrc_domain_get_crds(psc.mrc_domain);
  mrc_crds_set_type(crds, "multi_uniform");
  mrc_crds_set_param_int(crds, "sw", 2);
  mrc_crds_set_param_float3(crds, "h",  (float[3]) { psc.domain.length[0],
	psc.domain.length[1], psc.domain.length[2] });

  mrc_domain_set_from_options(psc.mrc_domain);
  mrc_domain_setup(psc.mrc_domain);
  mrc_domain_view(psc.mrc_domain);

  // set up index bounds,
  // sanity checks for the decomposed domain
  int gdims[3];
  mrc_domain_get_global_dims(psc.mrc_domain, gdims);
  struct mrc_patch *patches = mrc_domain_get_patches(psc.mrc_domain, &psc.nr_patches);
  psc.patch = calloc(psc.nr_patches, sizeof(*psc.patch));
  foreach_patch(p) {
    struct psc_patch *patch = &psc.patch[p];
    for (int d = 0; d < 3; d++) {
      patch->ldims[d] = patches[p].ldims[d];
      patch->off[d] = patches[p].off[d];
      patch->xb[d]  = patches[p].off[d] * psc.dx[d];
      
      int min_size = 1;
      if (patch->off[d] == 0 && // left-most patch in this dir
	  (psc.domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	   psc.domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
	min_size += psc.pml.size;
      }
      if (patch->off[d] + patch->ldims[d] == gdims[d] && // right-most patch in this dir
	  (psc.domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	   psc.domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
	min_size += psc.pml.size;
      }
      assert(psc.patch[p].ldims[d] >= min_size);
    }
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

  int np_total = 0;
  psc.particles.p = calloc(psc.nr_patches, sizeof(*psc.particles.p));
  foreach_patch(p) {
    int ilo[3], ihi[3];
    pml_find_bounds(p, ilo, ihi);

    int np = 0;
    if (_psc_case->Case->ops->init_npt) {
      for (int kind = 0; kind < 2; kind++) {
	for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	  for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	    for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	      double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
	      struct psc_particle_npt npt = { // init to all zero
	      };
	      psc_case_init_npt(_psc_case->Case, kind, xx, &npt);
	      
	      int n_in_cell = get_n_in_cell(npt.n);
	      np += n_in_cell;
	    }
	  }
	}
      }
    }

    particles_base_alloc(&psc.particles.p[p], np);
    np_total += np;
  }
  // calculate global particle label offset for unique numbering
  *particle_label_offset = 0; // necessary on proc 0
  MPI_Exscan(&np_total, particle_label_offset, 1, MPI_INT, MPI_SUM,
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

  foreach_patch(p) {
    int ilo[3], ihi[3];
    pml_find_bounds(p, ilo, ihi);
    particles_base_t *pp = &psc.particles.p[p];

    int i = 0;
    if (_psc_case->Case->ops->init_npt) {
      for (int kind = 0; kind < 2; kind++) {
	for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	  for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	    for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	      double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
	      struct psc_particle_npt npt = { // init to all zero
	      };
	      psc_case_init_npt(_psc_case->Case, kind, xx, &npt);
	      
	      int n_in_cell = get_n_in_cell(npt.n);
	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		particle_base_t *p = particles_base_get_one(pp, i++);
		
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
    pp->n_part = i;
  }
}
