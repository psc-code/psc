#include "psc.h"
#include "psc_case_private.h"
#include "psc_particles_as_c.h"

#include <mrc_params.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <time.h>

// FIXME, same code exists, essentially, in psc.c, so this should go away
// (eventually, with the cases)

static inline int
get_n_in_cell(struct psc *psc, struct psc_particle_npt *npt)
{
  if (psc->prm.const_num_particles_per_cell) {
    return psc->prm.nicell;
  }
  if (npt->particles_per_cell) {
    return npt->n * npt->particles_per_cell + .5;
  }
  return npt->n / psc->coeff.cori + .5;
}

// particles must not be placed in the pml regions
// check if pml bnds are set and restrict avoid particle placement inside pml regions

static void
pml_find_bounds(struct psc *psc, int p, int ilo[3], int ihi[3])
{
  struct psc_patch *patch = &psc->patch[p];
  for (int d = 0; d < 3; d++) {
    ilo[d] = 0;
    if (patch->off[d] == 0 && // left-most proc in this dir
	(psc->domain.bnd_fld_lo[d] == BND_FLD_UPML || 
	 psc->domain.bnd_fld_lo[d] == BND_FLD_TIME)) {
      ilo[d] += psc->pml.size+1;
    }
    ihi[d] = patch->ldims[d];
    if (ihi[d] + patch->off[d] == psc->domain.gdims[d] && // right-most proc in this dir
	(psc->domain.bnd_fld_hi[d] == BND_FLD_UPML || 
	 psc->domain.bnd_fld_hi[d] == BND_FLD_TIME)) {
      ihi[d] -= psc->pml.size+1;
    }
  }
}

int
psc_case_calc_nr_particles_in_patch(struct psc_case *_case, int p)
{
  struct psc *psc = _case->psc;
  int ilo[3], ihi[3];
  pml_find_bounds(psc, p, ilo, ihi);

  int np = 0;
  for (int kind = 0; kind < _case->nr_kinds; kind++) {
    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	  double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
	  struct psc_particle_npt npt = {}; // init to all zero
	  psc_case_init_npt(_case, kind, xx, &npt);
	  
	  int n_in_cell = get_n_in_cell(psc, &npt);
	  np += n_in_cell;
	}
      }
    }
  }
  return np;
}

void
psc_case_init_partition(struct psc_case *_case, int *nr_particles_by_patch,
			int *particle_label_offset)
{
  struct psc *psc = _case->psc;

  int np_total = 0;
  psc_foreach_patch(psc, p) {
    int np = psc_case_calc_nr_particles_in_patch(_case, p);
    nr_particles_by_patch[p] = np;
    np_total += np;
  }

  // calculate global particle label offset for unique numbering
  *particle_label_offset = 0; // necessary on proc 0
  MPI_Exscan(&np_total, particle_label_offset, 1, MPI_INT, MPI_SUM,
	     MPI_COMM_WORLD);
}

static inline int
sgn(int x)
{
  return x < 0 ? -1 : 0;
}

void psc_case_init_particles_patch(struct psc_case *_case, int p,
				   int particle_label_offset)
{
  struct psc *psc = _case->psc;
  double beta = psc->coeff.beta;
  
  int ilo[3], ihi[3];
  pml_find_bounds(psc, p, ilo, ihi);

  struct psc_particles *prts_base = psc_mparticles_get_patch(psc->particles, p);
  struct psc_particles *prts = psc_particles_get_as(prts_base, "c", MP_DONT_COPY);
  
  int i = 0;
  for (int kind = 0; kind < _case->nr_kinds; kind++) {
    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	  double xx[3] = { CRDX(p, jx), CRDY(p, jy), CRDZ(p, jz) };
	  struct psc_particle_npt npt = {}; // init to all zero
	  psc_case_init_npt(_case, kind, xx, &npt);
	  
	  int n_in_cell = get_n_in_cell(psc, &npt);
	  for (int cnt = 0; cnt < n_in_cell; cnt++) {
	    particle_t *prt = particles_get_one(prts, i++);
	    
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
	    
	    //Add a (negligible) delta to offset particles towards the center of the patch
	    //This prevents particles at the boundaries to be discarded as "out of domain"
	    const double eps = psc->prm.initial_particle_shift;
	    double offx = -sgn(jx - (psc->patch[p].ldims[0] / 2)) * eps * psc->dx[0];
	    double offy = -sgn(jy - (psc->patch[p].ldims[1] / 2)) * eps * psc->dx[1];
	    double offz = -sgn(jz - (psc->patch[p].ldims[2] / 2)) * eps * psc->dx[2];
	      
	    prt->xi = xx[0] + offx;
	    prt->yi = xx[1] + offy;
	    prt->zi = xx[2] + offz;
	    prt->pxi = px;
	    prt->pyi = py;
	    prt->pzi = pz;
	    prt->qni = npt.q;
	    prt->mni = npt.m;
	    //prt->lni = particle_label_offset + 1;
	    if (psc->prm.fortran_particle_weight_hack) {
	      prt->wni = npt.n;
	    } else {
	      prt->wni = npt.n / (n_in_cell * psc->coeff.cori);
	    }

	  }
	} // jx
      }
    }
  } // kind
  prts->n_part = i;
  psc_particles_put_as(prts, prts_base, 0);
}

void
psc_case_init_particles(struct psc_case *_case, int *nr_particles_by_patch,
			int particle_label_offset)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (_case->seed_by_time) {
    srandom(10*rank + time(NULL));
  } else {
    srandom(rank);
  }

  psc_foreach_patch(_case->psc, p) {
    psc_case_init_particles_patch(_case, p, particle_label_offset);
    assert(psc_mparticles_nr_particles_by_patch(_case->psc->particles, p) == nr_particles_by_patch[p]);
  }
}

