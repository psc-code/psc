
#include <mrc_profile.h>

#include <cstring>

#include "psc_balance.h"
#include "ddc_particles.hxx"
#include "bnd_particles.hxx"

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

// ======================================================================
// psc_bnd_particles

template<typename MP>
struct psc_bnd_particles_sub : BndParticlesBase
{
  using Mparticles = MP;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename Mparticles::real_t;
  using mparticles_t = PscMparticles<Mparticles>;
  using ddcp_t = ddc_particles<mparticles_t>;
  using ddcp_patch = typename ddcp_t::patch;

  // ----------------------------------------------------------------------
  // ctor

  psc_bnd_particles_sub(struct mrc_domain *domain, const Grid_t& grid)
    : ddcp{},
      balance_generation_cnt_{-1}
  {
    reset();
  }

  // ----------------------------------------------------------------------
  // dtor

  ~psc_bnd_particles_sub()
  {
    delete ddcp;
  }

  // ----------------------------------------------------------------------
  // reset

  void reset()
  {
    delete ddcp;
    ddcp = new ddc_particles<mparticles_t>(ppsc->mrc_domain_);
    balance_generation_cnt_ = psc_balance_generation_cnt;
  }

  // ----------------------------------------------------------------------
  // operator() (exchange on correct particle type)
  
  void operator()(Mparticles& mprts)
  {
    if (psc_balance_generation_cnt > balance_generation_cnt_) {
      reset();
    }
    for (int p = 0; p < mprts.n_patches(); p++) {
      ddcp_patch *dpatch = &ddcp->patches[p];
      dpatch->m_buf = &mprts[p].get_buf();
      dpatch->m_begin = 0;
    }

    process_and_exchange(mprts);
    
    //struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, "c", JXI, JXI + 3);
    //psc_bnd_particles_open_boundary(bnd, particles, mflds);
    //psc_mfields_put_as(mflds, psc->flds, JXI, JXI + 3);
  }

  // ----------------------------------------------------------------------
  // process_and_exchange

  void process_and_exchange(Mparticles& mprts)
  {
    static int pr_B, pr_C;
    if (!pr_B) {
      pr_B = prof_register("xchg_prep", 1., 0, 0);
      pr_C = prof_register("xchg_comm", 1., 0, 0);
    }
  
    prof_restart(pr_time_step_no_comm);
    prof_start(pr_B);
#pragma omp parallel for
    for (int p = 0; p < mprts.n_patches(); p++) {
      psc_balance_comp_time_by_patch[p] -= MPI_Wtime();
      process_patch(mprts, p);
      psc_balance_comp_time_by_patch[p] += MPI_Wtime();
    }
    prof_stop(pr_B);
    prof_stop(pr_time_step_no_comm);
    
    prof_start(pr_C);
    ddcp->comm();
    prof_stop(pr_C);

    //mprts->check();
  }
  
  // ----------------------------------------------------------------------
  // exchange_particles

  void exchange_particles(PscMparticlesBase mprts_base) override
  {
    mparticles_t mprts = mprts_base.get_as<mparticles_t>();
    (*this)(*mprts.sub());
    mprts.put_as(mprts_base);
  }

protected:
  void process_patch(Mparticles& mprts, int p);

protected:
  ddcp_t* ddcp;
  int balance_generation_cnt_;
};

// ----------------------------------------------------------------------
// psc_bnd_particles_sub::process_patch

template<typename MP>
void psc_bnd_particles_sub<MP>::process_patch(Mparticles& mprts, int p)
{
  struct psc *psc = ppsc;

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.

  const Grid_t::Patch& gpatch = psc->grid().patches[p];
  const int *b_mx = mprts[p].get_b_mx();
  real_t xm[3];
  for (int d = 0; d < 3; d++ ) {
    xm[d] = gpatch.xe[d] - gpatch.xb[d];
  }
  
  typename ddc_particles<mparticles_t>::patch *dpatch = &ddcp->patches[p];
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    dpatch->nei[dir1].send_buf.resize(0);
  }

  unsigned int n_begin = dpatch->m_begin;
  unsigned int n_end = dpatch->m_buf->size();
  unsigned int head = n_begin;

  for (int n = n_begin; n < n_end; n++) {
    particle_t *prt = &(*dpatch->m_buf)[n];
    real_t *xi = &prt->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
    real_t *pxi = &prt->pxi;
    
    Int3 b_pos = mprts[p].blockPosition(xi);
    
    if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] && // OPT, could be optimized with casts to unsigned
	b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      // fast path
      // particle is still inside patch: move into right position
      mprts[p].validCellIndex(*prt);
      (*dpatch->m_buf)[head++] = *prt;
      continue;
    }

    // slow path
    // handle particles which (seemingly) left the patch
    // (may end up in the same patch, anyway, though)
    bool drop = false;
    int dir[3];
    for (int d = 0; d < 3; d++) {
      if (b_pos[d] < 0) {
	if (!psc_at_boundary_lo(ppsc, p, d) || psc->domain_.bnd_part_lo[d] == BND_PART_PERIODIC) {
	  xi[d] += xm[d];
	  dir[d] = -1;
	  int bi = mprts[p].blockPosition(xi[d], d);
	  if (bi >= b_mx[d]) {
	    xi[d] = 0.;
	    dir[d] = 0;
	  }
	} else {
	  switch (psc->domain_.bnd_part_lo[d]) {
	  case BND_PART_REFLECTING:
	    xi[d] =  -xi[d];
	    pxi[d] = -pxi[d];
	    dir[d] = 0;
	    break;
	  case BND_PART_ABSORBING:
	    drop = true;
	    break;
	  default:
	    assert(0);
	  }
	}
      } else if (b_pos[d] >= b_mx[d]) {
	if (!psc_at_boundary_hi(ppsc, p, d) ||
	    psc->domain_.bnd_part_hi[d] == BND_PART_PERIODIC) {
	  xi[d] -= xm[d];
	  dir[d] = +1;
	  int bi = mprts[p].blockPosition(xi[d], d);
	  if (bi < 0) {
	    xi[d] = 0.;
	  }
	} else {
	  switch (psc->domain_.bnd_part_hi[d]) {
	  case BND_PART_REFLECTING: {
	    xi[d] = 2.f * xm[d] - xi[d];
	    pxi[d] = -pxi[d];
	    dir[d] = 0;
	    int bi = mprts[p].blockPosition(xi[d], d);
	    if (bi >= b_mx[d]) {
	      xi[d] *= (1. - 1e-6);
	    }
	    break;
	  }
	  case BND_PART_ABSORBING:
	    drop = true;
	    break;
	  default:
	    assert(0);
	  }
	}
      } else {
	// computational bnd
	dir[d] = 0;
      }
      if (!drop) {
	if (xi[d] < 0.f && xi[d] > -1e-6f) {
	  mprintf("d %d xi %g\n", d, xi[d]);
	  xi[d] = 0.f;
	}
	assert(xi[d] >= 0.f);
	assert(xi[d] <= xm[d]);
      }
    }
    if (!drop) {
      if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	mprts[p].validCellIndex(*prt);
	(*dpatch->m_buf)[head++] = *prt;
      } else {
	typename ddc_particles<mparticles_t>::dnei *nei = &dpatch->nei[mrc_ddc_dir2idx(dir)];
	nei->send_buf.push_back(*prt);
      }
    }
  }
  dpatch->m_buf->resize(head);
}

