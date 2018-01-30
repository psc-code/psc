
#include <mrc_profile.h>

#include <cstring>

#include "ddc_particles.hxx"

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

// ----------------------------------------------------------------------
// at_lo/hi_boundary

static inline bool
at_lo_boundary(int p, int d)
{
  return ppsc->patch[p].off[d] == 0;
}

static inline bool
at_hi_boundary(int p, int d)
{
  return ppsc->patch[p].off[d] + ppsc->patch[p].ldims[d] == ppsc->domain.gdims[d];
}

// ======================================================================
// psc_bnd_particles

template<typename MP>
struct psc_bnd_particles_sub
{
  using mparticles_t = MP;
  using particle_t = typename mparticles_t::particle_t;
  using real_t = typename mparticles_t::real_t;
  using ddcp_t = ddc_particles<mparticles_t>;
  using ddcp_patch = typename ddcp_t::patch;

  // ----------------------------------------------------------------------
  // interface to psc_bnd_particles_ops

  static void destroy(struct psc_bnd_particles *bnd);
  static void setup(struct psc_bnd_particles *bnd);
  static void reset(struct psc_bnd_particles *bnd);
  static void exchange_particles(struct psc_bnd_particles *bnd,
				 struct psc_mparticles *mprts_base);

protected:
  // ----------------------------------------------------------------------
  // ctor

  psc_bnd_particles_sub(struct mrc_domain *domain)
    : ddcp(nullptr)
  {
    ddcp = new ddc_particles<mparticles_t>(domain);
  }

  // ----------------------------------------------------------------------
  // dtor

  ~psc_bnd_particles_sub()
  {
    delete ddcp;
  }

  // ----------------------------------------------------------------------
  // reset

  void reset(struct mrc_domain *domain)
  {
    delete ddcp;
    ddcp = new ddc_particles<mparticles_t>(domain);
  }

  void process_patch(mparticles_t mprts, int p);
  void process_and_exchange(mparticles_t mprts);
  void exchange_particles(mparticles_t mprts);

protected:
  ddcp_t* ddcp;
};

// ----------------------------------------------------------------------
// psc_bnd_particles_sub::process_patch

template<typename MP>
void psc_bnd_particles_sub<MP>::process_patch(mparticles_t mprts, int p)
{
  struct psc *psc = ppsc;

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.

  const Grid_t::Patch& gpatch = psc->grid.patches[p];
  const real_t *b_dxi = mprts[p].get_b_dxi();
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
    
    int b_pos[3];
    mprts[p].get_block_pos(xi, b_pos);
    
    if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] && // OPT, could be optimized with casts to unsigned
	b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      // fast path
      // particle is still inside patch: move into right position
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
	if (!at_lo_boundary(p, d) || psc->domain.bnd_part_lo[d] == BND_PART_PERIODIC) {
	  xi[d] += xm[d];
	  dir[d] = -1;
	  int bi = fint(xi[d] * b_dxi[d]);
	  if (bi >= b_mx[d]) {
	    xi[d] = 0.;
	    dir[d] = 0;
	  }
	} else {
	  switch (psc->domain.bnd_part_lo[d]) {
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
	if (!at_hi_boundary(p, d) ||
	    psc->domain.bnd_part_hi[d] == BND_PART_PERIODIC) {
	  xi[d] -= xm[d];
	  dir[d] = +1;
	  int bi = fint(xi[d] * b_dxi[d]);
	  if (bi < 0) {
	    xi[d] = 0.;
	  }
	} else {
	  switch (psc->domain.bnd_part_hi[d]) {
	  case BND_PART_REFLECTING: {
	    xi[d] = 2.f * xm[d] - xi[d];
	    pxi[d] = -pxi[d];
	    dir[d] = 0;
	    int bi = fint(xi[d] * b_dxi[d]);
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
	(*dpatch->m_buf)[head++] = *prt;
      } else {
	typename ddc_particles<mparticles_t>::dnei *nei = &dpatch->nei[mrc_ddc_dir2idx(dir)];
	nei->send_buf.push_back(*prt);
      }
    }
  }
  dpatch->m_buf->resize(head);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub::process_and_exchange

template<typename MP>
void psc_bnd_particles_sub<MP>::process_and_exchange(mparticles_t mprts)
{
  // FIXME we should make sure (assert) we don't quietly drop particle which left
  // in the invariant direction

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
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub::exchange_particles

template<typename MP>
void psc_bnd_particles_sub<MP>::exchange_particles(mparticles_t mprts)
{
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
// psc_bnd_particles_sub::destroy

template<typename MP>
void psc_bnd_particles_sub<MP>::destroy(struct psc_bnd_particles *bnd)
{
  auto sub = static_cast<psc_bnd_particles_sub<MP>*>(bnd->obj.subctx);
  
  sub->~psc_bnd_particles_sub<MP>();
}
  
// ----------------------------------------------------------------------
// psc_bnd_particles_sub::setup

template<typename MP>
void psc_bnd_particles_sub<MP>::setup(struct psc_bnd_particles *bnd)
{
  auto sub = static_cast<psc_bnd_particles_sub<MP>*>(bnd->obj.subctx);

  new(sub) psc_bnd_particles_sub<MP>(bnd->psc->mrc_domain);
  //psc_bnd_particles_open_setup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_reset

template<typename MP>
void psc_bnd_particles_sub<MP>::reset(struct psc_bnd_particles *bnd)
{
  auto sub = static_cast<psc_bnd_particles_sub<MP>*>(bnd->obj.subctx);

  sub->reset(bnd->psc->mrc_domain);
  //psc_bnd_particles_open_reset(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles

template<typename MP>
void psc_bnd_particles_sub<MP>::exchange_particles(struct psc_bnd_particles *bnd,
						      struct psc_mparticles *mprts_base)
{
  auto sub = static_cast<psc_bnd_particles_sub<MP>*>(bnd->obj.subctx);
  mparticles_t mprts = mprts_base->get_as<MP>();
  
  sub->exchange_particles(mprts);

  mprts.put_as(mprts_base);
}
