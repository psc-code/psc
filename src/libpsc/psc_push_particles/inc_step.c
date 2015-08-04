
// ======================================================================
// EXT_PREPARE_SORT
//
// if enabled, calculate the new block position for each particle once
// it's moved, and also append particles that left the block to an extra
// list at the end of all local particles (hopefully there's enough room...)

#ifdef EXT_PREPARE_SORT

#ifdef PSC_PARTICLES_AS_SINGLE
static struct psc_particles_single *prts_sub;
#pragma omp threadprivate(prts_sub)
#endif

static inline void
ext_prepare_sort_before(struct psc_particles *prts)
{
  prts_sub = psc_particles_single(prts);
  memset(prts_sub->b_cnt, 0,
	 (prts_sub->nr_blocks + 1) * sizeof(*prts_sub->b_cnt));
}

static inline void
ext_prepare_sort(struct psc_particles *prts, int n, particle_t *prt,
		 int *b_pos)
{
  /* FIXME, only if blocksize == 1! */
  int *b_mx = prts_sub->b_mx;
  if (b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
      b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
    prts_sub->b_idx[n] = b_pos[2] * b_mx[1] + b_pos[1];
  } else { /* out of bounds */
    prts_sub->b_idx[n] = prts_sub->nr_blocks;
    assert(prts_sub->b_cnt[prts_sub->nr_blocks] < prts_sub->n_alloced);
    /* append to back */
    *particles_get_one(prts, prts->n_part + prts_sub->b_cnt[prts_sub->nr_blocks]) = *prt;
  }
  prts_sub->b_cnt[prts_sub->b_idx[n]]++;
}

#else

static inline void
ext_prepare_sort_before(struct psc_particles *prts)
{
}

static inline void
ext_prepare_sort(struct psc_particles *prts, int n, particle_t *prt,
		 int *b_pos)
{
}

#endif

// ======================================================================

// ----------------------------------------------------------------------
// push_one_

static inline void
push_one_(particle_t *prt, struct psc_fields *flds, struct psc_particles *prts, int n)
{
  // field interpolation
  int lg[3], lh[3];
  particle_real_t og[3], oh[3], xm[3];
#if PSC_PARTICLES_AS_CUDA2
  find_idx_off_pos_1st_rel(&prt->xi4.x, lg, og, xm, 0.f, prm.dxi); // FIXME passing xi hack
  find_idx_off_1st_rel(&prt->xi4.x, lh, oh, -.5f, prm.dxi);
#else  
  find_idx_off_pos_1st_rel(&prt->xi, lg, og, xm, 0.f, prm.dxi); // FIXME passing xi hack
  find_idx_off_1st_rel(&prt->xi, lh, oh, -.5f, prm.dxi);
#endif

  particle_real_t exq, eyq, ezq, hxq, hyq, hzq;
  INTERPOLATE_1ST(flds, exq, eyq, ezq, hxq, hyq, hzq);
  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
#if PSC_PARTICLES_AS_CUDA2
  int kind = cuda_float_as_int(prt->xi4.w);
#else
  int kind = prt->kind;
#endif
  particle_real_t dq = prm.dq_kind[kind];
  push_pxi(prt, exq, eyq, ezq, hxq, hyq, hzq, dq);

  particle_real_t vxi[3];
  calc_vxi(vxi, prt);
#if CALC_J == CALC_J_1VB_2D
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
  push_xi(prt, vxi, .5f * prm.dt);
  
  // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt
  calc_j_oop(flds, prt, vxi);
  
  // x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)
  push_xi(prt, vxi, .5f * prm.dt);
#else
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  push_xi(prt, vxi, prm.dt);
#endif
  
  int lf[3];
  particle_real_t of[3], xp[3];
#if PSC_PARTICLES_AS_CUDA2
  find_idx_off_pos_1st_rel(&prt->xi4.x, lf, of, xp, 0.f, prm.dxi);
#else
  find_idx_off_pos_1st_rel(&prt->xi, lf, of, xp, 0.f, prm.dxi);
#endif
  ext_prepare_sort(prts, n, prt, lf);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
  calc_j(flds, xm, xp, lf, lg, prt, vxi);
}

// ----------------------------------------------------------------------
// push_one

#if PSC_PARTICLES_AS_CUDA2

static inline void
push_one(struct psc_mparticles *mprts, struct psc_mfields *mflds, int n, int p)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  struct psc_fields *flds = psc_mfields_get_patch(mflds, p);

  particle_cuda2_t prt;
  _LOAD_PARTICLE_POS(prt, mprts_sub->h_xi4, n);
  _LOAD_PARTICLE_MOM(prt, mprts_sub->h_pxi4, n);
  
  push_one_(&prt, flds, NULL, n);

  _STORE_PARTICLE_POS(prt, mprts_sub->h_xi4, n);
  _STORE_PARTICLE_MOM(prt, mprts_sub->h_pxi4, n);
}

#else

static void
push_one(struct psc_fields *flds, struct psc_particles *prts, int n)
{
  particle_t *prt = particles_get_one(prts, n);

  push_one_(prt, flds, prts, n);
}

#endif
