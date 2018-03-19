
#pragma once

// ======================================================================
// SetupParticles

template<typename MP>
struct SetupParticles
{
  using Mparticles = MP;
  using particle_t = typename MP::particle_t;
  using mparticles_t = PscMparticles<Mparticles>;
  
  // ----------------------------------------------------------------------
  // get_n_in_cell
  //
  // helper function for partition / particle setup
  
  static int get_n_in_cell(struct psc *psc, struct psc_particle_npt *npt)
  {
    if (psc->prm.const_num_particles_per_cell) {
      return psc->prm.nicell;
    }
    if (npt->particles_per_cell) {
      return npt->n * npt->particles_per_cell + .5;
    }
    if (psc->prm.fractional_n_particles_per_cell) {
      int n_prts = npt->n / psc->coeff.cori;
      float rmndr = npt->n / psc->coeff.cori - n_prts;
      float ran = random() / ((float) RAND_MAX + 1);
      if (ran < rmndr) {
	n_prts++;
      }
      return n_prts;
    }
    return npt->n / psc->coeff.cori + .5;
  }

  // ----------------------------------------------------------------------
  // setup_particle
  
  static void setup_particle(struct psc *psc, particle_t *prt, struct psc_particle_npt *npt,
			     int p, double xx[3])
  {
    auto& kinds = psc->grid().kinds;
    double beta = psc->coeff.beta;
    
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
	      
    double pxi = npt->p[0] +
      sqrtf(-2.f*npt->T[0]/npt->m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
    double pyi = npt->p[1] +
      sqrtf(-2.f*npt->T[1]/npt->m*sqr(beta)*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4);
    double pzi = npt->p[2] +
      sqrtf(-2.f*npt->T[2]/npt->m*sqr(beta)*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6);

    if (psc->prm.initial_momentum_gamma_correction) {
      double gam;
      if (sqr(pxi) + sqr(pyi) + sqr(pzi) < 1.) {
	gam = 1. / sqrt(1. - sqr(pxi) - sqr(pyi) - sqr(pzi));
	pxi *= gam;
	pyi *= gam;
	pzi *= gam;
      }
    }
  
    assert(npt->kind >= 0 && npt->kind < kinds.size());
    prt->kind_ = npt->kind;
    assert(npt->q == kinds[prt->kind_].q);
    assert(npt->m == kinds[prt->kind_].m);
    /* prt->qni = kinds[prt->kind].q; */
    /* prt->mni = kinds[prt->kind].m; */
    prt->xi = xx[0] - psc->grid().patches[p].xb[0];
    prt->yi = xx[1] - psc->grid().patches[p].xb[1];
    prt->zi = xx[2] - psc->grid().patches[p].xb[2];
    prt->pxi = pxi * cos(psc->prm.theta_xz) + pzi * sin(psc->prm.theta_xz);
    prt->pyi = pyi;
    prt->pzi = - pxi * sin(psc->prm.theta_xz) + pzi * cos(psc->prm.theta_xz);
  }	      

  // ----------------------------------------------------------------------
  // setup_particles

  static void setup_particles(struct psc *psc, std::vector<uint>& n_prts_by_patch)
  {
    auto mprts_base = PscMparticlesBase{psc->particles};

    if (psc_ops(psc)->setup_particles) {
      mprts_base->reserve_all(n_prts_by_patch.data());
      psc_ops(psc)->setup_particles(psc, n_prts_by_patch, false);
      return;
    }

    if (psc_ops(psc)->init_npt) {
      setup_particles(mprts_base, psc, n_prts_by_patch);
    }
  }

  static void setup_particles(PscMparticlesBase mprts, psc* psc, std::vector<uint>& n_prts_by_patch)
  {
    setup_particles(mprts, psc, n_prts_by_patch, [&](int kind, double crd[3], psc_particle_npt& npt) {
	psc_ops(psc)->init_npt(psc, kind, crd, &npt);
      });
  }

  template<typename FUNC>
  static void setup_particles(PscMparticlesBase mprts_base, psc* psc, std::vector<uint>& n_prts_by_patch,
			      FUNC func)
  {
    auto mprts = mprts_base.get_as<mparticles_t>(MP_DONT_COPY | MP_DONT_RESIZE);
    setup_particles(*mprts.sub(), psc, n_prts_by_patch, func);
    mprts.put_as(mprts_base);
  }

  template<typename FUNC>
  static void setup_particles(Mparticles& mprts, psc* psc, std::vector<uint>& n_prts_by_patch,
			      FUNC func)
  {
    const auto& grid = mprts.grid();
    const auto& kinds = grid.kinds;
    
    mprts.reserve_all(n_prts_by_patch.data());

    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto ilo = Int3{}, ihi = grid.ldims;
  
      int nr_pop = psc->prm.nr_populations;
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	    double xx[3] = { .5 * (CRDX(p, jx) + CRDX(p, jx+1)),
			     .5 * (CRDY(p, jy) + CRDY(p, jy+1)),
			     .5 * (CRDZ(p, jz) + CRDZ(p, jz+1)) };
	    // FIXME, the issue really is that (2nd order) particle pushers
	    // don't handle the invariant dim right
	    if (grid.gdims[0] == 1) xx[0] = CRDX(p, jx);
	    if (grid.gdims[1] == 1) xx[1] = CRDY(p, jy);
	    if (grid.gdims[2] == 1) xx[2] = CRDZ(p, jz);
	  
	    int n_q_in_cell = 0;
	    for (int kind = 0; kind < nr_pop; kind++) {
	      struct psc_particle_npt npt = {};
	      if (kind < kinds.size()) {
		npt.kind = kind;
		npt.q    = kinds[kind].q;
		npt.m    = kinds[kind].m;
	      };
	      func(kind, xx, npt);

	      int n_in_cell;
	      if (kind != psc->prm.neutralizing_population) {
		n_in_cell = get_n_in_cell(psc, &npt);
		n_q_in_cell += npt.q * n_in_cell;
	      } else {
		// FIXME, should handle the case where not the last population is neutralizing
		assert(psc->prm.neutralizing_population == nr_pop - 1);
		n_in_cell = -n_q_in_cell / npt.q;
	      }
	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		particle_t prt;
		setup_particle(psc, &prt, &npt, p, xx);
		//p->lni = particle_label_offset + 1;
		if (psc->prm.fractional_n_particles_per_cell) {
		  prt.qni_wni_ = kinds[prt.kind_].q;
		} else {
		  prt.qni_wni_ = kinds[prt.kind_].q * npt.n / (n_in_cell * psc->coeff.cori);
		}
		mprts[p].push_back(prt);
	      }
	    }
	  }
	}
      }
      if (!psc->prm.fractional_n_particles_per_cell) {
	assert(mprts[p].size() == n_prts_by_patch[p]);
      }
    }
  }

  // ----------------------------------------------------------------------
  // setup_partition

  static std::vector<uint> setup_partition(struct psc *psc)
  {
    auto& kinds = psc->grid().kinds;
    std::vector<uint> n_prts_by_patch(psc->n_patches());

    if (psc_ops(psc)->setup_particles) {
      psc_ops(psc)->setup_particles(psc, n_prts_by_patch, true);
      return n_prts_by_patch;
    }
    if (!psc_ops(psc)->init_npt) {
      return n_prts_by_patch;
    }

    if (psc->prm.nr_populations < 0) {
      psc->prm.nr_populations = kinds.size();
    }
    if (psc->prm.neutralizing_population < 0) {
      psc->prm.neutralizing_population = psc->prm.nr_populations - 1;
    }

    const auto& grid = psc->grid();
    
    psc_foreach_patch(psc, p) {
      auto ilo = Int3{}, ihi = grid.ldims;

      int np = 0;
      int nr_pop = psc->prm.nr_populations;
      for (int kind = 0; kind < nr_pop; kind++) {
	for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	  for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	    for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	      double xx[3] = { .5 * (CRDX(p, jx) + CRDX(p, jx+1)),
			       .5 * (CRDY(p, jy) + CRDY(p, jy+1)),
			       .5 * (CRDZ(p, jz) + CRDZ(p, jz+1)) };
	      // FIXME, the issue really is that (2nd order) particle pushers
	      // don't handle the invariant dim right
	      if (grid.gdims[0] == 1) xx[0] = CRDX(p, jx);
	      if (grid.gdims[1] == 1) xx[1] = CRDY(p, jy);
	      if (grid.gdims[2] == 1) xx[2] = CRDZ(p, jz);

	      struct psc_particle_npt npt = {};
	      if (kind < kinds.size()) {
		npt.kind = kind;
		npt.q    = kinds[kind].q;
		npt.m    = kinds[kind].m;
	      };
	      psc_ops(psc)->init_npt(psc, kind, xx, &npt);

	      int n_in_cell = get_n_in_cell(psc, &npt);
	      if (psc->prm.fractional_n_particles_per_cell) {
		n_in_cell++; // we may get an extra particle
	      }
	      np += n_in_cell;
	    }
	  }
	}
      }
      n_prts_by_patch[p] = np;
    }
    return n_prts_by_patch;
  }
};

