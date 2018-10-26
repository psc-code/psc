
#pragma once

// ======================================================================
// SetupParticles

template<typename MP>
struct SetupParticles
{
  using Mparticles = MP;
  using particle_t = typename MP::particle_t;
  
  // ----------------------------------------------------------------------
  // get_n_in_cell
  //
  // helper function for partition / particle setup
  
  int get_n_in_cell(const Grid_t& grid, psc_particle_npt *npt)
  {
    if (fractional_n_particles_per_cell) {
      int n_prts = npt->n / grid.norm.cori;
      float rmndr = npt->n / grid.norm.cori - n_prts;
      float ran = random() / ((float) RAND_MAX + 1);
      if (ran < rmndr) {
	n_prts++;
      }
      return n_prts;
    }
    return npt->n / grid.norm.cori + .5;
  }

  // ----------------------------------------------------------------------
  // setup_particle
  
  void setup_particle(const Grid_t& grid, particle_inject *prt, psc_particle_npt *npt,
		      int p, double xx[3])
  {
    auto& kinds = grid.kinds;
    double beta = grid.norm.beta;
    
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

    if (initial_momentum_gamma_correction) {
      double gam;
      if (sqr(pxi) + sqr(pyi) + sqr(pzi) < 1.) {
	gam = 1. / sqrt(1. - sqr(pxi) - sqr(pyi) - sqr(pzi));
	pxi *= gam;
	pyi *= gam;
	pzi *= gam;
      }
    }
  
    assert(npt->kind >= 0 && npt->kind < kinds.size());
    prt->kind = npt->kind;
    assert(npt->q == kinds[prt->kind].q);
    assert(npt->m == kinds[prt->kind].m);
    /* prt->qni = kinds[prt->kind].q; */
    /* prt->mni = kinds[prt->kind].m; */
    prt->x[0] = xx[0];
    prt->x[1] = xx[1];
    prt->x[2] = xx[2];
    prt->u[0] = pxi;
    prt->u[1] = pyi;
    prt->u[2] = pzi;
  }	      

  // ----------------------------------------------------------------------
  // setup_particles

  template<typename FUNC>
  void setup_particles(Mparticles& mprts, std::vector<uint>& n_prts_by_patch,
		       FUNC func)
  {
    const auto& grid = mprts.grid();
    const auto& kinds = grid.kinds;
    
    mprts.reserve_all(n_prts_by_patch.data());

    for (int p = 0; p < mprts.n_patches(); ++p) {
      Int3 ilo = {}, ihi = grid.ldims;
      auto injector = mprts[p].injector();
  
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	    double xx[3] = {grid.patches[p].x_cc(jx),
			    grid.patches[p].y_cc(jy),
			    grid.patches[p].z_cc(jz)};
	    // FIXME, the issue really is that (2nd order) particle pushers
	    // don't handle the invariant dim right
	    if (grid.isInvar(0) == 1) xx[0] = grid.patches[p].x_nc(jx);
	    if (grid.isInvar(1) == 1) xx[1] = grid.patches[p].y_nc(jy);
	    if (grid.isInvar(2) == 1) xx[2] = grid.patches[p].z_nc(jz);
	  
	    int n_q_in_cell = 0;
	    for (int kind = 0; kind < kinds.size(); kind++) {
	      struct psc_particle_npt npt = {};
	      if (kind < kinds.size()) {
		npt.kind = kind;
		npt.q    = kinds[kind].q;
		npt.m    = kinds[kind].m;
	      };
	      func(kind, xx, npt);

	      int n_in_cell;
	      if (kind != neutralizing_population) {
		n_in_cell = get_n_in_cell(grid, &npt);
		n_q_in_cell += npt.q * n_in_cell;
	      } else {
		// FIXME, should handle the case where not the last population is neutralizing
		assert(neutralizing_population == kinds.size() - 1);
		n_in_cell = -n_q_in_cell / npt.q;
	      }
	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		using real_t = typename particle_t::real_t;
		int kind = npt.kind;
		real_t wni;
		if (fractional_n_particles_per_cell) {
		  wni = 1.;
		} else {
		  wni = npt.n / (n_in_cell * grid.norm.cori);
		}
		particle_inject prt{{}, {}, wni, kind};
		setup_particle(grid, &prt, &npt, p, xx);
		//p->lni = particle_label_offset + 1;
		injector(prt);
	      }
	    }
	  }
	}
      }
      if (!fractional_n_particles_per_cell) {
	assert(mprts[p].size() == n_prts_by_patch[p]);
      }
    }
  }

  // ----------------------------------------------------------------------
  // setup_partition

  template<typename FUNC>
  std::vector<uint> setup_partition(const Grid_t& grid, FUNC func)
  {
    const auto& kinds = grid.kinds;
    std::vector<uint> n_prts_by_patch(grid.n_patches());
    
    for (int p = 0; p < grid.n_patches(); ++p) {
      auto ilo = Int3{}, ihi = grid.ldims;
  
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	    double xx[3] = {grid.patches[p].x_cc(jx), grid.patches[p].y_cc(jy), grid.patches[p].z_cc(jz)};
	    // FIXME, the issue really is that (2nd order) particle pushers
	    // don't handle the invariant dim right
	    if (grid.isInvar(0) == 1) xx[0] = grid.patches[p].x_nc(jx);
	    if (grid.isInvar(1) == 1) xx[1] = grid.patches[p].y_nc(jy);
	    if (grid.isInvar(2) == 1) xx[2] = grid.patches[p].z_nc(jz);
	  
	    int n_q_in_cell = 0;
	    for (int kind = 0; kind < kinds.size(); kind++) {
	      struct psc_particle_npt npt = {};
	      if (kind < kinds.size()) {
		npt.kind = kind;
		npt.q    = kinds[kind].q;
		npt.m    = kinds[kind].m;
	      };
	      func(kind, xx, npt);

	      int n_in_cell;
	      if (kind != neutralizing_population) {
		n_in_cell = get_n_in_cell(grid, &npt);
		n_q_in_cell += npt.q * n_in_cell;
	      } else {
		// FIXME, should handle the case where not the last population is neutralizing
		assert(neutralizing_population == kinds.size() - 1);
		n_in_cell = -n_q_in_cell / npt.q;
	      }
	      n_prts_by_patch[p] += n_in_cell;
	    }
	  }
	}
      }
    }

    return n_prts_by_patch;
  }

  // the initial number of particles in a cell for this population will be st so that it achieves neutrality  
  int neutralizing_population = { -1 };
  bool fractional_n_particles_per_cell = { false };
  bool initial_momentum_gamma_correction = { false };
};

