
#pragma once

#include <psc_method.h>
#include <mrc_profile.h>
#include <psc_method.h>
#include <psc_diag.h>

#include <particles.hxx>

#include <push_particles.hxx>
#include <checks.hxx>
#include <output_particles.hxx>
#include <output_fields_c.hxx>

#ifdef VPIC
#include "../libpsc/vpic/vpic_iface.h"
void psc_method_vpic_print_status(struct psc_method *method); // FIXME
#endif

// ======================================================================
// PscParams

struct PscParams
{
  double cfl = { .75 };            // CFL number used to determine time step
  int nmax;                        // Number of timesteps to run
  double wallclock_limit = { 0. }; // Maximum wallclock time to run
  bool write_checkpoint = { false };
  int write_checkpoint_every_step = { 0 };

  bool detailed_profiling; // output profiling info for each process separately
  int stats_every;         // output timing and other info every so many steps
};
  
// ======================================================================
// Psc

template<typename PscConfig>
struct Psc
{
  using Mparticles_t = typename PscConfig::Mparticles_t;
  using MfieldsState = typename PscConfig::MfieldsState;
  using Balance_t = typename PscConfig::Balance_t;
  using Sort_t = typename PscConfig::Sort_t;
  using Collision_t = typename PscConfig::Collision_t;
  using PushParticles_t = typename PscConfig::PushParticles_t;
  using PushFields_t = typename PscConfig::PushFields_t;
  using Bnd_t = typename PscConfig::Bnd_t;
  using BndFields_t = typename PscConfig::BndFields_t;
  using BndParticles_t = typename PscConfig::BndParticles_t;
  using Checks_t = typename PscConfig::Checks_t;
  using Marder_t = typename PscConfig::Marder_t;

#ifdef VPIC
  using MaterialList = typename MfieldsState::MaterialList;
  using Material = typename MaterialList::Material;
#endif
  
  // ----------------------------------------------------------------------
  // ctor

  Psc()
  {
    time_start_ = MPI_Wtime();
    
    psc_ = psc_create(MPI_COMM_WORLD);
    psc_set_from_options(psc_);
  }

  // ----------------------------------------------------------------------
  // define_grid

  void define_grid(Grid_t::Domain& domain, GridBc& bc, Grid_t::Kinds& kinds,
		   double dt, Grid_t::NormalizationParams& norm_params)
  {
    auto coeff = Grid_t::Normalization{norm_params};
    grid_ = psc_setup_domain(psc_, domain, bc, kinds, coeff, dt);

#ifdef VPIC
    vgrid_ = Grid::create();
    vgrid_->setup(domain.dx, dt, norm_params.cc, norm_params.eps0);

    // define the grid
    define_periodic_grid(domain.corner, domain.corner + domain.length,
			 domain.gdims, domain.np);
    

    // set field boundary conditions
    for (int p = 0; p < grid().n_patches(); p++) {
      assert(p == 0);
      for (int d = 0; d < 3; d++) {
	bool lo = psc_at_boundary_lo(psc_, p, d);
	bool hi = psc_at_boundary_hi(psc_, p, d);

	if (lo && bc.fld_lo[d] != BND_FLD_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = -1;
	  set_domain_field_bc(bnd, bc.fld_lo[d]);
	}
	
	if (hi && bc.fld_hi[d] != BND_FLD_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = 1;
	  set_domain_field_bc(bnd, bc.fld_hi[d]);
	}
      }
    }

    // set particle boundary conditions
    for (int p = 0; p < grid().n_patches(); p++) {
      assert(p == 0);
      for (int d = 0; d < 3; d++) {
	bool lo = psc_at_boundary_lo(psc_, p, d);
	bool hi = psc_at_boundary_hi(psc_, p, d);

	if (lo && bc.prt_lo[d] != BND_PRT_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = -1;
	  set_domain_particle_bc(bnd, bc.prt_lo[d]);
	}
	
	if (hi && bc.prt_hi[d] != BND_PRT_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = 1;
	  set_domain_particle_bc(bnd, bc.prt_hi[d]);
	}
      }
    }

    grid_setup_communication();
#endif
  }
  
  // ----------------------------------------------------------------------
  // define_field_array

  void define_field_array(double damp = 0.)
  {
#ifdef VPIC
    // FIXME, mv assert innto MfieldsState ctor
    assert(!material_list_.empty());

    psc_->ibn[0] = psc_->ibn[1] = psc_->ibn[2] = 1;
#endif

#ifdef VPIC
    mflds_.reset(new MfieldsState{grid(), vgrid_, material_list_, damp});
#else
    mflds_.reset(new MfieldsState{grid()});
#endif

#ifdef VPIC
    hydro_.reset(new MfieldsHydro{grid(), vgrid_});
    interpolator_.reset(new MfieldsInterpolator{vgrid_});
    accumulator_.reset(new MfieldsAccumulator{vgrid_});
#endif
  }
  
  // ----------------------------------------------------------------------
  // init

  void init()
  {
    sort_.reset(new Sort_t{});
    pushp_.reset(new PushParticles_t{});
    pushf_.reset(new PushFields_t{});
    bnd_.reset(new Bnd_t{psc_->grid(), psc_->mrc_domain_, psc_->ibn});
    bndf_.reset(new BndFields_t{});
    bndp_.reset(new BndParticles_t{psc_->mrc_domain_, psc_->grid()});

    psc_setup_member_objs(psc_);
    initialize_stats();
  }

  // ----------------------------------------------------------------------
  // dtor

  ~Psc()
  {
    psc_destroy(psc_);
  }
  
  // ----------------------------------------------------------------------
  // initialize_stats
  
  void initialize_stats()
  {
    st_nr_particles = psc_stats_register("nr particles");
    st_time_step = psc_stats_register("time entire step");

    // generic stats categories
    st_time_particle = psc_stats_register("time particle update");
    st_time_field = psc_stats_register("time field update");
    st_time_comm = psc_stats_register("time communication");
    st_time_output = psc_stats_register("time output");
  }
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    psc_view(psc_);
    mprts_->view();

#ifdef VPIC
    initialize_vpic();
#else
    initialize_default();
#endif

    // initial output / stats
    mpi_printf(psc_comm(psc_), "Performing initial diagnostics.\n");
    diagnostics();
    print_status();

    mpi_printf(psc_comm(psc_), "Initialization complete.\n");
  }

  // ----------------------------------------------------------------------
  // integrate

  void integrate()
  {
    static int pr;
    if (!pr) {
      pr = prof_register("psc_step", 1., 0, 0);
    }

    mpi_printf(psc_comm(psc_), "*** Advancing\n");
    double elapsed = MPI_Wtime();

    bool first_iteration = true;
    while (psc_->timestep < p_.nmax) {
      prof_start(pr);
      psc_stats_start(st_time_step);

      if (!first_iteration &&
	  p_.write_checkpoint_every_step > 0 &&
	  psc_->timestep % p_.write_checkpoint_every_step == 0) {
	psc_write_checkpoint(psc_);
      }
      first_iteration = false;

      mpi_printf(psc_comm(psc_), "**** Step %d / %d, Code Time %g, Wall Time %g\n", psc_->timestep + 1,
		 p_.nmax, psc_->timestep * grid().dt, MPI_Wtime() - time_start_);

      prof_start(pr_time_step_no_comm);
      prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
    
      psc_->timestep++; // FIXME, too hacky
#ifdef VPIC
      if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
	vgrid_->step++;
	assert(vgrid_->step == psc_->timestep);
      }
#endif
      
      diagnostics();

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts_->get_n_prts();

      if (psc_->timestep % p_.stats_every == 0) {
	print_status();
      }

      if (p_.wallclock_limit > 0.) {
	double wallclock_elapsed = MPI_Wtime() - time_start_;
	double wallclock_elapsed_max;
	MPI_Allreduce(&wallclock_elapsed, &wallclock_elapsed_max, 1, MPI_DOUBLE, MPI_MAX,
		      MPI_COMM_WORLD);
      
	if (wallclock_elapsed_max > p_.wallclock_limit) {
	  mpi_printf(MPI_COMM_WORLD, "WARNING: Max wallclock time elapsed!\n");
	  break;
	}
      }
    }

    if (p_.write_checkpoint) {
      psc_write_checkpoint(psc_);
    }

    // FIXME, merge with existing handling of wallclock time
    elapsed = MPI_Wtime() - elapsed;

    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    mpi_printf(psc_comm(psc_), "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
	       elapsed, w, d, h, m, s );
  }

  virtual void step() = 0;

  // ----------------------------------------------------------------------
  // define_periodic_grid
  
  void define_periodic_grid(double xl[3], double xh[3], const int gdims[3], const int np[3])
  {
#ifdef VPIC
    // SimulationMixin::setTopology(np[0], np[1], np[2]); FIXME, needed for vpic_simulation,
    // I believe only because this info is written out in diagnostics_run
    vgrid_->partition_periodic_box(xl, xh, gdims, np);
#endif
  }
  
  // ----------------------------------------------------------------------
  // set_domain_field_bc
  
  void set_domain_field_bc(Int3 bnd, int bc)
  {
#ifdef VPIC
    int boundary = BOUNDARY(bnd[0], bnd[1], bnd[2]);
    int fbc;
    switch (bc) {
    case BND_FLD_CONDUCTING_WALL: fbc = Grid::pec_fields   ; break;
    case BND_FLD_ABSORBING:       fbc = Grid::absorb_fields; break;
    default: assert(0);
    }
    vgrid_->set_fbc(boundary, fbc);
#endif
  }

  // ----------------------------------------------------------------------
  // set_domain_particle_bc
  
  void set_domain_particle_bc(Int3 bnd, int bc)
  {
#ifdef VPIC
    int boundary = BOUNDARY(bnd[0], bnd[1], bnd[2]);
    int pbc;
    switch (bc) {
    case BND_PRT_REFLECTING: pbc = Grid::reflect_particles; break;
    case BND_PRT_ABSORBING:  pbc = Grid::absorb_particles ; break;
    default: assert(0);
    }
    vgrid_->set_pbc(boundary, pbc);
#endif
  }

  // ----------------------------------------------------------------------
  // define_material

#ifdef VPIC
  Material* define_material(const char *name,
			    double eps, double mu=1.,
			    double sigma=0., double zeta=0.)
  {
    auto m = material_list_.create(name,
		       eps,   eps,   eps,
		       mu,    mu,    mu,
		       sigma, sigma, sigma,
		       zeta,  zeta,  zeta);
    return material_list_.append(m);
  }
#endif

  void grid_setup_communication()
  {
#ifdef VPIC
    assert(vgrid_->nx && vgrid_->ny && vgrid_->ny);
    
    // Pre-size communications buffers. This is done to get most memory
    // allocation over with before the simulation starts running
    // FIXME, this isn't a great place. First, we shouldn't call mp
    // functions (semi-)directly. 2nd, whether we need these buffers depends
    // on b.c., which aren't yet known.
    
    // FIXME, this really isn't a good place to do this, as it requires layer breaking knowledge of
    // which communication will need the largest buffers...
    int nx1 = vgrid_->nx+1, ny1 = vgrid_->ny+1, nz1 = vgrid_->nz+1;
    vgrid_->mp_size_recv_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
    
    vgrid_->mp_size_send_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
#endif
  }

  // ----------------------------------------------------------------------
  // courant length
  
  double courant_length(const Grid_t::Domain& domain)
  {
    double inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!domain.isInvar(d)) {
	inv_sum += 1. / sqr(domain.dx[d]);
      }
    }
    if (!inv_sum) { // simulation has 0 dimensions (happens in some test?)
      inv_sum = 1.;
    }
    return sqrt(1. / inv_sum);
  }

  // ----------------------------------------------------------------------
  // create_diagnotics

  void create_diagnostics(int interval)
  {
#ifdef VPIC
    diag_mixin_.diagnostics_init(interval);
#endif
  }
  
  // ----------------------------------------------------------------------
  // setup_diagnostics

  void setup_diagnostics()
  {
#ifdef VPIC
    diag_mixin_.diagnostics_setup();
#endif
  }

private:

  // ----------------------------------------------------------------------
  // print_profiling

  void print_profiling()
  {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!p_.detailed_profiling) {
      prof_print_mpi(MPI_COMM_WORLD);
    } else {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      for (int i = 0; i < size; i++) {
	if (i == rank) {
	  mprintf("profile\n");
	  prof_print();
	}
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
  
#ifndef VPIC
  // ----------------------------------------------------------------------
  // initialize_default
  
  void initialize_default()
  {
    //pushp_.stagger(mprts, mflds); FIXME, vpic does it

    checks_->gauss(*mprts_, *mflds_);
  }
#endif

#ifdef VPIC

  // ----------------------------------------------------------------------
  // initialize_vpic
  
  void initialize_vpic()
  {
    // FIXME, just change the uses
    auto psc = psc_;
    auto& mprts = *mprts_;
    auto& mflds = mflds_;
    // Do some consistency checks on user initialized fields
    
    mpi_printf(psc_comm(psc), "Checking interdomain synchronization\n");
    double err;
    TIC err = CleanDivOps::synchronize_tang_e_norm_b(*mflds_); TOC(synchronize_tang_e_norm_b, 1);
    mpi_printf(psc_comm(psc), "Error = %g (arb units)\n", err);
    
    mpi_printf(psc_comm(psc), "Checking magnetic field divergence\n");
    TIC CleanDivOps::compute_div_b_err(*mflds_); TOC(compute_div_b_err, 1);
    TIC err = CleanDivOps::compute_rms_div_b_err(*mflds_); TOC(compute_rms_div_b_err, 1);
    mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
    TIC CleanDivOps::clean_div_b(*mflds_); TOC(clean_div_b, 1);
    
    // Load fields not initialized by the user
    
    mpi_printf(psc_comm(psc), "Initializing radiation damping fields\n");
    TIC AccumulateOps::compute_curl_b(*mflds_); TOC(compute_curl_b, 1);
    
    mpi_printf(psc_comm(psc), "Initializing bound charge density\n");
    TIC CleanDivOps::clear_rhof(*mflds_); TOC(clear_rhof, 1);
    ParticlesOps::accumulate_rho_p(*mprts_, *mflds_);
    CleanDivOps::synchronize_rho(*mflds_);
    TIC AccumulateOps::compute_rhob(*mflds_); TOC(compute_rhob, 1);
    
    // Internal sanity checks
    
    mpi_printf(psc_comm(psc), "Checking electric field divergence\n");
    TIC CleanDivOps::compute_div_e_err(*mflds_); TOC(compute_div_e_err, 1);
    TIC err = CleanDivOps::compute_rms_div_e_err(*mflds_); TOC(compute_rms_div_e_err, 1);
    mpi_printf(psc_comm(psc), "RMS error = %e (charge/volume)\n", err);
    TIC CleanDivOps::clean_div_e(*mflds_); TOC(clean_div_e, 1);
    
    mpi_printf(psc_comm(psc), "Rechecking interdomain synchronization\n");
    TIC err = CleanDivOps::synchronize_tang_e_norm_b(*mflds_); TOC(synchronize_tang_e_norm_b, 1);
    mpi_printf(psc_comm(psc), "Error = %e (arb units)\n", err);
    
    mpi_printf(psc_comm(psc), "Uncentering particles\n");
    if (!mprts_->empty()) {
      TIC InterpolatorOps::load(*interpolator_, *mflds_); TOC(load_interpolator, 1);
      
      for (auto& sp : mprts) {
	TIC ParticlesOps::uncenter_p(&sp, *interpolator_); TOC(uncenter_p, 1);
      }
    }
  }

#endif
  
  // ----------------------------------------------------------------------
  // diagnostics

  virtual void diagnostics()
  {
#ifdef VPIC
#if 0
    TIC user_diagnostics(); TOC(user_diagnostics, 1);
#endif
    if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
      diag_mixin_.diagnostics_run(*mprts_, *mflds_, *interpolator_, *hydro_, grid().domain.np);
    }
#else
    // FIXME
    psc_diag_run(psc_->diag, psc_, *mprts_, *mflds_);
    // FIXME
    (*outf_)(*mflds_, *mprts_);
#endif
    PscOutputParticlesBase{psc_->output_particles}.run(*mprts_);
  }

  // ----------------------------------------------------------------------
  // print_status

  void print_status()
  {
#ifdef VPIC
    if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
      psc_method_vpic_print_status(psc_->method);
    }
#endif
    psc_stats_log(psc_->timestep);
    print_profiling();
  }

public:
  const Grid_t& grid() { return *grid_; }

protected:
  double time_start_;

  PscParams p_;
  const Grid_t* grid_;
#ifdef VPIC
  Grid* vgrid_;
#endif
  psc* psc_;

#ifdef VPIC
  MaterialList material_list_;
#endif
  std::unique_ptr<MfieldsState> mflds_;
#ifdef VPIC
  std::unique_ptr<MfieldsHydro> hydro_;
  std::unique_ptr<MfieldsInterpolator> interpolator_;
  std::unique_ptr<MfieldsAccumulator> accumulator_;
  ParticleBcList particle_bc_list_;
#endif
  std::unique_ptr<Mparticles_t> mprts_;

  std::unique_ptr<Balance_t> balance_;
  std::unique_ptr<Sort_t> sort_;
  std::unique_ptr<Collision_t> collision_;
  std::unique_ptr<PushParticles_t> pushp_;
  std::unique_ptr<PushFields_t> pushf_;
  std::unique_ptr<Bnd_t> bnd_;
  std::unique_ptr<BndFields_t> bndf_;
  std::unique_ptr<BndParticles_t> bndp_;
  std::unique_ptr<Checks_t> checks_;
  std::unique_ptr<Marder_t> marder_;
  std::unique_ptr<OutputFieldsC> outf_;

#ifdef VPIC
  DiagMixin diag_mixin_;
#endif

  // FIXME, maybe should be private
  // need to make sure derived class sets these (? -- or just leave them off by default)
  int balance_interval;
  int sort_interval;
  int marder_interval;

  int num_comm_round = {3};
  
  int st_nr_particles;
  int st_time_step;
};
