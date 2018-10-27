
#pragma once

#include "psc.h"
#include "psc_fields_single.h"
#include "psc_particles_single.h"
#include "../libpsc/psc_push_particles/push_config.hxx"
#include "../libpsc/psc_push_particles/push_part_common.c"
#include "../libpsc/psc_push_particles/1vb/psc_push_particles_1vb.h"
#include "../libpsc/psc_push_particles/1vb.c"
#include "bnd_particles_impl.hxx"
#include "../libpsc/psc_checks/checks_impl.hxx"
#include "psc_push_fields_impl.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "setup_fields.hxx"
#include "setup_particles.hxx"

#ifdef USE_VPIC
#include "../libpsc/vpic/push_particles_vpic.hxx"
#include "../libpsc/vpic/bnd_particles_vpic.hxx"
#endif

#ifdef USE_CUDA
#include "../libpsc/cuda/push_particles_cuda_impl.hxx"
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#include "../libpsc/cuda/bnd_particles_cuda_impl.hxx"
#include "../libpsc/cuda/checks_cuda_impl.hxx"
#include "../libpsc/cuda/push_fields_cuda_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_2_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_3_impl.hxx"
#include "../libpsc/cuda/fields_item_moments_1st_cuda.hxx"
#endif

// ======================================================================
// Rng hackiness

#include "../vpic/PscRng.h"

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

// ======================================================================
// TestConfig

template<typename DIM, typename _Mfields, typename PUSHP, typename ORDER,
	 typename CHECKS = Checks_<typename PUSHP::Mparticles, typename PUSHP::MfieldsState, _Mfields, ORDER>,
	 typename BNDP = BndParticles_<typename PUSHP::Mparticles>,
	 typename PUSHF = PushFields<typename PUSHP::MfieldsState>,
	 typename BND = Bnd_<typename PUSHP::MfieldsState>,
	 typename MOMENT_N = ItemMomentLoopPatches<Moment_n_1st<typename PUSHP::Mparticles, _Mfields>>>
struct TestConfig
{
  using dim = DIM;
  using order = ORDER;
  using PushParticles = PUSHP;
  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using Mfields = _Mfields;
  using Checks = CHECKS;
  using BndParticles = BNDP;
  using PushFields = PUSHF;
  using Bnd = BND;
  using Moment_n = MOMENT_N;
};

using TestConfig2ndDouble = TestConfig<dim_xyz, MfieldsC,
				       PushParticles__<Config2ndDouble<dim_xyz>>,
				       checks_order_2nd>;
using TestConfig2ndDoubleYZ = TestConfig<dim_yz, MfieldsC,
					 PushParticles__<Config2ndDouble<dim_yz>>,
					 checks_order_2nd>;
using TestConfig2ndSingle = TestConfig<dim_xyz, MfieldsSingle,
				       PushParticles__<Config2nd<MparticlesSingle, MfieldsStateSingle, dim_xyz>>,
				       checks_order_2nd>;
using TestConfig1vbec3dSingle = TestConfig<dim_xyz, MfieldsSingle,
					   PushParticles1vb<Config1vbecSplit<MparticlesSingle, MfieldsStateSingle, dim_xyz>>,
					   checks_order_1st>;
using TestConfig1vbec3dSingleYZ = TestConfig<dim_yz, MfieldsSingle,
					     PushParticles1vb<Config1vbecSplit<MparticlesSingle, MfieldsStateSingle, dim_yz>>,
					     checks_order_1st>;
using TestConfig1vbec3dSingleXZ = TestConfig<dim_xz, MfieldsSingle,
					     PushParticles1vb<Config1vbecSplit<MparticlesSingle, MfieldsStateSingle, dim_xz>>,
					     checks_order_1st>;

#ifdef USE_VPIC
using TestConfigVpic = TestConfig<dim_xyz,
				  MfieldsSingle, // FIXME, this is not real nice, but might work...
				  PushParticlesVpic,
				  checks_order_1st,
				  Checks_<MparticlesVpic, MfieldsStateVpic, void, checks_order_1st>,
				  BndParticlesVpic>;
#endif

#ifdef USE_CUDA
using TestConfig1vbec3dCuda = TestConfig<dim_xyz, MfieldsCuda,
					 PushParticlesCuda<CudaConfig1vbec3dGmem<dim_xyz, BS144>>,
					 checks_order_1st,
					 ChecksCuda<MparticlesCuda<BS144>>,
					 BndParticlesCuda<MparticlesCuda<BS144>, dim_xyz>,
					 PushFieldsCuda,
					 BndCuda2<MfieldsStateCuda>>;
using TestConfig1vbec3dCuda444 = TestConfig<dim_xyz, MfieldsCuda,
					    PushParticlesCuda<CudaConfig1vbec3dGmem<dim_xyz, BS444>>,
					    checks_order_1st,
					    ChecksCuda<MparticlesCuda<BS444>>,
					    BndParticlesCuda<MparticlesCuda<BS444>, dim_xyz>,
					    PushFieldsCuda,
					    BndCuda3<MfieldsStateCuda>,
					    Moment_n_1st_cuda<MparticlesCuda<BS444>, dim_xyz>>;
using TestConfig1vbec3dCudaYZ = TestConfig<dim_yz, MfieldsCuda,
					   PushParticlesCuda<CudaConfig1vbec3d<dim_yz, BS144>>,
					   checks_order_1st,
					   ChecksCuda<MparticlesCuda<BS144>>,
					   BndParticlesCuda<MparticlesCuda<BS144>, dim_yz>,
					   PushFieldsCuda,
					   BndCuda<MfieldsStateCuda>,
					   Moment_n_1st_cuda<MparticlesCuda<BS144>, dim_yz>>;
#endif

// ======================================================================
// CurrentReference

struct CurrentReference
{
  int m;
  Int3 pos;
  double val;
};

// ======================================================================
// PushParticlesTest

template<typename T>
struct PushParticlesTest : ::testing::Test
{
  using dim = typename T::dim;
  using Mparticles = typename T::Mparticles;
  using MfieldsState = typename T::MfieldsState;
  using PushParticles = typename T::PushParticles;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename Mparticles::real_t;

  const real_t eps = 1e-5;
  const double L = 160;

  const real_t fnqx = .05, fnqy = .05, fnqz = .05;
  const real_t dx = 10., dy = 10., dz = 10.;
  
  Int3 ibn = { 2, 2, 2 };
  
  ~PushParticlesTest()
  {
    delete mprts;
    delete mflds;
  }

  void make_psc(const Grid_t::Kinds& kinds)
  {
    Int3 gdims = {16, 16, 16};
    if (dim::InvarX::value) { gdims[0] = 1; ibn[0] = 0; }
    if (dim::InvarY::value) { gdims[1] = 1; ibn[1] = 0; }
    if (dim::InvarZ::value) { gdims[2] = 1; ibn[2] = 0; }

    auto grid_domain = Grid_t::Domain{gdims, {L, L, L}};
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};
    
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 200;
    auto coeff = Grid_t::Normalization{norm_params};

    grid_ = Grid_t::psc_make_grid(grid_domain, grid_bc, kinds, coeff, 1., {});
  }
  
  const Grid_t& grid()
  {
    assert(grid_);
    return *grid_;
  }

  template<typename FUNC>
  void runSingleParticleTest(FUNC init_fields, particle_inject prt0, particle_inject prt1,
			     std::vector<CurrentReference> curr_ref = {})
  {
    auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
    make_psc(kinds);

    // init fields
    mflds = new MfieldsState{grid()};
    SetupFields<MfieldsState>::set(*mflds, init_fields);

    // init particle
    mprts = new Mparticles{grid()};
    {
      auto injector = (*mprts)[0].injector();
      injector(prt0);
    }
    
    //mprts->dump("mprts.dump");
  
    // do one step
    PushParticles pushp_;
    pushp_.push_mprts(*mprts, *mflds);

    // check against reference
    auto prt = *(*mprts)[0].get().begin();
    EXPECT_NEAR(prt.u()[0], prt1.u[0], eps);
    EXPECT_NEAR(prt.u()[1], prt1.u[1], eps);
    EXPECT_NEAR(prt.u()[2], prt1.u[2], eps);
    EXPECT_NEAR(prt.w(), prt1.w, eps);
    EXPECT_NEAR(prt.position()[0], prt1.x[0], eps);
    EXPECT_NEAR(prt.position()[1], prt1.x[1], eps);
    EXPECT_NEAR(prt.position()[2], prt1.x[2], eps);
    EXPECT_EQ(prt.kind(), prt1.kind);

    if (!curr_ref.empty()) {
      checkCurrent(curr_ref);
    }

    check_after_push(*mprts);
  }

  template<typename Mparticles>
  void check_after_push(Mparticles& mprts)
  {
  }

#ifdef USE_CUDA
  void check_after_push(MparticlesCuda<BS444>& mprts)
  {
    EXPECT_TRUE(mprts.check_after_push());
  }
#endif
  
  void checkCurrent(std::vector<CurrentReference>& curr_ref)
  {
    auto mflds_ref = MfieldsState{grid()};
    auto flds_ref = mflds_ref[0];
    for (auto& ref : curr_ref) {
      if (dim::InvarX::value) { ref.pos[0] = 0; }
      if (dim::InvarY::value) { ref.pos[1] = 0; }
      if (dim::InvarZ::value) { ref.pos[2] = 0; }
      flds_ref(ref.m, ref.pos[0], ref.pos[1], ref.pos[2]) += ref.val;
    }

    auto flds = (*mflds)[0];
    this->grid().Foreach_3d(2, 2, [&](int i, int j, int k) {
	for (int m = JXI; m <= JZI; m++) {
	  auto val = flds(m, i,j,k);
	  auto val_ref = flds_ref(m, i,j,k);
	  EXPECT_NEAR(val, val_ref, eps) << "ijk " << i << " " << j << " " << k << " m " << m;
	}
      });
  }
  
  template<typename particle_t>
  Vec3<double> push_x(const particle_t& prt0, particle_t& prt1)
  {
    Vec3<double> xi1 = { prt0.x[0] + vx(prt1),
			 prt0.x[1] + vy(prt1),
			 prt0.x[2] + vz(prt1) };
    
    if (!dim::InvarX::value) prt1.x[0] = xi1[0];
    if (!dim::InvarY::value) prt1.x[1] = xi1[1];
    if (!dim::InvarZ::value) prt1.x[2] = xi1[2];

    return xi1;
  }

  Grid_t* grid_ = {};
  Mparticles* mprts = {};
  MfieldsState* mflds = {};
};

