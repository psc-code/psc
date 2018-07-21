
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

#ifdef USE_CUDA
#include "../libpsc/cuda/push_particles_cuda_impl.hxx"
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#include "../libpsc/cuda/setup_particles_cuda.hxx"
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
// Getter etc

template<class Mparticles>
struct Getter
{
  Getter(Mparticles& mprts)
    : mprts_(mprts) {}

  typename Mparticles::patch_t& operator[](int p) { return mprts_[p]; }

private:
  Mparticles& mprts_;
};

#ifndef USE_CUDA

template<class Mparticles>
Getter<Mparticles> make_getter(Mparticles& mprts)
{
  return Getter<Mparticles>(mprts);
}

#endif

#ifdef USE_CUDA

template<class Mparticles>
struct GetterCuda
{
  GetterCuda(Mparticles& mprts)
    : mprts_(mprts.template get_as<MparticlesSingle>()) {}

  typename MparticlesSingle::patch_t& operator[](int p) { return mprts_[p]; }

private:
  MparticlesSingle& mprts_;
};

template<class Mparticles,
	 typename = typename std::enable_if<!(std::is_same<Mparticles, MparticlesCuda<BS144>>::value ||
					      std::is_same<Mparticles, MparticlesCuda<BS444>>::value)>::type>
Getter<Mparticles> make_getter(Mparticles& mprts)
{
  return Getter<Mparticles>(mprts);
}

template<class Mparticles,
	 typename = typename std::enable_if<std::is_same<Mparticles, MparticlesCuda<BS144>>::value ||
					    std::is_same<Mparticles, MparticlesCuda<BS444>>::value>::type>
GetterCuda<Mparticles> make_getter(Mparticles& mprts)
{
  return GetterCuda<Mparticles>(mprts);
}

#endif

// ======================================================================
// TestConfig

template<typename DIM, typename PUSHP, typename ORDER,
	 typename CHECKS = Checks_<typename PUSHP::Mparticles, typename PUSHP::Mfields, ORDER>,
	 typename BNDP = BndParticles_<typename PUSHP::Mparticles>,
	 typename PUSHF = PushFields<typename PUSHP::Mfields>,
	 typename BND = Bnd_<typename PUSHP::Mfields>,
	 typename MOMENT_N = ItemMomentLoopPatches<Moment_n_1st<typename PUSHP::Mparticles, typename PUSHP::Mfields>>>
struct TestConfig
{
  using dim = DIM;
  using order = ORDER;
  using PushParticles = PUSHP;
  using Mparticles = typename PushParticles::Mparticles;
  using Mfields = typename PushParticles::Mfields;
  using Checks = CHECKS;
  using BndParticles = BNDP;
  using PushFields = PUSHF;
  using Bnd = BND;
  using Moment_n = MOMENT_N;
};

using TestConfig2ndDouble = TestConfig<dim_xyz,
				       PushParticles__<Config2ndDouble<dim_xyz>>,
				       checks_order_2nd>;
using TestConfig2ndDoubleYZ = TestConfig<dim_yz,
					 PushParticles__<Config2ndDouble<dim_yz>>,
					 checks_order_2nd>;
using TestConfig2ndSingle = TestConfig<dim_xyz,
				       PushParticles__<Config2nd<MparticlesSingle, MfieldsSingle, dim_xyz>>,
				       checks_order_2nd>;
using TestConfig1vbec3dSingle = TestConfig<dim_xyz,
					   PushParticles1vb<Config1vbecSplit<MparticlesSingle, MfieldsSingle, dim_xyz>>,
					   checks_order_1st>;
using TestConfig1vbec3dSingleYZ = TestConfig<dim_yz,
					     PushParticles1vb<Config1vbecSplit<MparticlesSingle, MfieldsSingle, dim_yz>>,
					     checks_order_1st>;

#ifdef USE_CUDA
using TestConfig1vbec3dCuda = TestConfig<dim_xyz,
					 PushParticlesCuda<CudaConfig1vbec3dGmem<dim_xyz, BS144>>,
					 checks_order_1st,
					 ChecksCuda<MparticlesCuda<BS144>>,
					 BndParticlesCuda<MparticlesCuda<BS144>, dim_xyz>,
					 PushFieldsCuda,
					 BndCuda2<MfieldsCuda>>;
using TestConfig1vbec3dCuda444 = TestConfig<dim_xyz,
					    PushParticlesCuda<CudaConfig1vbec3dGmem<dim_xyz, BS444>>,
					    checks_order_1st,
					    ChecksCuda<MparticlesCuda<BS444>>,
					    BndParticlesCuda<MparticlesCuda<BS444>, dim_xyz>,
					    PushFieldsCuda,
					    BndCuda3<MfieldsCuda>,
					    Moment_n_1st_cuda<BS444, dim_xyz>>;
using TestConfig1vbec3dCudaYZ = TestConfig<dim_yz,
					   PushParticlesCuda<CudaConfig1vbec3d<dim_yz, BS144>>,
					   checks_order_1st,
					   ChecksCuda<MparticlesCuda<BS144>>,
					   BndParticlesCuda<MparticlesCuda<BS144>, dim_yz>,
					   PushFieldsCuda,
					   BndCuda,
					   Moment_n_1st_cuda<BS144, dim_yz>>;
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
  using Mfields = typename T::Mfields;
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
    ppsc = NULL; // FIXME, should use psc_destroy(ppsc), or really, get rid of ppsc...
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
    
    auto psc = psc_create(MPI_COMM_WORLD); // to create ppsc, mostly
    psc_default_dimensionless(psc);
    psc->coeff_ = psc_setup_coeff(psc);

    psc_setup_domain(psc, grid_domain, grid_bc, kinds, psc->coeff_, 1.);
  }
  
  const Grid_t& grid()
  {
    assert(ppsc);
    return ppsc->grid();
  }

  template<typename FUNC>
  void runSingleParticleTest(FUNC init_fields, particle_t prt0, particle_t prt1,
			     std::vector<CurrentReference> curr_ref = {})
  {
    auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
    make_psc(kinds);

    // init fields
    mflds = new Mfields{grid(), NR_FIELDS, ibn};
    SetupFields<Mfields>::set(*mflds, init_fields);

    // init particle
    mprts = new Mparticles{grid()};
    auto n_prts_by_patch = std::vector<uint>{1};
    SetupParticles<Mparticles>::setup_particles(*mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
	return prt0;
      });

    //mprts->dump("mprts.dump");
  
    // do one step
    PushParticles pushp_;
    pushp_.push_mprts(*mprts, *mflds);

    // check against reference
    for (auto& prt : make_getter(*mprts)[0]) {
      EXPECT_NEAR(prt.pxi, prt1.pxi, eps);
      EXPECT_NEAR(prt.pyi, prt1.pyi, eps);
      EXPECT_NEAR(prt.pzi, prt1.pzi, eps);
      EXPECT_NEAR(prt.qni_wni_, prt1.qni_wni_, eps);
      EXPECT_NEAR(prt.xi, prt1.xi, eps);
      EXPECT_NEAR(prt.yi, prt1.yi, eps);
      EXPECT_NEAR(prt.zi, prt1.zi, eps);
    }

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
    auto mflds_ref = Mfields{grid(), NR_FIELDS, ibn};
    auto flds_ref = mflds_ref[0];
    for (auto& ref : curr_ref) {
      if (dim::InvarX::value) { ref.pos[0] = 0; }
      flds_ref(ref.m, ref.pos[0], ref.pos[1], ref.pos[2]) = ref.val;
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
  
  Vec3<double> push_x(const particle_t& prt0, particle_t& prt1)
  {
    Vec3<double> xi1 = { prt0.xi + vx(prt1),
			 prt0.yi + vy(prt1),
			 prt0.zi + vz(prt1) };
    
    if (!dim::InvarX::value) prt1.xi = xi1[0];
    if (!dim::InvarY::value) prt1.yi = xi1[1];
    if (!dim::InvarZ::value) prt1.zi = xi1[2];

    return xi1;
  }

  Mparticles* mprts = {};
  Mfields* mflds = {};
};

