
#include <gtest/gtest.h>

#include "fields3d.hxx"
#include "psc_fields_c.h"
#include "psc_fields_single.h"
#ifdef USE_CUDA
#include "psc_fields_cuda.h"
#include "psc_fields_cuda.inl"
#endif
#include "setup_fields.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#endif

#include "OutputFieldsDefault.h"
#include "writer_mrc.hxx"
#ifdef PSC_HAVE_ADIOS2
#include "writer_adios2.hxx"
#endif

#include "kg/io.h"
#include "fields3d.inl"

#include "psc.h" // FIXME, just for EX etc

#include <gtensor/reductions.h>

static Grid_t make_grid()
{
  auto domain =
    Grid_t::Domain{{8, 4, 2}, {80., 40., 20.}, {-40., -20., 0.}, {2, 1, 1}};
  auto bc = psc::grid::BC{};
  auto kinds = Grid_t::Kinds({Grid_t::Kind(1., 1., "test_species")});
  auto norm = Grid_t::Normalization{};
  double dt = .1;
  return Grid_t{domain, bc, kinds, norm, dt};
}

template <typename T>
class MfieldsTest : public ::testing::Test
{};

#ifdef USE_CUDA
using MfieldsTestTypes = ::testing::Types<MfieldsSingle, MfieldsC, MfieldsCuda>;
#else
using MfieldsTestTypes = ::testing::Types<MfieldsSingle, MfieldsC>;
#endif

TYPED_TEST_SUITE(MfieldsTest, MfieldsTestTypes);

#ifdef PSC_HAVE_ADIOS2

TYPED_TEST(MfieldsTest, WriteRead)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, {}};

  setupFields(mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  auto io = kg::io::IOAdios2{};
  {
    auto writer = io.open("test.bp", kg::io::Mode::Write);
    writer.put("mflds", mflds);
    writer.close();
  }

  auto mflds2 = Mfields{grid, NR_FIELDS, {}};
  {
    auto reader = io.open("test.bp", kg::io::Mode::Read);
    reader.get("mflds", mflds2);
    reader.close();
  }

  EXPECT_EQ(gt::norm_linf(view_interior(mflds.gt(), mflds.ibn()) -
                          view_interior(mflds2.gt(), mflds2.ibn())),
            0);
}

TYPED_TEST(MfieldsTest, WriteWithGhostsRead)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, {2, 2, 2}};

  setupFields(mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  auto io = kg::io::IOAdios2{};

  {
    auto writer = io.open("test.bp", kg::io::Mode::Write);
    writer.put("mflds", mflds);
    writer.close();
  }

  auto mflds2 = Mfields{grid, NR_FIELDS, {}};
  {
    auto reader = io.open("test.bp", kg::io::Mode::Read);
    reader.get("mflds", mflds2);
    reader.close();
  }

  EXPECT_EQ(gt::norm_linf(view_interior(mflds.gt(), mflds.ibn()) -
                          view_interior(mflds2.gt(), mflds2.ibn())),
            0);
}

TYPED_TEST(MfieldsTest, WriteReadWithGhosts)
{
  using Mfields = TypeParam;

  auto grid = make_grid();
  auto mflds = Mfields{grid, NR_FIELDS, {}};

  setupFields(mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  auto io = kg::io::IOAdios2{};

  {
    auto writer = io.open("test.bp", kg::io::Mode::Write);
    writer.put("mflds", mflds);
    writer.close();
  }

  auto mflds2 = Mfields{grid, NR_FIELDS, {2, 2, 2}};
  {
    auto reader = io.open("test.bp", kg::io::Mode::Read);
    reader.get("mflds", mflds2);
    reader.close();
  }

  EXPECT_EQ(gt::norm_linf(view_interior(mflds.gt(), mflds.ibn()) -
                          view_interior(mflds2.gt(), mflds2.ibn())),
            0);
}

#endif

TYPED_TEST(MfieldsTest, OutputFieldsMRC)
{
  using Mfields = TypeParam;
  using real_t = typename Mfields::real_t;
  using Mparticles = MparticlesSimple<ParticleSimple<real_t>>;

  auto grid = make_grid();
  grid.ibn = {2, 2, 2};
  auto mflds = Mfields{grid, NR_FIELDS, {2, 2, 2}};
  auto mprts = Mparticles{grid};

  setupFields(mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  OutputFieldsParams outf_params{};
  OutputFieldsItemParams outf_item_params{};
  outf_item_params.pfield_interval = 1;
  outf_item_params.tfield_interval = 0;
  outf_item_params.tfield_average_every = 40;
  outf_params.fields = outf_item_params;
  outf_params.moments = outf_item_params;
  OutputFields<Mfields, Mparticles, dim_xyz, WriterMRC> outf{grid, outf_params};

  outf(mflds, mprts);
}

#ifdef PSC_HAVE_ADIOS2

TYPED_TEST(MfieldsTest, OutputFieldsADIOS2)
{
  using Mfields = TypeParam;
  using real_t = typename Mfields::real_t;
  using Mparticles = MparticlesSimple<ParticleSimple<real_t>>;

  auto grid = make_grid();
  grid.ibn = {2, 2, 2};
  auto mflds = Mfields{grid, NR_FIELDS, {2, 2, 2}};
  auto mprts = Mparticles{grid};

  setupFields(mflds, [](int m, double crd[3]) {
    return m + crd[0] + 100 * crd[1] + 10000 * crd[2];
  });

  OutputFieldsParams outf_params{};
  OutputFieldsItemParams outf_item_params{};
  outf_item_params.pfield_interval = 1;
  outf_item_params.tfield_interval = 0;
  outf_item_params.tfield_average_every = 40;
  outf_params.fields = outf_item_params;
  outf_params.moments = outf_item_params;
  OutputFields<Mfields, Mparticles, dim_xyz, WriterADIOS2> outf{grid,
                                                                outf_params};

  outf(mflds, mprts);
}

#endif

// ======================================================================
// main

int main(int argc, char** argv)
{
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
