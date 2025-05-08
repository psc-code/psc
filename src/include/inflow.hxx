#pragma once

#include "pushp.hxx"
#include "setup_particles.hxx"

template <typename MP, typename DIM>
class Inflow
{
public:
  using Mparticles = MP;
  using real_t = typename Mparticles::real_t;
  using Dim = DIM;

  // offset_in_cell: () -> double
  Inflow(const Grid_t& grid, psc_particle_npt npt,
         double (*offset_in_cell_dist)())
    : setup_particles_(grid),
      advance_(grid.dt),
      // FIXME np.p is a std::function and is called many times; better to use a
      // lambda
      np_{npt.kind, npt.n, setup_particles_.createMaxwellian(npt)},
      offset_in_cell_dist_(offset_in_cell_dist)
  {
    assert(INJECT_DIM_IDX_ >= 0);
  }

  auto get_advanced_prt(Double3 pos, real_t wni)
  {
    auto prt = setup_particles_.setupParticle(np_, pos, wni);

    auto v = advance_.calc_v(prt.u);
    advance_.push_x(prt.x, v, 1.0);

    return prt;
  }

  template <typename Injector>
  void inject_into_boundary_cell(const Grid_t& grid, Injector& injector,
                                 Int3 boundary_cell_global_idx)
  {
    assert(boundary_cell_global_idx[INJECT_DIM_IDX_] == 0);
    boundary_cell_global_idx[INJECT_DIM_IDX_] = -1;

    int n_in_cell = setup_particles_.get_n_in_cell(np_);
    double wni = setup_particles_.getWeight(np_, n_in_cell);

    for (int cnt = 0; cnt < n_in_cell; cnt++) {
      Double3 offset = {offset_in_cell_dist_(), offset_in_cell_dist_(),
                        offset_in_cell_dist_()};
      auto pos = (Double3(boundary_cell_global_idx) + offset) * grid.domain.dx +
                 grid.domain.corner;
      auto prt = get_advanced_prt(pos, wni);

      if (prt.x[INJECT_DIM_IDX_] < grid.domain.corner[INJECT_DIM_IDX_]) {
        continue;
      }

      injector(prt);
    }
  }

  template <typename Injector>
  void inject_into_boundary_patch(const Grid_t& grid, Injector& injector,
                                  const Grid_t::Patch& boundary_patch)
  {
    Int3 ilo = boundary_patch.off;
    Int3 ihi = ilo + grid.ldims;

    assert(ilo[INJECT_DIM_IDX_] == 0);

    int dim1 = std::min((INJECT_DIM_IDX_ + 1) % 3, (INJECT_DIM_IDX_ + 2) % 3);
    int dim2 = std::max((INJECT_DIM_IDX_ + 1) % 3, (INJECT_DIM_IDX_ + 2) % 3);

    for (Int3 cell_idx = ilo; cell_idx[dim1] < ihi[dim1]; cell_idx[dim1]++) {
      for (cell_idx[dim2] = ilo[dim2]; cell_idx[dim2] < ihi[dim2];
           cell_idx[dim2]++) {
        inject_into_boundary_cell(grid, injector, cell_idx);
      }
    }
  }

  template <typename Injectors>
  void inject_into_boundary(const Grid_t& grid, Injectors& injectors_by_patch)
  {
    for (int patch_idx = 0; patch_idx < grid.n_patches(); patch_idx++) {
      const auto& patch = grid.patches[patch_idx];

      if (patch.off[INJECT_DIM_IDX_] != 0) {
        continue;
      }

      auto& injector = injectors_by_patch[patch_idx];
      inject_into_boundary_patch(grid, injector, patch);
    }
  }

  template <typename MfieldsState>
  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    const Grid_t& grid = mprts.grid();
    auto injectors_by_patch = mprts.injector();

    inject_into_boundary(grid, injectors_by_patch);

    // TODO update j
  }

  AdvanceParticle<real_t, Dim> advance_;
  SetupParticles<Mparticles> setup_particles_;
  psc_particle_np np_;
  double (*offset_in_cell_dist_)();
  const static int INJECT_DIM_IDX_ = (!Dim::InvarX::value   ? 0
                                      : !Dim::InvarY::value ? 1
                                      : !Dim::InvarZ::value ? 2
                                                            : -1);
};
