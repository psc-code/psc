
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "fields3d.hxx"
#include "fields_traits.hxx"

#include <psc_method.h>

#include "../libpsc/vpic/vpic_iface.h" // FIXME path


struct fields_vpic_t : fields3d<float, LayoutAOS>
{
  using Base = fields3d<float, LayoutAOS>;

  using Base::Base;
};

struct MfieldsStateVpic
{
  using real_t = float;
  using fields_t = fields_vpic_t;
  
  enum {
    EX = 0, EY = 1, EZ = 2, DIV_E_ERR = 3,
    BX = 4, BY = 5, BZ = 6, DIV_B_ERR = 7,
    TCAX = 8, TCAY = 9, TCAZ = 10, RHOB = 11,
    JFX = 12, JFY = 13, JFZ = 14, RHOF = 15,
    N_COMP = 20,
  };

  MfieldsStateVpic(const Grid_t& grid, const MaterialList& material_list, double damp = 0.)
    : grid_{grid}
  {
    assert(grid.n_patches() == 1);

    Simulation* sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    
    vmflds_fields_ = FieldArray::create(sim->vgrid_, material_list, damp);
  }

  const Grid_t& grid() const { return grid_; }
  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }
  Int3 ibn() const { return {1,1,1}; }
  
  fields_vpic_t operator[](int p)
  {
    assert(p == 0);
    int ib[3], im[3];
    float* data = vmflds().getData(ib, im);
    return fields_vpic_t(ib, im, N_COMP, data);
  }

  FieldArray& vmflds() { return *vmflds_fields_; }

  // static const Convert convert_to_, convert_from_;
  // const Convert& convert_to() override { return convert_to_; }
  // const Convert& convert_from() override { return convert_from_; }

private:
  FieldArray* vmflds_fields_;
  const Grid_t& grid_;
};

template<>
struct Mfields_traits<MfieldsStateVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
