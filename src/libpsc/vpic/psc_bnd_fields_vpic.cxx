
#include "bnd_fields.hxx"

struct BndFieldsVpic : BndFieldsBase
{
  void fill_ghosts_E(MfieldsBase& mflds_base) override {}
  void fill_ghosts_H(MfieldsBase& mflds_base) override {}
  void add_ghosts_J(MfieldsBase& mflds_base) override {}
};

// ----------------------------------------------------------------------
// psc_bnd_fields: subclass "vpic"

struct psc_bnd_fields_ops_vpic : psc_bnd_fields_ops {
  using Wrapper = PscBndFieldsWrapper<BndFieldsVpic>;
  psc_bnd_fields_ops_vpic() {
    name                  = "vpic";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_bnd_fields_vpic_ops;


