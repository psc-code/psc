

#include "fields_item_jeh.hxx"
#include "fields_item_fields.hxx"

static FieldsItemFieldsOps<Item_j_nc> psc_output_fields_item_j_nc_ops;
static FieldsItemFieldsOps<Item_j_cc> psc_output_fields_item_j_ops;
static FieldsItemFieldsOps<Item_j_ec> psc_output_fields_item_j_ec_ops;

static FieldsItemFieldsOps<Item_e_nc> psc_output_fields_item_e_nc_ops;
static FieldsItemFieldsOps<Item_e_cc> psc_output_fields_item_e_ops;
static FieldsItemFieldsOps<Item_e_ec> psc_output_fields_item_e_ec_ops;

static FieldsItemFieldsOps<Item_h_nc> psc_output_fields_item_h_nc_ops;
static FieldsItemFieldsOps<Item_h_cc> psc_output_fields_item_h_ops;
static FieldsItemFieldsOps<Item_h_fc> psc_output_fields_item_h_fc_ops;

static FieldsItemFieldsOps<Item_poyn> psc_output_fields_item_poyn_ops;
static FieldsItemFieldsOps<Item_e2> psc_output_fields_item_e2_ops;
static FieldsItemFieldsOps<Item_h2> psc_output_fields_item_h2_ops;

static FieldsItemFieldsOps<Item_jdote> psc_output_fields_item_jdote_ops;
static FieldsItemFieldsOps<Item_divb> psc_output_fields_item_divb_ops;
static FieldsItemFieldsOps<Item_divj<MfieldsStateDouble, MfieldsC>> psc_output_fields_item_divj_ops;

// FIXME, ugly way of making sure the linker doesn't drop all our self-registering
// FieldsItems

void registerFieldsItemMoments1st();
void registerFieldsItemMoments1stNcSingle();
void registerFieldsItemMoments1stNcDouble();

void registerFieldsItemFields()
{
  registerFieldsItemMoments1st();
  registerFieldsItemMoments1stNcSingle();
  registerFieldsItemMoments1stNcDouble();
}

