
#include "output_fields.hxx"
#include "output_fields_c.hxx"

#include "psc_output_fields_private.h"

// ----------------------------------------------------------------------

using PscOutputFields_t = PscOutputFields<OutputFieldsC>;

// ----------------------------------------------------------------------
// OutputFieldsC_destroy

static void
OutputFieldsC_destroy(struct psc_output_fields *out)
{
  PscOutputFields_t outf{out};

  outf->~OutputFieldsC();
}

// ----------------------------------------------------------------------
// OutputFieldsC_setup

static void
OutputFieldsC_setup(struct psc_output_fields *out)
{
  PscOutputFields_t outf{out};

  OutputFieldsCParams prm = *outf.sub();
  new(outf.sub()) OutputFieldsC{psc_output_fields_comm(out), prm};
}

// ----------------------------------------------------------------------
// OutputFieldsC_run

static void
OutputFieldsC_run(struct psc_output_fields *out,
			MfieldsBase& mflds, MparticlesBase& mprts)
{
  PscOutputFields_t outf{out};

  (*outf.sub())(mflds, mprts);
}

// ======================================================================
// psc_output_fields: subclass "c"

#define VAR(x) (void *)offsetof(struct OutputFieldsC, x)

// FIXME pfield_out_[xyz]_{min,max} aren't for pfield only, better init to 0,
// use INT3

static struct param OutputFieldsC_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "output_fields"      , VAR(output_fields)        , PARAM_STRING("j,e,h")   },
  { "pfield_first"       , VAR(pfield_first)         , PARAM_INT(0)            },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)           },
  { "tfield_first"       , VAR(tfield_first)         , PARAM_INT(0)            },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)           },
  { "tfield_length"      , VAR(tfield_length)        , PARAM_INT(10)           },
  { "tfield_every"       , VAR(tfield_every)         , PARAM_INT(1)            },
  { "pfield_out_x_min"   , VAR(rn[0])                , PARAM_INT(0)            },  
  { "pfield_out_x_max"   , VAR(rx[0])                , PARAM_INT(1000000000)  },     // a big number to change it later to domain.ihi or command line number
  { "pfield_out_y_min"   , VAR(rn[1])                , PARAM_INT(0)           }, 
  { "pfield_out_y_max"   , VAR(rx[1])                , PARAM_INT(1000000000)  },
  { "pfield_out_z_min"   , VAR(rn[2])                , PARAM_INT(0)            }, 
  { "pfield_out_z_max"   , VAR(rx[2])                , PARAM_INT(1000000000)  },
  {},
};
#undef VAR

struct psc_output_fields_ops_c : psc_output_fields_ops {
  psc_output_fields_ops_c() {
    name                  = "c";
    size                  = sizeof(struct OutputFieldsC);
    param_descr           = OutputFieldsC_descr;
    setup                 = OutputFieldsC_setup;
    destroy               = OutputFieldsC_destroy;
    run                   = OutputFieldsC_run;
  }
} psc_output_fields_c_ops;
