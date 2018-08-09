
#include <mrc_obj.h>

#pragma once

#include "psc.h"
#include "fields3d.hxx"

MRC_CLASS_DECLARE(psc_output_fields, struct psc_output_fields);

BEGIN_C_DECLS

void psc_output_fields_run(struct psc_output_fields *output_fields,
			   MfieldsBase& mflds, MparticlesBase& mprts);

END_C_DECLS
