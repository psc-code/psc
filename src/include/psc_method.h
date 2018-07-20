
#ifndef PSC_METHOD_H
#define PSC_METHOD_H

#include "psc.h"
#include "particles.hxx"
#include "fields3d.hxx"

MRC_CLASS_DECLARE(psc_method, struct psc_method);

void psc_method_output(struct psc_method *method, struct psc *psc,
		       int stats_every,
		       MfieldsBase& mflds, MparticlesBase& mprts);

#endif

