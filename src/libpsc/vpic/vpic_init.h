
#ifndef VPIC_INIT_H
#define VPIC_INIT_H

#include "vpic.h"

void user_init(vpic_simulation *simulation, const vpic_params *vpic_prm,
	       const vpic_harris_params *vpic_harris_prm);

void vpic_simulation_diagnostics(vpic_simulation *simulation);

#endif
