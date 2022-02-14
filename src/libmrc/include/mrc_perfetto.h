
#ifndef MRC_PERFETTO_H
#define MRC_PERFETTO_H

#include <mrc.h>

BEGIN_C_DECLS;

void perfetto_initialize();
void perfetto_finalize();
void perfetto_event_begin(const char* s);
void perfetto_event_end();

END_C_DECLS;

#endif
