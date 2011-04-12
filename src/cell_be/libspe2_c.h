
#ifndef LIBSPE2_H
#define LIBSPE2_H

#define SPE_MBOX_ANY_NONBLOCKING (2)
#define SPE_DEFAULT_ENTRY (3)

typedef int (*spe_program_handle_t)(unsigned long long, unsigned long long, unsigned long long );
typedef int spe_context_ptr_t;

int spe_in_mbox_write(int speid, unsigned int *cmd, int l, int flags);
int spe_out_mbox_status(int speid);
int spe_out_mbox_read(int speid, unsigned int *buf, int cnt);

int spe_context_create(int i, void *p);
int spe_context_run(int speid, unsigned int *entry, unsigned int runflags, void *argp,
		    void *p1, void *p2);
int spe_program_load(int speid, spe_program_handle_t *handle);

#endif
