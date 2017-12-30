
#ifndef MRC_H
#define MRC_H

#include <mrc_common.h>
#include <stdbool.h>

#ifdef __cplusplus
# define BEGIN_C_DECLS extern "C" {
# define END_C_DECLS }
#else
# define BEGIN_C_DECLS /* empty */
# define END_C_DECLS /* empty */
#endif

BEGIN_C_DECLS

// ======================================================================
// global flags
//
// these are used to change fundamental libmrc behavior
// they can be combined using logical |, &, and set with
// the two functions below.

enum {
  // Normally, options for an mrc_obj should be prefixed with
  // the name of the object.
  // If an option is used without prefix, libmrc will print a warning
  MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING = 1,

  // Normally libmrc captures "--help" and shows help for all command
  // line options that would be recognized as it goes.
  // This flag disables the behavior and leaves it to the application
  // to handle --help.
  MRC_FLAG_IGNORE_OPTION_HELP                 = 2,
};

void mrc_set_flags(unsigned long flags);
void mrc_clear_flags(unsigned long flags);

void libmrc_finalize(void);

// private to libmrc

extern unsigned long mrc_flags;

END_C_DECLS

#endif
