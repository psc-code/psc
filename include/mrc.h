
#ifndef MRC_H
#define MRC_H

#include <mrc_common.h>
#include <stdbool.h>

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

#endif
