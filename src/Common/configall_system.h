/*
 * This header file is included by the *Config.h in the individual
 * COIN packages when the code is compiled in a setting that doesn't
 * use the configure script (i.e., HAVE_CONFIG_H is not defined).
 * This header file includes the system and compile dependent header
 * file defining macros that depend on what compiler is used.
 */

#ifdef _MSC_VER
# include "configall_system_msc.h"
#else
# error "Trying to use configall_system for unknown compiler."
#endif
