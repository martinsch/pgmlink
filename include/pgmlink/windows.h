#ifndef PGMLINK_WINDOWS_H
#define PGMLINK_WINDOWS_H

// prevent the global namespace to become polluted with
// badly named Windows macros

#if defined(_WIN32)
# define VC_EXTRALEAN
# ifndef NOMINMAX
#  define NOMINMAX
#  define _PGMLINK_UNDEFINE_NOMINMAX
# endif
# include <windows.h>
# ifdef _PGMLINK_UNDEFINE_NOMINMAX
#  undef NOMINMAX
#  undef _PGMLINK_UNDEFINE_NOMINMAX
# endif
# ifdef DIFFERENCE
#  undef DIFFERENCE
# endif
# ifdef IN
#  undef IN
# endif
# ifdef OUT
#  undef OUT
# endif
#endif

#endif // PGMLINK_WINDOWS_H
