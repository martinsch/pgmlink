#ifndef PGMLINK_WINDOWS_H
#define PGMLINK_WINDOWS_H

// prevent the global namespace to become polluted with
// badly named Windows macros

#if defined(_WIN32)
# ifdef IN
#  undef IN
# endif
# ifdef OUT
#  undef OUT
# endif
#endif

#endif // PGMLINK_WINDOWS_H
