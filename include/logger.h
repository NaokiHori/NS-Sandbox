#if !defined(LOGGER_H)
#define LOGGER_H

#include <stdio.h> // fprintf, stderr

#define LOGGER_FAILURE(...) \
  do { \
    fprintf( \
        stderr, \
        "[FAILED] (%s:%d): %s\n", \
        __FILE__, __LINE__, __VA_ARGS__ \
    ); \
  } while(0)

#endif // LOGGER_H
