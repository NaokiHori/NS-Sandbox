#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h> // size_t
#include <stdbool.h> // true, false

#define X_PERIODIC false
#define Y_PERIODIC true

#define NX 384
#define NY 128

#define LX 3.
#define LY 1.

#define DX (LX / NX)
#define DY (LY / NY)

extern const size_t ux_imin;
extern const size_t uy_jmin;

#endif // DOMAIN_H
