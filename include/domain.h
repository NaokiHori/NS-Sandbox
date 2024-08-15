#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h> // size_t
#include <stdbool.h> // true, false

#define X_PERIODIC true
#define Y_PERIODIC false

#define NX 128
#define NY 384

#define LX 1.
#define LY 3.

#define DX (LX / NX)
#define DY (LY / NY)

extern const size_t ux_imin;
extern const size_t uy_jmin;

#endif // DOMAIN_H
