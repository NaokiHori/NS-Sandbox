#if !defined(DOMAIN_H)
#define DOMAIN_H

#if !defined(NX)
#error "NX is not defined"
#endif

#if !defined(NY)
#error "NY is not defined"
#endif

#if !defined(LX)
#error "LX is not defined"
#endif

#if !defined(LY)
#error "LY is not defined"
#endif

#include <stddef.h> // size_t
#include <stdbool.h> // true, false

#define X_PERIODIC true
#define Y_PERIODIC false

#define DX (LX / NX)
#define DY (LY / NY)

extern const size_t ux_imin;
extern const size_t uy_jmin;

#endif // DOMAIN_H
