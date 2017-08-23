
#if PSC_FIELDS_AS_SINGLE

#define PFX(x) psc_fields_single_ ## x
#define MPFX(x) psc_mfields_single_ ## x

#elif PSC_FIELDS_AS_C

#define PFX(x) psc_fields_c_ ## x
#define MPFX(x) psc_mfields_c_ ## x

#elif PSC_FIELDS_AS_FORTRAN

#define PFX(x) psc_fields_fortran_ ## x
#define MPFX(x) psc_mfields_fortran_ ## x

#endif
