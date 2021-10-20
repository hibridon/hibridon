#ifdef HAS_CONCATENATION_OPERATOR
#define CONCATENATE(a,b) a ## b
#else
#define IDENTITY(x) x
#define CONCATENATE(a,b) IDENTITY(a)IDENTITY(b)
#endif

! grovec : GROwable VECtor
! grovec provides a way to store a sparse array of size n in the form of blocks that are allocated when needed
module mod_grovec
  implicit none
  private

#define GROVEC_CLASS_NAME dgrovec_type
#define GROVEC_ELEMENT_TYPE real(8)
#include "grovec_type.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

#define GROVEC_CLASS_NAME igrovec_type
#define GROVEC_ELEMENT_TYPE integer
#include "grovec_type.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

contains

#define GROVEC_CLASS_NAME dgrovec_type
#define GROVEC_ELEMENT_TYPE real(8)
#include "grovec_imp.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

#define GROVEC_CLASS_NAME igrovec_type
#define GROVEC_ELEMENT_TYPE integer
#include "grovec_imp.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

end module mod_grovec
