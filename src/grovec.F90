#ifdef HAS_CONCATENATION_OPERATOR
#define CONCATENATE(a,b) a ## b
#else
#define IDENTITY(x) x
#define CONCATENATE(a,b) IDENTITY(a)IDENTITY(b)
#endif

module mod_grovec
  implicit none
  private

#define GROVEC_CLASS_NAME dgrovec
#define GROVEC_ELEMENT_TYPE real(8)
#include "grovec_type.F90"

#define GROVEC_CLASS_NAME igrovec
#define GROVEC_ELEMENT_TYPE integer
#include "grovec_type.F90"

contains

#define GROVEC_CLASS_NAME dgrovec
#define GROVEC_ELEMENT_TYPE real(8)
#include "grovec_imp.F90"

#define GROVEC_CLASS_NAME igrovec
#define GROVEC_ELEMENT_TYPE integer
#include "grovec_imp.F90"

end module mod_grovec

program test_grovec
use mod_grovec, only: dgrovec, igrovec
type(dgrovec) :: dg1
type(igrovec) :: ig1
dg1 = dgrovec(block_size=1024*1024, num_blocks=1024)
call dg1%set_element(int(1, 8), 3.d0)
call dg1%append(4.d0)
call dg1%set_element(int(1024, 8), 1024.d0)
call dg1%set_element(int(1024*1024+1, 8), 1025.d0)
write (6,*) 'dg1(1)=', dg1%get_element(1)
write (6,*) 'dg1(2)=', dg1%get_element(2)
write (6,*) 'dg1(1024)=', dg1%get_element(1024)
write (6,*) 'dg1(1024*1024+1)=', dg1%get_element(1024*1024+1)
ig1 = igrovec(block_size=1024, num_blocks=16)
call ig1%append(42)
write (6,*) 'ig1(1)=', ig1%get_element(1)
end program test_grovec