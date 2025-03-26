!comdeck syusr
#include "hiutil.inc.F90"

subroutine syusr (irpot, readpt, iread)
!  dummy syusr subroutine
use mod_hipot, only: loapot
implicit none
integer, intent(inout) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
character*(*), intent(in) :: fname
UNUSED(irpot)
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  return
endif
entry ptrusr (fname, readpt)
UNUSED(fname)
entry savusr (readpt)
entry chkusr
return
end
