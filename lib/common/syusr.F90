!comdeck syusr
#include "hiutil.inc.F90"

subroutine syusr (irpot, readpt, iread)
!  dummy syusr subroutine
implicit none
interface
  subroutine loapot(iunit, filnam)
    integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
    character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)    
  end subroutine
end interface
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
