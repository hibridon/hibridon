!comdeck ground
#include "hiutil.inc.F90"
subroutine ground(wf, r, nch, nphoto, mxphot)
implicit none
real(8), intent(out) :: wf(nch*nphoto) ! array of dimension nch*nphoto, containing, on return,
! ground state wavefunction in each of nch components
! nphoto is number of difference ground state wavefunctions
real(8), intent(out) :: r  ! value of separation coordinate
integer, intent(in) :: nch  ! total number of channels (row dimension of q)
integer, intent(in) :: nphoto ! number of different wavefunctions calculated
! column index of q vector
integer, intent(in) :: mxphot  ! maximum size of q vector (mxphot .ge. nch*nphoto)

real(8), intent(in) :: yymin
integer, intent(in) :: nny

UNUSED(wf)
UNUSED(r)
UNUSED(nch)
UNUSED(nphoto)
UNUSED(mxphot)

entry wfintern(wf,yymin,nch,nny)

UNUSED(wf)
UNUSED(yymin)
UNUSED(nch)
UNUSED(nny)

return
end
