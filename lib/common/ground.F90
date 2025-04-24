!comdeck ground
#include "hiutil.inc.F90"
subroutine ground(wf, r, nch, nphoto, mxphot)
implicit none
real(8), intent(out) :: wf(nch*nphoto) ! array of dimension nch*nphoto, containing, on return,
! ground state wavefunction in each of nch components
! nphoto is number of difference ground state wavefunctions
real(8), intent(in) :: r  ! value of separation coordinate
integer, intent(in) :: nch  ! total number of channels (row dimension of q)
integer, intent(in) :: nphoto ! number of different wavefunctions calculated
! column index of q vector
integer, intent(in) :: mxphot  ! maximum size of q vector (mxphot .ge. nch*nphoto)

UNUSED(wf)
wf(1) = 0.0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(r)
UNUSED(nch)
UNUSED(nphoto)
UNUSED(mxphot)
end subroutine

subroutine wfintern(wf, yymin, nch, nphoto, nny, ifull)
implicit none
real(8), intent(out) :: wf(nch*nphoto) ! array of dimension nch*nphoto, containing, on return,
! ground state wavefunction in each of nch components
! nphoto is number of difference ground state wavefunctions
real(8), intent(in) :: yymin
integer, intent(in) :: nch  ! total number of channels (row dimension of q)
integer, intent(in) :: nphoto ! number of different wavefunctions calculated
! column index of q vector
integer, intent(in) :: nny
logical, intent(in) :: ifull

UNUSED(wf)
wf(1) = 0.0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(yymin)
UNUSED(nch)
UNUSED(nny)
UNUSED(ifull)

return
end subroutine
