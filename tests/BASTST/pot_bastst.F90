! Dummy potential for testing purposes

#include "common/syusr.F90"
#include "common/bausr.F90"
#include "common/ground.F90"

!  -----------------------------------------------------------------------
! subroutine called by testpot to interactively provide potential values
!  -----------------------------------------------------------------------
subroutine driver
use mod_covvl, only: vvl
use mod_parpot, only: potnam=>pot_name
implicit none
    write(*,*) 'This is a dummy potential for testing purposes'
    return
end subroutine driver

!  -----------------------------------------------------------------------
! subroutine to initialize the potential
!  -----------------------------------------------------------------------
subroutine loapot(iunit, filnam)
use mod_parpot, only: potnam=>pot_name
implicit none
    integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
    character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)
    potnam = 'DUMMY POTENTIAL FOR TESTING PURPOSES'
    return
end subroutine loapot

!  -----------------------------------------------------------------------
!  calculates the r-dependent coefficients in the collision of Ar with N2
!
!  on return:
!  vv0
!    contains the isotropic term (n=0) in the potential.
!  vvl 
!    the coefficients for each angular term in the coupling potential
!    [ vvl(i) for i = 1, nlam ] are returned in module mod_covvl
!    vvl(1) contains the anisotropic (n=2) term in the potential
!  -----------------------------------------------------------------------
subroutine pot (vv0, r)
use mod_covvl, only: vvl
use mod_conlam, only: nlam
    implicit none
    real(8), intent(out) :: vv0
    real(8), intent(in) :: r
    vv0  = 0d0
    vvl = 0d0
    return
end subroutine pot
