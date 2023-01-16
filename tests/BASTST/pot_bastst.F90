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
use mod_selb, only: ibasty
use mod_parpot, only: potnam=>pot_name
implicit none

    integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
    character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)
    potnam = 'DUMMY POTENTIAL FOR TESTING PURPOSES'

    call init_pot_parameters(ibasty)

    return
end subroutine loapot



subroutine init_pot_parameters(ibasty)
use mod_parbas, only: ntv, lammin, lammax, mproj, ivcol, ivrow
use mod_conlam, only: nlam
implicit none
integer, intent(in) :: ibasty
    ivcol = 0
    ivrow = 0
    select case (ibasty)
        case (1)
            ntv(1) = 1
            lammin(1) = 2
            lammax(1) = 2
        case (2)
            ntv(1) = 1
            lammin(1) = 1
            lammax(1) = 10
            mproj(1) = 0
        case (3)
            ntv(1) = 1
            lammin(1) = 1  ; lammin(2) = 2
            lammax(1) = 10 ; lammax(2) = 10
            mproj(1) = 0   ; mproj(2) = 2
    
        case default
            write(*,*) 'ERROR: unknown ibasty for BASTST tests'
            !call exit(1)
    end select
end subroutine init_pot_parameters


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
