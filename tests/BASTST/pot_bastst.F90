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
use mod_chiral, only: lms_chiral => lms
use mod_asymln, only: lms_asymln => lms
use mod_1sg1sg, only: lms_1sg1sg => lms
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
            ntv(1) = 1 ; ntv(2) = 1
            lammin(1) = 1  ; lammin(2) = 2
            lammax(1) = 10 ; lammax(2) = 10
            mproj(1) = 0   ; mproj(2) = 2
        case (4)
            ntv = 1
            lammin(1) = 0 
            lammax(1) = 10 
            mproj(1) = 0  
        case (28)
            nlam = 3
            allocate(lms_1sg1sg(nlam))
            lms_1sg1sg(1)%l1 = 0 ; lms_1sg1sg(1)%l2 = 0 ;lms_1sg1sg(1)%ltot = 0
            lms_1sg1sg(2)%l1 = 0 ; lms_1sg1sg(2)%l2 = 2 ;lms_1sg1sg(2)%ltot = 2
            lms_1sg1sg(3)%l1 = 1 ; lms_1sg1sg(3)%l2 = 0 ;lms_1sg1sg(3)%ltot = 1
        case (29)
            nlam = 5
            allocate(lms_chiral(nlam))
            lms_chiral(1)%l1 = 0 ;lms_chiral(1)%m1 = 0
            lms_chiral(2)%l1 = 1 ;lms_chiral(2)%m1 = 1
            lms_chiral(3)%l1 = 2 ;lms_chiral(3)%m1 = 2
            lms_chiral(4)%l1 = 3 ;lms_chiral(4)%m1 = 0
            lms_chiral(5)%l1 = 4 ;lms_chiral(5)%m1 = 1
        case (30)
            nlam = 3
            allocate(lms_asymln(nlam))
            lms_asymln(1)%l1 = 0 ;lms_asymln(1)%m1 = 0 ; lms_asymln(1)%l2 = 0 ;lms_asymln(1)%ltot = 0
            lms_asymln(2)%l1 = 1 ;lms_asymln(2)%m1 = 1 ; lms_asymln(2)%l2 = 1 ;lms_asymln(2)%ltot = 2
            lms_asymln(3)%l1 = 2 ;lms_asymln(3)%m1 = 2 ; lms_asymln(3)%l2 = 2 ;lms_asymln(3)%ltot = 4
        
        case default
            write(*,*) 'ERROR: unknown ibasty for BASTST tests'
            call exit(1)
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
