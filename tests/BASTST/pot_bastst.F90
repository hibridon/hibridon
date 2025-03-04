! Dummy potential for testing purposes
! the aim of this dummy potential is to provide a potential
! that works with each basis supported by hibridon
! important: this potentiel produces meaningless physical results!
#include "unused.h"
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
    UNUSED_DUMMY(iunit)
    UNUSED_DUMMY(filnam)
    potnam = 'DUMMY POTENTIAL FOR TESTING PURPOSES'

    call init_pot_parameters(ibasty)

    return
end subroutine loapot



subroutine init_pot_parameters(ibasty)
use mod_parbas, only: ntv, lammin, lammax, mproj, ivcol, ivrow
use mod_conlam, only: nlam, nlammx
use mod_cosysi, only: ispar
use mod_chiral, only: lms_chiral => lms
use mod_asymln, only: lms_asymln => lms
use mod_1sg1sg, only: lms_1sg1sg => lms
use mod_2pi1sg, only: lms_2pi1sg => lms
use mod_stp1sg, only: lms_stp1sg => lms
use mod_hiba19_sgpi1, only: numvib_sgpi1 => numvib, ivibpi_sgpi1 => ivibpi
implicit none
integer, intent(in) :: ibasty
    ivcol = 0
    ivrow = 0
    ntv = 1
    nlam = 1
    select case (ibasty)
        case (1)
            lammin(1) = 2
            lammax(1) = 2
        case (2)
            lammin(1) = 1
            lammax(1) = 10
            mproj(1) = 0
        case (3)
            lammin(1) = 1  ; lammin(2) = 2
            lammax(1) = 10 ; lammax(2) = 10
            mproj(1) = 0   ; mproj(2) = 2
        case (4)
            lammin(1) = 0 
            lammax(1) = 10 
            mproj(1) = 0  
        case (5)
            lammin(1) = 0 
            lammax(1) = 10 
            mproj(1) = 0
        case (7)
            lammin(1) = 1
            lammax(1) = 4
            mproj(1) = 0  
        case (8)
            lammin(1) = 2
            lammax(1) = 4
            mproj(1) = 0
        case (9)
            nlam = 45
            lammin(1) = 0
            lammax(1) = 6
            mproj(1) = 0
        case (10)
            lammin(1) = 1
            lammax(1) = 13
        case (11)
            lammin(1) = 1 ; lammin(2) = 4
            lammax(1) = 1 ; lammax(2) = 4
            mproj(1)  = 0 ; mproj(2)  = 4
        case (12)
            lammin(1)=1
            lammax(1)=9
            mproj(1)=0
        case (13)
            lammin(1)=1
            lammax(1)=9
            mproj(1)=0
        case (14)
            lammin(1) = 1 ; lammin(2) = 2
            lammax(1) = 1 ; lammax(2) = 2
            mproj(1)  = 0 ; mproj(2)  = 2
        case (15)
            nlam = 10
            lammin(1) = 1
            lammax(1) = 10
            mproj(1) = 0
        case (16)
            ispar(1) = 4
            nlam = 12
            nlammx = 12
            mproj(1) = 0
            mproj(2) = 1
            mproj(3) = 2
            mproj(4) = 3
            lammin(1) = 0
            lammin(2) = 1
            lammin(3) = 2
            lammin(4) = 3
            lammax(1) = 6
            lammax(2) = 5
            lammax(3) = 6
            lammax(4) = 5
        case (17)
            ispar(1) = 4
            nlam = 12
            nlammx = 12
            mproj(1) = 0
            mproj(2) = 1
            mproj(3) = 2
            mproj(4) = 3
            lammin(1) = 0
            lammin(2) = 1
            lammin(3) = 2
            lammin(4) = 3
            lammax(1) = 6
            lammax(2) = 5
            lammax(3) = 6
            lammax(4) = 5
        case (18)
            ispar(1) = 4
            nlam = 20
            nlammx = 20
            mproj(1) = 0
            mproj(2) = 3
            mproj(3) = 6
            mproj(4) = 9
            lammin(1) = 0
            lammin(2) = 3
            lammin(3) = 6
            lammin(4) = 9
            lammax(1) = 8
            lammax(2) = 9
            lammax(3) = 8
            lammax(4) = 9
        case (19)
            nlam = 4
            numvib_sgpi1 = 1
            ivibpi_sgpi1(1) = 0
            lammin = 1
            lammax = 1
        case(20)
            nlam = 3
            allocate(lms_2pi1sg(3))
            lms_2pi1sg(1)%l1 = 0; lms_2pi1sg(1)%l2 = 0 ; lms_2pi1sg(1)%l = 0 ; lms_2pi1sg(1)%is_diag = .true.
            lms_2pi1sg(2)%l1 = 0; lms_2pi1sg(2)%l2 = 2 ; lms_2pi1sg(2)%l = 2 ; lms_2pi1sg(2)%is_diag = .true.
            lms_2pi1sg(3)%l1 = 0; lms_2pi1sg(3)%l2 = 4 ; lms_2pi1sg(3)%l = 4 ; lms_2pi1sg(3)%is_diag = .true.
        case (21)
            nlam = 3
            allocate(lms_stp1sg(3))
            lms_stp1sg(1)%l1 = 0; lms_stp1sg(1)%mu1 = 0 ; lms_stp1sg(1)%l2 = 0 ; lms_stp1sg(1)%l = 0
            lms_stp1sg(2)%l1 = 0; lms_stp1sg(2)%mu1 = 0 ; lms_stp1sg(2)%l2 = 2 ; lms_stp1sg(2)%l = 2
            lms_stp1sg(3)%l1 = 0; lms_stp1sg(3)%mu1 = 0 ; lms_stp1sg(3)%l2 = 4 ; lms_stp1sg(3)%l = 4
        case (22)
            lammin(1) = 1
            lammax(1) = 12
            mproj(1) = 0
        case (23)
            lammin(1) = 1
            lammax(1) = 6
            mproj(1) = 0
        case (24)
            nlam = 4
            ispar(1) = 4
            lammin(1) = 3 ; lammax(1) = 3
            lammin(2) = 4 ; lammax(2) = 4
            lammin(3) = 6 ; lammax(3) = 6
            lammin(4) = 7 ; lammax(4) = 7
        case (25)
            nlam = 3
            allocate(lms_1sg1sg(nlam))
            lms_1sg1sg%l1 = 0 ; lms_1sg1sg%l2 = 0 ;lms_1sg1sg%ltot = 0
            lms_1sg1sg%l1 = 0 ; lms_1sg1sg%l2 = 2 ;lms_1sg1sg%ltot = 2
            lms_1sg1sg%l1 = 0 ; lms_1sg1sg%l2 = 4 ;lms_1sg1sg%ltot = 4
        case (26)
            nlam = 3
            allocate(lms_1sg1sg(nlam))
            lms_1sg1sg%l1 = 0 ; lms_1sg1sg%l2 = 0 ;lms_1sg1sg%ltot = 0
            lms_1sg1sg%l1 = 0 ; lms_1sg1sg%l2 = 2 ;lms_1sg1sg%ltot = 2
            lms_1sg1sg%l1 = 0 ; lms_1sg1sg%l2 = 4 ;lms_1sg1sg%ltot = 4
        case (27)
            nlam = 24
            lammin(1) = 0 ; lammin(2) = 2 ; lammin(3) = 4 ; lammin(4) = 6
            lammax(1) = 8 ; lammax(2) = 8 ; lammax(3) = 8 ; lammax(4) = 8
            mproj(1)  = 0 ; mproj(2)  = 2 ; mproj(3)  = 4 ; mproj(4)  = 6
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
use mod_selb, only: ibasty
    implicit none
    real(8), intent(out) :: vv0
    real(8), intent(in) :: r
    vv0  = 0d0
    vvl = 0d0

    ! set non-zero values of potential for atom-atom col
    select case (ibasty)
        case (22)
        vvl(6) = 2d-4; vvl(9) = 1d-4
        case (23)
        vvl(5) = 1d-4
    end select


    return
end subroutine pot
