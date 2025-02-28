!system: N2-Ar using canonical Pattengill, LaBudde, Bernstein potential
! reference:  m.d. pattengill, r. a. labudde, r.b. bernstein,
!  and c.f. curtiss, j. chem. phys. 55, 5517 (1971)]

! declare unused user subroutines syusr, bausr and ground
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
    integer :: ios
    real(8) :: r, vv0
    write(6, *) potnam
    do while(.true.)
        write (6, *) ' r (bohr)'
        read (5, *, iostat=ios) r
        if (ios /= 0) exit
        call pot(vv0, r)  ! values of the potential are returned in vv0 and vvl
        write (6, "(' vsum', /, 7(1pe16.8))") vv0, vvl
    end do
end subroutine driver

!  -----------------------------------------------------------------------
! subroutine to initialize the potential
!  -----------------------------------------------------------------------
subroutine loapot(iunit, filnam)
use mod_parbas, only: ntv, ivcol, ivrow, lammin, lammax, mproj
use mod_parpot, only: potnam=>pot_name
implicit none
    integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
    character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)
    integer, parameter :: l1i = 1  ! lambda1 index
    potnam = 'PATTENGILL-LABUDDE-BERNSTEIN AR-N2'
    lammin(l1i) = 2 ; lammax(l1i) = 2  ! lambda1's range is [2, 2]
    mproj(l1i) = 0
    ntv(l1i) = 1
    ivcol(l1i, l1i) = 0 ; ivrow(l1i, l1i) = 0
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
    real(8) :: eps, a6, a12, r0, eps_au, rat, w(6), rr0 
    data eps, a6, a12, r0 / 83.05d0, 0.13d0, 0.5d0, 3.929d0/  ! from reference paper above
    w = (/ 1.0d0, 0.2d0, 0.1d0, 0.0d0, 0.025d0, 0.02d0/)  ! weights

    rr0 = r0 / 0.52917715d0  ! convert distance r0 to bohr
    eps_au = eps / 219474.6d0  ! convert energy eps to hartree

    rat = rr0 / r
    vv0   = eps_au * (rat**12 - 2.d0 * rat**6)
    vvl(1) = eps_au * (a12 * rat**12 - 2.d0 * a6 * rat**6)
    vvl(1:nlam) = vvl(1) * w(1:nlam)
    return
end subroutine pot
