!     ----------------------------------------------------------------- 
!     pot routine for the RCCSD(T)-F12a/aVTZ NH - H2 PES 
!     computed by Y. Kalugina
!
!     points fit with LAmax = 10, LBmax = 4
!
!     to be used with the 3Sigma-1Sigma basis routine (in hiba3sg1sg.f)
!     data file used:  pot_nhh2_10_4.dat
!
!     pot matrix elements are computed assuming the angular dependence defined
!     in the appendix of green, jcp 62, 2271 (1975)
!
!     revised from pot_c2hh2.f
!     current revision:  27-jan-2021 (p. dagdigian)
!     ------------------------------------------------------------------
!     Dummy subroutines for user-defined bases.
#include "common/syusr.F90"
#include "common/bausr.F90"
#include "common/ground.F90"
!     ----------------------------------------------------------------- 
!     Module containing the shared arrays for this pot routine
!
module mod_nhh2
implicit none
integer :: nr1
real(8), dimension(:), allocatable :: rr1
real(8), dimension(:, :), allocatable :: coef1, spl_b1, &
  spl_c1, spl_d1
real(8) econv, xmconv, sq4pi
parameter (econv=219474.6315343234d0, &
  xmconv=0.0005485799094979479d0, sq4pi=3.544907701811032d0)
end module mod_nhh2
!     ----------------------------------------------------------------- 
!     loapot subroutine loads pot parameters for NH - H2 interaction
!     ----------------------------------------------------------------- 
subroutine loapot (iunit, file_name)
use mod_1sg1sg
use mod_nhh2
use mod_conlam, only: nlam
use mod_covvl, only: vvl
use mod_parpot, only: pot_name, pot_label
implicit none
!     common/parbas is replaced by module ba1sg1sg to allow more
!     parameters be passed between the pot routine and the basis routine
character*(*) :: file_name
character(255) :: file_path
integer :: iunit, ir, iv, nv
!     A call to this subroutine with a string containing a space will be
!     made at the time hibridon loads. Input file is not available at
!     the time.
if (file_name .eq. " ") return
pot_name = 'NH - H2 PES'
call datfln(trim(file_name), file_path)
open (unit=iunit, file=file_path, status="old")
! 
read (iunit, *) nr1
if (allocated(rr1)) deallocate(rr1)
allocate(rr1(nr1))
read (iunit, *) rr1
read (iunit, *) nv
nlam = nv
if (allocated(lms)) deallocate(lms)
allocate(lms(nv))
if (allocated(spl_b1)) deallocate(spl_b1)
allocate(spl_b1(nr1, nv))
if (allocated(spl_c1)) deallocate(spl_c1)
allocate(spl_c1(nr1, nv))
if (allocated(spl_d1)) deallocate(spl_d1)
allocate(spl_d1(nr1, nv))
if (allocated(coef1)) deallocate(coef1)
allocate(coef1(nr1, nv))
!     
do iv = 1, nv
   read (iunit, *) lms(iv)%l1, lms(iv)%l2, lms(iv)%ltot
   read (iunit, *) (coef1(ir, iv), ir = 1, nr1)
end do
!******
!  fix normalization of coeffs and convert to a.u.
!      coef1 = coef1 * sq4pi / econv
!******
!  convert to a.u.
coef1 = coef1 / econv
!     
close(unit=iunit)
!     Spline parameters prepared here
do iv = 1, nv
   call spline(nr1, rr1, coef1(1, iv), spl_b1(1, iv), &
        spl_c1(1, iv), spl_d1(1, iv))
end do
return
end subroutine loapot
!     ------------------------------------------------------------------
!     Main program for makepot
subroutine driver
use mod_1sg1sg
use mod_nhh2
use mod_conlam, only: nlam
use mod_covvl, only: vvl
implicit none
character(40), parameter :: data_file_name= &
  'pot_nhh2_10_4.dat'
real(8) :: r, vv0
integer :: i, nv
call loapot(10, data_file_name)
10 print *, 'R (bohr), Ctrl+D to exit:'
read (5, *, end=99) r
call pot(vv0, r)

!******
!      write (6, 20) (lms(i)%l1, lms(i)%l2, lms(i)%ltot,
!     $     vvl(i) * econv / sq4pi, i=1, nlam)
!******

write (6, 20) (lms(i)%l1, lms(i)%l2, lms(i)%ltot, &
     vvl(i) * econv, i=1, nlam)
20 format (3(3(i2, 1x), 1x, 1pe16.8, 2x))
goto 10
! 99   return
99 stop
end subroutine driver
!     ------------------------------------------------------------------
subroutine pot(vv0, r_inp)
use mod_1sg1sg
use mod_nhh2
use mod_conlam, only: nlam
use mod_covvl, only: vvl
implicit none
double precision vv0, r_inp, r
double precision seval
integer iv
vv0 = 0.d0
if (r_inp .lt. 3.75d0) then
   r = 3.75d0
else
   r = r_inp
end if
do iv = 1, nlam
   vvl(iv) = seval(nr1, r, rr1, coef1(1, iv), spl_b1(1, iv), &
        spl_c1(1, iv), spl_d1(1, iv))
end do
return
end subroutine pot
!     ------------------------------------------------------------------
subroutine datfln(filenm, fullnm)
character (len=*) :: filenm, fullnm
fullnm = 'potdata/' // trim(filenm)
return
end subroutine datfln
!     ----------------------------------eof-----------------------------
