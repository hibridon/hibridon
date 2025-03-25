!     ----------------------------------------------------------------- 
!     pot routine for the CN-H2 PES computed by Kalugina et al.,
!     jcp 139, 074301 (2013)
!
!     to be used with the 1Sigma-1Sigma or 2Sigma-1Sigma basis routine
!     data file used:  pot_cnh2_lique.dat
!
!     pot matrix elements are computed assuming the angular dependence defined
!     in the appendix of green, jcp 62, 2271 (1975)
!
!     revised from pot_ohh2.f, written by qianli ma
!     current revision:  8-jun-2017 (p. dagdigian)
!     ------------------------------------------------------------------
!     Dummy subroutines for user-defined bases.
#include "common/syusr.F90"
#include "common/bausr.F90"
#include "common/ground.F90"
!     ----------------------------------------------------------------- 
!     Module containing the shared arrays for this pot routine
!
module mod_cnh2
implicit none
integer :: nr1
real(8), dimension(:), allocatable :: rr1
real(8), dimension(:, :), allocatable :: coef1, spl_b1, &
  spl_c1, spl_d1
real(8) econv, xmconv, sq4pi
parameter (econv=219474.6315343234d0, &
  xmconv=0.0005485799094979479d0, sq4pi=3.544907701811032d0)
end module mod_cnh2
!     ----------------------------------------------------------------- 
!     loapot subroutine loads pot parameters for CN-H2 interaction
!     ----------------------------------------------------------------- 
subroutine loapot (iunit, filnam)
use mod_1sg1sg
use mod_cnh2
use mod_conlam, only: nlam
use mod_parpot, only: pot_name
implicit none
integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)    
!     common/parbas is replaced by module ba1sg1sg to allow more
!     parameters be passed between the pot routine and the basis routine
character(255) :: file_path
integer :: ir, iv, nv
!     A call to this subroutine with a string containing a space will be
!     made at the time hibridon loads. Input file is not available at
!     the time.
if (filnam .eq. " ") return
pot_name = 'KALUGINA ET AL. CN-H2 RCCSD(T) PES'
call datfln(trim(filnam), file_path)
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
!  fix normalizaotion of coeffs and convert to a.u.
coef1 = coef1 * sq4pi / econv
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
use mod_cnh2
use mod_conlam, only: nlam
use mod_covvl, only: vvl
implicit none
character(40), parameter :: data_file_name='pot_cnh2_lique.dat'
real(8) :: r, vv0
integer :: i
call loapot(10, data_file_name)
10 print *, 'R (bohr), Ctrl+D to exit:'
read (5, *, end=99) r
call pot(vv0, r)
write (6, 20) (lms(i)%l1, lms(i)%l2, lms(i)%ltot, &
     vvl(i) * econv / sq4pi, i=1, nlam)
20 format (3(3(i2, 1x), 1x, 1pe16.8, 2x))
goto 10
! 99   return
99 stop
end subroutine driver
!     ------------------------------------------------------------------
subroutine pot(vv0, r_inp)
use mod_1sg1sg
use mod_cnh2
use mod_conlam, only: nlam
use mod_covvl, only: vvl
implicit none
double precision vv0, r_inp, r
double precision seval
integer iv
vv0 = 0.d0
if (r_inp .lt. 4.5d0) then
   r = 4.5d0
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
