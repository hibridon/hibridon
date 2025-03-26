!     ------------------------------------------------------------------
!     Pot routine for the 2013 Ma, Dagdigian, Klos, and Alexander
!     ab initio OH--H2 PES at MRCI/aug-cc-pVQZ level.
!     
!     Author: Qianli Ma
!     
!     To be used in combination of the 2Pi--1Sigma basis.
!     
!     Data file used: pot_ohh2.dat
!     ------------------------------------------------------------------
!     Dummy subroutines for user-defined bases.
#include "common/syusr.F90"
#include "common/bausr.F90"
#include "common/ground.F90"
!     ------------------------------------------------------------------
!     Module containing the shared arrays for this pot routine.
module pot_ohh2
use mod_2pi1sg
implicit none
integer :: nr
real(8), dimension(:), allocatable :: rr
real(8), dimension(:, :), allocatable :: coef, spl_b, spl_c, spl_d
real(8) econv, xmconv
parameter (econv=219474.6315343234d0, &
     xmconv=0.0005485799094979479d0)
end module pot_ohh2
!     ------------------------------------------------------------------
!     Main program for makepot
subroutine driver
use pot_ohh2
use mod_covvl, only: vvl
use mod_hipot, only: loapot, pot
implicit none
character(40), parameter :: data_file_name= &
    'pot_ohh2_mrci_core1_highsym.dat'
real(8) :: r, vv0
integer :: i
call loapot(10, data_file_name)
10 print *, 'R (bohr), Ctrl+D to exit:'
read (5, *, end=99) r
call pot(vv0, r)
write (6, 20) (lms(i)%l1, lms(i)%l2, lms(i)%l, lms(i)%is_diag, &
     vvl(i) * econv, i=1, nv)
20 format (3(3(i2, 1x), l1, 1x, 1pe16.8, 2x))
goto 10
! 99   return
99 stop
end subroutine driver
!     ------------------------------------------------------------------
!     Load the data file of the potential
subroutine loapot(iunit, filnam)
use pot_ohh2
use mod_parpot, only: pot_name, pot_label
implicit none
integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)    
!     common/parbas is replaced by module ba2pi1sg to allow more
!     parameters be passed between the pot routine and the basis routine
character(255) :: file_path
integer :: ir, iv
!     A call to this subroutine with a string containing a space will be
!     made at the time hibridon loads. Input file is not available at
!     the time.
if (filnam .eq. " ") return
!     
pot_name = 'MA-DAGDIGIAN-KLOS-ALEXANDER OH--H2 MRCI PES'
call datfln(trim(filnam), file_path)
open (unit=iunit, file=file_path, status="old")
!     
read (iunit, *) nr
if (allocated(rr)) deallocate(rr)
allocate(rr(nr))
read (iunit, *) rr
read (iunit, *) nv
if (allocated(lms)) deallocate(lms)
allocate(lms(nv))
if (allocated(spl_b)) deallocate(spl_b)
allocate(spl_b(nr, nv))
if (allocated(spl_c)) deallocate(spl_c)
allocate(spl_c(nr, nv))
if (allocated(spl_d)) deallocate(spl_d)
allocate(spl_d(nr, nv))
if (allocated(coef)) deallocate(coef)
allocate(coef(nr, nv))
!     
do iv = 1, nv
   read (iunit, *) lms(iv)%l1, lms(iv)%l2, lms(iv)%l, &
        lms(iv)%is_diag
   read (iunit, *) (coef(ir, iv), ir = 1, nr)
end do
coef = coef / econv
!     
close(unit=iunit)
!     Spline parameters prepared here
do iv = 1, nv
   call spline(nr, rr, coef(1, iv), spl_b(1, iv), spl_c(1, iv), &
        spl_d(1, iv))
end do
return
end subroutine loapot
!     ------------------------------------------------------------------
subroutine pot(vv0, r_inp)
use pot_ohh2
use mod_covvl, only: vvl
implicit none
double precision vv0, r_inp, r
double precision seval
integer iv
vv0 = 0d0
if (r_inp .lt. 3.5d0) then
   r = 3.5d0
else
   r = r_inp
end if
do iv = 1, nv
   vvl(iv) = seval(nr, r, rr, coef(1, iv), spl_b(1, iv), &
        spl_c(1, iv), spl_d(1, iv))
end do
return
end subroutine pot
!     ------------------------------------------------------------------

subroutine datfln(filenm, fullnm)
character (len=*) :: filenm, fullnm
fullnm = 'potdata/' // trim(filenm)
return
end subroutine datfln
