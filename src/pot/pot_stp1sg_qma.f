C     ------------------------------------------------------------------
C     Pot routine for the collision between a symmetric top and a 1Sigma
C     molecule.
C     
C     Author: Qianli Ma
C
C     The data file used is to be defined in the input file.
C     ------------------------------------------------------------------
C     Dummy subroutines for user-defined bases.
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
C     ------------------------------------------------------------------
C     Module containing the shared arrays for this pot routine.
      module pot_stp1sg_qma
      use mod_stp1sg
      implicit none
      integer :: nr
      real(8), dimension(:), allocatable :: rr
      real(8), dimension(:, :), allocatable :: coef, spl_b, spl_c, spl_d
      real(8) econv, xmconv
      parameter (econv=219474.6315343234d0,
     $     xmconv=0.0005485799094979479d0)
      end module pot_stp1sg_qma
C     ------------------------------------------------------------------
C     Main program for makepot
      subroutine driver
      use pot_stp1sg_qma
      implicit none
      common /covvl/ vvl
      real(8), dimension(1) :: vvl
      real(8) :: r, vv0
      integer :: i
      character(40) :: filenm
      print *, 'Please input the name of the potdata file:'
      read (5, *, end=99) filenm
      call loapot(10, filenm)
      do while (.true.)
         print *, 'R (bohr), Ctrl+D to exit:'
         read (5, *, end=99) r
         call pot(vv0, r)
         write (6, 20) (lms(i)%l1, lms(i)%l2, lms(i)%l, lms(i)%mu1,
     $        vvl(i) * econv, i=1, nv)
      end do
 20   format (3(4(i2, 1x), 1x, 1pe16.8, 2x))
* 99   return
 99   continue
*      end subroutine driver
      end
C     ------------------------------------------------------------------
C     Load the data file of the potential
      subroutine loapot(iunit, file_name)
      use pot_stp1sg_qma
      implicit none
C     common/parbas is replaced by module bastp1sg to allow more
C     parameters be passed between the pot routine and the basis routine
      common /coptnm/ pot_name, pot_label
      character(48) :: pot_name, pot_label
      character*(*) :: file_name
      character(255) :: file_path
      integer :: iunit, ir, iv
C     A call to this subroutine with a string containing a space will be
C     made at the time hibridon loads. Input file is not available at
C     the time.
      if (file_name .eq. " ") return
C     
      pot_name = 'SYMMETRIC TOP--1SIGMA GENERAL POT ROUTINE'
      call datfln(trim(file_name), file_path)
      open (unit=iunit, file=file_path, status="old")
C     
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
C     
      do iv = 1, nv
         read (iunit, *) lms(iv)%l1, lms(iv)%mu1, lms(iv)%l2, lms(iv)%l
         read (iunit, *) (coef(ir, iv), ir = 1, nr)
      end do
      coef = coef / econv
C     
      close(unit=iunit)
C     Spline parameters prepared here
      do iv = 1, nv
         call spline(nr, rr, coef(1, iv), spl_b(1, iv), spl_c(1, iv),
     $        spl_d(1, iv))
      end do
      return
      end subroutine loapot
C     ------------------------------------------------------------------
      subroutine pot(vv0, r_inp)
      use pot_stp1sg_qma
      implicit none
      common /covvl/ vvl
      real(8), dimension(1) :: vvl
      double precision vv0, r_inp, r
      double precision seval
      integer iv
      vv0 = 0d0
      r = r_inp
      if (r .lt. rr(1)) r = rr(1)
      do iv = 1, nv
         vvl(iv) = seval(nr, r, rr, coef(1, iv), spl_b(1, iv),
     $        spl_c(1, iv), spl_d(1, iv))
      end do
      return
      end subroutine pot
C     ------------------------------------------------------------------
      subroutine datfln(filenm, fullnm)
      character (len=*) :: filenm, fullnm
      fullnm = 'potdata/' // trim(filenm)
      return
      end subroutine datfln
