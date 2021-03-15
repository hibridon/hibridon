c     ----------------------------------------------------------------- 
c     pot routine for the HCl-H2 PES computed by Lanza et al.,
c     jcp 140, 064316 (2014)
c
c     to be used with the 1Sigma-1Sigma basis routine
c     data file used:  pot_hclh2_lique.dat
c
c     pot matrix elements are computed assuming the angular dependence defined
c     in the appendix of green, jcp 62, 2271 (1975)
c
c     revised from pot_ohh2.f, written by qianli ma
c     current revision:  1-jun-2017 (p. dagdigian)
c     ------------------------------------------------------------------
c     Dummy subroutines for user-defined bases.
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
c     ----------------------------------------------------------------- 
C     Module containing the shared arrays for this pot routine
c
      module mod_hclh2
      implicit none
      integer :: nr1
      real(8), dimension(:), allocatable :: rr1
      real(8), dimension(:, :), allocatable :: coef1, spl_b1, 
     $  spl_c1, spl_d1
      real(8) econv, xmconv, sq4pi
      parameter (econv=219474.6315343234d0,
     $  xmconv=0.0005485799094979479d0, sq4pi=3.544907701811032d0)
      end module mod_hclh2
c     ----------------------------------------------------------------- 
c     loapot subroutine loads pot parameters for HCl-H2 interaction
c     ----------------------------------------------------------------- 
      subroutine loapot (iunit, file_name)
      use mod_1sg1sg
      use mod_hclh2
      implicit none
C     common/parbas is replaced by module ba2pi1sg to allow more
C     parameters be passed between the pot routine and the basis routine
      common /coptnm/ pot_name, pot_label
      common /conlam/ nlam
      character(48) :: pot_name, pot_label
      character*(*) :: file_name
      character(255) :: file_path
      integer :: iunit, ir, iv, nv, ifrsts, junks, nlam, nlammx, lamnum
      real(8), dimension(1) :: vvl
      common /covvl/ ifrsts, junks, vvl
C     A call to this subroutine with a string containing a space will be
C     made at the time hibridon loads. Input file is not available at
C     the time.
      if (file_name .eq. " ") return
      pot_name = 'LANZA ET AL. HCl-H2 CCSD(T) PES'
      call datfln(trim(file_name), file_path)
      open (unit=iunit, file=file_path, status="old")
C 
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
C     
      do iv = 1, nv
         read (iunit, *) lms(iv)%l1, lms(iv)%l2, lms(iv)%ltot
         read (iunit, *) (coef1(ir, iv), ir = 1, nr1)
      end do
c  fix normalizaotion of coeffs and convert to a.u.
      coef1 = coef1 * sq4pi / econv
C     
      close(unit=iunit)
C     Spline parameters prepared here
      do iv = 1, nv
         call spline(nr1, rr1, coef1(1, iv), spl_b1(1, iv), 
     $        spl_c1(1, iv), spl_d1(1, iv))
      end do
      return
      end subroutine loapot
C     ------------------------------------------------------------------
C     Main program for makepot
      subroutine driver
      use mod_1sg1sg
      use mod_hclh2
      implicit none
      common /conlam/ nlam
      common /covvl/ vvl
       real(8), dimension(1) :: vvl
      character(40), parameter :: data_file_name='pot_hclh2_lique.dat'
      real(8) :: r, vv0
      integer :: i, nv, nlam
      call loapot(10, data_file_name)
 10   print *, 'R (bohr), Ctrl+D to exit:'
      read (5, *, end=99) r
      call pot(vv0, r)
      write (6, 20) (lms(i)%l1, lms(i)%l2, lms(i)%ltot,
     $     vvl(i) * econv / sq4pi, i=1, nlam)
 20   format (3(3(i2, 1x), 1x, 1pe16.8, 2x))
      goto 10
C 99   return
 99   stop
      end subroutine driver
C     ------------------------------------------------------------------
      subroutine pot(vv0, r_inp)
      use mod_1sg1sg
      use mod_hclh2
      implicit none
      common /conlam/ nlam
      common /covvl/ vvl
      real(8), dimension(1) :: vvl
      double precision vv0, r_inp, r
      double precision seval
      integer iv, nlam
      vv0 = 0.d0
      if (r_inp .lt. 3.5d0) then
         r = 3.5d0
      else
         r = r_inp
      end if
      do iv = 1, nlam
         vvl(iv) = seval(nr1, r, rr1, coef1(1, iv), spl_b1(1, iv),
     $        spl_c1(1, iv), spl_d1(1, iv))
      end do
      return
      end subroutine pot
C     ------------------------------------------------------------------
      subroutine datfln(filenm, fullnm)
      character (len=*) :: filenm, fullnm
      fullnm = 'potdata/' // trim(filenm)
      return
      end subroutine datfln
c     ----------------------------------eof-----------------------------
