c     ----------------------------------------------------------------- 
c     pot routine for the  propylene oxide - helium PES
c
c     The PES was computed at the CCSD(T)-F12/CBS(3,4) level.
c     The original fit is based on the interpolated moving least 
c     squares (IMLS) method by Richard Dawes.
c
c     Fit is a least-squares based on a random sampling of geometries 
c     (30 000 points by distance), as usually. The number of basis 
c     functions is 1681 (!). Due to the large anisotropy the fit is good 
c     only for interaction energies lower than ~ 500 cm-1. This means that 
c     rate coefficients can be computed for temperatures lower than ~ 70K.
c
c     to be used with the hibaastp2.f basis routine 
c
c     data file used:  /hibxx/bin/progs/potdata/pot_propoxid.dat
c     if hibridon is run from another directory, this directory should
c     include /potdata/pot_propoxid.dat
c
c     see notes by C. Rist and P. Dagdigian (jan-2019)
c
c     written by p. dagdigian
c     current revision:  12-feb-2019 (p. dagdigian)
c     ------------------------------------------------------------------
c     Dummy subroutines for user-defined bases.
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------
c     This module contains (explictly) the number of terms and their
c     indices in the expansion of the PES.  Its contents should be
c     filled in the pot routine.
c
c     This module replaces lammin, lammax, mproj in hibridon
c     commented ourt here since it is in hibaastp2.f
c
*      module mod_chiral
*      implicit none
c
*      type lm_type
*      integer :: l1, m1
*      end type lm_type
c
*      type(lm_type), dimension(:), allocatable :: lms
*      end module mod_chiral
c     ----------------------------------------------------------------- 
C     Module containing the shared arrays for this pot routine
c
      module mod_propoxid
      implicit none
      integer :: nr1
      real(8), dimension(:), allocatable :: rr1
      real(8), dimension(:, :), allocatable :: coef1, spl_b1, 
     $  spl_c1, spl_d1
      real(8) econv, xmconv, sq4pi
      parameter (econv=219474.6315343234d0,
     $  xmconv=0.0005485799094979479d0, sq4pi=3.544907701811032d0)
      end module mod_propoxid
c     ----------------------------------------------------------------- 
c     loapot subroutine loads pot parameters for propylene oxide - helium 
c     interaction
c     ----------------------------------------------------------------- 
      subroutine loapot (iunit, file_name)
      use mod_chiral
      use mod_propoxid
      implicit none
C     common/parbas is replaced by module ba1sg1sg to allow more
C     parameters be passed between the pot routine and the basis routine
      common /coptnm/ pot_name, pot_label
      common /conlam/ nlam
      character(48) :: pot_name, pot_label
      character*(*) :: file_name
      character(255) :: file_path
      integer :: iunit, ir, iv, nv, ifrsts, junks, nlam, nlammx, lamnum
      real(8), dimension(1) :: vvl
      common /covvl/ ifrsts, junks, vvl
C  A call to this subroutine with a string containing a space will be
C  made at the time hibridon loads. Input file is not available at
C  the time.
      if (file_name .eq. " ") return
      pot_name = 'Propylene oxide - He PES'
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
         read (iunit, *) lms(iv)%l1, lms(iv)%m1
         read (iunit, *) (coef1(ir, iv), ir = 1, nr1)
      end do
c  convert to a.u.
      coef1 = coef1 / econv
C     
      close(unit=iunit)
C  spline parameters prepared here
      do iv = 1, nv
         call spline(nr1, rr1, coef1(1, iv), spl_b1(1, iv), 
     $        spl_c1(1, iv), spl_d1(1, iv))
      end do
      return
      end subroutine loapot
C     ------------------------------------------------------------------
C     Main program for makepot
      subroutine driver
      use mod_chiral
      use mod_propoxid
      implicit none
      common /conlam/ nlam
      real(8), dimension(1) :: vvl
      common /covvl/ vvl
      character(40), parameter :: data_file_name='pot_propoxid.dat'
      real(8) :: r, vv0
      integer :: i, nv, nlam
      call loapot(10, data_file_name)
 10   print *, 'R (bohr), Ctrl+D to exit:'
      read (5, *, end=99) r
      call pot(vv0, r)
      write (6, 20) (lms(i)%l1, lms(i)%m1,
     $     vvl(i) * econv, i=1, nlam)
 20   format (3(2(i3, 1x), 1x, 1pe16.8, 2x))
      goto 10
 99   stop
      end subroutine driver
C     ------------------------------------------------------------------
      subroutine pot(vv0, r_inp)
      use mod_chiral
      use mod_propoxid
      implicit none
      common /conlam/ nlam
      real(8), dimension(1) :: vvl
      common /covvl/ vvl
      double precision vv0, r_inp, r, c6, vvx, switch
      double precision seval
      integer iv, nlam
      vv0 = 0.d0
      if (r_inp .lt. 4.5d0) then
         r = 4.5d0
         write(6,111)
111      format('r is too small.  reset to 4.5 bohr')
      else
         r = r_inp
      end if
c  switching function for long range
      switch = 0.5d0 * (tanh(0.5d0 * (r - 18.d0)) + 1.d0)
c  r^-6 fit of isotropic term [vvl(1)] at R = 20 bohr
      c6 = -8.8039e7 / econv
      do iv = 1, nlam
         vvl(iv) = seval(nr1, r, rr1, coef1(1, iv), spl_b1(1, iv),
     $        spl_c1(1, iv), spl_d1(1, iv))
         vvx = (1.d0 - switch) * vvl(iv)
         if (iv .eq. 1) then
c  merge isotropic term with asymptotic term
              vvl(iv) = vvx + switch * c6 / (r**6)
         else
              vvl(iv) = vvx
         end if
      end do
c
      return
      end subroutine pot
C     ------------------------------------------------------------------
      subroutine datfln(filenm, fullnm)
      character (len=*) :: filenm, fullnm
      fullnm = 'potdata/' // trim(filenm)
      return
      end subroutine datfln
c     ----------------------------------eof-----------------------------
