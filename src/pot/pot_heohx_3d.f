c     pot_heohx_3d.f
c
c     Pot routine for OH(X 2Pi) + He, ab initio 3D PES's
c     PES calculated by K. Gubbels and coworkers
c
c     Reference:  K. B. Gubbels, Q. Ma, M. H. Alexander, P. J. Dagdigian, 
c     D. Tanis, G. C. Gronenboom, A. van der Avoird, and 
c     S. Y. T. van de Meerakker, J. Chem. Phys. 136, 144308 (2012).
c
c     author for hibridon compatible routine: Qianli Ma
c
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
c         pot_heohx_3d_data
c
c     Included header file (placed in hibxx/src/pot):
c     pot_heohx_3d_common.f
c
c     Dummy subroutines for using built-in basis routine
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
c
c
c     Begin of pot routine
c     *****************************************************************
      subroutine driver
      implicit double precision (a-h, o-z)
      include "pot_heohx_3d_common.f"
c
c     Main function for `makepot`, print vvl's interactively
c
      double precision r, vv0
      integer i
c
      write(6,*) 'He--OH(X) 3D potential'
 1    write(6,*) 'R (bohr), Ctrl+D to exit:'
      read (5, *, end=99) r
      call pot(vv0, r)
*  convert from atomic units for printout
      econv=219474.6d0
      write(6,115) vv0*econv,(vvl(i)*econv,i=1,CSLMAX)
115   format('vsum: ',14(1pe15.8))
      write(6,215) (vvl(i)*econv,i=CSLMAX+1,2*CSLMAX-1)
215   format('vdif: ',14(1pe15.8))
      goto 1
 99   end
*----------------------------------------------------------------------
c     end subroutine driver
c
c
c     *****************************************************************
      subroutine loapot(iunit, filnam)
      implicit double precision (a-h, o-z)
      include "common/parbas"
      include "common/parpot"
      include "pot_heohx_3d_common.f"
c
      integer iunit
      character*(*) filnam
c
c     Assign pot constants
c
      potnam = 'He--OH(X) RCCSDT/AV5Z 3D PES'
      lammin(1) = 1
      lammax(1) = CSLMAX
      lammin(2) = 2
      lammax(2) = CSLMAX
      lamnum(1) = CSLMAX
      lamnum(2) = CSLMAX - 1
      nlam = CSNVVL
      nlammx = nlam
      mproj(1) = 0
      mproj(2) = 2
      ntv(1) = 1
      ivcol(1, 1) = 0
      ivrow(1, 1) = 0
      return
      end
c     end subroutine loapot
c
c
c     *****************************************************************
      subroutine pot(vv0, r)
      implicit double precision (a-h, o-z)
      include "pot_heohx_3d_common.f"
c
      double precision vv0, r
c
      double precision rraw(CSNR), vraw(CSNR, CSNVVL + 1)
      double precision b(CSNR, CSNVVL + 1), c(CSNR, CSNVVL + 1),
     $     d(CSNR, CSNVVL + 1)
      integer i, ir, l, ilarr0, ilarr2
      double precision tmpr, tmpl, tmpv0, tmpv2
c
      double precision seval
c
      logical isfst
      save isfst
      data isfst /.true./
c
      if (isfst) then
         open (unit=10, file='potdata/pot_heohx_3d_data')
         do ir = 1, CSNR
            do l = 0, CSLMAX
               read (10, *) tmpr, tmpl, tmpv0, tmpv2
               ilarr0 = 1 + l
               ilarr2 = CSLMAX + l
               rraw(ir) = tmpr
               vraw(ir, ilarr0) = tmpv0
               if (l .ge. 2) vraw(ir, ilarr2) = tmpv2
            enddo
         enddo
c     Prepare for spline
         do i = 1, CSNVVL + 1
            call spline(CSNR, rraw, vraw(1, i), b(1, i), c(1, i),
     $           d(1, i))
         enddo
         close (10)
         isfst = .false.
      endif
c
      vv0 = seval(CSNR, r, rraw, vraw(1, 1), b(1, 1), c(1, 1), d(1, 1))
      do i = 2, CSNVVL + 1
         vvl(i - 1) = seval(CSNR, r, rraw, vraw(1, i), b(1, i),
     $        c(1, i), d(1, i))
      enddo
c
      return
      end
c     end subroutine pot
