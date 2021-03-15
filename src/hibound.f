c     hibound library, subroutines included:
c
c     1. bound        Bound state program
c     2. vmat_bound   Auxiliary subroutine for bound
c     3. gauger       Nodes and weights for Gauss-Hermite quadrature
c
c     ------------------------------------------------------------------
      subroutine bound(nch, nmax)
c     Bound state program, to be called by propag in hibrid3.f
c
c     Author: Susan K Gregurick, 3-aug-92
c     Revision: Moonbong Yang, 27-Sep-94
c     Completely rewritten: mha, 24-apr-1997
c     Revision: mha 11-may-1997
c     Revision: 10-jun-2006
c
c     Major revision: Q. Ma, 27-jun-2013 (dynamically allocated memory;
c     more efficient alogorithm; allow using Gauss-Hermite quadrature;
c     reduced memory usage; limit the number of positive eigenvalue in
c     output; bug fix on the dimension of w array; and other minor
c     improvements)
c
c     Revision: Q. Ma, 08-may-2014 (Query scratch array size before
c     diagonalizing matrices)
c
c     This subroutine may require large amount of memory. For larger
c     systems, it is important to link to the lapack library supporting
c     64-bit integers.
c
c     Definition of variables in call list:
c     nch     Number of coupled equations
c     nmax    Leading dimension of the w matrix used in Hibridon
c     ------------------------------------------------------------------
      implicit none
      integer, intent(in) :: nch, nmax
      include "common/parpot"
      common /corpar/ r1, r2, c, spac, delr, hsimp, eigmin, tolai, xmu
      real(8) :: r1, r2, c, spac, delr, hsimp, eigmin, tolai, xmu
      common /cofile/ input, output, jobnam, savfil
      character(40) :: jobnam, input, output, savfil
      common /colpar/ lpar(26), wavefl
      logical :: lpar, wavefl
      common /coconv/ econv, xmconv
      common /coered/ ered, rmu
      real(8) :: econv, xmconv, ered, rmu
      common /coiout/ niout
      integer :: niout
c
      real(8), parameter :: pi=dacos(-1d0)
c     maxpos is the maximum number of positive energies to print.
c     maxchn is the maximum number of channels to print in analyzing
c     wave functions when wavefl is set true.
c     maxvec is the maximum number of wave vectors to print.
      integer, parameter :: maxpos=20, maxchn=10, maxvec=20, ione=1
      character(40) :: evalfil
      integer :: vmax, ndim, i, j, k, irow, icol, nr, ir,
     $     nbound, lwork, ifail, info, m, lenx, liwork, irow1, irow2,
     $     icol1, icol2, ngh, i1, j1, il, iu
      real(8) :: del, alph, hexac, rnow, wt, xij, xfact, ratio
      integer, dimension(1) :: iworks
c
      real(8), dimension(:, :), allocatable :: cchn, schn, tchn, s, h,
     $     wr, z
      real(8), dimension(:), allocatable :: w, work, xgh, wgh, rgh,
     $     wgh1, g
      integer, dimension(:), allocatable :: iwork
c
      vmax = int((r2 - r1) / spac + 1.d-8)
      del = (r2 - r1) / dble(vmax)
      alph = (c / del) ** 2
      vmax = vmax + 1
      ndim = nch * vmax
c     Open file for eigenvalues
      call gennam (evalfil, jobnam, 0, 'evl', lenx)
      call openf(1, evalfil, 'sf', 0)
      call version(1)
      write (1, 1) potnam, label
 1    format (' ** POTNAM: ', (a), /, ' ** LABEL: ', (a))
      write (6, 5) r1, r2, spac, del, vmax, ndim, c, alph
      write (9, 5) r1, r2, spac, del, vmax, ndim, c, alph
      write (1, 5) r1, r2, spac, del, vmax, ndim, c, alph
 5    format (/,' ** BOUND STATE CALCULATION:  R1 =',
     :     f5.2,';  R2 =', f5.2,/,
     :     6x,'SPAC-TRY =', f6.3, ', SPAC-EXACT =', f9.6,
     :     ', NO. GAUSSIANS =', i3, ', SIZE OF TOTAL BASIS =', i5, /,
     :     6x, 'C =', f5.2, ', ALPH =', g12.5)
c
c     Calculate explicitly C, S, T matrices for one coupled channel
c     Refer to Eq. (7) -- (11) of Hamilton and Light
      allocate(cchn(vmax, vmax))
      allocate(schn(vmax, vmax))
      allocate(tchn(vmax, vmax))
      do i = 1, vmax
         do j = i, vmax
            cchn(i, j) = 0.5d0 * ((i - j) * c) ** 2
            cchn(j, i) = cchn(i, j)
         end do
      end do
      schn = exp(-cchn)
      tchn = schn * 0.5d0 * alph * (1 - 2 * cchn) / rmu
c     User cchn to store schn for diagonalization
      cchn = schn
      allocate(w(vmax))
      lwork = 8 * vmax
      allocate(work(lwork))
      allocate(iwork(5 * vmax))
      call dsyevx('N', 'I', 'L', vmax, cchn, vmax, 0d0, 0d0, 1, 1, 0d0,
     $     m, w, 0d0, 1, work, lwork, iwork, ifail, info)
      deallocate(iwork)
      deallocate(work)
      write (1, 19) w(1)
      write (6, 19) w(1)
      write (9, 19) w(1)
 19   format (' ** MINIMUM EIGENVALUE OF S MATRIX =', 1pe11.4)
      if (w(1) .lt. eigmin) then
         write (6, 25) eigmin
         write (9, 25) eigmin
 25      format ('    THIS IS .LT.',1pe8.1,
     :        '; INCREASE PARAMETER C OR DECREASE EIGMIN')
         deallocate(cchn)
         deallocate(schn)
         deallocate(tchn)
         deallocate(w)
         return
      endif
      deallocate(w)
      deallocate(cchn)
c     Allocate for the full (dim nch * vmax) S, H matrix.
c     Each vmax x vmax block of this matrix includes all n channels
c     and one distributed gaussian function.
      allocate(s(ndim, ndim))
      s = 0d0
      allocate(h(ndim, ndim))
      h = 0d0
      allocate(wr(nmax, nmax))
c     Choose quadrature method: 100 G-H quadrature; 200 Simpson's rule
c     Only blocks (row i, column j) with j=i or j=i+1 are integrated
c     out; other blocks will be filled staring from 300.
c     This choice minimize (i - j) ** 2, which improve the numerical
c     stability of "ratio" defined below.
      if (hsimp .gt. 0d0) goto 200
c
c     Gauss-Hermite quadrature integration
 100  ngh = int(dabs(hsimp))
      if (ngh .le. 1) ngh = 1
      write (1, 110) ngh
      write (6, 110) ngh
      write (9, 110) ngh
 110  format (' ** ', i3, ' POINT(S) GAUSS-HERMITE QUADRATURE')
      xfact = dsqrt(0.5d0 / alph)
      allocate(xgh(ngh))
      allocate(wgh(ngh))
      call gauher(ngh, xgh, wgh)
      wgh = wgh / dsqrt(pi)
      allocate(rgh(ngh))
      allocate(wgh1(ngh))
      do i = 1, vmax
         irow1 = (i - 1) * nch + 1
         irow2 = i * nch
         do j = i, i + 1
            if (j .gt. vmax) cycle
            wgh1 = wgh * dexp(-2d0 * alph
     $           * (dble(i - j) * del / 2d0) ** 2)
            icol1 = (j - 1) * nch + 1
            icol2 = j * nch
            xij = r1 + dble(i + j - 2) / 2d0 * del
            rgh = xfact * xgh + xij
            do ir = 1, ngh
               rnow = rgh(ir)
               call vmat_bound(wr, rnow, nch, nmax)
               wr = wr * wgh1(ir)
               h(irow1:irow2, icol1:icol2) = wr(:nch, :nch)
     $              + h(irow1:irow2, icol1:icol2)
            end do
         end do
      end do
      deallocate(xgh)
      deallocate(wgh)
      deallocate(wgh1)
      deallocate(rgh)
      goto 300
c
c     Simpson's 3/8 rule integration
 200  nr = (r2 - r1 + 2d0 * delr) / hsimp
      if (nr .lt. 2) nr = 2
      hexac = (r2 - r1 + 2d0 * delr) / nr
      write (1, 210) delr, hsimp, nr + 1, hexac
      write (6, 210) delr, hsimp, nr + 1, hexac
      write (9, 210) delr, hsimp, nr + 1, hexac
 210  format (' ** SIMPSON''S 3/8 RULE INTEGRATION OF POTENTIAL',
     $     /, 6x, 'DEL-R =', f5.2, ', H-START =', f6.3,
     $     ', N-POINTS =', i4, ', H-EXACT =', f9.6)
      allocate(g(vmax))
      do ir = 1, nr + 1
         rnow = (ir - 1) * hexac + r1 - delr
         call vmat_bound(wr, rnow, nch, nmax)
         if (ir .eq. 1 .or. ir .eq. nr) then
            wt = hexac * 0.375d0
         else if (mod(ir, 3) .eq. 1) then
            wt = hexac * 0.75d0
         else
            wt = hexac * 1.125d0
         end if
         wr = wr * wt
         do i = 1, vmax
            g(i) = dexp(-alph * (rnow - r1 - (i - 1) * del) ** 2)
         end do
         g = g * dsqrt(dsqrt(2 * alph / pi))
         do i = 1, vmax
            do j = i, i + 1
               if (j .gt. vmax) cycle
               irow1 = (i - 1) * nch + 1
               irow2 = i * nch
               icol1 = (j - 1) * nch + 1
               icol2 = j * nch
               h(irow1:irow2, icol1:icol2) = wr(:nch, :nch) * g(i)
     $              * g(j) + h(irow1:irow2, icol1:icol2)
            end do
         end do
      end do
      deallocate(g)
c
 300  deallocate(wr)
c     Expand H matrix using symmetry properties
      do i = 1, vmax
         irow = (i - 1) * nch + 1
         do j = 1, vmax
            if ((i .eq. j) .or. (j - i .eq. 1)) cycle
            icol = (j - 1) * nch + 1
            i1 = (i + j) / 2
            j1 = i1 + mod(i + j, 2)
            ratio = dexp(-0.5d0 * c ** 2
     $           * (dble(i - j) ** 2 - dble(i1 - j1) ** 2))
            irow1 = (i1 - 1) * nch + 1
            icol1 = (j1 - 1) * nch + 1
            h(irow:irow+nch-1, icol:icol+nch-1) = ratio
     $           * h(irow1:irow1+nch-1, icol1:icol1+nch-1)
         end do
      end do
c     Expand S matrix and add T matrix to H matrix
      do i = 1, vmax
         irow = (i - 1) * nch + 1
         do j = 1, vmax
            icol = (j - 1) * nch + 1
            do k = 0, nch - 1
               s(irow + k, icol + k) = schn(i, j)
               h(irow + k, icol + k) = h(irow + k, icol + k)
     $              + tchn(i, j)
            end do
         end do
      end do
      deallocate(schn)
      deallocate(tchn)
c
      allocate(w(ndim))
      if (wavefl) then
         if (tolai .ge. 0d0) then
            write (6, 401)
 401        format (' ** CALL DSYGVD TO CALCULATE ALL EIGENVALUES ',
     $           'AND EIGENVECTORS.')
            allocate(work(1))
            call dsygvd(1, 'V', 'L', ndim, h, ndim, s, ndim, w,
     $           work, -1, iworks, -1, info)
            lwork = work(1)
            liwork = iworks(1)
            write (1, *) 'REQUIRED ARRAY LENGTH', lwork, liwork
            deallocate(work)
            allocate(work(lwork))
            allocate(iwork(liwork))
            call dsygvd(1, 'V', 'L', ndim, h, ndim, s, ndim, w,
     $           work, lwork, iwork, liwork, info)
            m = ndim
         else
            il = 1
            iu = int(abs(tolai))
            if (iu .lt. 1) iu = 1
            write (6, 411) iu
 411        format (' ** CALL DSYGVX TO CALCULATE ', i3,
     $           ' LOWEST EIGENVALUES AND CORRESPONDING EIGENVECTORS.')
            allocate(work(1))
            call dsygvx(1, 'V', 'I', 'L', ndim, h, ndim, s, ndim,
     $           0d0, 0d0, il, iu, 0d0, m, w, z, ndim,
     $           work, -1, iworks, ifail, info)
            lwork = work(1)
            liwork = 5 * ndim
            write (1, *) 'REQUIRED ARRAY LENGTH', lwork, liwork
            deallocate(work)
            allocate(work(lwork))
            allocate(iwork(liwork))
            allocate(z(ndim, iu))
            call dsygvx(1, 'V', 'I', 'L', ndim, h, ndim, s, ndim,
     $           0d0, 0d0, il, iu, 0d0, m, w, z, ndim,
     $           work, lwork, iwork, ifail, info)
            h(:, :m) = z(:, :m)
            deallocate(z)
         end if
      else
         if (tolai .ge. 0d0) then
            write (6, 421)
 421        format (' ** CALL DSYGVD TO CALCULATE ALL EIGENVALUES.')
            allocate(work(1))
            call dsygvd(1, 'N', 'L', ndim, h, ndim, s, ndim, w,
     $           work, -1, iworks, -1, info)
            lwork = work(1)
            liwork = iworks(1)
            write (1, *) 'REQUIRED ARRAY LENGTH', lwork, liwork
            deallocate(work)
            allocate(work(lwork))
            allocate(iwork(liwork))
            call dsygvd(1, 'N', 'L', ndim, h, ndim, s, ndim, w,
     $           work, lwork, iwork, liwork, info)
            m = ndim
         else
            il = 1
            iu = int(abs(tolai))
            write (6, 431) iu
 431        format (' ** CALL DSYGVX TO CALCULATE ', i3,
     $           ' LOWEST EIGENVALUES.')
            if (iu .lt. 1) iu = 1
            allocate(work(1))
            call dsygvx(1, 'N', 'I', 'L', ndim, h, ndim, s, ndim,
     $           0d0, 0d0, il, iu, 0d0, m, w, 0d0, ione,
     $           work, -1, iworks, ifail, info)
            lwork = work(1)
            liwork = 5 * ndim
            write (1, *) 'REQUIRED ARRAY LENGTH', lwork, liwork
            deallocate(work)
            allocate(work(lwork))
            allocate(iwork(liwork))
            call dsygvx(1, 'N', 'I', 'L', ndim, h, ndim, s, ndim,
     $           0d0, 0d0, il, iu, 0d0, m, w, 0d0, ione,
     $           work, lwork, iwork, ifail, info)
        end if
      end if
      deallocate(iwork)
      deallocate(work)
      if (info .ne. 0) then
         write (0, *) ' *** bound: dsygvd or dsygvx fails', info
         stop
      end if
      deallocate(s)
c
      write (6, 55) evalfil, label
 55   format(' ** NEGATIVE AND LOWEST POSITIVE EIGENVALUES ',
     $     'SAVED IN FILE ', a, /, 6x, 'LABEL:', a, /,
     $     ' ** NEGATIVE EIGENVALUES ARE (CM-1):')
      write (1, 57)
 57   format(' ** EIGENVALUES (CM-1)')
      nbound = 0
      do i = 1, m
         write (1, 65) w(i) * econv
         if (w(i) .le. epsilon(0d0)) then
            nbound = nbound + 1
            write (9, 65) w(i) * econv
            write (6, 65) w(i) * econv
         end if
         if (i - nbound .gt. maxpos) exit
      end do
 65   format(2x, f11.3)
c
      if (wavefl) then
         write (1, 75) nch, vmax, min(m, maxvec)
 75      format (' NCH =', i3, '; N-GAUSSIANS = ', i3, '; NSTATE =', i3)
         write (1, 80) alph, (r1 + (i - 1) * del, i = 1, vmax)
 80      format (/,' GAUSSIAN FUNCTIONS - ALPH =',1pg12.5,
     $        '; CENTERED AT:', /, 0p40f8.4)
         write (1, 85)
 85      format (' WAVEFUNCTION EXPANSION COEFFICIENTS:',
     $        ' ROWS ARE CHANNELS, COLUMNS ARE GAUSSIANS')
         if (nbound .le. 0) then
            write (6, 90)
            write (9, 90)
            write (1, 90)
 90         format (' ** NO BOUND STATES (E < 0)')
            goto 500
         end if
         do j = 1, min(m, maxvec)
            write (1, 95) j
 95         format (' STATE =', i2)
            do i = 1, nch
               write (1, 97) h(i::nch, j)
 97            format (40(1pg13.5))
            end do
         end do
      end if
c
c     Quick analysis of wave vectors
      if (wavefl) then
         call wpr_bound(1, ndim, nch, vmax, h,
     $        min(m, maxvec), min(maxchn, nch), r1, r2, alph)
      end if
c
 500  close(1)
      deallocate(h)
      deallocate(w)
      return
      end subroutine bound
c     ------------------------------------------------------------------
      subroutine wpr_bound(iunit, n, nch, ng, z, m, mch, r1, r2, alph)
      implicit none
      integer, intent(in) :: iunit, n, nch, ng, m, mch
      real(8), dimension(n, m), intent(in) :: z
      real(8), intent(in) :: r1, r2, alph
      real(8), dimension(:), allocatable :: chnwt, r
      real(8), dimension(:, :), allocatable :: g, w
      real(8) :: chwtmx, del
      integer :: i, j, ich, mchnow, ichnow, ir, ig, irrow, igrow
c
c     Compute distributed Gaussian
      allocate(g(ng, ng))
      del = (r2 - r1) / dble(ng)
      do ir = 1, ng
         do ig = 1, ng
         g(ir, ig) = sqrt(sqrt(2d0 * alph / dacos(-1d0)))
     $        * exp(-alph * (del * (ir - ig)) ** 2)
         end do
      end do
c     Convert coefficients of distributed Gaussian to actuall wave
c     function
      allocate(w(n, m))
      w = 0d0
      do ir = 1, ng
         irrow = (ir - 1) * nch + 1
         do ig = 1, ng
            igrow = (ig - 1) * nch + 1
            w(irrow:irrow+nch-1, :) = w(irrow:irrow+nch-1, :)
     $           + g(ir, ig) * z(igrow:igrow+nch-1, :)
         end do
      end do
      deallocate(g)
c
      allocate(chnwt(nch))
      do j = 1, m
         write (iunit, *)
         write (iunit, 10) j
 10      format (' ** ANALYZING WAVE VECTOR NO. ', i3)
         chnwt = 0d0
         do i = 1, nch
            chnwt(i) = sum(w(i::nch, j) ** 2)
         end do
         do mchnow = 1, mch
            chwtmx = 0d0
            do ich = 1, nch
               if (chnwt(ich) .gt. chwtmx) then
                  chwtmx = chnwt(ich)
                  ichnow = ich
               end if
            end do
            write (iunit, 11) ichnow, chwtmx
 11         format ('   * CHANNEL NO. ', i4, ' WEIGHT ', f12.6,
     $           ' WAVE FUNCTION VS R:')
            chnwt(ichnow) = 0d0
            write (iunit, 12) w(ichnow::nch, j)
 12         format (10e12.4)
         end do
      end do
      deallocate(chnwt)
      deallocate(w)
      return
      end subroutine wpr_bound
c     ------------------------------------------------------------------
      subroutine vmat_bound(wr, r, nch, nmax)
      implicit none
      integer, intent(in) :: nch, nmax
      real(8), intent(in) :: r
      real(8), dimension(nmax, nmax), intent(out) :: wr
      common /coered/ ered, rmu
      real(8) :: ered, rmu
      real(8) :: xirmu
      integer :: i, j, ig, jg, ilast, ntop
c
      xirmu = 0.5d0 / rmu
      call potmat(wr, r, nch, nmax)
      do i = 1, nch
         do j = i + 1, nch
            wr(i, j) = wr(j, i)
         end do
      end do
      wr = wr * xirmu
c     Add back collision energy substracted in subroutine potmat
      do i = 1, nch
         wr(i, i) = wr(i, i) + ered
      end do
      return
      end subroutine vmat_bound
c     ------------------------------------------------------------------
      subroutine gauher(n, x, w)
c     Generate nodes and weights for the Gauss-Hermite quadrature
c     From Numerical Recipes 3rd Ed, p. 185
c     Translated from C++ to fortran by Q. Ma
      implicit none
      integer, intent(in) :: n
      real(8), dimension(0:n-1), intent(out) :: x, w
      real(8), parameter :: eps=1d-14, maxit=10,
     $     pim4=dsqrt(dsqrt(1.0/dacos(-1d0)))
      integer :: i, its, j, m
      real(8) :: p1, p2, p3, pp, z, z1
c
      m = (n + 1) / 2
      do i = 0, m - 1
         select case (i)
            case (0)
               z = dsqrt(dble(2 * n + 1)) -
     $              1.85575d0 * dble(2 * n + 1) ** (-0.16667d0)
            case (1)
               z = z - 1.14d0 * dble(n) ** 0.426d0 / z
            case (2)
               z = 1.86d0 * z - 0.86d0 * x(0)
            case (3)
               z = 1.91d0 * z - 0.91d0 * x(1)
            case default
               z = 2.0d0 * z - x(i - 2)
         end select
         do its = 1, maxit
            p1 = pim4
            p2 = 0
            do j = 0, n - 1
               p3 = p2
               p2 = p1
               p1 = z * dsqrt(2d0 / dble(j + 1)) * p2
     $              - dsqrt(dble(j) / dble(j + 1)) * p3
            end do
            pp = dsqrt(dble(2 * n)) * p2
            z1 = z
            z = z1 - p1 / pp
            if (dabs(z - z1) .le. eps) exit
         end do
         x(i) = z
         x(n - 1 - i) = -z
         w(i) = 2.0 / pp ** 2
         w(n - 1 - i) = w(i)
      end do
      return
      end subroutine gauher
