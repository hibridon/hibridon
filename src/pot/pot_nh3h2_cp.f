c     pot_nh3h2_cp.f
c     authors: Qianli Ma
c
c     Pot routine for a collision between a symmetric top and a linear
c     molecule.
c
c     PES to be read from data file specified in the input file.
c
c     This basis routine should have ibasty = 100, so hibridon will
c     treat it as a mol-mol problem with j12.
c
c     DUMMY SUBROUTINES FOR ALL USER-DEFINED BASIS/POT.
      include "common/ground"
c
c     ------------------------------------------------------------------
c     THE FOLLOWING SOUBROUTINE WILL BE THE MAIN FUNCTION FOR MAKEPOT.
c     NOTE THAT VVL IS IN HARTREES.
      subroutine driver
      implicit none
c
      include 'pot_nh3h2_qma_common.f'
      character*40 filenm
      double precision r, vv0
      integer i
c
      print *, 'Please input the name of potdata file:'
      read (5, *, end=99) filenm
      call loapot(10, filenm)
 10   print *, 'R (bohr), Ctrl+D to exit:'
      read (5, *, end=99) r
      call pot(vv0, r)
      write (6, 20) (vvl(i) * econv, i=1, nv)
 20   format (10(1pe16.8))
      goto 10
 99   end
c     ------------------------------------------------------------------
c     DATA FILE, IF REQUIRED, CAN BE LOADED WITH THIS SOUBROUTINE.
      subroutine loapot(iunit, filnam)
      implicit none
c
      include 'pot_nh3h2_qma_common.f'
      include 'common/parpot'
c
      character*(*) filnam
      integer iunit, ir, iv
      character*255 datfl
c
c     WHEN HIBRIDON LOADS, A STRING CONTAINING ONLY ONE SPACE WILL BE
c     PASSED TO THIS SUBROUTINE.  MAKE SURE NOTHING WILL BE DONE.
      if (filnam .eq. ' ') return
c
c     LOAD DATA TO A COMMON BLOCK ONLY USED IN THIS FILE, AND DO PRE-
c     PROCESSING IF NECESSARY.
      call datfln(filnam, datfl)
      open (unit=iunit, file=datfl, status='old')
      read (iunit, *) nr
      if (nr .gt. MAX_NR) call raise('too many R in data file.')
      read (iunit, *) (rr(ir), ir=1, nr)
      read (iunit, *) nv
      if (nv .gt. MAX_NV) call raise('too many vlm terms in data file.')
      do iv = 1, nv
         read (iunit, *) lb1(iv), mu1(iv), lb2(iv), mu2(iv)
         read (iunit, *) (v_pot(ir, iv), ir=1, nr)
      end do
      close (unit=iunit)
      do iv = 1, nv
         call dscal(nr, 1d0 / econv, v_pot(1, iv), 1)
      end do
c
c     spline parameters are prepared here; pot routine only call the
c     accompanying function seval.
      do iv = 1, nv
         call spline(nr, rr, v_pot(1, iv), spl_b(1, iv), spl_c(1, iv),
     $        spl_d(1, iv))
      end do
      return
      end
c     ------------------------------------------------------------------
c     A POT ROUTINE RETURNS VVL ARRAY (COVVL BLOCK) FOR A GIVEN R.
c     ALWAYS SET VV0 = 0 TO AVOID CONFUSION.  THE ISOTROPIC TERM CAN BE
c     EASILY TREAT AS ONE TERM IN THE POTENTIAL EXPANSION.  NOTE THAT
c     VVL SHOULD BE IN HARTREES.
      subroutine pot(vv0, r)
      implicit none
c
      include 'pot_nh3h2_qma_common.f'
      double precision vv0, r
      double precision seval
      integer iv
c
      vv0 = 0d0
      do iv = 1, nv
         vvl(iv) = seval(nr, r, rr, v_pot(1, iv),
     $        spl_b(1, iv), spl_c(1, iv), spl_d(1, iv))
      end do
      return
      end
c     ------------------------------------------------------------------
c     THIS SUBROUTINE GOVERNS THE INPUT/OUTPUT OF THE BASIS ROUTINE.
c     ONLY IREAD IS USED: RETURN DIRECTLY IF ZERO.
      subroutine syusr(irpot, readpt, iread)
      implicit none
c
      include 'pot_nh3h2_qma_common.f'
      integer irpot, iread
      logical readpt
      character*(*) fname
c     NUMBER OF BASIS-SPECIFIC VARIABLES, MODIFY ACCORDINGLY.
      integer icod, ircod, lencod
      parameter (icod=5, ircod=5, lencod=icod+ircod+3)
      common /cosys/ scod(lencod)
      character*8 scod
c     INTEGER VARIABLES.  LEAVE THE FIRST TWO AS IT IS.
      common /cosysi/ nscode, isicod, nterm, ipotsy, iop, ninv, jmax,
     $     ipotsy2, j2max, j2min
      integer nscode, isicod, nterm, ipotsy, iop, ninv, jmax, ipotsy2,
     $     j2max, j2min
c     REAL VARIABLES.  LEAVE THE FIRST TWO AS IT IS.
      common /cosysr/ isrcod, junkr, brot, crot, delta, emax, drot
      integer isrcod, junkr
      double precision brot, crot, delta, emax, drot
      character*40 potfil
      save potfil
c     DEFINE THE NAMES HERE
      scod(1)='NTERM'
      scod(2)='IPOTSY'
      scod(3)='IOP'
      scod(4)='NINV'
      scod(5)='JMAX'
      scod(6)='IPOTSY2'
      scod(7)='J2MAX'
      scod(8)='J2MIN'
      scod(9)='BROT'
      scod(10)='CROT'
      scod(11)='DELTA'
      scod(12)='EMAX'
      scod(13)='DROT'
      nscode = lencod
      isicod = icod
      isrcod = ircod
c     KEEP THE FOLLOWING LINE
      if (iread .eq. 0) return
c     READ THE LAST FEW LINES OF THE INPUT FILE
      read (8, *, err=80) ipotsy, iop, ninv, ipotsy2
      read (8, *, err=80) jmax
      read (8, *, err=80) j2min, j2max
      read (8, *, err=80) brot, crot, delta, emax
      read (8, *, err=80) drot
      read (8, *, err=80) potfil
      call loapot(10, potfil)
      close (8)
      return
 80   call raise('error read from input file.')
      return
c     ------------------------------------------------------------------
      entry ptrusr(fname, readpt)
      return
c     ------------------------------------------------------------------
      entry savusr(readpt)
c     WRITE THE LAST FEW LINES OF THE INPUT FILE.
      write (8, 220) ipotsy, iop, ninv, ipotsy2
 220  format (4i4, 14x,'   ipotsy, iop, ninv, ipotsy2')
      write (8, 230) jmax
 230  format (i4, 26x, '   jmax')
      write (8,231) j2min, j2max
 231  format (2i4, 22x,'   j2min, j2max')
      write (8, 250) brot, crot, delta, emax
 250  format (3f8.4, f8.2, ' brot, crot, delta, emax' )
      write (8, 251) drot
 251  format (f12.6, 18x,'   drot')
      write (8, *) potfil
      return
      end
c     ------------------------------------------------------------------
c     THE BASIS ROUTINE.  JHOLD, EHOLD, ISHOLD, NLEVEL, NLEVOP ARE TO
c     BE RETURNED FOR LEVEL INFORMATION; J, L, IS, IEPS, EINT, CENT, N,
c     NTOP, AND POSSIBLY J12 ARE TO BE RETURNED FOR CHANNEL INFORMATION;
c     ALSO CALCULATE NON-ZERO COUPLING MATRIX ELEMENTS AND RETURN IN V2,
c     IV2, LAMNUM.
      subroutine bausr(j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     $     k, ieps, jtemp, ktemp, rcut, jtot, flaghf, flagsu, csflag,
     $     clist, bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
      implicit none
c
      include 'pot_nh3h2_qma_common.f'
      double precision rcut
      integer jtot, nu, numin, jlpar, nmax
      logical flaghf, flagsu, csflag, clist, bastst, ihomo
      integer j(*), k(*), l(*), is(*), ieps(*), jhold(*), ishold(*)
      double precision ehold(*)
      integer nlevel, nlevop, n, ntop, jtemp(*), ktemp(*)
c
      common /cosysi/ nscode, isicod, nterm, ipotsy, iop, ninv, jmax,
     $     ipotsy2, j2max, j2min
      integer nscode, isicod, nterm, ipotsy, iop, ninv, jmax, ipotsy2,
     $     j2max, j2min
      common /cosysr/ isrcod, junkr, brot, crot, delta, emax, drot
      integer isrcod, junkr
      double precision brot, crot, delta, emax, drot
      common /conlam/ nlam, nlammx, lamnum
      integer nlam, nlammx, lamnum(MAX_NV)
      common /coered/ ered, rmu
      double precision ered, rmu
      common /cov2/ nv2max, junkv, v2
      integer nv2max, junkv
      double precision v2(1)
      common /coiv2/ iv2
      integer iv2(1)
      common /cocent/ cent
      double precision cent(1)
      common /coeint/ eint
      double precision eint(1)
      common /coj12/ j12
      integer j12(1)
      common /coamat/ ietmp
      integer ietmp(1)
      common /coipar/ junkip, iprint
      integer junkip(9), iprint
c
      integer nlist, ki, ji, numeps, iep, iepsil, isym, ji2, i1, i2,
     $     jsave, ksave, isave, njk, ipar, j12max, j12min, ji12,
     $     lmax, lmin, li, lpar, j2par, i, lamsum, ilam, iv, icol,
     $     irow, jc, jr, j2c, j2r, kc, kr, j12c, j12r, lc, lr, inum
      double precision roteng, esave, vee, vstp1sg
c
      if (flaghf) call raise('FLAGHF = .TRUE. FOR SINGLET SYSTEM')
      if (flagsu) call raise('FLAGSU = .TRUE. FOR MOL-MOL COLLISION')
      if (csflag) call raise('CS calculation not implemented')
      nlist = 0
      do ki = 0, jmax
         do ji = ki, jmax
            numeps = 2
            if (ki .eq. 0) numeps = 1
            do iep = 1, numeps
               if (iop .eq. -1 .and. mod(ki, ipotsy) .ne. 0) cycle
               if (iop .eq. 1 .and. mod(ki, ipotsy) .eq. 0) cycle
               roteng = brot * ji * (ji + 1) + (crot - brot) * ki ** 2
               if (roteng .gt. emax) cycle
               iepsil = 1 - 2 * (iep - 1)
               isym = -iepsil * (-1) ** ji
               if (isym .eq. -1) roteng = roteng + delta
               if ((ninv .ne. 2) .and. (ninv .ne. isym)) cycle
               do ji2 = j2min, j2max, ipotsy2
                  nlist = nlist + 1
                  jtemp(nlist) = 10 * ji + ji2
                  ktemp(nlist) = ki
                  ietmp(nlist) = iepsil
                  ehold(nlist) = (roteng + drot * ji2 * (ji2 + 1))
     $                 / econv
               end do
            end do
         end do
      end do
      do i1 = 1, nlist - 1
         esave = ehold(i1)
         do i2 = i1 + 1, nlist
            if (ehold(i2) .lt. esave) then
               esave = ehold(i2)
               ehold(i2) = ehold(i1)
               ehold(i1) = esave
               jsave = jtemp(i2)
               jtemp(i2) = jtemp(i1)
               jtemp(i1) = jsave
               ksave = ktemp(i2)
               ktemp(i2) = ktemp(i1)
               ktemp(i1) = ksave
               isave = ietmp(i2)
               ietmp(i2) = ietmp(i1)
               ietmp(i1) = isave
            end if
         end do
      end do
      do i = 1, nlist
         jhold(i) = jtemp(i)
         ishold(i) = -ietmp(i) * (-1) ** (jhold(i) / 10) * ktemp(i)
      end do
      if (bastst) call pr_lev_nh3h2(nlist, jtemp, ktemp, ietmp, ehold)
c
*     now set up channel and level list for scattering calculation
      n = 0
      do nlevel = 1, nlist
         ki = ktemp(nlevel)
         ji = jtemp(nlevel) / 10
         ji2 = mod(jtemp(nlevel), 10)
         iepsil = ietmp(nlevel)
         ipar = -iepsil * (-1) ** (ji + ki)
         do ji12 = iabs(ji - ji2), ji + ji2
            do li = iabs(jtot - ji12), jtot + ji12
               lpar = (-1) ** (li - jtot + ji2)
               if (ipar * lpar .ne. jlpar) cycle
               n = n + 1
               if (n .gt. nmax) call raise('too many channels.')
               j(n) = jtemp(nlevel)
               is(n) = -iepsil * (-1) ** ji * ki
               ieps(n) = iepsil
               j12(n) = ji12
               eint(n) = ehold(nlevel)
               l(n) = li
               cent(n) = li * (li + 1)
               if (bastst .and. iprint .ge. 2)
     $              write (6, 280) n, j(n), is(n), ieps(n), j12(n),
     $              l(n), eint(n) * econv
 280           format (6i6, f12.3)
            end do
         end do
      end do
      nlevel = nlist
*  also determine number of levels which are open
      nlevop = 0
      do i = 1, nlist
         if (ehold(i) .le. ered) nlevop = nlevop + 1
      end do
*  return if no channels
      if (n .eq. 0) return
      if (nu .eq. numin) then
         ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
         if (mod(ntop, 2) .eq. 0 .and. ntop .lt. nmax)
     :        ntop = ntop + 1
      else
         if (n .gt. ntop) call raise('nch greater than ntop.')
      end if
c
c     Calculate coupling matrix elements
      nlam = nv
      i = 0
      lamsum = 0
      ilam = 0
      do iv = 1, nv
         ilam = ilam + 1
         inum = 0
         do icol = 1, n
            do irow = icol, n
               jc = j(icol) / 10
               jr = j(irow) / 10
               j2c = mod(j(icol), 10)
               j2r = mod(j(irow), 10)
               kc = iabs(is(icol))
               kr = iabs(is(irow))
               j12c = j12(icol)
               j12r = j12(irow)
               lc = l(icol)
               lr = l(irow)
               vee = vstp1sg(jr, lr, jc, lc, j2r, j2c, j12r, j12c,
     $              jtot, kr, kc, lb1(iv), mu1(iv), lb2(iv), mu2(iv),
     $              ieps(irow), ieps(icol))
               if (vee .eq. 0d0) cycle
               if (i .eq. nv2max) call raise(
     $              'too many non-zero V-matrix terms.')
               i = i + 1
               inum = inum + 1
               v2(i) = vee
               iv2(i) = ntop * (icol - 1) + irow
               if (bastst .and. iprint .ge. 2)
     $              write (6, 345) ilam, lb1(iv), icol, irow, i,
     $              iv2(i), vee
 345           format (i4, 2i7, 2i6, i6, f17.10)
            end do
         end do
         lamnum(ilam) = inum
         lamsum = lamsum + inum
         if (bastst) write (6, 347) ilam, lb1(iv), mu1(iv), lb2(iv),
     $         mu2(iv), lamnum(ilam)
 347     format ('ILAM=', i3, ' LAM=', i3, ' MU=', i3,
     $        ' LAM2=', i3, ' MU2=', i3, ' LAMNUM(ILAM) = ', i6)
      end do
      if (bastst .and. iprint .ge. 2) then
         write (6, 350) i
 350     format ('number of non-zero terms: ', 6i)
         write (6, 351) (v2(i1), i1=1, i)
 351     format (10(f7.2), 1x)
         write (6, 352) (iv2(i1), i1=1, i)
 352     format (10(i7), 1x)
      end if
      return
      end
c     ------------------------------------------------------------------
      subroutine pr_lev_nh3h2(n, js, ks, iepss, es)
      implicit none
c
      integer n, js(*), ks(*), iepss(*)
      double precision es(*)
      integer i, j1, j2, isym
      double precision ecm, econv
      parameter (econv=219474.6315343234)
      write (6, 125)
 125  format (/, 10x,
     $     'SORTED LEVEL LIST', /, '   N   J   K  EPS INV J2   ',
     $     'EINT(CM-1)')
      do i = 1, n
         j2 = mod(js(i), 10)
         j1 = js(i) / 10
         isym = -iepss(i) * (-1) ** j1
         ecm = es(i) * econv
         write (6, 126) i, j1, ks(i), iepss(i), isym, j2, ecm
 126     format (6i4, f10.3)
      end do
      return
      end
c     ------------------------------------------------------------------
      subroutine raise(mesg)
      implicit none
      character*(*) mesg
      write (0, *) 'hibridon: error: ', mesg
      stop
      end
C     ------------------------------------------------------------------
      real(8) function vstp1sg(j1p, lp, j1, l, j2p, j2, j12p, j12,
     $     jtot, kp, k, lam1, mu1, lam2, lam, epsp, eps)
      implicit none
      integer :: j1p, kp, epsp, j2p, j12p, lp, j1, k, eps, j2, j12, l,
     $     jtot, lam1, mu1, lam2, lam
      integer :: iphase, tj1p, tj1, tlam1, tkp, tmu1, tk
      real(8) :: pref, threej, sixj, ninej, tf3jm0, tf3j, f6j, f9j
      real(8), parameter :: machep = epsilon(0d0)
C     
      vstp1sg = 0d0
      iphase = epsp * eps * (-1) ** (j1p + j1 + lam2 + lam + mu1)
      if (iphase .eq. -1) return
      threej = tf3jm0(2 * j2p, 2 * lam2, 2 * j2)
      if (dabs(threej) .lt. machep) return
      threej = threej * tf3jm0(2 * lp, 2 * lam, 2 * l)
      if (dabs(threej) .lt. machep) return
C     
      tj1p = 2 * j1p
      tlam1 = 2 * lam1
      tj1 = 2 * j1
      tkp = 2 * kp
      tmu1 = 2 * mu1
      tk = 2 * k
      threej = threej * (tf3j(tj1p, tlam1, tj1, -tkp, tmu1, tk)
     $     + epsp * eps * tf3j(tj1p, tlam1, tj1, tkp, tmu1, -tk)
     $     + epsp * tf3j(tj1p, tlam1, tj1, tkp, tmu1, tk)
     $     + eps * tf3j(tj1p, tlam1, tj1, -tkp, tmu1, -tk))
      if (dabs(threej) .lt. machep) return
C     
      sixj = f6j(j12, l, jtot, lp, j12p, lam)
      if (dabs(sixj) .lt. machep) return
      ninej = f9j(j1, j2, j12, j1p, j2p, j12p, lam1, lam2, lam)
      if (dabs(ninej) .lt. machep) return
C     
      iphase = (-1) ** (jtot + lam1 - lam2 + j1 - j2 + j12p - l - lp + k
     $     + mu1)
      if (mu1 .eq. 0) then
         pref = 0.5d0
      else
         pref = 1d0
      end if
      if (kp .eq. 0) pref = pref * dsqrt(0.5d0)
      if (k .eq. 0) pref = pref * dsqrt(0.5d0)
      pref = pref * dsqrt(dble((tj1p + 1) * (2 * j2p + 1)
     $     * (2 * j12p + 1)
     $     * (2 * lp + 1) * (tj1 + 1) * (2 * j2 + 1) * (2 * j12 + 1)
     $     * (2 * l + 1) * (2 * lam + 1)))
C     
      vstp1sg = pref * dble(iphase) * threej * sixj * ninej
      return
      end function vstp1sg
C
      subroutine datfln(filenm, fullnm)
      character (len=*) :: filenm, fullnm
      fullnm = 'potdata/' // trim(filenm)
      return
      end subroutine datfln
