c     Basis subroutine for the collision of a symmetric top with
c     a singlet-Sigma molecule
c
c     system dependent parameters:
c       ipotsy - cylindrical symmetry of the potential.  presently
c          only ipotsy=3 supported (e.g. NH3, CD3)
c       iop  - bitwise flag to set rotational basis (see the .html help
c       file for details)
c          CH3:  ortho levels, iop=2; para levels, iop=4 or 8
c          CD3:  A1 levels, iop=1; A1 levels, iop=2; E levels, iop=4 or 8
c          NH3:  ortho levels, iop=3; para levels, iop=12
c          ND3:  A1 levels, 17; A2 levels, 34; E levels, iop=12 (added ny p.j.d.)
c       j1max - maximum rotational momentum in channel expansion of sym. top
c       e1max - maximum energy of state to be included in sym. top expansion
c       j2min, j2max - minimum and maximum rotational angular momentum of
c          the linear molecule
c       brot, crot - rotational constants of the symmetric top
c       drot - rotational constant of the linear molecule
c       delta - inversion splitting of the symmetric top (assumed to
c          be the same for all levels
c       see the .html help file for the angular expansion of the potential used
c
c     Author: Qianli Ma
c     written:  Jul-2013
c
c     extended to include the 3 nuclear symmetry permutations of ND3
c     (20-nov-2014, p. dagdigian)
c
c     current revision:  24-nov-2104 (p. dagdigian)
c     ------------------------------------------------------------------
c     This module contains (explictly) the number of terms and their
c     indices in the expansion of the PES.  Its contents should be
c     filled in the pot routine.
c     This module replaces lammin, lammax, mproj in hibridon.
      module mod_stp1sg
      implicit none
c
      type lm_type
      integer :: l1, l2, l, mu1
      end type lm_type
c
      type lev_type
      integer :: j1, k1, eps1, j2
      real(8) :: e
      end type lev_type
c
      type chn_type
      integer :: ilev, j12, l
      end type chn_type
c
      integer :: nv
      type(lm_type), dimension(:), allocatable :: lms
      type(lev_type), dimension(:), allocatable :: levs
      type(chn_type), dimension(:), allocatable :: chns
      end module mod_stp1sg
c     ------------------------------------------------------------------
      subroutine bastp1sg(jchn, lchn, ischn, jlev, elev, islev, nlev,
     $     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, bastst,
     $     ihomo, nu, numin, jlpar, twomol, nchn, nmax, ntop)
      use mod_stp1sg
      implicit none
c
c     The following arrays store the parameters of channels and levels.
c     Note that j12 is stored in a common block
      integer :: jchn(*), lchn(*), ischn(*), jlev(*), islev(*),
     $     j12chn(1)
      real(8) :: elev(*), echn(1), cchn(1)
      common /coj12/ j12chn
      common /cocent/ cchn
      common /coeint/ echn
c     The following parameters store the number of channels, opened
c     levels and levels
      integer :: nchn, nlevop, nlev
c     Various parameters used through out hibridon
      integer :: jtot, jlpar, nu, numin
      logical :: flaghf, flagsu, csflag, clist, bastst, ihomo, twomol
c     Limits of channels
      integer :: nmax, ntop
c
      real(8) :: rcut
c
      common /cosysi/ nscode, isicod, ipotsy, iop, ipotsy2, j1max,
     $     j2min, j2max
      integer :: nscode, isicod, ipotsy, iop, ipotsy2, j1max,
     $     j2min, j2max
      common /cosysr/ isrcod, junkr, brot, crot, delta, e1max, drot
      integer :: isrcod, junkr
      real(8) :: brot, crot, delta, e1max, drot
      common /conlam/ nlam, nlammx, lamnum
      integer :: nlam, nlammx, lamnum(1)
      common /coered/ ered, rmu
      real(8) :: ered, rmu
      common /cov2/ nv2max, junkv, v2
      integer :: nv2max, junkv
      real(8), dimension(1) :: v2
      common /coiv2/ iv2
      integer, dimension(1) :: iv2
      common /coipar/ junkip, iprint
      integer, dimension(9) :: junkip
      integer :: iprint
c
      integer :: i, ilev, iv, irow, icol, inum, i1, i2, lamsum
      real(8) :: vee, vstp1sg
      real(8), parameter :: machep=epsilon(0d0)
c
      if (flaghf)
     $     call raise_2pi1sg('FLAGHF = .TRUE. FOR SINGLET SYSTEM')
      if (ihomo) call raise_2pi1sg('IHOMO NOT USED IN THIS BASIS,'
     $     //' PLEASE SET IT FALSE')
      if (flagsu)
     $     call raise_2pi1sg('FLAGSU = .TRUE. FOR MOL-MOL COLLISION')
      if (csflag)
     $     call raise_2pi1sg('CS CALCULATION NOT IMPLEMENTED')
      if (nv .gt. nlammx)
     $     call raise_2pi1sg('NLAMMX TOO SMALL FOR THE POTENTIAL.')
c$$$      if (rcut .ge. 0d0 .and. .not. bastst)
c$$$     $     call raise_2pi1sg('RCUT NOT IMPLEMENTED, PLEASE SET IT '
c$$$     $     // 'NEGATIVE.')
      if (.not. twomol)
     $     call raise_2pi1sg('TWOMOL IS FALSE FOR MOL-MOL COLLISION.')
c
      call genlev_stp1sg(ipotsy, iop, ipotsy2, j1max, j2min, j2max,
     $     brot, crot, delta, drot, e1max)
      call sortlev_stp1sg()
      if (bastst) call prtlev_stp1sg()
      call genchn_stp1sg(jtot, jlpar)
      if (bastst .and. iprint .ge. 1) call prtchn_stp1sg()
c
c     Copy level and channel info to hibridon arrays
      nlev = size(levs)
      nchn = size(chns)
      ntop = max(nchn, nlev)
      if (nchn .gt. nmax)
     $     call raise_2pi1sg('TOO MANY CHANNELS.')
      if (nchn .eq. 0) return
      do ilev = 1, nlev
         jlev(ilev) = levs(ilev)%j1 * 10 + levs(ilev)%j2
         islev(ilev) = -levs(ilev)%eps1 * (-1) ** levs(ilev)%j1
     $        * levs(ilev)%k1
         elev(ilev) = levs(ilev)%e
      end do
      do i = 1, nchn
         jchn(i) = jlev(chns(i)%ilev)
         ischn(i) = islev(chns(i)%ilev)
         j12chn(i) = chns(i)%j12
         lchn(i) = chns(i)%l
         cchn(i) = dble((lchn(i) * (lchn(i) + 1)))
         echn(i) = elev(chns(i)%ilev)
      end do
c
c     Count number of asymtotically open levels
      nlevop = 0
      do ilev = 1, nlev
         if (elev(ilev) .le. ered) then
            nlevop = nlevop + 1
         else
            exit
         end if
      end do
c
c     Calculate coupling matrix elements
      nlam = nv
      i = 0
      lamsum = 0
      do iv = 1, nv
         inum = 0
         do icol = 1, nchn
            i1 = chns(icol)%ilev
            do irow = icol, nchn
               i2 = chns(irow)%ilev
               vee = vstp1sg(levs(i1)%j1, chns(icol)%l, levs(i2)%j1,
     $              chns(irow)%l, levs(i1)%j2, levs(i2)%j2,
     $              chns(icol)%j12, chns(irow)%j12, jtot,
     $              levs(i1)%k1, levs(i2)%k1, lms(iv)%l1, lms(iv)%mu1,
     $              lms(iv)%l2, lms(iv)%l, levs(i1)%eps1, levs(i2)%eps1)
               if (dabs(vee) .lt. machep) cycle
               if (i .eq. nv2max) call raise_2pi1sg(
     $              'TOO MANY NON-ZERO V-MATRIX TERMS.')
               i = i + 1
               inum = inum + 1
               v2(i) = vee
               iv2(i) = ntop * (icol - 1) + irow
               if (bastst .and. iprint .ge. 2) write (6, 346) iv,
     $              lms(iv)%l1, lms(iv)%l2, lms(iv)%l, lms(iv)%mu1,
     $              icol, irow, i, iv2(i), vee
 346           format(i4, 4i3, 2i4, 2i5, f17.10)
            end do
         end do
         lamnum(iv) = inum
         lamsum = lamsum + inum
         if (bastst) write (6, 347) iv, lms(iv)%l1, lms(iv)%l2,
     $        lms(iv)%l, lms(iv)%mu1, lamnum(iv)
 347     format (' ILAM=', i3, '  L1=', i3, '  L2=', i3,
     $        '  L=', i3, '  MU=', i3, '  LAMNUM(ILAM) = ', i6)
      end do
      if (bastst) write (6, 350) i
 350  format ('TOTAL NUMBER OF NON-ZERO V2 TERMS: ', i6)
      return
      end subroutine bastp1sg
c     ------------------------------------------------------------------
      subroutine genchn_stp1sg(jtot, jlpar)
      use mod_stp1sg
      implicit none
      integer, intent(in) :: jtot, jlpar
      integer :: nchn, ichn, ilev, j12, l
c     First pass (nchn = -1): Count the number of channels
c     Second pass (nchn is the number of channels): Save channel info
      nchn = -1
      do while (.true.)
         ichn = 0
         if (nchn .ne. -1) then
            if (allocated(chns)) deallocate(chns)
            allocate(chns(nchn))
         end if
         do ilev = 1, size(levs)
            do j12 = iabs(levs(ilev)%j1 - levs(ilev)%j2),
     $           levs(ilev)%j1 + levs(ilev)%j2
               do l = iabs(jtot - j12), jtot + j12
                  if (jlpar .ne. -levs(ilev)%eps1 *
     $                 (-1) ** (levs(ilev)%k1 + levs(ilev)%j1
     $                 + levs(ilev)%j2 + l - jtot)) cycle
                  ichn = ichn + 1
                  if (nchn .eq. -1) cycle
                  chns(ichn)%ilev = ilev
                  chns(ichn)%j12 = j12
                  chns(ichn)%l = l
               end do
            end do
         end do
         if (nchn .eq. -1) then
            nchn = ichn
         else
            exit
         end if
      end do
      return
      end subroutine genchn_stp1sg
c     ------------------------------------------------------------------
      subroutine sortlev_stp1sg()
      use mod_stp1sg
      implicit none
      type(lev_type) :: temp
      type(lev_type), dimension(:), allocatable :: levs1
      integer :: i, j
      real(8) :: esave
      common /coconv/ econv
      real(8) :: econv
      do i = 1, size(levs) - 1
         esave = levs(i)%e
         do j = i + 1, size(levs)
            if (levs(j)%e .lt. esave) then
               esave = levs(j)%e
               temp = levs(j)
               levs(j) = levs(i)
               levs(i) = temp
            end if
         end do
      end do
      return
      end subroutine sortlev_stp1sg
c     ------------------------------------------------------------------
      subroutine genlev_stp1sg(ipotsy, iop, ipotsy2, j1max, j2min,
     $     j2max, brot, crot, delta, drot, e1max)
      use mod_stp1sg
      implicit none
      integer, intent(in) :: ipotsy, iop, ipotsy2, j1max, j2min, j2max
      real(8), intent(in) :: brot, crot, delta, drot, e1max
      integer :: j1, k1, eps1, j2, ilev, nlev, igroup
      real(8) :: roteng
      common /coconv/ econv
      real(8) :: econv
c
c     Only C3v/D3h symmetric top implemented
      if (ipotsy .ne. 3) call raise_2pi1sg('ONLY C3v/D3h SYMMETRIC '
     $     //'TOP IMPLEMENTED')
c     First pass (nlev .eq. -1): Count the number of levels
c     Second pass (nlev is the number of channels): Save level info
      nlev = -1
      do while (.true.)
         if (nlev .ne. -1) then
            if (allocated(levs)) deallocate(levs)
            allocate(levs(nlev))
         end if
         ilev = 0
         do j1 = 0, j1max
            do k1 = 0, j1
               do eps1 = -1, 1, 2
                  if (k1 .eq. 0 .and. eps1 .eq. -1) cycle
                  if (mod(k1, 3) .eq. 0) then
                     igroup = 0
                  else
                     igroup = 2
                  end if
                  if ((k1.ne.0) .and. (mod(k1,3).eq.0) .and.
     $               ((iop.eq.17) .or. (iop.eq.34))) goto 25
                  if (eps1 * (-1) ** j1 .eq. 1) then
                     igroup = igroup + 1
                     if (k1.eq.0 .and. iop.eq.17) goto 25
                  else
                     if (k1.eq.0 .and. iop.eq.34) goto 25
                  end if
                  if (.not. btest(iop, igroup)) cycle
25                roteng = brot * j1 * (j1 + 1)
     $                 + (crot - brot) * k1 ** 2
                  if (eps1 * (-1) ** j1 .eq. 1) then
                     if (iop.eq.17 .and. k1.eq.0) then
                        roteng = roteng
                     else
                        roteng = roteng + delta
                     end if
                  else
                     if (iop.eq.17 .and. k1.eq.0) then
                        roteng = roteng + delta
                     end if
                  end if
                  if (roteng .gt. e1max) cycle
                  if (nlev .eq. -1) then
                     ilev = ilev + (j2max - j2min) / ipotsy2 + 1
                     cycle
                  end if
                  do j2 = j2min, j2max, ipotsy2
                     ilev = ilev + 1
                     levs(ilev)%j1 = j1
                     levs(ilev)%k1 = k1
                     levs(ilev)%eps1 = eps1
                     levs(ilev)%j2 = j2
                     levs(ilev)%e = (roteng + drot * j2 * (j2 + 1))
     $                    / econv
                  end do
               end do
            end do
         end do
         if (nlev .eq. -1) then
            nlev = ilev
         else
            exit
         end if
      end do
      return
      end subroutine genlev_stp1sg
c     ------------------------------------------------------------------
      subroutine prtlev_stp1sg()
      use mod_stp1sg
      implicit none
      integer :: i
      common /coconv/ econv
      real(8) :: econv
      write (6, 10)
 10   format (/' ** CC SYMMETRIC TOP-LINEAR MOLECULE'
     $   //10x, 'SORTED LEVEL LIST'/
     $   '    N   J1   K1  EPS1  J2   ENT(CM-1)')
      do i = 1, size(levs)
         write (6, 20) i, levs(i)%j1, levs(i)%k1, levs(i)%eps1,
     $        levs(i)%j2, levs(i)%e * econv
      end do
20    format(5i5,f11.3)
      return
      end subroutine prtlev_stp1sg
c     ------------------------------------------------------------------
      subroutine prtchn_stp1sg()
      use mod_stp1sg
      implicit none
      integer :: i, ilev
      common /coconv/ econv
      real(8) :: econv
      write (6, 10)
 10   format (/, ' ** CHANNEL LIST')
      do i = 1, size(chns)
         ilev = chns(i)%ilev
         write (6, 280) i, levs(ilev)%j1, levs(ilev)%k1,
     $        levs(ilev)%eps1, levs(ilev)%j2, chns(i)%j12,
     $        chns(i)%l, levs(ilev)%e * econv
 280     format (i4, 1x, 'J1=', i2, ' K1=', i2, ' EPS1=', i2,
     $        2x, 'J2=', i2, 2x, 'J12=', i2, 1x, 'L=', i3, 2x,
     $        'E=', f9.4)
      end do
      write (6, *)
      return
      end subroutine prtchn_stp1sg
c     ------------------------------------------------------------------
      real(8) function vstp1sg(j1p, lp, j1, l, j2p, j2, j12p, j12,
     $     jtot, kp, k, lam1, mu1, lam2, lam, epsp, eps)
      implicit none
      integer :: j1p, kp, epsp, j2p, j12p, lp, j1, k, eps, j2, j12, l,
     $     jtot, lam1, mu1, lam2, lam
      integer :: iphase, tj1p, tj1, tlam1, tkp, tmu1, tk
      real(8) :: pref, threej, sixj, ninej, tf3jm0, tf3j, f6j, f9j
      real(8), parameter :: machep = epsilon(0d0)
c
      vstp1sg = 0d0
      iphase = epsp * eps * (-1) ** (j1p + j1 + lam2 + lam + mu1)
      if (iphase .eq. -1) return
      threej = tf3jm0(2 * j2p, 2 * lam2, 2 * j2)
      if (dabs(threej) .lt. machep) return
      threej = threej * tf3jm0(2 * lp, 2 * lam, 2 * l)
      if (dabs(threej) .lt. machep) return
c
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
c
      sixj = f6j(j12, l, jtot, lp, j12p, lam)
      if (dabs(sixj) .lt. machep) return
      ninej = f9j(j1, j2, j12, j1p, j2p, j12p, lam1, lam2, lam)
      if (dabs(ninej) .lt. machep) return
c
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
c
      vstp1sg = pref * dble(iphase) * threej * sixj * ninej
      return
      end function vstp1sg
c     ------------------------------------------------------------------
c     THIS SUBROUTINE GOVERNS THE INPUT/OUTPUT OF THE BASIS ROUTINE.
c     ONLY IREAD IS USED: RETURN DIRECTLY IF ZERO.
      subroutine systp1sg(irpot, readpt, iread)
      implicit none
c
      integer irpot, iread
      logical readpt
      character*(*) fname
c     NUMBER OF BASIS-SPECIFIC VARIABLES, MODIFY ACCORDINGLY.
      integer icod, ircod, lencod
      parameter (icod=6, ircod=5, lencod=icod+ircod)
      common /cosys/ scod(lencod)
      character*8 scod
c     INTEGER VARIABLES.  LEAVE THE FIRST TWO AS IT IS.
      common /cosysi/ nscode, isicod, ipotsy, iop, ipotsy2, j1max,
     $     j2min, j2max
      integer :: nscode, isicod, ipotsy, iop, ipotsy2, j1max,
     $     j2min, j2max
c     REAL VARIABLES.  LEAVE THE FIRST TWO AS IT IS.
      common /cosysr/ isrcod, junkr, brot, crot, delta, e1max, drot
      integer :: isrcod, junkr
      real(8) :: brot, crot, delta, e1max, drot
      character*40 potfil
      save potfil
c     DEFINE THE NAMES HERE
      scod(1)='IPOTSY'
      scod(2)='IOP'
      scod(3)='IPOTSY2'
      scod(4)='J1MAX'
      scod(5)='J2MIN'
      scod(6)='J2MAX'
      scod(7)='BROT'
      scod(8)='CROT'
      scod(9)='DELTA'
      scod(10)='E1MAX'
      scod(11)='DROT'
      nscode = lencod
      isicod = icod
      isrcod = ircod
c     KEEP THE FOLLOWING LINE
      if (iread .eq. 0) return
c     READ THE LAST FEW LINES OF THE INPUT FILE
      read (8, *, err=80) ipotsy, iop
      read (8, *, err=80) j1max, e1max
      read (8, *, err=80) j2min, j2max, ipotsy2
      read (8, *, err=80) brot, crot, delta
      read (8, *, err=80) drot
      read (8, *, err=80) potfil
      call loapot(10, potfil)
      close (8)
      return
 80   call raise_2pi1sg('error read from input file.')
      return
c     ------------------------------------------------------------------
      entry ptrstp1sg(fname, readpt)
      return
c     ------------------------------------------------------------------
      entry savstp1sg(readpt)
c     WRITE THE LAST FEW LINES OF THE INPUT FILE.
      write (8, 220) ipotsy, iop
 220  format (2i4, 25x,'ipotsy, iop')
      write (8, 230) j1max, e1max
 230  format (i4, f8.2, 21x, 'j1max, e1max')
      write (8, 231) j2min, j2max, ipotsy2
 231  format (3i4, 21x,'j2min, j2max, ipotsy2')
      write (8, 250) brot, crot, delta
 250  format (3f10.4, 3x, 'brot, crot, delta, emax' )
      write (8, 251) drot
 251  format (f10.4, 23x,'drot')
      write (8, *) potfil
      return
      end
c     ------------------------------------------------------------------
