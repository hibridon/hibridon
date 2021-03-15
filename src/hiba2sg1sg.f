* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of a 2sigma+ diatomic with a closed-shell (1sigma+) diatomic.
*  the second molecule can be homonuclear.
*
*  this subroutine is an extension of ba1sg1sg, which treats collisions
*  ot two unlike 1sigma molecules
*
*  pot matrix elements are computed assuming the angular dependence defined
*  in the appendix of green, jcp 62, 2271 (1975)
*
*  author:  paul dagdigian
*  current revision date:  20-jun-2017 by pjd
* --------------------------------------------------------------------
c     This module contains (explictly) the number of terms and their
c     indices in the expansion of the PES.  Its contents should be
c     filled in the pot routine.
c
c     This module replaces lammin, lammax, mproj in hibridon and
c     is compiled with ba1sg1sg
c
c      module mod_1sg1sg
c      implicit none
c
c      type lm_type
c      integer :: l1, l2, ltot
c      end type lm_type
c
c      type(lm_type), dimension(:), allocatable :: lms
c      end module mod_1sg1sg
* --------------------------------------------------------------------
      subroutine ba2s1sg (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains combined rotational quantum numbers for each
*              channel.  in terms of the rotational quantum numbers of each
*              molecule we have:  j = 10 * j1 + j2
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       additional quantum index of each channel:  on return this is
*              set equal to j12 (the vector resultant of j1 and j2) times nsym
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level, irrespective of projection degeneracy and the
*              degeneracy associated with different values of j12
*              note that jhold = 10 * j1hold + j2hold
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return this is set equal to nsym
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    sc1:      scratch vector (not used here)
*    sc2:      scratch vector (not used here)
*    sc3:      scratch vector (not used here)
*    sc4:      scratch vector (not used here)
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*    jtot:     total angular momentum
*              in cc calculation jtot+1/2 is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons (this variable
*              should always be false in the present case)
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              the homonuclear option is not specifically implemented here
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              parity=(-1)**jlpar (by definition parity=(-1)**(j1+j2+l)
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*  variables in common block /cosysr/
*    isrcod:   number of real system dependent variables
*    b1rot:    rotational constant in cm-1 of molecuole 1
*    d1rot:    centrifugal distortion constant in cm-1 of molecule 1
*    gamma:    spin-rotation constant of molecule 1
*    b2rot:    rotational constant in cm-1 of molecuole 2
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   numbr of integer system dependent variables
*    n1max:    maximum rotational angular momentum n of molecule 1
*    j2min:    minimum rotational angular momentum of molecule 2
*    j2max:    maximum rotational angular momentum of molecule 2
*    ipotsy2:  symmetry of potential.  set to 2 for homonuclear
*              molecule 2, set to 1 for heteronuclear molecule 2
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*               nlam is set equal to nterm; see above
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
*  subroutines called:
*   vlmlml:    returns molecule-molecule angular coupling coefficient for
*              particular choice of channel index
* --------------------------------------------------------------------
      use mod_1sg1sg
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"
      common /coipar/ iiipar(9), iprint
      common /cosysi/ nscode, isicod, n1max, j2min, j2max,
     :  ipotsy2
      common /cosysr/ isrcod, junkr, b1rot, d1rot, gamma, b2rot
      common /coselb/ ibasty
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coj12/  j12(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      dimension j(1), l(1), is(1), jhold(1), ehold(1),
     :    sc1(1), sc2(1), sc3(1), sc4(1), ishold(1)
      data zero, ione, itwo, ithree / 0.d0, 1, 2, 3 /
*  check for consistency in the values of flaghf and csflag
      if (.not. flaghf) then
        write (9, 5)
5       format (' *** FLAGHF = .FALSE. FOR DOUBLET SYSTEM; ABORT ***' )
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (flagsu) then
        write (6, 10)
        write (9, 10)
10      format ('  *** FLAGSU = .TRUE. FOR',
     :         ' MOLECULE-MOLECULE COLLISION; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (csflag) then
        write (6, 25)
        write (9, 25)
25      format (' *** CSFLAG=.TRUE.; NOT ALLOWED IN 2SIG-1SIG BASIS;',
     :          ' ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*  check that ipotsy2 is either 1 or 2
      if (ipotsy2.lt.1 .or. ipotsy2.gt.2) then
        write (6, 27) ipotsy2
        write (9, 27) ipotsy2
27      format (' *** IPOTSY2 .NE. 1 .OR.2;',
     :          ' ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (bastst) then
        write (6,  90) rmu * xmconv, ered * econv, jtot, jlpar
        write (9,  90) rmu * xmconv, ered * econv, jtot, jlpar
90      format(/2x, ' **  CC  2SIGMA-1SIGMA ** RMU=',f10.6,
     :     '  E=', f9. 3, 2x, 'JTOT=', i4,
     :     '  JLPAR=', i2)
      end if
      write (9, 95) rcut
95    format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f7.2)
*  assign quantum numbers and energies for rotational levels
      nlevel = 0
      do n1 = 0, n1max
        do iep1 = 2, 1, -1
          if (n1.gt.0 .or. (n1.eq.0 .and. iep1.eq.1)) then
            do j2 = j2min, j2max, ipotsy2
              nlevel = nlevel + 1
              if (iep1 .eq. 1) then
                j1  = n1
              else
                j1 = n1 - 1
              end if
              jhold(nlevel) = 10 * j1 + j2
              ieps1 = (1 - iep1) * 2 + 1
              ishold(nlevel) = ieps1
              nn1 = n1 * (n1 + 1.d0)
c
c  DON'T subtract energy of lowest H2 level
              ehold(nlevel) = b1rot * nn1 - d1rot * nn1**2
     :            + b2rot * j2 * (j2 + 1.d0)
              if (ifs1.eq.1) then
                ehold(nlevel) = ehold(nlevel) + 0.5d0 * gamma * n1
              else
                ehold(nlevel) = ehold(nlevel)
     :             - 0.5d0 * gamma * (n1 + 1.d0)
              end if
              ehold(nlevel) = ehold(nlevel) / econv
            end do
          end if
        end do
      end do
*
*  form list of all energetically open rotational levels included in the
*  calculations and their energies
      nlevop = 0
      do 80 i = 1, nlevel - 1
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        else
          do 75 ii = i + 1, nlevel
            if (ehold(ii) .le. ered) then
              nlevop = nlevop + 1
              ikeep = jhold(i)
              jhold(i) = jhold(ii)
              jhold(ii) = ikeep
              ikeep = ishold(i)
              ishold(i) = ishold(ii)
              ishold(ii) = ikeep
              ekeep = ehold(i)
              ehold(i) = ehold(ii)
              ehold(ii) = ekeep
              goto 80
            end if
75        continue
        end if
80    continue
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
*
*  set up CC scattering channel basis
      n = 0
      xjtot = jtot
      do i = 1, nlevel
        j1 = jhold(i)/10
        j2 = mod(jhold(i),10)
        j12mn = abs(j1 + 0.5d0 - j2)
        j12mx = j1 + j2
        j12min = j12mn
        j12max = j12mx
        do j12i = j12min, j12max
          lmin = abs(jtot - j12i)
          lmax = jtot + j12i + 1
          do li = lmin, lmax
c  parity of 2sigma molecule is iep1 * (-1) ^ (j - 1/2)
            ix = ishold(i) * (- 1) ** (j1 + j2 - jtot + li)
            if (ix .eq. jlpar) then
              n = n + 1
              if (n .le. nmax) then
                j(n) = jhold(i)
                is(n) = ishold(i)
                l(n) = li
                j12(n) = j12i
                eint(n) = ehold(i)
                cent(n) = li * (li + 1.d0)
              else
                write (9, 230) n, nmax
                write (6, 230) n, nmax
230             format(/' *** NCHANNELS=', i5,' .GT. ',
     :              'MAX DIMENSION OF ',i5,'; ABORT')
                if (bastst) then
                  return
                else
                  call exit
                end if
              end if
            end if
          end do
        end do
      end do
*
*  order channel basis functions in order of increasing energy
*  sort algorithm from Numerical Recipes
      if (n .gt. 1) then
        do 212 jj = 2, n
          ekeep = eint(jj)
          jkeep = j(jj)
          ikeep = is(jj)
          lkeep = l(jj)
          j12kp = j12(jj)
          ckeep = cent(jj)
          do 211 ii = jj-1, 1, -1
            if (eint(ii).le.ekeep) goto 210
            eint(ii+1) = eint(ii)
            j(ii+1) = j(ii)
            is(ii+1) = is(ii)
            l(ii+1) = l(ii)
            j12(ii+1) = j12(ii)
            cent(ii+1) = cent(ii)
211       continue
          ii = 0
210       eint(ii+1) = ekeep
          j(ii+1) = jkeep
          is(ii+1) = ikeep
          l(ii+1) = lkeep
          j12(ii+1) = j12kp
          cent(ii+1) = ckeep
212     continue
      end if
*
*  now check to see if any of the open channels are closed at r = rcut
      if (rcut .gt. 0.d0 .and. .not.boundc) then
        emin = 1.e+7
        do 120 i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is open asymptotically
            if (jtot * (jtot + 1) / (2.d0 * rmu * rcut ** 2)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r  = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest rotational energy for which this condition is met
            end if
          end if
120     continue
*  now eliminate all channels with eint .gt. emin if any of the channels
*  are open asymptoptically but closed at r  = rcut
        if (emin .lt. ered) then
          nn = 0
          do 130 i = 1, n
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              nn = nn + 1
              j(nn) = j(i)
              eint(nn) = eint(i)
              is(nn) = is(i)
              l(nn) = l(i)
              j12(nn) = j12(i)
              cent(nn) = cent(i)
            end if
130       continue
*  reset number of channels
          n = nn
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
      if (nn .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure that this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) then
          ntop = ntop + 1
        else
          if (n .gt. ntop) then
            write(6,160) nn, n, ntop
            write(9,160) nn, n, ntop
160         format (' *** NCH = ',i3,' AT NU = ',i2,' .GT. NTOP = ',i3,
     :          '; ABORT ***',/,
     :      '      CHECK RCUT')
            call exit
          end if
        end if
      end if
*
*  now list channels if requested
      if (clist) then
        if (bastst) write (6, 200)
        write(9, 200)
200     format(/'   N     J1     EPS1    J2     J12      L   ',
     :      '   EINT(CM-1)')
        do i = 1, n
          ecm = eint(i) * econv
          jj1 = j(i)/10
          fj1 = jj1 + 0.5d0
          j2 = mod(j(i),10)
          fj12 = j12(i) + 0.5d0
          write (6, 220) i, fj1, is(i), j2, fj12, l(i), ecm
          write (9, 220) i, fj1, is(i), j2, fj12, l(i), ecm
220       format (i4, f8.1, 2i7, f8.1, i7, f13.3)
        end do
      end if
c
c     Calculate coupling matrix elements
      if (bastst .and. iprint.eq.2) then
        write (6, 285)
        write (9, 285)
285     format (/' ILAM  L1   L2  LTOT   ICOL IROW      I',
     :    '      IV2         VEE')
      end if
      i = 0
      lamsum = 0
      do 400 ilam = 1, nlam
c     ilam denotes a particular L1,L2,L term
        ll1 = lms(ilam)%l1
        ll2 = lms(ilam)%l2
        lltot = lms(ilam)%ltot
        inum = 0
        do 355 icol = 1, n
          j1c = j(icol)/10
          iepsc = is(icol)
          j2c = mod(j(icol),10)
          j12c = j12(icol)
          lc = l(icol)
          do 350 irow = icol, n
            ij = ntop * (icol - 1) + irow
            j1r = j(irow)/10
            iepsr = is(irow)
            j2r = mod(j(irow),10)
            j12r = j12(irow)
            lr = l(irow)
c     initialize potential to zero
            vee = zero
            call v2sgsg(j1r,iepsr,j2r,j12r,lr,j1c,iepsc,
     :          j2c,j12c,lc,jtot,ll1,ll2,lltot,vee)
            if (vee .ne. zero) then
              i = i + 1
              if (i .le. nv2max) then
                inum = inum + 1
                v2(i) = vee
                iv2(i) = ij
                if (bastst .and. iprint.ge.2) then
                  write (6, 290) ilam, ll1, ll2, lltot,
     :                icol, irow, i, iv2(i), vee
                  write (9, 290) ilam, ll1, ll2, lltot,
     :                icol, irow, i, iv2(i), vee
290               format (i4, 3i5, 2x, 2i5, 2i8, e20.7)
                end if
              end if
            end if
350       continue
355     continue
        if (i .le. nv2max) lamnum(ilam) = inum
        if (bastst) then
          write (6, 370) ilam, lamnum(ilam)
          write (9, 370) ilam, lamnum(ilam)
370       format ('ILAM=',i4,' LAMNUM(ILAM) = ',i7)
        end if
        lamsum = lamsum + lamnum(ilam)
400   continue
      if (i .gt. nv2max) then
        write (6, 410) i, nv2max
        write (6, 410) i, nv2max
410     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :      ' .GT. NV2MAX=',i6,'; ABORT ***')
      end if
      if (clist .and. bastst) then
        write (6, 420) lamsum
        write (9, 420) lamsum
420     format (' *** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :      i8)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine v2sgsg(j1r,iepsr,j2r,j12r,lr,j1c,iepsc,j2c,j12c,
     :          lc,jtot,l1,l2,ltot,vee)
*  subroutine to calculate v-lambda matrices for close-coupled
*  treatment of collisions of a 2sigma and a 1sigma molecule
*  electronic state with a (different) molecule in a 1sigma state
*
*  matrix elements adapted from 1sigma-1sigma pot matrix elements
*  originally derived by s. green, jcp 62, 2271 (1973)
*
*  subroutines called:
*     xf3j, xf6j, xf9j
*
*  written  by p. dagdigian
*  current revision date:  19-jun-2017 (p.dagdigian)
* --------------------------------------------------------------------
      use mod_1sg1sg
      implicit double precision (a-h,o-z)
      data zero, half, one, two/ 0.d0, 0.5d0, 1.d0, 2.d0/,
     :  sq4pi3 / 44.546623974653656d0 /
      vee = zero
      xj1r = j1r + half
      xj2r = j2r
      xj12r = j12r + half
      xlr = lr
      xj1c = j1c + half
      xj2c = j2c
      xj12c = j12c + half
      xlc = lc
      xjtot = jtot + half
      xl1 = l1
      xl2 = l2
      xltot = ltot
      iph = xj1r + xj1c + xl1
      phase = iepsr * iepsc * (-1.d0) ** iph
      if (phase .eq. 1.d0) return

c      iph = xjtot + xl1 - xl2 + xj1c - xj2c + xj12r
c     :    - xlc - xlr - 0.5d0

      iph = xjtot + xj1c + xj2c + xj12r - 0.5d0

      phase = 1.d0
      if (iph .ne. 2*(iph/2)) phase = -1.d0
      facj = (two*xltot+one) * sqrt((two*xl1+one) * (two*xl2+one)
     :    * (two*xj1r + one) * (two*xj2r + one) * (two*xj12r + one)
     :    * (two*xlr + one)
     :    * (two*xj1c + one) * (two*xj2c + one) * (two*xj12c + one)
     :    * (two*xlc + one) )
      cg1 = xf3j (xlr, xltot, xlc, zero, zero, zero)
      if (cg1 .eq. zero) return
      cg2 = xf3j (xj1r, xl1, xj1c, -half, zero, half)
      if (cg2 .eq. zero) return
      cg3 = xf3j (xj2r, xl2, xj2c, zero, zero, zero)
      if (cg3 .eq. zero) return
      f6j = xf6j (xj12c, xlc, xjtot, xlr, xj12r, xltot)
      if (f6j .eq. zero) return
      f9j = xf9j (xj1c, xj2c, xj12c, xj1r, xj2r, xj12r,
     :  xl1, xl2, xltot)
      if (f9j .eq. zero) return
      vee = phase * facj * cg1 * cg2 * cg3 * f6j * f9j / sq4pi3
      return
      end
* --------------------------------eof---------------------------------

