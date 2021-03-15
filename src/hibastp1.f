* ----------------------------------------------------------------------
      subroutine bastp1(j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, k, ieps, jtemp, ktemp, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of a symmetric top molecule having no inversion douling with a
*  structureless atom or with an uncorrugated surface (basis type = 18)
*
*  NOTE:  this subroutine is a modification of bastp.  the parameter
*  ninv (appropriate to molecules like NH3 that have inversion
*  doubling) has been dropped.  input parameters changed also.
*
*  this version now has the capability of setting up calculations
*  for all 3 nuclear spin modifications of CD3.
*
*  author:  paul dagdigian
*  revision:  13-mar-2013 by p.j.dagdigian (capable of setting up 3
*     nuclear spin modifications)
*  revision:  4-jun-2013 by q. ma (fix a bug in counting anisotropic
*     terms)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index (eps * k) for
*              each channel
*  Note that we have adopted the following convention for the symmetry
*  index "is" so that on return the symmetric top molecular levels can be
*  uniquely identified by the two labels "j" and "is":
*           is = +/- k
*       where k is the projection quantum number
*       the plus and minus sign refers to levels with eps = +1 and -1,
*       respectively.
*  the symmetric top wave functions can be written as
*                   [|j,k> + eps * |j,-k>]
*  for k=0, where only eps=+1 states exist, the index "is" equals zero
*  The k > 0 levels of a symmetric top are nondegenerate.  For
*  symmetrical molecules (ihomo = T), the oath levels (k not a multiple
*  of 3, iop = -1) come in degenerate pairs (eps = +/-1), which are not
*  collisional connected.  We include only one set of these in the
*  basis [(1,1), (2,-2), (2,-1), (3,2), (3,1), (4,-4), ...].  The para levels
*  (k a multiple of 3, iop = 1) are nondegenerate.
*
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each energetically
*              distinct level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    k:        on return contains projection quantum number for each channel
*    ieps:     on return contains epsilon label for each channel
*    jtemp:    scratch array used to create channel list
*    ktemp:    scratch array used to create channel list
*    ietmp:    scratch array used to create channel list
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              to turn of this test, set rcut negative
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*              in cc calculation jtot is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true., then the molecule posseses interchange symmetry
*              (e.g. CH3), so that only the ortho or para levels will be
*              included depending on the value of the parameter iop in common
*              /cosysi/ (see below)
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(parity+l-jtot)=jlpar
*              where parity designates the parity of the molecular state
*              [by definition parity = ieps (-1)^(j + k) for symmetric
*              top wave functions]
*              in cs calculation jlpar is set equal to 1 in calling program
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     A=B rotational constant for prolate top
*    crot:     C rotational constant for prolate top
*    emax:     the maximum rotational energy (in cm-1) for a channel to be
*              included in the basis
*  variables in common block /cosysi/
*    nscod:    total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    ipotsy:   cylindrical symmetry of potential.  If ihomo = .true., only
*              terms with mu equal to an integral multiple of ipotsy
*              can be included in the potential.  Example:  for CH3, ipotsy = 3.
*    iop:      nuclear spin modification label for molecular levels. If ihomo=.true.
*              then for molecules with equivalent nuclei of spin = 1/2 (e.e. CH3) only
*              levels (denoted para) with k a multiple of 3, and only even j levels
*              for k = 0, will be included if iop=1; only levels (denoted ortho)
*              with k not a multiple of 3 will be included if iop=-1.
*                For a molecule with nuclear spins = 1 (e.g. CD3), set iop = -1
*              for rotational levels (with A1 nuclear spin symmetry) with k a
*              multiple of 3 and the even j levels in the k = 0 stack;
*              set iop = 2 for rotational levels (with A2 nuclear spin symmetry)
*              with k a multiple of 3 and the odd j levels in the k = 0 stack;
*              set iop= 1 for rotational levels (with E nuclear spin symmetry)
*              with k not equal to a multiple of 3 (e.g. k = 1, 2, 4, 5, 7, 8,…).
*    jmax:     the maximum rotational angular momentum for the symmetric top
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*               the zero of energy is assumed to be the j=0, k=0 level
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               stored in packed column form that is (1,1), (2,1), (3,1) ...
*               (n,1), (2,2), (3,2) ... (n,2), etc.
*  variable in common block /coiv2/
*   iv2:        row+column index of v2 matrix for each non-zero element
*  variable in common block /coconv/
*   econv:        conversion factor from cm-1 to hartrees
*   xmconv:       converson factor from amu to atomic units
*  subroutines called:
*   vlmstp:    returns angular coupling coefficient for particular
*              choice of channel index (this subroutine is in file hibastp.f
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical flaghf, csflag, clist, flagsu, ihomo, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop, jmax
      common /coipar/ iiipar(9), iprint
      common /cosysr/ isrcod, junkr, brot, crot, emax
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      common /coamat/ ietmp(1)
      dimension j(1), l(1), jhold(1), ehold(1), is(1), k(1),
     :          ieps(1), jtemp(1), ishold(1), ktemp(1)
      zero = 0.d0
      two = 2.d0
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (6, 5)
        write (9, 5)
5       format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
        stop
      end if
      if (flagsu .and. .not. csflag) then
        write (6, 6)
        write (9, 6)
6       format
     :   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
        stop
      end if
      do 25  i = 1, nterm
        if (ihomo .and. mod(mproj(i), ipotsy) .ne. 0 ) then
          write (20, 20) i, mproj(i), ipotsy
          write (6, 20) i, mproj(i), ipotsy
20        format(/' *** MPROJ(',i2,')=',i2,
     :     ' .NE. INTEGER MULTIPLE OF', i2, ' FOR IHOMO=T; ABORT ***')
          stop
        end if
25    continue
      xjtot = jtot
      nsum = 0
      do 35  i = 1, nterm
        if (mproj(i) .gt. lammin(i) ) then
          write (6, 30) i, mproj(i), lammin(i)
          write (9, 30) i, mproj(i), lammin(i)
30        format (' *** MPROJ=',i2,' > LAMMIN=',i2,
     :            ' FOR TERM',i2,'; ABORT ***')
          stop
        end if
        idl = 1
        if (ihomo) idl = 2
        nsum = nsum + (lammax(i) - lammin(i)) / idl + 1
35    continue
      if (nlammx .lt. nsum) then
        write (6, 40) nsum, nlammx
        write (9, 40) nsum, nlammx
40      format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2,
     :          ' .GT. NLAMMX=', i2,'; ABORT')
        stop
      end if
      if (nsum .ne. nlam) then
        write (6, 45) nsum, nlam
        write (9, 45) nsum, nlam
45      format (' *** NLAM IN BASIS=', i3,' .NE. NLAM FROM SYSDAT=',
     :           i3, '; ABORT ***')
        stop
      end if
      if (abs(iop).ne.1 .and. iop.ne.2) then
        write (6,145) iop
        write (9,145) iop
145     format (' *** INVALID VALUE OF IOP (',i2,');',
     :    ' ABORT ***')
        stop
      end if
      if (bastst) write (6, 46) nsum
      write (9, 46) nsum
46    format (' ** TOTAL NUMBER OF ANISTROPIC TERMS IN POTENTIAL =',
     :        i3)
      nlam = nsum
      if (clist) then
        if (flagsu) then
          if (ihomo) then
            if (bastst)
     :      write (6,55) rmu * xmconv, brot, crot, ipotsy, iop,
     :                   ered * econv, jtot, nu
            write (9,55) rmu * xmconv, brot, crot, ipotsy, iop,
     :                   ered * econv, jtot, nu
55        format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=', f7.3,/,
     :            '     POT-SYM=', i2,'  O/P=',i2,
     :           '  E=', f7.2, '       LBAR=', i5, '  NU=', i3)
          else
            write (6,60) rmu * xmconv, brot, crot, ipotsy,
     :                   ered * econv, jtot, nu
            write (9,60) rmu * xmconv, brot, crot, ipotsy,
     :                   ered * econv, jtot, nu
60        format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=', f7.3,/,
     :            '     POT-SYM=', i2,
     :           '     E=', f7.2, '       LBAR=', i5, '  NU=', i3)
          end if
        else
          if (csflag) then
            if (ihomo) then
              if (bastst)
     :        write (6,65) rmu * xmconv, brot, crot, ipotsy, iop,
     :                   ered * econv, jtot, nu
              write (9,65) rmu * xmconv, brot, crot, ipotsy, iop,
     :                   ered * econv, jtot, nu
65          format(/,' **  CS SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,/,
     :            '     POT-SYM=', i2,'  O/P=',i2,
     :             '  E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
            else
              if (bastst)
     :        write (6,75) rmu * xmconv, brot, crot, ipotsy,
     :                   ered * econv, jtot, nu
              write (9,75) rmu * xmconv, brot, crot, ipotsy,
     :                   ered * econv, jtot, nu
75          format(/,' **  CS SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,/,
     :            '     POT-SYM=', i2,
     :             '     E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
            end if
          else
            if (ihomo) then
              if (bastst)
     :        write (6,80) rmu * xmconv, brot, crot, ipotsy, iop,
     :                     ered * econv, jtot, jlpar
              write (9,80) rmu * xmconv, brot, crot, ipotsy, iop,
     :                     ered * econv, jtot, jlpar
80          format(/,' **  CC SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,/,
     :            '     POT-SYM=', i2,'  O/P=',i2,
     :             '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
            else
              if (bastst)
     :        write (6,85) rmu * xmconv, brot, crot,
     :                     ered * econv, jtot, jlpar
              write (9,85) rmu * xmconv, brot, crot,
     :                     ered * econv, jtot, jlpar
85          format(/,' **  CC SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,/,
     :             '     E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
            end if
          end if
        end if
        if (.not. flagsu) write (9,90) rcut
90      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      end if
*
*  set up rotational basis
*
*  first set up list of j,k levels within the requested
*  nuclear spin modification (iop value)
      nlist = 0
      do 105  ki = 0, jmax
        do 105  ji = ki, jmax
          numeps = 2
          if (ki.eq.0) numeps = 1
          do 105 iep = 1, numeps
*  check that level has correct nuclear permutation symmetry
            if (ihomo) then
*  levels with iop = -1 (called ortho levels for CH3) and iop = 2
*  levels have k a multiple of 3.
*  all other levels have iop = 1 (called para levels for CH3)
              if (iop.eq.-1 .and.
     :          mod(ki,ipotsy).ne.0) goto 105
              if (iop.eq.2 .and.
     :          mod(ki,ipotsy).ne.0) goto 105
              if (iop.eq.1 .and. mod(ki,ipotsy).eq.0) goto 105
*  for ki = 0, the even j levels are ortho (iop = -1) and the odd j
*  levels exist only for iop = 2.
              if (iop.eq.-1 .and. mod(ji,2).ne.0 .and.
     :          ki.eq.0) goto 105
              if (iop.eq.2 .and. mod(ji,2).eq.0 .and.
     :          ki.eq.0) goto 105
*  for iop = -1, levels with ki > 0 and a multiple of ipotsy,
*  only the level with (-1) ** (ji+iep) exists
*  for iop = 2 levels with ki > 0 and a multiple of ipotsy,
*  only the level with -(-1) ** (ji+iep) exists
              if (iop.eq.-1 .and. ki.gt.0 .and. mod(ki,ipotsy).eq.0
     :          .and. (-1)**(ji+iep).eq.1) goto 105
              if (iop.eq.2 .and. ki.gt.0 .and. mod(ki,ipotsy).eq.0
     :          .and. (-1)**(ji+iep).eq.-1) goto 105
*  include only one of the degenerate pair of levels for iop = -1
              ies = 1 - 2*(iep - 1)
              if (iop.eq.1 .and. ies*(-1)**ji.eq.1) goto 105
            end if
*  check to make sure rotational energy is less than emax
            roteng = brot * ji * (ji + 1) + (crot - brot) * ki * ki
            if (roteng .gt. emax) go to 105
            nlist = nlist + 1
            ehold(nlist)= roteng / econv
            jtemp(nlist) = ji
            ktemp(nlist) = ki
            ietmp(nlist) = 1 - 2*(iep - 1)
105   continue
*  now sort this list in terms of increasing energy
      if (nlist .gt. 1) then
        do 120 i1 = 1, nlist - 1
          esave = ehold(i1)
          do 115 i2 = i1 + 1, nlist
            if (ehold(i2) .lt. esave) then
*  state i2 has a lower energy than state i1, switch them
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
115       continue
120     continue
      end if
c
c     Determine the number of levels that are open asymtotically
c
      nlevop = 0
      do while (nlevop .lt. nlist .and. ehold(nlevop + 1) .le. ered)
         nlevop = nlevop + 1
      end do
c
*     print this list if bastst = .true.
      if (bastst) then
         write (6, 130)
130     format (/,10x,
     :    'SORTED LEVEL LIST',/,'   N   J   K  EPS   EINT(CM-1)')
        do  140 i = 1, nlist
          ecm = ehold(i) * econv
          write (6, 135) i, jtemp(i), ktemp(i), ietmp(i), ecm
135       format (4i4,f10.3)
140     continue
      end if
*  now set up channel and level list for scattering calculation
      n = 0
      nlevel = 0
      do 170  njk = 1, nlist
        ki = ktemp(njk)
        ji = jtemp(njk)
        nlevel = nlevel + 1
        jhold(nlevel) = ji
        ishold(nlevel) = ietmp(njk) * ki
*  here for cs calculations; include levels only if j is at least
*  equal to coupled states projection index
        if (csflag) then
*  here for cs calculations; include state only if j at least equal to coupled
*  states projection index
          if (ji .ge. nu) then
            n = n + 1
            if (n .gt. nmax) then
              write (6, 150) n, nmax
              write (9, 150) n, nmax
150           format(/' *** NCHANNELS=', i5,
     :            ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
              stop
            end if
            j(n) = ji
            k(n) = ki
            ieps(n) = ietmp(njk)
            is(n) = ietmp(njk) * ki
            eint(n) = ehold(njk)
            l(n) = jtot
            cent(n) = jtot * (jtot + 1)
          end if
        else if (.not. csflag) then
*  here for cc calculations.  first calculate range of orbital angular
*  momentum quantum numbers allowed for this state
*
*  determine parity of molecular level.  expression for parity
*  taken from eq. 5 of m.h.alexander, jcp 77, 1855 (1982)
          ipar = ietmp(njk) * (-1) ** (ji + ki)
          lmax = jtot + ji
          lmin = iabs (jtot - ji)
          do 155  li = lmin, lmax
            lpar = (-1) ** (li - jtot)
*  check to see if this channel has the desired parity, if so, include it
            if (ipar * lpar .eq. jlpar) then
              n = n + 1
              if (n .gt. nmax) then
                write (6, 150) n, nmax
                write (9, 150) n, nmax
                stop
              end if
              is(n) = ietmp(njk) * ki
              k(n) = ki
              ieps(n) = ietmp(njk)
              j(n) = ji
              eint(n) = ehold(njk)
              l(n) = li
              cent(n) = li * (li + 1)
            end if
155       continue
        end if
170   continue
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions or for rcut < 0
      if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
        emin = 1. e+7
        do 290  i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is open asymptotically
            if ( jtot * (jtot + 1) / (two * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this condition is met
            end if
          end if
290     continue
*  now eliminate all channels with eint .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn = 0
          do 300 i = 1, n
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              nn = nn + 1
              eint(nn) = eint(i)
              j(nn) = j(i)
              ieps(nn) = ieps(i)
              is(nn) = is(i)
              cent(nn) = cent(i)
              k(nn) = k(i)
              l(nn) = l(i)
            end if
300       continue
*  reset number of channels
          n = nn
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
      if (nu .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      else
        if (n.gt.ntop) then
          write (6, 303) nu, n, ntop
          write (9, 303) nu, n, ntop
303       format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3,
     :            '; ABORT **',/,
     :    '     CHECK RCUT')
          call exit
        endif
      end if
*  now list channels if requested
      if (clist) then
        if (.not.csflag) then
          if (bastst) write (6, 305)
          write (9,305)
305       format
     :     (/'   N   J  EPS  K   L    EINT(CM-1)')
        else
          if (bastst) write (6, 310) nu
          write (9,310) nu
310       format
     :     (/'   N   J  EPS  K   L    EINT(CM-1) ** NU = ',i2)
        end if
        do 330  i = 1, n
          ecm = eint(i) * econv
          if (bastst)
     :      write (6, 320) i, j(i), ieps(i), k(i), l(i), ecm
            write (9, 320) i, j(i), ieps(i), k(i), l(i), ecm
320         format (5i4, f10.3)
330     continue
      end if
*  now calculate coupling matrix elements
*
*  increment lb by 2 for symmetrical molecules
      idell = 1
      if (ihomo) idell = 2
*  i counts v2 elements
*  inum counts v2 elements for given lambda
*  ilam counts numver of v2 matrices
*  ij is address of given v2 element in present v2 matrix
      i = 0
      if (bastst.and. iprint.ge. 2) then
        write (6, 340)
        write (9, 340)
340     format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
      end if
      lamsum = 0
      ilam = 0
      do 400 iterm = 1, nterm
        lbmin = lammin(iterm)
*  if bastst = .true., then get the matrix elements of the lb=0 term
*  in the potential
        if (bastst .and. iterm .eq. 1) lbmin = 0
        do 390 lb = lbmin, lammax(iterm), idell
*  ilam is the index for the next term in the potential matrix
*  lb is the actual value of lambda
          ilam = ilam + 1
          mu = mproj(iterm)
          inum = 0
          ij = 0
          do 355  icol = 1, n
            do 350  irow = icol, n
              ij = ntop * (icol - 1) +irow
              lrow = l(irow)
              if (csflag) lrow = nu
              call vlmstp (j(irow), lrow, j(icol), l(icol), jtot,
     :                     k(irow), k(icol), lb, mu, ieps(irow),
     :                     ieps(icol), vee, csflag)
              if (vee .ne. zero) then
                i = i + 1
                if (i .le. nv2max) then
                  inum = inum + 1
                  v2(i) = vee
                  iv2(i) = ij
                  if (bastst.and. iprint.ge.2) then
                    write (6, 345) ilam, lb, icol, irow, i, iv2(i),
     :                             vee
                    write (9, 345) ilam, lb, icol, irow, i, iv2(i),
     :                             vee
345                 format (i4, 2i7, 2i6, i6, g17.8)
                  end if
                end if
              end if
350         continue
355       continue
          if (i .le. nv2max) lamnum(ilam) = inum
          if (bastst) then
            write (6, 370) ilam, lb, mu, lamnum(ilam)
            write (9, 370) ilam, lb, mu, lamnum(ilam)
370         format ('ILAM=',i3,' LAM=',i3,' MU=',i3,
     :         ' LAMNUM(ILAM) = ',i6)
          end if
          lamsum = lamsum + lamnum(ilam)
390     continue
400   continue
      if ( i.gt. nv2max) then
         write (6, 450) i, nv2max
         write (9, 450) i, nv2max
450      format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
         stop
      end if
      if (clist .and. bastst) then
        write (6, 460) lamsum
        write (9, 460) lamsum
460     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i6)
      end if


*      write(6,*) (lamnum(i),i=1,nlam)



      return
      end
