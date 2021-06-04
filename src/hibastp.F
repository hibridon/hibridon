* ----------------------------------------------------------------------
      subroutine bastp (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, k, ieps, jtemp, ktemp, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of a symmetric top molecule possessing inversion doubling (e.g. NH3)
*  with a structureless atom or with an uncorrugated surface
*  (basis type = 6)
*
*  authors:  millard alexander
*  revision:  11-aug-2011 by p.j.dagdigian
*    rotational channel basis now chosen by j < jmax and ehold < emax
*    as in bastp1 basis routine.  also changed printout of basis in bastst
*    also, added inversion splitting, labeled as delta
*  revision:  4-jun-2013 by q. ma (fix a bug in counting anisotropic
*     terms)
*  revision:  22-oct-2013 by p.j.dagdigian (capable of setting up 3
*     nuclear spin modifications or nuclear spin = 1)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry of inversion vibration for
*              each channel times k
*  Note that we have adopted the following convention for the symmetry
*  index "is" so that on return the symmetric top molecular levels can be
*  uniquely identified by the two labels "j" and "is":
*           is = +/- k
*       where k is the projection quantum number
*       the plus sign refers to the "plus" inversion levels and the minus
*       sign to the "minus" inversion levels
*       for k=0, where only eps=+1 states exist, the index "is" is equal to zero
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
*              (e.g. NH3), so that only the ortho or para levels will be
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
*              [by definition parity=isym (-1)^k - see Eq. (A21) of
*              S. Green, J. Chem. Phys. 73, 2740 (1980) ]
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
*    delta:    inversion splitting
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
*              can be included in the potential.  Example:  for NH3, ipotsy = 3.
*    iop:      nuclear spin modification label for molecular levels. If ihomo=.true.
*              then for molecules with equivalent nuclei of spin = 1/2 (e.e. NH3) only
*              para levels with k a multiple of 3, and only symmectric inverson levels
*              with k = 0, will be included if iop=1; only ortho levels
*              with k not a multiple of 3 will be included if iop=-1.
*                For a molecule with nuclear spins = 1 (e.g. ND3), set iop as follows:
*              iop =  2 for rotational levels with A1 nuclear spin symmetry
*              iop = -1 for rotational levels with A2 nuclear spin symmetry
*              iop =  1 for rotational levels with E nuclear spin symmetry
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
*              choice of channel index
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical flaghf, csflag, clist, flagsu, ihomo, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop, jmax
      common /coipar/ iiipar(9), iprint
      common /cosysr/ isrcod, junkr, brot, crot, delta, emax
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
5      format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
        stop
      end if
      if (flagsu .and. .not. csflag) then
        write (6, 6)
        write (9, 6)
6     format
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
        nsum = nsum + lammax(i) - lammin(i) + 1
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
      if (bastst) write (6, 46) nsum
      write (9, 46) nsum
46    format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS IN POTENTIAL =',
     :        i3)
      nlam = nsum
      if (clist) then
        if (flagsu) then
          if (ihomo) then
            if (bastst)
     :      write (6,55) rmu * xmconv, brot, crot, delta, ipotsy, iop,
     :                   ered * econv, jtot, nu
            write (9,55) rmu * xmconv, brot, crot, delta, ipotsy, iop,
     :                   ered * econv, jtot, nu
55        format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=', f7.3,
     :           '  DELTA=',f7.3,/,'     POT-SYM=', i2,'  O/P=',i2,
     :           '  E=', f7.2, '       LBAR=', i5, '  NU=', i3)
          else
            write (6,60) rmu * xmconv, brot, crot, delta, ipotsy,
     :                   ered * econv, jtot, nu
            write (9,60) rmu * xmconv, brot, crot, delta, ipotsy,
     :                   ered * econv, jtot, nu
60        format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=', f7.3,
     :           '  DELTA=',f7.3,/,'     POT-SYM=', i2,
     :           '     E=', f7.2, '       LBAR=', i5, '  NU=', i3)
          end if
        else
          if (csflag) then
            if (ihomo) then
              if (bastst)
     :        write (6,65) rmu * xmconv, brot, crot, delta, ipotsy, iop,
     :                   ered * econv, jtot, nu
              write (9,65) rmu * xmconv, brot, crot, delta, ipotsy, iop,
     :                   ered * econv, jtot, nu
65          format(/,' **  CS SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,
     :           '  DELTA=',f7.3,/,'     POT-SYM=', i2,'  O/P=',i2,
     :           '  E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
            else
              if (bastst)
     :        write (6,75) rmu * xmconv, brot, crot, delta, ipotsy,
     :                   ered * econv, jtot, nu
              write (9,75) rmu * xmconv, brot, crot, delta, ipotsy,
     :                   ered * econv, jtot, nu
75          format(/,' **  CS SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,
     :           '  DELTA=',f7.3/,'     POT-SYM=', i2,
     :           '     E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
            end if
          else
            if (ihomo) then
              if (bastst)
     :        write (6,80) rmu * xmconv, brot, crot, delta, ipotsy, iop,
     :                     ered * econv, jtot, jlpar
              write (9,80) rmu * xmconv, brot, crot, delta, ipotsy, iop,
     :                     ered * econv, jtot, jlpar
80          format(/,' **  CC SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,
     :           '  DELTA=',f7.3/,'     POT-SYM=', i2,'  O/P=',i2,
     :           '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
            else
              if (bastst)
     :        write (6,85) rmu * xmconv, brot, crot, delta,
     :                     ered * econv, jtot, jlpar
              write (9,85) rmu * xmconv, brot, crot, delta,
     :                     ered * econv, jtot, jlpar
85          format(/,' **  CC SYMMETRIC TOP **',
     :      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,
     :           '  DELTA=',f7.3/,
     :           '     E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
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
*  levels with iop = -1 (called ortho levels for NH3) and
*  iop =2 (called para levels for NH3) have ki a multiple of 3
*  all other levels have iop = 1 (called para levels for NH3)
              if (iop.eq.-1 .and. mod(ki,ipotsy).ne.0) goto 105
              if (iop.eq.2 .and. mod(ki,ipotsy).ne.0) goto 105
              if (iop.eq.1 .and. mod(ki,ipotsy).eq.0) goto 105
            end if
*  check to make sure rotational energy is less than emax
            roteng = brot * ji * (ji + 1) + (crot - brot) * ki * ki
            if (roteng .gt. emax) go to 105
            nlist = nlist + 1
            jtemp(nlist) = ji
            ktemp(nlist) = ki
            ietmp(nlist) = 1 - 2*(iep - 1)
*  add inversion splitting to - inversion levels
            isym = -ietmp(nlist) * (-1) ** jtemp(nlist)
            if (iop.eq.2) isym = -isym
            if (isym .eq. -1) roteng = roteng + delta
            ehold(nlist)= roteng / econv
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
*  print this list if bastst = .true.
      if (bastst) then
        write (6, 130)
130     format (/,10x,
     :    'SORTED LEVEL LIST',/,'   N   J   K  EPS INV   EINT(CM-1)')
        do  140 i = 1, nlist
          isym = -ietmp(i) * (-1) ** jtemp(i)
          if (iop.eq.2) isym = -isym
          ecm = ehold(i) * econv
          write (6, 135) i, jtemp(i), ktemp(i), ietmp(i), isym, ecm
135       format (5i4,f10.3)
140     continue
      end if
*  now set up channel and level list for scattering calculation
      n = 0
      nlevel = 0
      do 170  njk = 1, nlist
        ji = jtemp(njk)
        ki = ktemp(njk)
        nlevel = nlevel + 1
        jhold(nlevel) = ji
*  ishold set equal to k times inversion symmetry
        ishold(nlevel) = -ietmp(njk) * (-1.d0)**ji * ki
*  originally, ishold set equal to k times epsilon
*        ishold(nlevel) = ietmp(njk) * ki
*
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
*  is set equal to k times inversion symmetry
            is(n) = -ietmp(njk) * (-1.d0)**ji * ki
*  originally, ishold set equal to k times epsilon
*            is(n) = ietmp(njk) * ki
            eint(n) = ehold(njk)
            l(n) = jtot
            cent(n) = jtot * (jtot + 1)
          end if
        else if (.not. csflag) then
*
*  here for cc calculations.  first calculate range of orbital angular
*  momentum quantum numbers allowed for this state
*
*  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem. Phys.
*  73, 2740 (1980) ]
            isym = -ietmp(njk) * (-1) ** ji
            if (iop.eq.2) isym = -isym
            ipar = isym * (-1) ** ki
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
                j(n) = ji
                k(n) = ki
                ieps(n) = ietmp(njk)
*  is set equal to k times inversion symmetry
                is(n) = -ietmp(njk) * (-1.d0)**ji * ki
*  originally, ishold set equal to k times epsilon
*                is(n) = ishold(nlevel)
                eint(n) = ehold(nlevel)
                l(n) = li
                cent(n) = li * (li + 1)
              end if
155         continue
          end if
160     continue
170   continue
*  also determine number of levels which are open
      nlevop = 0
      do 250  i = 1, nlevel
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        end if
250   continue
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
        end if
      end if
*  now list channels if requested
      if (clist) then
        if (.not.csflag) then
          if (bastst) write (6, 305)
          write (9,305)
305       format
     :     (/'   N   J   K  EPS INV  L    EINT(CM-1)')
        else
          if (bastst) write (6, 310) nu
          write (9,310) nu
310       format
     :     (/'   N   J   K  EPS INV  L    EINT(CM-1) ** NU = ',i2)
        end if
        do 330  i = 1, n
          if (k(i) .ne. 0) isym = is(i) / k(i)
          ecm = eint(i) * econv
          if (bastst)
     :      write (6, 320) i, j(i), k(i), ieps(i), isym, l(i), ecm
            write (9, 320) i, j(i), k(i), ieps(i), isym, l(i), ecm
320         format (6i4, f10.3)
330     continue
      end if
*  now calculate coupling matrix elements
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
        do 390 lb = lbmin, lammax(iterm)
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
                  if (bastst .and. iprint.ge.2) then
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
460     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS ',
     :           i6)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine vlmstp (jp, lp, j, l, jtot, kp, k, lambda, mu,
     :                   iepsp, ieps, v, csflag)
*  subroutine to calculate v-lambda matrices for close-coupled or coupled-stat
*  treatment of collisions of a symmetric top with an atom
*  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
*  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
*  the expressions for the full matrix elements are given in eqs. (46-48) and
*  (50) of the same article.
*  note, that in this article the bra indices (left side of matrix elements) u
*  are primed, while in the conventions of the present subroutine the bra
*  indices are primed and the ket indices (right side of matrix elements),
*  unprimed.
*
*  author:  millard alexander
*  current revision date:  27-mar-90
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    kp:       k quantum number of bra
*    k:        k quantum number of ket
*    iepsp:    symmetry index of bra
*    ieps:     symmetry index of ket
*    lambda:   order of legendre term in expansion of potential
*    mu:       absolute value of index of legendre term in expansion of
*              potential
*    v:        on return, contains desired matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index)
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum and l and lp correspond
*                to the orbital angular momenta
*  subroutines called:
*     xf3j, xf6j, prmstp
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag
      data one, zero / 1.0d0, 0.0d0 /
      v = zero
      iphase = ieps * iepsp * (-1) ** (jp + j + lambda + mu)
      if (iphase .eq. -1) return
      kdif = k - kp
      if (iabs(kdif) .eq. mu) then
        omeg = one
        if (csflag) then
*  the signed value of mu in the cs matrix elements is given by the
*  3j symbol (j' lambda j / -k' mu k) so that mu = k' - k = - kdif
          musign = - kdif
*  the multiplicative factor is given by Eq. (52) of S. Green, j. chem. phys.
*  64, 3463 (1976)
          if (kdif .gt. 0) omeg =  (-1) ** mu
        else if (.not. csflag) then
*  the signed value of mu in the cc matrix elements is given by the
*  3j symbol (j' j lambda / k' -k mu) so that mu = k - k'= kdif
          musign = kdif
*  the multiplicative factor is given by Eq. (48) of S. Green, j. chem. phys.
*  64, 3463 (1976)
          if (kdif .lt. 0)  omeg =  (-1) ** mu
        end if
*  contribution from (jp, kp, lp / Y(lambda, mu) / j, k, l), that is
*  the first term in Eq. (46) of S. Green, j. chem. phys. 64, 3463 (1976)
        call prmstp (jp, lp, j, l, jtot, kp, k, lambda, musign,
     :               vprm, csflag)
        v = v + omeg * vprm
      end if
      if (k + kp .eq. mu) then
*  n.b. for k = 0 and/or kp = 0, we recompute the same primitive matrix
*  element (here we follow MOLSCAT, although this might be somewhat inefficient
*  this is the second term in Eq. (46) of S. Green, j. chem. phys. 64, 3463
        if (.not.csflag) then
*  cc contribution from (jp, -kp, lp / Y(lambda, mu) / j, k, l)
          call prmstp (jp, lp, j, l, jtot, -kp, k, lambda, mu,
     :                 vprm, csflag)
          v = v + vprm * iepsp
        else if (csflag) then
*  cs contribution from (jp, kp, lp / Y(lambda, mu) / j, -k, l)
          call prmstp (jp, lp, j, l, jtot, kp, -k, lambda, mu,
     :                 vprm, csflag)
          v = v + vprm * ieps
        end if
      end if
      return
      end
* ----------------------------------------------------------------------
      subroutine prmstp (jp, lp, j, l, jtot, kp, k, lambda, mu,
     :                   vprm, csflag)
*  subroutine to calculate primitive v-lambda matrix elements for close-coupled
*  treatment of collisions of a symmetric top with an atom
*  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
*  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
*  note, that in this article the bra indices are unprimed and the ket indices
*  primed, while in the conventions of the present subroutine the bra indices
*  are primed and the ket indices, unprimed.
*  author:  millard alexander
*  current revision date:  27-mar-90
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    kp:       k quantum number of bra
*    k:        k quantum number of ket
*    lambda:   order of legendre term in expansion of potential
*    mu:       index of legendre term in expansion of potential
*    vrpm:     on return, contains primitive matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index)
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum
*  subroutines called:
*     xf3j, xf6j
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag
      data half, one, two, zero, four / 0.5d0, 1.d0, 2.d0, 0.d0, 4.d0/
      data pi /3.1415926536d0 /
      vprm = zero
      xjp = jp
      xj = j
      xkp = kp
      xk = k
      xjtot = jtot
      if (csflag) then
        nu = lp
        xnu = nu
      end if
      xlp = float(lp)
      xl = float(l)
      xlamda = float(lambda)
      xmu = float(mu)
      xnorm = (two * xjp + one) * (two * xj + one) *
     :        (two * lambda + one) / (four * pi)
*  special normalization for k and/or kp = 0
*  see Eq. (46) of S. Green, j. chem. phys. 64, 3463 (1976)
      if (k .eq. 0) xnorm = xnorm * half
      if (kp .eq. 0) xnorm = xnorm * half
*  the desired matrix element is constructed in x
      if (.not. csflag) then
*  here for cc matrix elements
        x = xf3j (xlp, xlamda, xl, zero, zero, zero)
        if (x .eq. zero) return
        x = x * xf3j (xjp, xj, xlamda, xkp, - xk, xmu)
        if (x .eq. zero) return
        x = x * xf6j (xj, xl, xjtot, xlp, xjp, xlamda)
        if (x .eq. zero) return
        iphase = jp + j + kp - jtot
        xnorm = xnorm * (two * lp + one) * (two * l + one)
      else if (csflag) then
*  here for cs matrix elements
        iphase = - k - nu
        x = xf3j (xjp, xlamda, xj, -xnu, zero, xnu)
        if (x .eq. zero) return
        x = x * xf3j (xjp, xlamda, xj, -xkp, xmu, xk)
      end if
      vprm = ( (-1) ** iphase) * x * sqrt(xnorm)
      return
      end
