* --------------------------------------------------------------------
      subroutine ba2sg (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, nrot, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      implicit double precision (a-h,o-z)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of a 2sigma molecule in a hund's case(a) basis with a
*  structureless atom or with an uncorrugated surface
*  author:  millard alexander
*  fixed is label to take care of +/- symmetry of the sigma state
*  current revision date: 9-dec-2011 by p.j.dagdigian
*  --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index eps of each channel
*              eps = +1 or -1
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
*    nrot:     n quantum number of each level
*    sc2-sc4:  scratch arrays (only sc2 is used here)
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*              in cc calculation jtot+1/2 is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true., then homonuclear molecule
*              if .false., then heteronuclear molecule
*              if the molecule is homonuclear (ihomo = .true.), either
*              only the s or the a rotational levels are included
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(parity+l-jtot)=jlpar
*              where parity is +1 for e levels and -1 for f levels
*              with the standard definition that e levels are those for
*              which eps*(-1)**(j-1/2-s) = 1
*              and f levels are those for which eps(-1)**(j-1/2-s) = -1
*              here s=0 for sigma-plus levels and s=1 for sigma-minus
*              levels
*              in cs calculation jlpar is set equal to 1 in calling program
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*    note!!!   if flaghf = .true., then the true values of the rotational
*    quantum numbers, the total angular momentum, and the coupled-states
*    projection index are equal to the values stored in j, jtot, and nu
*    plus 1/2
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brotsg:   rotational constant in cm-1
*    gsr:      2sigma state spin-rotation constant in cm-1
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables (integer
*              plus real)
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    nrmax:    the maximum case (b) rotational angular momenta
*    npar:     number of symmetry doublets included (npar=2 will ensure
*              both spin doublets)
*    isym:     if isym=+1, then the electronic symmetry is sigma-plus
*              if isym=-1, then the electronic symmetry is sigma-minus
*    igu:      if igu=+1, then the inversion symmetry is gerade
*              if igu=-1, then the inversion symmetry is ungerade
*    isa:      s/a symmetry index, if the molecule is homonuclear (ihomo=t)
*              then, if isa=+1 then only the s-levels are included in the
*              basis, if isa=-1, then only the a-levels are included
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*               zero of energy is taken to be n=0
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
*   vlm2sg:    returns angular coupling coefficient for particular
*              choice of channel index
* --------------------------------------------------------------------
      logical clist, csflag, flaghf, flagsu, ihomo, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysr/ isrcod, junkr, brotsg, gsr, drotsg, hrotsg
      common /cosysi/ nscode, isicod, nterm, nrmax, npar, isym, igu,
     :                isa
      common /coipar/ iiipar(9), iprint
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      dimension j(1), l(1), jhold(1), ehold(1), is(1), nrot(1),
     :          ieps(2), ishold(1), sc2(1)
      data ieps / -1, 1 /
      half = 0.5d0
      zero = 0.d0
      xjtot = jtot + half
      xnu = nu + half
*  check for consistency in the values of flaghf and csflag
      if (.not. flaghf) then
        write (6, 7)
        write (9, 7)
7       format (' *** FLAGHF = .FALSE. FOR DOUBLET SYSTEM; ABORT ***')
        stop
      end if
      if (flagsu .and. .not. csflag) then
        write (6, 8)
        write (8, 8)
8       format
     :   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
        stop
      end if
*  check that isym equals +1 or -1
      if (abs(isym).ne.1) then
        write (6, 112)
        write (6, 112)
112     format (' *** ISYM MUST EQUAL +1 OR -1; ABORT ***')
        call exit
      end if
*  check for consistency in values of nterm, lammin, lammax, mproj
      if (nterm .gt. 1) then
        write (6, 9) nterm
        write (9, 9) nterm
9       format (' *** NTERM=',i2,
     :   ' .GT. 1 FOR 2-SIGMA BASIS; ABORT ***')
        stop
      end if
      nsum = 0
      do 13  i = 1, nterm
        if (mproj(i) .gt. lammin(i) ) then
          write (6, 11) i, mproj(i), lammin(i)
          write (9, 11) i, mproj(i), lammin(i)
11        format (' *** MPROJ=',i2,' > LAMMIN=',i2,
     :            ' FOR TERM',i2,'; ABORT ***')
          stop
        end if
        if (ihomo) then
          if (mod(lammax(i)-lammin(i),2) .ne. 0) then
            write (6, 12) i, lammin(i), lammax(i)
            write (9, 12) i, lammin(i), lammax(i)
12          format
     :         (' *** IHOMO=T BUT ODD NO. OF TERMS FOR I=',i2,
     :        /,'     LAMMIN=',i2,' LAMMAX=',i2,'; ABORT ***')
            stop
          end if
            nsum = nsum + (lammax(i) - lammin(i) )/2+1
        else
          nsum = nsum + lammax(i) - lammin(i)+1
        end if
13    continue
      if (nlammx .lt. nsum) then
        write (6, 14) nsum, nlammx
        write (9, 14) nsum, nlammx
14      format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2,
     :          ' .GT. NLAMMX=', i2,'; ABORT')
        stop
      end if
      if (bastst) write (6, 15) nsum
      if (bastst) write (9, 15) nsum
      if (nsum.ne.nlam) write (6, 16) nsum, nlam
15    format (' ** TOTAL NUMBER OF ANISTROPIC TERMS IN POTENTIAL =',
     :        i3)
16    format (' ** TOTAL NUMBER OF ANISTROPIC TERMS IN POTENTIAL =',
     :        i3,' BUT NLAM=',i3,'  NLAM WILL BE REDEFINED')
      nlam = nsum
      if (clist) then
        if (flagsu) then
          if (ihomo) then
            write (9,20) rmu * xmconv, brotsg, gsr, isym, igu, isa,
     :                   ered * econv, jtot, xnu
            if (bastst)
     :      write (6,20) rmu * xmconv, brotsg, gsr, isym, igu, isa,
     :                   ered * econv, jtot, xnu
20          format(/,' **  2SIGMA UNCORRUGATED SURFACE **',
     :      '     RMU=', f9.4,'  BROT=', f7.3,'  G-SR=',g10.3,
     :      '  +/-=', i2,' g/u=',i2,
     :    '  s/a=', i2,/,'     E=', f7.2, '       LBAR=', i5,
     :    '  NU=', f5.1)
          else
            write (9,24) rmu * xmconv, brotsg, gsr, isym, ered * econv,
     :                 jtot, xnu
            if (bastst)
     :      write (6,24) rmu * xmconv, brotsg, gsr, isym, ered * econv,
     :                 jtot, xnu
24          format(/,' **  2SIGMA UNCORRUGATED SURFACE **',
     :      '     RMU=', f9.4,'  BROT=', f7.3,'  G-SR=',g10.3,
     :      '  +/-=', i2,
     :       /,'     E=', f7.2, '       LBAR=', i5,
     :    '  NU=', f5.1)
          end if
        else
          if (ihomo) then
            if (csflag) then
              write (9,25) rmu * xmconv, brotsg, gsr, isym, igu, isa,
     :                     ered * econv, jtot, xnu
              if (bastst)
     :        write (6,25) rmu * xmconv, brotsg, gsr, isym, igu, isa,
     :                     ered * econv, jtot, xnu
25            format(/,' ** CS 2SIGMA ** RMU=', f8.4,
     :             '  BROT=', f7.3,'  G-SR=',g10.3,
     :         ' +/-=', i2,' g/u=', i2,
     :        ' s/a=', i2,/,'     E=', f7.2,'  LBAR=', i5, 2x,
     :        '  NU=', f5.1)
            else
              if (bastst)
     :        write (6,30) rmu * xmconv, brotsg, gsr, isym, igu, isa,
     :                     ered * econv, xjtot, jlpar
              write (9,30) rmu * xmconv, brotsg, gsr, isym, igu, isa,
     :                     ered * econv, xjtot, jlpar
30            format(/,' **  CC 2SIGMA ** RMU=', f8.4,
     :             '  BROT=', f7.3,'  G-SR=',g10.3,
     :            ' +/-=', i2,' g/u=', i2,
     :     ' s/a=', i2,/,'     E=', f7.2,'  JTOT=', f5.1,
     :       2x,' JLPAR=', i2)
            end if
          else
            if (csflag) then
              write (9,31) rmu * xmconv, brotsg, gsr, isym,
     :                     ered * econv, jtot, xnu
              if (bastst)
     :        write (6,31) rmu * xmconv, brotsg, gsr, isym,
     :                     ered * econv, jtot, xnu
31            format(/,' **  CS 2SIGMA ** RMU=', f9.4,
     :             '  BROT=', f7.3,'  G-SR=',g10.3, '  +/-=', i2,
     :               /,'     E=', f7.2,'  LBAR=', i5, 2x,
     :     ' NU=', f5.1)
            else
              if (bastst)
     :        write (6,32) rmu * xmconv, brotsg, gsr, isym,
     :                     ered * econv, xjtot, jlpar
              write (9,32) rmu * xmconv, brotsg, gsr, isym,
     :                     ered * econv, xjtot, jlpar
32            format(/,' **  CC 2SIGMA ** RMU=', f9.4,
     :             '  BROT=', f7.3,'  G-SR=',g10.3, '  +/-=', i2,
     :        /,'     E=', f7.2,'  JTOT=', f5.1,
     :          2x,' JLPAR=', i2)
            end if
          end if
        end if
        if (.not. flagsu) write (9,35) rcut
35      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      end if
*  first set up list of all case (a) levels
      n = 0
*  the n=0, eps=+1, j=1/2 level is
*     s for sigma-g-plus or for sigma-u-minus
*     a for sigma-g-minus or for sigma-g-plus
      isplus = igu * isym * isa
      do 45 nnrot = 0, nrmax
        do 40 ip = 1, npar
*  include only the eps=+1 level for n=0
          if (nnrot .eq. 0 .and. ip .eq. 1) go to 40
*  now calculate j for each level
          if (ip .eq. 1) then
            ji = nnrot -1
          else
            ji = nnrot
          end if
*  actual half integer value of j is ji + 1/2
*  if homonuclear, include level only if allowed
*  see table i of alexander and corey, j. chem. phys. 84, 100 (1986)
          if (.not.ihomo .or. isa.eq.0. or.
     :       (ihomo .and. ieps(ip)*(-1)**ji.eq.isplus)) then
             n = n + 1
             is (n) = ieps(ip)
             j(n) = ji
             nrot(n) = nnrot
*  now assign energies for case (a) level and store in array eint
*    the matrix elements are given by a. j. kotlar, r. w. field,
*    and j. i. steinfeld, j. mol. spectr. 80, 86 (1980)
             x = float(j(n)) + 1.
             nn1 = x * (x - is(n))
             eint(n) = brotsg * nn1 - drotsg*nn1**2 + hrotsg*nn1**3
     :                       - half * (1 - is(n) * x) * gsr
*  next statement to take care of +/- symmetry of the sigma state
             is(n) = is(n) * isym
          end if
40      continue
45    continue
*  n now contains the number of levels
*  assume zero of energy is n=0
      emin = zero
c      emin = 1. e+7
c      do 65  i = 1, n
c        if (eint(i) .lt. emin) emin = eint(i)
c65    continue
*  form list of all energetically distinct rotational levels included in the
*  channel basis and their energies (with zero of energy set at lowest level)
      nlevel = 0
      do 70  i = 1, n
        nlevel = nlevel + 1
        ehold(nlevel) = (eint(i) - emin) / econv
        jhold(nlevel) = j(i)
        ishold(nlevel) = is(i)
70    continue
*  now sort this list to put closed levels at end
*  also determine number of levels which are open
      nlevop = 0
      do 80  i = 1, nlevel - 1
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
              go to 80
            end if
75        continue
        end if
80    continue
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
*  set up coupled-states channel basis (if desired)
      if (csflag) then
        nn = 0
        do 85  i = 1, n
*  include only channels with j at least equal to coupled states
*  projection index
          if (j(i) .ge. nu) then
            nn = nn + 1
            eint(nn) = (eint(i) - emin) / econv
            j(nn) = j(i)
            is(nn) = is(i)
            nrot(nn) = nrot(i)
            l(nn) = jtot
            if (.not. boundc) then
                 cent(nn) = jtot * (jtot + 1)
            else
                 xjtot=jtot+0.5d0
                 xj=j(nn)+0.5d0
                 xnu = nu+0.5d0
                 cent(nn)=xjtot*(xjtot+1)+xj*(xj+1)-2*xnu*xnu
            endif

          end if
85      continue
*  set number of coupled states channels
        n = nn
*  set up close-coupled channel basis (if desired)
      else if (.not. csflag) then
*  first move all indices of rotational levels to the top of the vectors
*  e.g. move j(n) to j(nmax), j(n-1) to j(nmax-1),  etc
        do 90  i = 1, n
*  move (n-i+1)th element to (nmax-i+1)th element
          inew = nmax - i + 1
          iold = n - i + 1
          eint(inew) = eint(iold)
          j(inew) = j(iold)
          is(inew) = is(iold)
          nrot(inew) = nrot(iold)
90      continue
        nn = 0
        do 100  i = 1, n
*  now take (nmax-n+i)th element, duplicate it as many times as is
*  required for rotational degeneray, and store the new elements in
*  nn, nn+1, ....
          ipoint = nmax - n + i
          ji = j(ipoint)
*  ipar is the symmetry of the molecular wavefunction
*  eps(-1)**(j-1/2) for sigma plus states
*  eps(-1)**(j-1/2-1) for sigma minus state
*  ipar=1 for e levels and -1 for f levels
          ipar = (-1) ** (ji + (isym - 1)/2 ) * is(ipoint)
          lmax = jtot + ji + 1
          lmin = iabs (jtot - ji)
          do 95  li = lmin, lmax
            ix = ipar * (-1) ** (li - jtot)
            if (ix .eq. jlpar) then
              nn = nn + 1
              eint(nn) = (eint(ipoint) - emin) / econv
              j(nn) = ji
              is(nn) = is(ipoint)
              nrot(nn) = nrot(ipoint)
              l(nn) = li
              cent(nn) = li * ( li + 1)
            end if
95        continue
100     continue
*  set number of close coupled channels
        n = nn
      end if
      if (n .gt. nmax) then
        write (6, 110) n, nmax
        write (9, 110) n, nmax
110     format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF',
     :         i4,' ABORT ***')
        stop
      end if
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions or for rcut < 0
      if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
        emin = 1. e+7
        do 120  i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is open asymptotically
            if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this condition is met
            end if
          end if
120     continue
*  now eliminate all channels with eint .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn = 0
          do 130 i = 1, n
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              nn = nn + 1
              eint(nn) = eint(i)
              j(nn) = j(i)
              is(nn) = is(i)
              nrot(nn) = nrot(i)
              cent(nn) = cent(i)
              l(nn) = l(i)
            end if
130       continue
*  reset number of channels
          n = nn
        end if
      end if
*  if no channels, return
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
          write (6, 160) nu, n, ntop
          write (9, 160) nu, n, ntop
160       format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3,
     :            '; ABORT **',/,
     :    '     CHECK RCUT')
          call exit
        endif
      end if
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
*  now list channels if requested
      if (clist) then
        if (.not.csflag) then
          if (bastst) write (6, 200)
          write (9,200)
200       format
     :     (/'   N   J  EPS   N    L   EINT(CM-1)')
        else
          if (bastst) write (6, 210) xnu
          write (9,210) xnu
210       format
     :     (/'   N   J  EPS   N    L   EINT(CM-1)',
     :       '    **  NU=', f4.1/)
        end if
        do 230  i = 1, n
          fj = j(i) + half
          ecm = eint(i) * econv
          if (bastst)
     :    write (6, 220) i, fj, is(i), nrot(i), l(i), ecm
          write (9, 220) i, fj, is(i), nrot(i), l(i), ecm
220       format (i4, f5.1, 2i4, i6, 1f10.3)
230     continue
      end if
*  now calculate coupling matrix elements
*  i counts v2 elements
*  inum counts v2 elements for given lambda
*  ilam counts numver of v2 matrices
*  ij is address of given v2 element in present v2 matrix
      i = 0
      if (bastst.and.iprint.ge.2) then
        write (6, 285)
        write (9, 285)
285     format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
      end if
      lamsum = 0
      do 400 ilam = 1, nlam
*  ilam is the angular expansion label
        lb = ilam
        inum = 0
        do 355  icol = 1, n
          do 350  irow = icol, n
            ij = ntop * (icol - 1) +irow
            lrow = l(irow)
            if (csflag) lrow = nu
*  always initialize potential to zero
            vee = zero
            call vlm2sg (j(irow), lrow, j(icol), l(icol), jtot,
     :                   lb, is(irow), is(icol), vee, csflag)
            if (vee .ne. zero) then
              i = i + 1
              if (i .le. nv2max) then
                inum = inum + 1
                v2(i) = vee
                iv2(i) = ij
                if (bastst.and.iprint.ge.2) then
                  write (6, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
                  write (9, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
290               format (i4, 2i7, 2i6, i6, g17.8)
                end if
              end if
            end if
350       continue
355     continue
        if (i .le. nv2max) lamnum(ilam) = inum
        if (bastst) then
          write (6, 370) ilam, lamnum(ilam)
          write (9, 370) ilam, lamnum(ilam)
370     format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
        end if
        lamsum = lamsum + lamnum(ilam)
400   continue
      if ( i.gt. nv2max) then
         write (6, 410) i, nv2max
         write (9, 410) i, nv2max
410      format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
         stop
      end if
      if (clist .and. bastst) then
        write (6, 420) lamsum
        write (9, 420) lamsum
420     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i4)
      end if
      if(bastst) write (9, 600) nlam, nlammx, (lamnum(i), i=1,nlam)
600   format (' nlam, nlammx, lamnum', 10i4)
      return
      end
*  -----------------------------------------------------------------------
      subroutine vlm2sg (jp, lp, j, l, jtot, lambda, iepsp, ieps,
     :                   v, csflag)
*  subroutine to calculate v-lambda matrices for close-coupled and
*  coupled-states treatments of collisions of a molecule in a 2sigma
*  electronic state
*  latest revision date:  21-feb-89
*  the cc matrix elements are given in eq. (15) of m.h. alexander,
*  j. chem. phys. 76, 3637 (1982)
*  the cs matrix elements are given in eq. (24) of m.h. alexander,
*  j. chem. phys. 76, 3637 (1982), with corrections corresponding
*  to eq. 14 of t. orlikowski and m.h. alexander, j. chem. phys. 79, 6006
*  (1983)
*  note that for cc collisions of a 2sigma molecule with a flat surface, the
*  coupling matrix elements are assumed to be identical to the cs matrix
*  elements here
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*              minus 1/2
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*              minus 1/2
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    lambda:   order of legendre term in expansion of potential
*    iepsp:    symmetry index of bra
*    ieps:     symmetry index of ket
*    v:        on return, contains matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index) minus 1/2
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum minus 1/2
*    for collisions of a 2sigma molecule with a surface, nu is equivalent
*    to m (the projection of j along the surface normal) minus 1/2
*  subroutines called:
*     xf3j, xf6j
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer jp, j, jtot, lp, l, lambda, iepsp, ieps, nu, iphase
      logical csflag
      half = 0.5d0
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      v = zero
      xjp = jp + half
      xj = j + half
      xjtot = jtot + half
      if (csflag) then
        nu = lp
        xnu = nu + half
      end if
      xlp = float(lp)
      xl = float(l)
      xlamda = float(lambda)
      iphase = ieps * iepsp * ((-1) ** (jp + j + lambda + 1))
      if (iphase .eq. 1) return
      if (csflag) then
*  note true phase factor is nu - omega, where nu is a half/integer and
*  omega=1/2.  here, however nu is nu-true + 1/2, so nu-true - omega = nu
        iphase = nu
        xnorm = (two * xjp + one) * (two * xj + one)
        xnu = nu + half
        xnum = - xnu
*  indices in denominator of 3j correspond
*  to eq. 14 of t. orlikowski and m.h. alexander, j. chem. phys. 79, 6006
*  (1983) rather than eq. (24) of m.h. alexander,
*  j. chem. phys. 76, 3637 (1982)
        x = xf3j (xjp, xlamda, xj, xnum, zero, xnu)
      else
        iphase = jp + j + 1 + jtot
        xnorm = (two * xjp + one) * (two * xj + one) * (two * lp + one)
     :        * (two * l + one)
        x = xf3j (xlp, xlamda, xl, zero, zero, zero)
        if  (x .eq. zero) return
        x = x * xf6j (xjp, xlp, xjtot, xl, xj, xlamda)
      end if
      if (x .eq. zero) return
      x = x * xf3j (xjp, xlamda, xj, -half, zero, half)
      iphase = (-1) ** iabs(iphase)
      v = iphase * x * sqrt(xnorm)
      return
      end
*  -----------------------------------------------------------------------
