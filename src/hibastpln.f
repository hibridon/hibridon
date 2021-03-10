       subroutine bastpln(j, l, is, jhold, ehold, ishold, nlevel,
     :                   nlevop, ktemp, jtemp,
     :                  ieps, jt12, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin,
     :                  jlpar, n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of a rigid symetric top with a rigid diatomic  molecule
*  authors:  millard alexander and peter vohralik
*  current revision date:  29-nov-1995 by c.r.
*  current revision date:  7-apr-2003 by mha
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains combined rotational quantum numbers for each
*              channel.  in terms of the rotational quantum numbers of each
*              molecule we have:  j = 10 * j1 + j2
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry of inversion vibration for
*              each channel
*  Note that we have adopted the following convention for the symmetry
*  index "is" so that on return the symmetric top molecular levels can
*  be uniquely identified by the two labels "j" and "is":
*           is = +/- k
*       where k is the projection quantum number
*       the plus sign refers to the "plus" inversion levels and the
*       minus sign to the "minus" inversion levels
*       for k=0, where only eps=+1 states exist, the index "is" is equal
*       to zero.
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level, irrespective of projection degeneracy and the
*              degeneracy associated with different values of j12
*              note that jhold = 10 * j1hold + j2hold
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each energetically
*              distinct level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    ktemp:    scratch array used to create channel list
*    ieps:     on return contains epsilon label for each channel
*    jt12:      on return contains j12 = j1 + j2 for each channel
*    jtemp:    scratch array used to create channel list.
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
*              if ihomo.true., then symmetric top posseses interchange
*              symmetry (e.g. NH3),
*              so that only the ortho or para levels will be
*              included depending on the value of parameters iop
*              in common /cosysi/ (see below)
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              jlpar (by definition jlpar=(-1)**(j2+l-jtot)*parity
*              jtot is added to be consistent with bastp routine but is
*              not needed here)
*              where parity is for the symmetric-top molecule
*              parity=isym*(-1)*k=ieps*(-1)**(k+j1)
*
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*  variables in common block /cosysr/
*    isrcod:   number of real system dependent variables
*    brot:     A=B rotational constant for prolate top
*    crot:     C rotational constant for prolate top
*    drot:     rotational constant for linear molecule
*  variables in common block /cosysi/
*    nscod:    total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   numbr of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    ipotsy:   cylindrical symmetry of potential.  If ihomo = .true.,
*              only
*              terms with mu1 equal to an integ multiple of ipotsy
*              can be included in the potential.  Example:  for NH3,
*              ipotsy = 3
*    iop:      ortho/para label for molecular states. If ihomo=.true.
*              then only para states will be included if iop=1
*              and only ortho states if iop=-1
*    ninv:     number of inversion doublets included
*              if ninv = +1, only + inversion levels included
*              if ninv = -1, only - inversion levels included
*              if ninv = 2, both inversion levels included
*    kmax:     the maximum projection quantum number included
*    jmax0:    the maximum rotational angular momenta for the k=0 stack
*    jmax1:    the maximum rotational angular momenta for the k=1 stack
*    jmax2:    the maximum rotational angular momenta for the k=2 stack
*    jmax3:    the maximum rotational angular momenta for the k=3 stack
*    jmax4:    the maximum rotational angular momenta for the k=4 stack
*    jmax5:    the maximum rotational angular momenta for the k=5 stack
*    jmax6:    the maximum rotational angular momenta for the k=6 stack
*    ipotsy2:  symmetry of potential. if linear molecule is homonuclear
*              then ipotsy=2 and only terms with lambda2  even can be
*              included in the potential,else ipotsy=1.
*    j2max:    the maximum rotational angular momentum for linear
*              molecule
*    j2min:    the minimum rotational angular momentum for linear
*              molecule
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
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix
*               elements stored in packed column form that is (1,1)
*                (2,1), (3,1), (n,1), (2,2), (3,2) ... (n,2), etc.
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
*  subroutines called:
*   vlmstpln:  returns molecule-molecule angular coupling coefficient for
*              particular choice of channel index
*  variable in common block /co2mol/
*   twomol     if .true. collision between symmetric top and linear
*              molecule, if .false. collision symmetric top-atom.
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst, twomol
      include "common/parbas"
      include "common/parbasl"
      common /coipar/ iiipar(9), iprint
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop, ninv,
     :                kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :                jmax6, ipotsy2, j2max, j2min
      common /cosysr/ isrcod, junkr, brot, crot, drot
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      common /co2mol/ twomol
      dimension j(1), l(1), jhold(1), ehold(1), is(1),
     :          ishold(1),
     :          nlami(45), jmx(0:6)
      dimension ieps(1), ktemp(1), jtemp(1), jt12(1)
      data ione, itwo, ithree / 1, 2, 3 /
      data one, zero, two / 1.0d0, 0.0d0, 2.d0 /
      twomol=.true.
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (9, 5)
5       format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***' )
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
        write (6, 15)
        write (9, 15)
15      format (' *** CSFLAG=.TRUE.; NOT IMPLEMENTED IN BASIS;',
     :          ' ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      jmx(0) = jmax0
      jmx(1) = jmax1
      jmx(2) = jmax2
      jmx(3) = jmax3
      jmx(4) = jmax4
      jmx(5) = jmax5
      jmx(6) = jmax6
      do  19  ki = 0, kmax
        if (ki .gt. jmx(ki) ) then
          write (6, 18) ki, jmx(ki)
          write (9, 18) ki, jmx(ki)
18        format (' *** K =',i3,' .GT. JMAX(K) =',i3,'; ABORT ***')
          stop
        end if
19    continue
      do 25  i = 1, nterm
        if (ihomo .and. mod(mproj(i), ipotsy) .ne. 0 ) then
          write (9, 20) i, mproj(i), ipotsy
          write (6, 20) i, mproj(i), ipotsy
20        format(/' *** MPROJ(',i2,')=',i2,
     :     ' .NE. INTEGER MULTIPLE OF', i2, ' FOR IHOMO=T; ABORT ***')
          stop
        end if
        if (twomol) then
          if (mod(lam2(i), ipotsy2) .ne. 0 ) then
            write (9, 21) i, lam2(i), ipotsy2
            write (6, 21) i, lam2(i), ipotsy2
21          format(/' *** LAM2(',i2,')=',i2,
     :      ' .NE. INTEGER MULTIPLE OF', i2,
     :         ' FOR HOMONUCLEAR MOLECULE ABORT ***')
            stop
          endif
        endif
25    continue

      nsum = 0
      do 35  i = 1, nterm
        if (mproj(i) .gt. lammin(i) ) then
          write (6, 30) i, mproj(i), lammin(i)
          write (9, 30) i, mproj(i), lammin(i)
30        format (' *** MPROJ=',i2,' > LAMMIN=',i2,
     :            ' FOR TERM',i2,'; ABORT ***')
          stop
        end if
        if (twomol) then
          if (m2proj(i).gt.lammin(i)) then
            write (6, 31) i, m2proj(i), lammin(i)
            write (9, 31) i, m2proj(i), lammin(i)
31          format (' *** M2PROJ=',i2,' > LAMMIN=',i2,
     :            ' FOR TERM',i2,'; ABORT ***')
            stop
          endif
          if (m2proj(i).gt.lam2(i) ) then
            write (6, 32) i, m2proj(i), lam2(i)
            write (9, 32) i, m2proj(i), lam2(i)
32          format (' *** M2PROJ=',i2,' > LAM2=',i2,
     :            ' FOR TERM',i2,'; ABORT ***')
            stop
          end if
        end if
        if(lammax(i) .lt. 0) goto 35
        nlami(i) = lammax(i) - lammin(i) + 1
        nsum = nsum + nlami(i)
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
46    format (' ** TOTAL NUMBER OF ANISTROPIC TERMS IN POTENTIAL =',
     :        i3)
      nlam = nsum
      if (bastst) then
        write (9, 50) nlam
        write (6, 50) nlam
50      format (' *** TOTAL NUMBER OF ANISTROPIC TERMS =',i3)
      end if
      if (clist) then
        if (flagsu) then
          write(6,51)
          write(9,51)
51        format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **',
     :             ' NON IMPLEMENTED HERE, ABORT**')
          stop
        elseif (csflag) then
          write(6,52)
          write(9,52)
52        format(/,' **  SYMMETRIC TOP - COUPLED STATES **',
     :             ' NON IMPLEMENTED HERE, ABORT**')
          stop
        elseif (twomol) then
          if (ihomo) then
            if (bastst)
     :      write (6,80) rmu * xmconv, brot, crot, drot, ipotsy, iop,
     :                     ipotsy2, ered * econv, jtot, jlpar
            write (9,80) rmu * xmconv, brot, crot, drot, ipotsy, iop,
     :                     ipotsy2, ered * econv, jtot, jlpar
80          format(/,' **  CC SYMMETRIC TOP-LINEAR MOLECULE **',
     :          /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,
     :            '  DROT=',f7.3,/,
     :            '     POT-SYM=', i2,'  O/P=',i2, ' HOMOLIN=', i2,
     :             '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
          else
            if (bastst)
     :      write (6,85) rmu * xmconv, brot, crot, drot,
     :                    ipotsy2, ered * econv, jtot, jlpar
            write (9,85) rmu * xmconv, brot, crot, drot,
     :                    ipotsy2, ered * econv, jtot, jlpar
85          format(/,' **  CC SYMMETRIC TOP **',
     :          /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3,
     :             '  DROT=',f7.3,/,'  HOMOLIN=', i2,
     :             '     E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
          end if
        else
           write(6,86)
           write(9,86)
86         format(/,' **  SYMMETRIC TOP - ATOM **',
     :            ' NON IMPLEMENTED HERE, ABORT**')
           stop
        end if
      end if
      if (.not. flagsu) write (9,90) rcut
90      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
*  first set up list of all j,k j2 states included in basis

      nlist = 0
      do 105  ki = 0, kmax
        do 105  ji = ki, jmx(ki)
          if (twomol) then
            do 103 ji2 = j2min, j2max, 2
              nlist = nlist+1
              ehold(nlist)= brot * ji * (ji + 1) + (crot - brot)*ki*ki
     :                      + drot * ji2 * ( ji2 + 1)
              jtemp(nlist) = 10*ji + ji2
              ktemp(nlist) = ki
103         continue
          else
            nlist = nlist + 1
            ehold(nlist)= brot * ji * (ji + 1) + (crot - brot)*ki*ki
            jtemp(nlist) = ji
            ktemp(nlist) = ki
          end if
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
            end if
115       continue
120     continue
      end if
*  print this list if bastst = .true.
      if (bastst) then
        if(twomol) then
          write (6, 125)
125       format (/,10x,
     :    'SORTED LEVEL LIST',/,'   N   J   K   J2   EINT(CM-1)')
          do  127 i = 1, nlist
            ji2=mod(jtemp(i),10)
            ji1=jtemp(i)/10
            write (6, 126) i, ji1, ktemp(i), ji2, ehold(i)
126         format (4i4,f10.3)
127       continue
        else
          write (6, 130)
130       format (/,10x,
     :    'SORTED LEVEL LIST',/,'   N   J   K    EINT(CM-1)')
          do  140 i = 1, nlist
            write (6, 135) i, jtemp(i), ktemp(i), ehold(i)
135         format (3i4,f10.3)
140       continue
        end if
      end if
*  now set up channel and level list
      n = 0
      nlevel = 0
      do 170  njk = 1, nlist
        ki = ktemp(njk)
        ji = jtemp(njk)
        if (twomol) then
          ji=jtemp(njk)/10
          ji2=mod(jtemp(njk),10)
        else
          ji2=0
        endif
        if (ihomo) then
*  check to see if ki corresponds to an allowed k stack if interchange
*  symmetry, this if statement is specifically for NH3
          if (iop .eq. -1 .and. mod(ki, ipotsy) .ne. 0) go to 170
          if (iop .eq. 1 .and. mod(ki, ipotsy) .eq. 0) go to 170
        end if
        iepmin = -1
        if (ki .eq. 0) iepmin = 1
*  now sum over inversion states except for k = 0 stack, where only
*  epsilon = + 1 exists
        do 160 iepsil = iepmin, 1, 2
*  determine the product of epsilon times (-1)^j.  this is equal to
*  the inversion symmetry of the level [Eq. (A20) of S. Green, J. Chem.
*  Phys. 73, 2740 (1980) ]
          isym = iepsil * (-1) ** ji
*  if only one inversion state is included in the basis, exclude this
*  state if it has the wrong inversion symmetry
          if (ninv .eq. 1 .and. isym .ne. 1) go to 160
          if (ninv .eq. -1 .and. isym .ne. -1) go to 160
*  here if this state is to be included
          nlevel = nlevel + 1
          if (twomol) then
            ehold(nlevel)=
     :       (brot * ji * (ji + 1) + (crot - brot) * ki * ki +
     :        drot * ji2* (ji2 + 1 )) / econv
            jhold(nlevel) = 10*ji + ji2
          else
            ehold(nlevel)=
     :       (brot * ji * (ji + 1) + (crot - brot) * ki * ki) / econv
            jhold(nlevel) = ji
          end if
*  the vector ihold will contain k times the inversion symmetry of the
*  state
          ishold(nlevel) = ki * isym
          if (csflag) then
*  here for cs calculations; include state only if j at least equal to
*  coupled states projection index
            if (ji .ge. nu) then
              n = n + 1
              if (n .gt. nmax) then
                write (6, 150) n, nmax
                write (9, 150) n, nmax
150             format(/' *** NCHANNELS=', i4,
     :             ' .GT. MAX DIMENSION OF',i4,' ABORT ***')
                stop
              end if
              is(n) = ishold(nlevel)
              ieps(n) = iepsil
              j(n) = ji
              eint(n) = ehold(nlevel)
              l(n) = jtot
              cent(n) = jtot * (jtot + 1)
            end if
            if (twomol) then
              write(6,151)
              write(9,151)
151           format(/,'*** COUPLED STATES NOT IMPLEMENTED FOR TWO'
     :               ,' MOLECULES IN BASIS***',/,'*** ABORT***')
              stop
            end if
          else if (.not. csflag) then
*  here for cc calculations.  first calculate range of orbital angular
*  momentum quantum numbers allowed for this state
*  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem.
*  Phys.
*  73, 2740 (1980) ]
            ipar = isym * (-1) ** ki
            j12max=ji + ji2
            j12min=iabs(ji-ji2)
            do 156 ji12=j12min, j12max
              lmax = jtot + ji12
              lmin = iabs (jtot - ji12)
              do 155  li = lmin, lmax
                lpar = (-1) ** (li - jtot)
                j2par = (-1)**ji2
*  check to see if this channel has the desired parity, if so, include
*  it
                if (ipar * lpar * j2par .eq. jlpar) then
                  n = n + 1
                  if (n .gt. nmax) then
                    write (6, 150) n, nmax
                    write (9, 150) n, nmax
                    stop
                  end if
                  is(n) = ishold(nlevel)
                  ieps(n) = iepsil
                  if (twomol) then
                    j(n) = 10*ji + ji2
                    jt12(n)=ji12
                  else
                    j(n) = ji
                  endif
                  eint(n) = ehold(nlevel)
                  l(n) = li
                  cent(n) = li * (li + 1)
                end if
155           continue
156         continue
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
* here if channel is open asymptotically
            if ( jtot * (jtot + 1) / (two * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this condition
*  is met
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
              if (twomol) then
                jt12(nn)=jt12(i)
              endif
              cent(nn) = cent(i)
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
     :     (/'   N   J  EPS INV  K  J2  J12  L    EINT(CM-1)')
        else
          if (bastst) write (6, 310) nu
          write (9,310) nu
310       format
     :     (/'   N   J  EPS INV  K   L    EINT(CM-1) ** NU = ',i2)
        end if
        do 330  i = 1, n
          ki = iabs(is(i))
          if (is(i) .ne. 0) inv = is(i) / ki
          ecm = eint(i) * econv
          if (twomol) then
            ji = j(i)/10
            ji2 = mod(j(i),10)
          else
            ji = j(i)
            ji2 = 0
          endif
          if (bastst)
     :      write (6, 320) i, ji, ieps(i), inv, ki, ji2, jt12(i),
     :                     l(i), ecm
            write (9, 320) i, ji, ieps(i), inv, ki, ji2, jt12(i),
     :                     l(i), ecm
320         format (8i4, f10.3)
330     continue
      end if

*  now calculate coupling matrix elements
*  i counts v2 elements
*  inum counts v2 elements for given lambda
*  ilam counts number of v2 matrices
*  ij is address of given v2 element in present v2 matrix
      i = 0
      if (bastst.and. iprint.ge. 2) then
        write (6, 340)
*        write (76,340)
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
          if (twomol) then
            lb2 = lam2(iterm)
            mu2 = m2proj( iterm)
          endif
          inum = 0
          ij = 0
          do 355  icol = 1, n
            do 350  irow = icol, n
              ij = ntop * (icol - 1) +irow
              lr = l(irow)
              lc = l(icol)
              kc = iabs (is(icol))
              kr = iabs (is(irow))

              if (csflag) lrow = nu
              if (twomol) then
                jc = j(icol)/10
                jr = j(irow)/10
                j2c = mod(j(icol), 10)
                j2r = mod(j(irow), 10)
                kc = iabs (is(icol))
                kr = iabs (is(irow))
                j12r = jt12(irow)
                j12c = jt12(icol)
* this call to vlmstp is a test for sym-atom limit
*                call vlmstp (jr, lr, jc, lc, jtot,
*     :                     kr, kc, lb, mu, ieps(irow),
*     :                     ieps(icol), vee, csflag)

                call vlmstpln (jr, lr, jc, lc, j2r, j2c,
     :                     j12r, j12c, jtot, kr, kc, lb, mu, lb2, mu2,
     :                     ieps(irow), ieps(icol), vee, csflag)

              else
                call vlmstp (j(irow), lr, j(icol), lc, jtot,
     :                     kr, kc, lb, mu, ieps(irow),
     :                     ieps(icol), vee, csflag)
              endif
              if (vee .ne. zero) then
                i = i + 1
                if (i .le. nv2max) then
                  inum = inum + 1
                  v2(i) = vee
                  iv2(i) = ij
                  if (bastst.and. iprint.ge.2) then
                    write (6, 345) ilam, lb, icol, irow, i, iv2(i),
     :                             vee
*                    write (76, 345) ilam, lb, icol, irow, i, iv2(i),
*     :                             vee
                    write (9, 345) ilam, lb, icol, irow, i, iv2(i),
     :                             vee
345                 format (i4, 2i7, 2i6, i6, f17.10)
                  end if
                end if
              end if
350         continue
355       continue
          if (i .le. nv2max) lamnum(ilam) = inum
          if (bastst) then
            write (6, 370) ilam, lamnum(ilam)
*            write (76,370) ilam, lamnum(ilam)
            write (9, 370) ilam, lamnum(ilam)
370         format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
          end if
          lamsum = lamsum + lamnum(ilam)
390     continue
400   continue
      if ( i.gt. nv2max) then
         write (6, 450) i, nv2max
*         write(76, 450) i, nv2max
         write (9, 450) i, nv2max
450      format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
         stop
      end if
      if (clist .and. bastst) then
        write (6, 460) lamsum
*        write(76, 460) lamsum
        write (9, 460) lamsum
460     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i4)
      end if
      return
      end

* --------------------------------------------------------------------
      subroutine vlmstpln (jp, lp, j, l, j2p, j2,
     :                     j12p, j12, jtot, kp, k, lambda, mu,
     :                     lambda2, mu2,
     :                     iepsp, ieps, v, csflag)

*  subroutine to calculate v-lambda matrices for close-coupled or
*  coupled-state treatment of collisions of a symmetric top with an atom
*  the primitive cc  elements are given in claire's thesis
*  note, that in this article the bra indices (left side of matrix
*  elements) u
*  are primed, while in the conventions of the present subroutine the
*  bra indices are primed and the ket indices (righ side of matrix
*  elements), unprimed.
*
*  author: claire rist
*  current revision date: 17-jan-1992
*
* -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element
*              (bra) for symmetric top
*    lp:       orbital angular momentum of left side of matrix element
*              (bra)
*    j:        rotational quantum number of right side of matrix element
*              (ket) for symmetric top
*    l:        orbital angular momentum of right side of matrix element
*              (ket)
*    j2p:      rotational quantum number of left side of matrix element
*              (bra) for linear molecule
*    j2:       rotational quantum number of left side of matrix element
*              (ket) for linear molecule
*    j12p:     rotational quantum number of left side of matrix element
*              (bra) j12p = jp + j2p
*    j12:      rotational quantum number of left side of matrix element
*              (ket) j12 = j + j2
*    jtot:     total angular momentum
*    kp:       k quantum number of bra
*    k:        k quantum number of ket
*    iepsp:    symmetry index of bra
*    ieps:     symmetry index of ket
*    lambda:   order of rotation matrix term in expansion of potential
*    mu:       absolute value of index of rotation matrix term in expansion of
*              expansion of potential
*    lambda2:  order of rotation matrix term in expansion of potential
*    mu2:      absolute value of index of rotation matrix term in
*              expansion of potential
*    v:        on return, contains desired matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index)
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum and l and lp correspond
*                to the orbital angular momenta
*  subroutines called:
*     xf3j, xf6j, xf9j, prmstp
*     f9j, cf9j (Claire's routine)
*
* -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag
      data one, zero, two  / 1.0d0, 0.0d0, 2.0d0 /
      v = zero
      iphase = ieps * iepsp * (-1) ** (jp + j + lambda + mu)
      if (iphase .eq. -1) return
      kdif = k - kp
      if (iabs(kdif) .eq. mu) then
        omeg = one
*  the signed value of mu in the cc matrix elements is given by the
*  3j symbol (j' lambda j/ k' mu -k) so that mu = k - k'= kdif
        musign = kdif
        if (kdif .lt. 0)  omeg =  (-1) ** mu
*  contribution from these 3j coeff (j' lambda j/k' mu -k)...
        call prmstpln (jp, lp, j2p, j12p, j, l, j2, j12,
     :               jtot, kp, k, lambda, musign, lambda2, mu2,
     :               vprm, csflag)
        v = v + omeg * vprm
      end if
      if ((k + kp) .eq. mu) then
*  n.b. for k = 0 and/or kp = 0, we recompute the same primitive matrix
*  element
*  cc contribution from  3j coeff (j' lambda j/ -k' mu, -k)
        call prmstpln (jp, lp, j2p, j12p, j, l, j2, j12,
     :               jtot, -kp, k, lambda, mu, lambda2, mu2,
     :                 vprm, csflag)
        v = v + vprm * iepsp
      end if
      s4pir8 = sqrt ( 4.d0 * acos(-1.d0) )
      s4pi = s4pir8
      return
      end
* ----------------------------------------------------------------------
      subroutine prmstpln (jp, lp, j2p, j12p, j, l, j2, j12,
     :               jtot, kp, k, lambda1, mu1,lambda2, mu2,
     :                   vprm, csflag)
*  subroutine to calculate primitive v-lambda matrix elements for
*  close-couple
*  treatment of collisions of a symmetric top with a linear molecule
*  author: claire rist
*  current revision date:  17-jan-1992
*
* -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element
*              (bra) for symmetric top
*    lp:       orbital angular momentum of left side of matrix element
*              (bra)
*    j2p:      rotational quantum number of left side of matrix element
*              (bra) for linear molecule
*    j12p:     total molecular rotational quantum number of left side of
*              matrix element (bra)
*    j:        rotational quantum number of right side of matrix element
*              (ket) for symmetric top
*    l:        orbital angular momentum of right side of matrix element
*              (ket)
*    j2:       rotational quantum number of right side of matrix element
*              (ket) for linear molecule
*    j12:      total molecular rotational quantum number of right side
*              of matrix element
*    jtot:     total angular momentum
*    kp:       k quantum number of bra
*    k:        k quantum number of ket
*    lambda1:  order of rotation matrix term in expansion of potential
*    mu1:      index of rotation matrix term in expansion of potential
*    lambda2:  order of rotation matrix term in expansion of potential
*    mu2:      index of rotation matrix term in expansion of potential
*    vrpmstp:  on return, contains primitive matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index)
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum
*  subroutines called:
*     xf3j, xf6j, xf9j
*     f9j, cf9j (Claire's routines)
*
* -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag
      data half, one, two, zero, four / 0.5d0, 1.d0, 2.d0, 0.d0, 4.d0/
      data pi /3.1415926536d0 /
      vprm = zero
      xjp = dble(jp)
      xj = dble(j)
      xkp = dble(kp)
      xk = dble(k)
      xj2p = dble(j2p)
      xj2 = dble(j2)
      xj12p = dble(j12p)
      xj12 = dble(j12)
      xjtot = dble(jtot)
      xlp = dble(lp)
      xl = dble(l)
      xlambda1 = dble(lambda1)
      xmu1 = dble(mu1)
      xlambda2 = dble(lambda2)
      xmu2 = dble(mu2)
* Last term in j12 added by Claire to reproduce the isotropic
* matrix elements
      xnorm = (two * xjp + one) * (two * xj + one) *
     :        (two * xj2p + one) * (two * xj2 + one) *
     :        (two * xlp + one) * (two * xl + one)*
     :        (two * j12p + one) * (two * j12 + one)
*  special normalization for k and/or kp = 0
      if (k .eq. 0) xnorm = xnorm * half
      if (kp .eq. 0) xnorm = xnorm * half
*  special normalization for mu1 and/or mu2 = 0
*      if (mu1 .eq. 0) xnorm = xnorm * half * half
      if (mu2 .eq. 0) xnorm = xnorm * half * half
*  the desired matrix element is constructed in x
*  here for cc matrix elements
      x = xf3j (xjp, xlambda1, xj, xkp, xmu1, - xk)
      x2 = xf3j (xj2p, xlambda2, xj2, zero, zero, zero)
      x = x * xf3j (xj2p, xlambda2, xj2, zero, zero, zero)
      if (x .eq. zero) return
      lmin = max(iabs(lambda1-lambda2), iabs(l-lp))
      lmax = min(lambda1+lambda2, l+lp)
      sum = zero
      do 500 lambda = lmin, lmax
        xlambda = dble(lambda)
        xs = xf3j (xlp, xlambda, xl, zero, zero, zero)
        if (xs .eq. zero) go to 500
        xs = xs * xf3j (xlambda1,xlambda2,xlambda,xmu2,-xmu2,zero)
        if (xs .eq. zero) go to 500
        xs = xs * xf6j (xj12, xl, xjtot, xlp, xj12p, xlambda)
        if (xs .eq. zero) go to 500
*        xs4 = 1/sqrt((2*xjp+1)*(2*xj+1)*(2*xlambda+1))
*        xs = xs * xs4
        xs = xs * f9j(lambda1, j, jp, lambda2, j2, j2p,
     :       lambda, j12, j12p)
        if (xs .eq. zero) go to 500
        xs = xs * (two * xlambda + one)
*        write(6,*)'lambda,  xs', lambda, xs
        sum = sum + xs
500   continue
      x = x * sum
*      iphase = jtot+k+mu2
      iphase = jtot + k + mu2 + jp + j2p + j12p
      vprm = ( (-1) ** iphase) * x * sqrt(xnorm)
* Claire factor ( 1+ eps*epsp*(-1)**(j+jp+lambda1+mu1))=2
      vprm = two * vprm
c     provisoire: check with nh3h2(j=0) calculations
c     multiplicationfactor:
*      s4pir8 = sqrt ( 4.d0 * acos(-1.d0) )
*      s4pi = s4pir8
*      vprm = vprm *sqrt(2*xlambda1+1) /s4pi
*      cnorm= two*sqrt(2*xlambda1+1)*sqrt(xnorm)/s4pi
      return
      end


      function f9j(i1,i2,i3,i4,i5,i6,i7,i8,i9)

      implicit real*8 (a-h,o-z)
      data zero /0.d0/
      kmin=max0(iabs(i6-i2),iabs(i8-i4))
      kmin=max0(kmin,iabs(i9-i1))
      kmax=min0((i6+i2),(i8+i4))
      kmax=min0(kmax,(i9+i1))
      sum=zero
      f9j=zero
      do 20 k=kmin,kmax
        b2=f6j(i1,i2,i3,i6,i9,k)
        if(b2.eq.zero) go to 20
        b3=f6j(i4,i7,i1,i9,k,i8)
        if(b3.eq.zero) go to 20
        b4=f6j(i5,i8,i2,k,i6,i4)
        if(b4.eq.zero) go  to 20
        sum=sum+(2*k+1)*b2*b3*b4
 20        continue
      f9j=sum
      return
      end


      function cf9j(x1,x2,x3,x4,x5,x6,x7,x8,x9)

      implicit real*8 (a-h,o-z)
      data two, zero, tol /2.d0, 0.d0, 1.d-15/
      i1 = nint(two*x1)
      i2 = nint(two*x2)
      i3 = nint(two*x3)
      i4 = nint(two*x4)
      i5 = nint(two*x5)
      i6 = nint(two*x6)
      i7 = nint(two*x7)
      i8 = nint(two*x8)
      i9 = nint(two*x9)
      xkmin=max(abs(x6-x2),abs(x8-x4))
      xkmin=max(xkmin,abs(x9-x1))
      xkmax=min((x6+x2),(x8+x4))
      xkmax=min(xkmax,(x9+x1))
      kmin = nint(xkmin*two)
      kmax = nint(xkmax*two)
      sum=zero
      cf9j=zero
      do 30 k=kmin,kmax,2
        xk = dble(k)/two
        b2=xf6j(x1,x2,x3,x6,x9,xk)
        if(b2.eq. zero) go to 30
        b3=xf6j(x4,x7,x1,x9,xk,x8)
        if(b3.eq. zero) go to 30
        b4=xf6j(x5,x8,x2,xk,x6,x4)
        if(b4.eq.zero) go  to 30
        sum=sum+(2*xk+1)*b2*b3*b4
 30        continue
      cf9j=sum
      return
      end







