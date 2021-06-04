* ----------------------------------------------------------------------
      subroutine basgpi1(j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, ivhold, c12, c32, csig,
     :                  rcut, jtot,flaghf, flagsu, csflag, clist,
     :                  bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
* ----------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collisional
*  transfer between a 2sigma state and one or more vibrational levels
*  in a 2pi state in collisions with a structureless atom or with an
*  uncorrugated surface.  it is assumed that there is no isolated-
*  molecule mixing of the sigma and pi states.
*
*  this is a modification of basgpi, whose authors re given below
*  authors:  original version by millard alexander and didier lemoine
*            rewritten and extended to vibration and sigma-pi coupling
*            by h.-j. werner
*  author of basgpi1:  p.j.dagdigian
*  current revision date:  12-dec-2012 by p.j.dagdigian (to make sign
*  of c21-c32 coefficients for Pi F2 levels match those in ba2pi basis
*  routine)
*  --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index of each channel
*  note that we have adopted the following convention for the symmetry
*  index "is" so that on return the doublet pi molecular levels can be
*  uniquely identified by the two labels "j" and "is":
*           for Fi = 1   is = (100+10*ivpi)*eps
*           for Fi = 2   is = (200+10*ivpi)*eps
*           for sigma    is = 300*eps
*        where  eps = +1 or -1 is the "true" case (a) symmetry index
*        and ivpi is the 2pi vibrational quantum number
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each energetically
*              distinct level, similar to is
*    ivhold:   on return, contains vibr. quantum numbers for each level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!
*    jtot:     total angular momentum
*              in cc calculation jtot+1/2 is the total angular momentum
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
*    ihomo:    if .true., then homonuclear molecule
*              only the s or a levels will be included depending on the
*              value of the parameter isa in common /cosysi/ (see below)
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
*              here s=1 for sigma-minus and s=0 otherwise
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
*    esg:      energy of 2sigma level in cm-1
*    bsg:      rotational constant for 2sigma state in cm-1
*    dsg:      centrifugal distortion constant for 2sigma state in cm-1
*    gsr:      2sigma state spin-rotation constant in cm-1
*    epi:      array of energies of 2pi levels in cm=1
*    bpi:      array of rotational constants for 2pi levels in cm-1
*    dpi:      array of centrifugal distortion constants for 2pi levels in cm-1
*    aso:      array if spin-orbit constants for 2pi levels in cm-1
*    p, q:     array if lambda-doubling constants for 2pi levels in cm-1
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential.  This equals 4 for coupling of the
*              sigma state with one pi level, 7 for coupling with 2 pi
*              levels, 10 for coupling with 3 pi levels, etc.
*              term must be consistent with the pot routine.
*    isym:     if isym=+1, then the electronic symmetry is sigma-plus
*              if isym=-1, then the electronic symmetry is sigma-minus
*    isa:      s/a symmetry index, if the molecule is homonuclear (ihomo=t)
*              then, if isa=+1 then only the s-levels are included in the
*    igusg:    if igu=+1, then the inversion symmetry is gerade
*              if igu=-1, then the inversion symmetry is ungerade
*    nmaxsg:   maximum case (b) rotational angular momenta for 2sigma state
*    nparsg:   number of symmetry doublets included (npar=2 will ensure
*              both spin doublets) for 2sigma state.  nparsg=1 for eps=+1 levels,
*              nparsg=-1 for eps=-1 levels.
*              basis, if isa=-1, then only the a-levels are included
*    igupi:    permutation inversion symmetry of 2pi state
*              igu=1 for gerade states, igu=-1 for ungerade states
*              for heteronuclear molecules igu should be +1
*    nparpi:   number of 2pi symmetry doublets included (nparpi=2 will ensure
*              both lambda doublets).  otherwise, set nparpi=+1 for e levels,
*              nparpi=-1 for f levels only.
*    numvpi:   number of 2pi vibrational levels (must be >= 1) (place holder)
*    jmax:     array of maximum rotational angular momenta for each 2pi level
*              with convention omega .le. j .le. jmax+0.5
*  variables in common block /covibp/
*    ivpi:     array of vbibrational quantum numbers of the 2pi state to be
*              included in the calculation.  must be consistent with order of
*              levels listed in common block cvibp, from pot routine.
*  variables in common block /covibp/
*    numvib:   number of 2pi vibrational levels for which  the pot routine
*              can provide Vpi, V1, V2 for coupling with the 2sigma state
*    ivibpi:   array of the vibrational quantum numbers of 2pi levels for
*              which the pot routine can provide the potentials
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
*  variables in common block /cosgpi/
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
*  variable in common block /coconv/
*    econv:     conversion factor from cm-1 to hartrees
*    xmconv:    converson factor from amu to atomic units
* --------------------------------------------------------------------
*  subroutines called:
*   iswap:      interchanges a pair of integer variables (in hibasgpi.f)
*   swap:       interchanges a pair of real variables (in hibasgpi.f)
*   vlm2sg:     computes coupling matrix elements between 2sigma levels
*               (in hiba2sg.f)
*   vlm2pi:     computes coupling matrix elements between case a 2pi levels
*               (in hiba2pi.f)
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag, clist, flaghf, flagsu, ihomo, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysi/ nscode, isicod, nterm, isym, isa, igusg,
     :  nmaxsg, nparsg, igupi, nparpi, numvpi, ispar(10)
      common /cosysr/ isrcod, junkr, esg, bsg, dsg, gsr, rspar(40)
      common /coipar/ iiipar(9), iprint
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /cotq1/ vec(3,3)
      common /coisc1/ ivec(1)
      common /coisc2/ nrot(1)
      common /coisc3/ ifi(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /covibp/ ivpi(5)
      common /covpot/ numvib,ivibpi(5)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      dimension j(1), l(1), is(1), jhold(1), ehold(1), ishold(1),
     :  ivhold(1), c12(1), c32(1), csig(1), ieps(2),
     :  epi(maxvib), bpi(maxvib), dpi(maxvib), aso(maxvib), p(maxvib),
     :  q(maxvib), jmax(maxvib)
      parameter (nvmax=10)
* this determines which eps level is first in channel list (arbitrary)
      data ieps / -1, 1 /
      data izero, ione, itwo / 0,     1,    2 /
      pi2 = 1.570796326794897d0
      zero = 0.d0
      tzero=1.d-12
      one = 1.d0
      two = 2.d0
      four = 4.d0
      half = 0.5d0
      quart = 0.25d0
      xjtot = jtot + half
      xnu = nu + half
*  check for consistency in the values of flaghf and csflag
      if (.not.flaghf) then
        write (6, 2)
        write (9, 2)
2       format (' *** FLAGHF = .FALSE. FOR DOUBLET SYSTEM; ABORT ***')
        call exit
      end if
      if (flagsu .and. .not. csflag) then
        write (6, 3)
        write (9, 3)
3       format
     :   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
        call exit
      end if
*  check that isym equals +1 or -1
      if (abs(isym).ne.1) then
        write (6, 12)
        write (6, 12)
12      format (' *** ISYM MUST EQUAL +1 OR -1; ABORT ***')
        call exit
      end if
      istep = 1
      if (ihomo) istep = 2
      do 6  i = 1, nterm
        if (ihomo) then
          if (mod(lammax(i)-lammin(i),2) .ne. 0) then
            write (6, 4) i, lammin(i), lammax(i)
            write (9, 4) i, lammin(i), lammax(i)
4           format (' *** IHOMO=T BUT ODD NO. OF TERMS FOR I=', i2/,
     :              '     LAMMIN=', i2, ' LAMMAX=', i2, '; ABORT ***')
            call exit
          end if
        end if
6     continue
      iicod = 0
      ircod = 0
      do 8 iv=1,numvpi
        jmax(iv) = ispar(iicod+1)
        iicod = iicod + 1
        epi(iv) = rspar(ircod+1)
        bpi(iv) = rspar(ircod+2)
        dpi(iv) = rspar(ircod+3)
        aso(iv) = rspar(ircod+4)
        p(iv) = rspar(ircod+5)
        q(iv) = rspar(ircod+6)
        ircod = ircod+6
8     continue
      if (clist) then
        if (flagsu) then
          write (9, 7) rmu * xmconv, ered * econv, jtot, xnu
          write (6, 7) rmu * xmconv, ered * econv, jtot, xnu
7         format (/,' ** 2SIG + 2PI INT. COUPLING UNCORR. SURFACE **',
     :            /,'    RMU=', f9.4,'  E=', f7.2, '  LBAR=', i5,
     :              '  NU=', f5.1)
        else
          write (9, 27) igusg,isym,isa,igupi,isa
          write (6, 27) igusg,isym,isa,igupi,isa
27        format (/,' **  2SIG + 2PI INT. COUPLING **'/
     :              '     2SIG:   g/u=', i2,'  +/- = ', i2,
     :              '  s/a = ',i2/
     :              '     2PI:    g/u=', i2,'  s/a = ',i2)
          if (csflag) then
            write (9, 30) rmu * xmconv, ered * econv, jtot, xnu
            write (6, 30) rmu * xmconv, ered * econv, jtot, xnu
30          format (/,' ** CS **    RMU=', f9.4,'  E=', f8.2,
     :              '  LBAR=', i5, '  NU=', f5.1)
          else
            write (9, 31) rmu * xmconv, ered * econv, xjtot, jlpar
            write (6, 31) rmu * xmconv, ered * econv, xjtot, jlpar
31          format (/,' ** CC **    RMU=', f9.4,'  E=', f8.2,
     :              '  JTOT=', f5.1, ' JLPAR=', i2)
          end if
          write(6,35)
          write(9,35)
35        format(/' 2SIG PARAMETERS:'/
     :       ' NMAX',5x,'ESIG',5x,'BROT',8x,'DROT',8x,'GSR')
          write(6,36) nmaxsg, esg, bsg, dsg, gsr
          write(9,36) nmaxsg, esg, bsg, dsg, gsr
36        format(1x,i4,f11.3,1x,3(1pe12.4))
          write(9,32)
          write(6,32)
32        format(/' 2PI PARAMETERS:'/
     :      '  V JMAX',5x,'EVIB',5x,'BROT',8x,'DROT',8x,
     :      'ASO',9x,'P',11x,'Q')
          do 34 jj=1,numvpi
            i=ivpi(jj)
            write (6, 33) i,jmax(jj),epi(jj),bpi(jj),dpi(jj),aso(jj),
     :        p(jj),q(jj)
            write (9, 33) i,jmax(jj),epi(jj),bpi(jj),dpi(jj),aso(jj),
     :        p(jj),q(jj)
33          format(1x,i2,i5,f9.3,1pe13.4,4(1pe12.4))
34        continue
        end if
        if (.not. flagsu) write (9, 392) rcut
392     format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :         f6.2)
      end if
*
*    every vibrational energy level is referred to the lowest
*    vibrational level
*
      n=0
*  first set up list of 2sigma rotational levels
      n = 0
*  the n=0, eps=+1, j=1/2 level is
*     s for sigma-g-plus or for sigma-u-minus
*     a for sigma-g-minus or for sigma-g-plus
      isplus = igusg * isym * isa
      do 145 nnrot = 0, nmaxsg
        npbot = 1
        nptop = 2
        if (nparsg.eq.1) npbot = 2
        if (nparsg.eq.-1) nptop = 1
        do 140 ip = npbot, nptop
*  include only the eps=+1 level for n=0
          if (nnrot .eq. 0 .and. ip .eq. 1) go to 140
*  now calculate j for each level
            if (ip .eq. 1) then
              ji = nnrot - 1
            else
              ji = nnrot
            end if
*  actual half integer value of j is ji + 1/2
*  if homonuclear, include level only if allowed
*  see table i of alexander and corey, j. chem. phys. 84, 100 (1986)
            if (.not.ihomo .or. isa.eq.0. or.
     :        (ihomo .and. ieps(ip)*(-1)**ji.eq.isplus)) then
              n = n + 1
              if (n .gt. nmax) then
                write (6, 141) n, nmax
                write (9, 141) n, nmax
141             format(/' *** NCHANNELS=', i4,
     :            ' .GT. MAX DIMENSION OF',i4,' ABORT ***')
               stop
             end if
             is(n) = ieps(ip)
             j(n) = ji
             nrot(n) = nnrot
*  now assign energies for case (a) level and store in array eint
*    the matrix elements are given by a. j. kotlar, r. w. field,
*    and j. i. steinfeld, j. mol. spectr. 80, 86 (1980)
             x = float(j(n)) + 1.
             nn1 = x * (x - is(n))
             eint(n) = esg + bsg * nn1 - dsg*nn1**2
     :            - half * (1 - is(n) * x) * gsr
*  parameters to designate level as belonging to 2sigma state
             is(n) = 300 * is(n) * isym
             c12(n) = 0.
             c32(n) = 0.
             csig(n) = 1.
             ifi(n) = 3
             ivhold(n) = 100
          end if
140     continue
145   continue
*
*  now add 2pi levels
*
*  first set up list of all case (a) levels for omega = 1/2
*  for homonuclear molecules in gerade electronic states the e levels
*  are s for j=1/2, a for j=3/2, etc. while the f levels are a for
*  j=1/2, s for j=3/2, etc.  this reverses for ungerade states.
*  see table i of alexander and corey, j. chem. phys. 84, 100 (1986)
      do 103 iv=1,numvpi
        n0 = n
        do 45 ji = 0, jmax(iv)
          npbot=1
          nptop=2
          if (nparpi .eq. 1) npbot=2
          if (nparpi .eq. -1) nptop=1
          do 45 ip = npbot,nptop
*  ipar = 1 for e levels, ipar=-1 for f levels
            ipar = ieps(ip) * (-1) ** ji
            if (.not. ihomo .or. isa.eq.0. or.
     :        (ihomo .and. (ipar*igupi .eq. isa) ) ) then
              n = n + 1
              if (n .gt. nmax) then
                write (6, 40) n, nmax
                write (9, 40) n, nmax
40              format(/' *** NCHANNELS=', i4,
     :           ' .GT. MAX DIMENSION OF',i4,' ABORT ***')
                stop
              end if
              is(n) = ieps(ip)
              j(n) = ji
            end if
45      continue
*  now assign omega values and energies for case (a) levels
        do 50 i = n0 + 1, n
*  initialize the array ifi to be 1 (corresponding to F1)
          ifi(i) = 1
*  now set up arrays of internal energies for the case (a) levels
*    eint contains the energies of the omega = 1/2 levels
*    cent contains the energies of the omega = 3/2 levels
*    the matrix elements are given by a. j. kotlar, r. w. field,
*    and j. i. steinfeld, j. mol. spectr. 80, 86 (1980)
*    x is j + 1/2
          x = j(i) + one
          xsq = x*x
          eint(i) = epi(iv) - half * aso(iv) + bpi(iv) * xsq
     :      + dpi(iv) * ((1 - xsq**2) - xsq)
     :      + half * (one - is(i) * x) * (p(iv)
     :      + q(iv) * (one - is(i) * x))
          cent(i) = epi(iv) + half * aso(iv) + bpi(iv) * (xsq - two)
     :     + dpi(iv) * ((1 - xsq) - (xsq - two)**2)
     :     + half * (xsq - one) * q(iv)
*  now calculate the mixing angle due to the j.s term in the hamiltonian
*  store this angle temporarily in the array c12 and store the 1/2 - 3/2
*  coupling matrix element temporarily in the array c32
          if (j(i) .eq. 0) then
*  no mixing if j = 1/2
            c32(i) = 0.
            c12(i) = 0.
          else if (j(i) .gt. 0) then
            c32(i) = - quart * sqrt(xsq - one)
     :        * (four * bpi(iv) - two * four * dpi(iv) * (xsq - one)
     :        + p(iv) + two * (one - is(i) * x) * q(iv))
            adum1=two*c32(i)
            adum2=eint(i) - cent(i)
            c12(i) = half * atan2(adum1,adum2) + pi2
*  below to indicate the this level is not 2sigma
            csig(i) = 0.
            ivhold(i) = ivpi(iv)
          end if
50      continue
*  n - n0 now contains the number of omega = 1/2 channels levels
*  check if j=1/2 is to be assigned to F1 or F2 manifold
        y = aso(iv)/bpi(iv)
        if (y .le. 2) then
*  reorder omega=1/2 channels so that j=1/2 comes at the end, appropriate
*  to F2
          nmn = 1
          if (j(n0+2) .eq. 0) nmn=2
* nmn is the number of j=0 levels
          eh1 = eint(n0+1)
          eh2 = eint(n0+2)
          centh1 = cent(n0+1)
          centh2 = cent(n0+2)
          c121 = c12(n0+1)
          c122 = c12(n0+2)
          c321 = c32(n0+1)
          c322 = c32(n0+2)
          csig1 = csig(n0+1)
          csig2 = csig(n0+2)
          ish1 = is(n0+1)
          ish2 = is(n0+2)
          ivh1 = ivhold(n0+1)
          ivh2 = ivhold(n0+2)
          do 53 i = n0+nmn+1, n
            j(i-nmn) = j(i)
            eint(i-nmn) = eint(i)
            cent(i-nmn) = cent(i)
            c12(i-nmn) = c12(i)
            c32(i-nmn) = c32(i)
            csig(i-nmn) = csig(i)
            is(i-nmn) = is(i)
            ivhold(i-nmn) = ivhold(i)
53        continue
          if (nmn .eq. 2) then
            eint(n-1) = eh1
            eint(n) = eh2
            cent(n-1) = centh1
            cent(n) = centh2
            c12(n-1) = c121
            c12(n) = c122
            c32(n-1) = c321
            c32(n) = c322
            csig(n-1) = csig1
            csig(n) = csig2
            j(n-1) = 0
            j(n) = 0
            ifi(n-1) = 2
            ifi(n) = 2
            is(n-1) = ish1
            is(n) = ish2
            ivhold(n-1) = ivh1
            ivhold(n) = ivh2
          else
            eint(n) = eh1
            cent(n) = centh1
            c12(n) = c121
            c32(n) = c321
            csig(n) = csig1
            j(n) = 0
            ifi(n) = 2
            is(n) = ish1
            ivhold(n) = ivh1
          endif
        endif
*  now add the omega = 3/2 levels
        nn = n
        do 60 i = n0 + 1, n
          if ( j(i) .eq. 0) then
*  here for j = 1/2, which can exist only for omega = 1/2 and hence
*  is always a pure case (a) state
            c12(i) = 1.
            c32(i) = 0.
            csig(i) = 0.
            ivhold(i) = ivpi(iv)
          else if ( j(i) .ge. 1) then
*  here for j .ge. 3/2
*  the index nn points to the nominal omega = 3/2 channels
            nn = nn + 1
            if (nn .gt. nmax) then
              write (6, 40) nn, nmax
              write (9, 40) nn, nmax
              stop
            end if
            j(nn) = j(i)
            l(nn) = l(i)
* assign F2 label to each level
            ifi(nn) = 2
            is(nn) = is(i)
            ivhold(nn) = ivhold(i)
*  now determine array of mixing coefficients:
*    the c12 array indicates the amount of omega = 1/2 component in each
*    mixed level and the c32 array indicates the amount of omega = 3/2
*    component
*  also save the j.s matrix element temporarily in the variable h12
*  and save the case(a) energies temporarily in the variables e12 and e32
            h12 = c32(i)
            e12 = eint(i)
            e32 = cent(i)
            sn = sin(c12(i))
            cs = cos(c12(i))
            c32(i) = sn
            c12(i) = cs
            c32(nn) = cs
            c12(nn) = - sn
*  below to designate that level is not in 2sigma state
            csig(i) = 0.
            ivhold(i) = ivpi(iv)
            csig(nn) = 0.
            ivhold(nn) = ivpi(iv)
*  now calculate the energies of the mixed states and store these in the
*  array eint
            eint(i) = e12 * c12(i)**2 + e32 * c32(i)**2
     :            + two * c12(i) * c32(i) * h12
            eint(nn) = e12 * c12(nn)**2 + e32 * c32(nn)**2
     :            + two * c12(nn) * c32(nn) * h12
          end if
60      continue
*  attach full label to "is"
        do 62 i=n0 + 1,nn
          is(i) = (100*ifi(i) + 10*ivpi(iv)) * is(i)
*  make c12 coeff for Fi = 1 levels negative, to match ba2pi basis routine
          if (j(i).ne.0) then
            if ((abs(c12(i)) .gt. abs(c32(i))) .and.
     :          (c12(i).gt.0.)) then
              c12(i) = -c12(i)
              c32(i) = -c32(i)
*  make c32 coeff for Fi = 2 levels positive
            else if ((abs(c32(i)) .gt. abs(c12(i)))
     :          .and. (c32(i).lt.0.)) then
              c12(i) = -c12(i)
              c32(i) = -c32(i)
            end if
          end if
62      continue
*  set n to equal total number of levels
        n = nn
*  find lowest energy
        emin = 1. e+7
        do 65  i = 1, n
          if (eint(i) .lt. emin) emin = eint(i)
65      continue
103   continue
*
*  form list of all energetically distinct rotational levels included in the
*  channel basis and their energies (with zero of energy set at lowest level)
      nlevel = 0
      do 70  i = 1, n
        nlevel = nlevel + 1
        ehold(nlevel) = (eint(i) - emin) / econv
        jhold(nlevel) = j(i)
        ishold(nlevel) = is(i)
70    continue
*
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
              call iswap(jhold(i), jhold(ii))
              call iswap(ishold(i), ishold(ii))
              call iswap(ivhold(i), ivhold(ii))
              call rswap(ehold(i), ehold(ii))
              go to 80
            end if
75        continue
        end if
80    continue
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
*
*  set up coupled-states channel basis (if desired)
      if (csflag) then
        nn = 0
        do 85  i = 1, n
*  include only channels with j at least equal to coupled states
*  projection index
          if (j(i) .ge. nu) then
            nn = nn + 1
            if (nn .gt. nmax) then
              write (6, 40) nn, nmax
              write (9, 40) nn, nmax
              stop
            end if
            eint(nn) = (eint(i) - emin) / econv
            j(nn) = j(i)
            c12(nn) = c12(i)
            c32(nn) = c32(i)
            csig(nn) = csig(i)
            is(nn) = is(i)
            ifi(nn) = ifi(i)
            l(nn) = jtot
            if (.not. boundc) then
                 cent(nn) = jtot * (jtot + 1)
            else
                 xjtot=jtot+0.5d0
                 xj=j(nn)+0.5d0
                 cent(nn)=xjtot*(xjtot+1)+xj*(xj+1)-2*xnu*xnu
            endif
          end if
85      continue
*  set number of coupled states channels
        n = nn
*
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
          c12(inew) = c12(iold)
          c32(inew) = c32(iold)
          csig(inew) = csig(iold)
          is(inew) = is(iold)
          ifi(inew) = ifi(iold)
90      continue
        nn = 0
        do 100  i = 1, n
*  now take (nmax-n+i)th element, duplicate it as many times as is
*  required for rotational degeneracy, and store the new elements in
*  nn, nn+1, ....
          ipoint = nmax - n + i
          ji = j(ipoint)
          lmax = jtot + ji + 1
          lmin = iabs (jtot - ji)
          isx = sign(ione,is(ipoint))
          do 95  li = lmin, lmax
            ix = (-1) ** (ji + li - jtot) * isx
            if (ix .eq. jlpar) then
              nn = nn + 1
              if (nn .gt. nmax) then
                write (6, 40) nn, nmax
                write (9, 40) nn, nmax
                stop
              end if
              eint(nn) = (eint(ipoint) - emin) / econv
              j(nn) = ji
              c12(nn) = c12(ipoint)
              c32(nn) = c32(ipoint)
              csig(nn) = csig(ipoint)
              is(nn) = is(ipoint)
              ifi(nn) = ifi(ipoint)
              l(nn) = li
              cent(nn) = li * ( li + 1)
            end if
95        continue
100     continue
*  set number of close coupled channels
        n = nn
      end if
*
*  now check to see if any of the open channels are closed at r = rcut
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
              c12(nn) = c12(i)
              c32(nn) = c32(i)
              csig(nn) = csig(i)
              is(nn) = is(i)
              cent(nn) = cent(i)
              ifi(nn) = ifi(i)
              l(nn) = l(i)
            end if
130       continue
*  reset number of channels
          n = nn
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
*
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
        end if
      end if
*  CHANGE HERE FOR TEST
      if (boundc) then
*         ntop = ntop - 2*numin
*  mha and jacek's modification; make sure ntop equals the number of channels
         ntop = n + 1
      end if
*  now list channels if requested
      if (clist) then
        if (.not.csflag) then
          if (bastst) write (6, 200)
          write (9,200)
200       format
     :     (/'   N   V   J    EPS   FI    L    EINT(CM-1)    C-1/2',
     :        '      C-3/2      CSIG')
        else
          if (bastst) write (6,210) xnu
          write (9,210) xnu
210       format
     :     (/'   N   V   J    EPS   FI    L    EINT(CM-1)   C-1/2',
     :       '      C-3/2      CSIG','  **  NU=', f4.1)
        end if
        do 230  i = 1, n
          fj = j(i) + half
          ie = sign(ione,is(i))
          iv = abs((is(i) - 100*(is(i)/100))/10)
          ecm = eint(i) * econv
          if (bastst)
     :    write (6, 220) i, iv, fj, ie, ifi(i), l(i), ecm,
     :                   c12(i), c32(i),csig(i)
          write (9, 220) i, iv, fj, ie, ifi(i), l(i), ecm,
     :                   c12(i), c32(i),csig(i)
220       format (i4, i4, f5.1, i6, i4, i6, f13.3, 3f11.5)
230     continue
      endif
*
*  now calculate coupling matrix elements
      if (bastst .and.iprint.gt.1) then
        write (6, 235)
        write (9, 235)
235     format (/' ILAM   LAMBDA   MU     I   ICOL  IROW',
     :           '  IV2       VEE')
      end if
*
*  the following only for safety.  check the vibrational levels requested
*  consistent with pot routine
      if (numvpi.gt.numvib) then
        write (6,237) numvpi
237     format(' *** NVIBPI =',i3,' TOO LARGE FOR NUMBER OF VIB LEVELS',
     :    ' SUPPORTED BY POT ROUTINE; ABORT ***')
        stop
      end if
*  checks that pot routine contains all required vib. levels in the same order
      ifind = 1
      do 250 iv=1,numvpi
        if (ivpi(iv).ne.ivibpi(iv)) then
          ifind = 0
          write (6,246) ivpi(iv)
246       format(' *** V =',i3,' IS IN INPUT FILE.  POT ROUTINE NOT',
     :      ' EXPECTING THIS VIBRATIONAL LEVEL')
        end if
250   continue
      if (ifind.eq.0) then
        write(6,247)
247     format(/' *** VIBRATIONAL LEVELS NOT CONSISTENT WITH POT',
     :    ' ROUTINE.  ABORT ***')
        stop
      end if
*
*  i counts v2 elements
*  ilam counts number of v2 matrices
*  inum counts v2 elements for given ilam
*  ij is address of given v2 element in present v2 matrix
*  nlam is the total number of anisotropic terms
      nlam = 0
      i = 0
      lamsum = 0
      istep = 1
      if (ihomo) istep = 2
      ilam = 0
      do 600 iterm=1,nterm
*  below for 2sigma potential Vsig
        if (iterm.eq.1) then
          ipot = 1
          goto 499
        end if
*  below for 2pi potential Vpi
        if (mod(iterm,3).eq.2) ipot = 2
*  below for 2pi potential V2
        if (mod(iterm,3).eq.0) ipot = 3
*  below for 2sigma-2pi coupling potential V1
        if (mod(iterm,3).eq.1 .and. iterm.ne.1) ipot = 4
499     nlami = (lammax(iterm) - lammin(iterm))/istep + 1
        do 500 il=1,nlami
          nlam = nlam + 1
          lb = lammin(iterm) + (il - 1) * istep
          mu = mproj(iterm)
          inum = 0
          do 450 icol = 1, n
          do 450 irow = icol, n
            ij = ntop * (icol - 1) + irow
            vee = zero
            isignr = sign(ione,is(irow))
            isignc = sign(ione,is(icol))
            iomr = abs(is(irow))/100
            iomc = abs(is(icol))/100
            lrow = l(irow)
            lcol = l(icol)
            ivibr = mod(abs(is(irow)),100)/10
            ivibc = mod(abs(is(icol)),100)/10
            if (csflag) then
              lrow = nu
              lcol = nu
            end if
            goto (460, 470, 480, 490), ipot
*  Vsig
460         if (iomr.eq.3 .and. iomc.eq.3) then
              call vlm2sg(j(irow), lrow, j(icol), lcol,
     :          jtot, lb, isignr, isignc, vee, csflag)
            end if
            goto 440
*  Vpi
470         if (iomr.lt.3 .and. iomc.lt.3 .and.
     :          ivibr.eq.ivibc .and. ivibr.eq.(iterm-2)/3) then
              call vlm2pi(j(irow), lrow, j(icol), lcol,
     :          jtot, izero, izero, lb, isignr, isignc,
     :          v1212, csflag)
              call vlm2pi(j(irow), lrow, j(icol), lcol,
     :          jtot, ione, ione, lb, isignr, isignc,
     :          v3232, csflag)
              vee = c12(irow) * c12(icol) * v1212
     :          + c32(irow) * c32(icol) * v3232
            end if
            goto 440
*  V2
480         if (iomr.lt.3 .and. iomc.lt.3 .and.
     :          ivibr.eq.ivibc .and. ivibr.eq.(iterm-3)/3) then
              call vlm2pi(j(irow), lrow, j(icol), lcol,
     :          jtot, izero, ione, lb, isignr, isignc,
     :          v1212, csflag)
              call vlm2pi(j(irow), lrow, j(icol), lcol,
     :          jtot, ione, izero, lb, isignr, isignc,
     :          v3232, csflag)
              vee = c12(irow) * c32(icol) * v1212
     :          + c32(irow) * c12(icol) * v3232
            end if
            goto 440
*  V1
490         if (iomr.eq.3 .and. iomc.lt.3 .and.
     :          ivibc.eq.(iterm-4)/3) then
              call vlmsp1 (j(irow), lrow, j(icol),
     :          lcol, jtot, lb, izero,
     :          isignr, isignc, v12, csflag)
              call vlmsp1 (j(irow), lrow, j(icol),
     :          lrow, jtot, lb, ione,
     :          isignr, isignc, v32, csflag)
              vee = c12(icol) * v12 + c32(icol) * v32
            end if
            if (iomc.eq.3 .and. iomr.lt.3 .and.
     :          ivibr.eq.(iterm-4)/3) then
              call vlmsp1 (j(icol), lcol, j(irow),
     :          lrow, jtot, lb, izero,
     :          isignc, isignr, v12, csflag)
              call vlmsp1 (j(icol), lcol, j(irow),
     :          lrow, jtot, lb, ione,
     :          isignc, isignr, v32, csflag)
              vee = c12(irow) * v12 + c32(irow) * v32
            end if
440         if (abs(vee).lt.1.d-15) goto 450
            i = i + 1
            inum = inum + 1
            if (i.gt.nv2max) goto 450
            v2(i) = vee
            iv2(i) = ij
            if (.not.bastst .or. iprint.le.1) goto 450
            write (6,495) nlam, lb, mu, i, icol, irow, iv2(i), vee
            write (9,495) nlam, lb, mu, i, icol, irow, iv2(i), vee
495         format (i4, 3i7, 2i6, i6, g17.8)
450       continue
          if (i.le.nv2max) lamnum(nlam) = inum
          if (bastst) then
            write (6, 360) iterm, nlam, lb, mu, lamnum(nlam)
            write (9, 360) iterm, nlam, lb, mu, lamnum(nlam)
360         format ('ITERM=',i3,' ILAM=', i3,' LAMBDA=',i3,
     :        ' MU=',i3,' LAMNUM(ILAM)=',i8)
          end if
          lamsum = lamsum + lamnum(il)
500     continue
600   continue
      if (nlam.gt.nlammx) then
        write(6,610) nlam,nlammx
        write(9,610) nlam,nlammx
610     format (' *** NLAM = ',i6,' .GT. NLAMMX=',i6,'; ABORT ***')
        call exit
      end if
      if (i.gt.nv2max) then
        write (6, 620) i, nv2max
        write (9, 620) i, nv2max
620     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
        call exit
      end if
      if (clist) then
        write (6, 630) i
        write (9, 630) i
630     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i9)
      end if
      if (bastst .and. iprint .ge. 1) then
        write (6, 660)
660     format(/1x,'   N   V   J   EPS   FI    L   EINT(CM-1)',
     :   '     C-1/2       C-3/2       C-SIG')
        do 680  i = 1, n
          fj = j(i) + half
          ie = sign(ione,is(i))
          iv = abs((is(i) - 100*(is(i)/100))/10)
          iom = abs(is(i)/100)
          ecm = eint(i) * econv
          write (6, 670) i,iv,fj,ie,iom,l(i), ecm,
     :        c12(i), c32(i), csig(i)
670       format (1x,2i4, f5.1, i5, i4, i6, f13.3, 3f12.6)
680     continue
      end if
      return
1000  write (6, 1010) n, nmax
      write (9, 1010) n, nmax
1010  format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF',
     :         i4,' ABORT ***')
      call exit
      end
* --------------------------------------------------------------------
      subroutine vlmsp1 (jp, lp, j, l, jtot, lambda, iomeg, iepsp,
     :                   ieps, v, csflag)
*  subroutine to calculate v-lambda matrices for close-coupled and
*  coupled-states treatments of collisions of 2sigma <-> 2pi transitions
*
*  V1 term
*
*  modification of vlm2sg -- written by p.j.dagdigian - 6-dec-2011
*  current revision date:  2-apr-2012
*
*  the cc matrix elements are given in eqs. (33,36,37) of m.h. alexander
*  and g.c. corey, j. chem. phys. 84, 100 (1986), and in eq.(7) of
*  dagdigian et al., j. chem. phys. 98, 8580 (1993)
*  note that for cc collisions of a 2sigma molecule with a flat surface, the
*  coupling matrix elements are assumed to be identical to the cs matrix
*  elements here
*  -----------------------------------------------------------------------
*  variables in call list:
*    it is assumed that the bra is a level in the 2sigma state and the ket
*    is a level in the 2pi state
*    jp:       rotational quantum number of left side of matrix element (bra)
*              minus 1/2 (2sigma state)
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*              minus 1/2 (2pi state)
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    lambda:   order of legendre term in expansion of potential
*    iomeg:    omega quantum number of ket (=0 for omega=1/2, =1 for omega=3/2,
*              for the 2pi state)
*    iepsp:    symmetry index of bra (i.e. for level in 2sigma state)
*    ieps:     symmetry index of ket (i.e. for level in 2pi state)
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
      data zero, half, one, thrtwo, two, halfm, onem
     :  /0.d0, 0.5d0, 1.d0, 1.5d0, 2.d0, -0.5d0, -1.d0/
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
      if (iomeg .eq. 0) then
*  omega = 1/2 for 2pi level
        x = x * ieps * xf3j (xjp, xlamda, xj, halfm, one, halfm)
      else
*  omega = 3/2 for 2pi level
*  factor of (-1) takes care of (-1) ** iomeg phase factor
        x = x * xf3j (xjp, xlamda, xj, halfm, onem, thrtwo) * (-1)
      end if
      iphase = (-1) ** iabs(iphase)
      v = iphase * x * sqrt(xnorm)
      return
      end
* --------------------------------------------------------------------
