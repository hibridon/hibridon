* --------------------------------------------------------------------
      subroutine ba3p2s (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist,
     :     bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine coupling potential (electrostatic + spin-orbit)
*  for collision involving the 3P state of atom in a p^4 electronic
*  configuration with an atom in a 2S state
*  this is a revision of the basis routine written originally by p.dagdigian
*  the latter is stored in /hibxx/src as hiba3p2s.f~origPJD
*  the revision is based no notes written by m. alexander, p. dagdigian, and
*  j. klos:
*     "The Calculation of (some) atomic and molecular spin-orbit coupling
*     matrix elements"
*     "The case (e) scattering basis for O(3P)-H(2S) collisions" and addendum
*     "MAtrix element of the potential in HIBRIDON"
*
*  the spin-orbit part of the W matrix computed as W(R) - W(R=inf)
*  here, the spin-orbit parameters are assumed to be independent of R
*
*  j12 array is in common block /coj12/ to pass to other subrs
*
*  author:  paul dagdigian
*  current revision date:  4-sep-2018 by p.dagdigian
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains electronic angular momentum quantum number
*              for each channel (j = 0, 1, 2 for the 3P spin-orbit levels
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains electronic spin of each channel (0 or 1)
*    jhold:    on return contains electronic angular momentum quantum number
*              for each level
*    ehold:    on return contains energy in hartrees of each level
*    ishold:   on return contains spin multiplicity of each level
*    nlevel:   on return contains number of energetically distinct
*              levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              levels used in channel basis which are open
*              asymptotically
*    isc1:     scratch vector of length at least nmax
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              NOTE:  rcut option not implemented here
*    jtot:     total angular momentum
*    flaghf:   if .true., then system with half-integer spin
*              if .false., then system with integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              if the molecule is homonuclear (ihomo = .true.), the
*              rotational levels included go from jmin to jmax in steps
*              of 2 and only even lambda terms in the anisotropy are
*              included
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(l-jtot)=jlpar
*              n.b. jlpar=+1 corresponds to f levels, jlpar=-1, to e levels
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
*  variables in common block /cosysr/
*    isrcod:   number of real system dependent variables (none here)
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different types of electronic coupling terms
*              this should be 1 here
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variable in common block /coj12/
*    j12:      array containing vector sum of ja + j (similar
*              situation with molecule-molecule scattering,
*              see hibastpln basis routine)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of case(a) interaction potentials and spin-orbit
*               matrix elements actually used
*    nlammx:    the maximum number of angular coupling terms
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
*  variable in common block /coconv/
*     econv:    conversion factor from cm-1 to hartrees
*     xmconv:   converson factor from amu to atomic units
*  subroutines called:
*   vlm13p:    returns angular coupling coefficient for particular
*              choice of channel index
*
*  subroutines called:
*  xf3j:     evaluates 3j symbol
*  xf9j:     evaluates 9j symbol
* ------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysi/ nscode, isicod, nterm, nstate
      common /cosysr/ isrcod
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coj12/  j12(1)
      common /coered/ ered, rmu
      common /coskip/ nskip, iskip
      common /covvl/  vvl(6)
*  econv is conversion factor from cm-1 to hartrees
*  xmconv is converson factor from amu to atomic units
      common /coconv/ econv, xmconv
*  arrays in argument list
      dimension j(1), l(1), is(1), jhold(1), ehold(1), ishold(1)
      data fjh / 0.5d0 /
*  check for consistency in the values of flaghf, csflag, and flagsu
      if (flaghf) then
        write (6, 5)
        write (9, 5)
5       format (' *** FLAGHF = .TRUE. FOR 3P/2S ATOM; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
c
      if (csflag) then
        write (6, 8)
        write (9, 8)
8       format('  *** CSFLAG = .TRUE. FOR 3P/2S ATOM; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
c
      if (flagsu) then
        write (6, 109)
        write (9, 109)
109     format('  *** FLAGSU = .TRUE. FOR 3P/2S ATOM; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
c
      if (bastst) then
        write (6, 14) nlam
        write (9, 14) nlam
14      format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
        write (6,20) npot,
     :    rmu * xmconv, ered * econv, jtot, jlpar
        write (9,20) npot,
     :    rmu * xmconv, ered * econv, jtot, jlpar
20      format(/,' **  CC 3P/2S ATOM ; NPOT =',i2,' ** RMU=', f9.4,
     :    ' E=', f10.2, /, '         JTOT=', i5, 2x,' JLPAR=',i2)
      end if
*
*  form list of spin-orbit levels and determine the energetically open levels
*  order of the j levels for p^4 orbital occupancy
*
*  get atomic spin-orbit parameter (apar at large R)
      r = 100.d0
      call pot(vv0 ,r)
      a  = vvl(5) * econv
      nlevel = 0
      do jj = 2, 0, -1
        nlevel = nlevel + 1
        jhold(nlevel) = jj
        ishold(nlevel) = 3
        eso = - 0.5d0 * a * (jj * (jj + 1.d0) - 4.d0)
*  set zero of energy as E(j=2)
        ehold(nlevel) = (eso + 0.5d0 * a * 2.d0) / econv
      end do
*
*  form list of energetically open levels and sort list
*  to put closed channels at end
      nlevop = 0
      if (nlevel .gt. 1) then
        do 180 i = 1, nlevel - 1
          if (ehold(i) .le. ered) then
            nlevop = nlevop + 1
          else
            do 175 ii = i + 1, nlevel
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
                goto 180
              endif
175         continue
          endif
180     continue
        if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
      else
*  here for only one level
        if (ehold(1) .le. ered) then
          nlevop = 1
        else
          nlevop = 0
        endif
      endif
      if (nlevop .le. 0) then
        write (9,185)
        write (6,185)
185     format('*** NO OPEN LEVELS IN BA3P2S; ABORT')	
        if (bastst) return
        call exit
      endif
*
*  set up CC channel basis
      n = 0
      do 120 jj = 1, nlevel
        ji = jhold(jj)
*   determine j12 values - vector coupling of j(O) + j(H) [=1/2]
        fj12mn = abs(ji - fjh)
        fj12mx = ji + fjh
        j12min = fj12mn
        j12max = fj12mx
        do 119 j12i = j12min, j12max
          lmin = abs(jtot-j12i)
          lmax = jtot + j12i + 1
          do 110 li=lmin,lmax
*  parity in symmetry-adapted basis
            ix = (-1) ** li
            if (ix .eq. jlpar) then
*  here for correct orbital angular momentum
              n = n + 1
              if (n .gt. nmax) go to 130
              l(n) = li
              cent(n) = li * (li + 1.)
              is(n) = 3
              j(n) = ji
              j12(n) = j12i
              eint(n) = ehold(jj)
            end if
110       continue
119     continue
120   continue
130   if (n .gt. nmax) then
        write (9, 140) n, nmax
        write (6, 140) n, nmax
140     format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF',
     :         i4,'; ABORT')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
      ntop = max(n, nlevop)
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
*        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
*     :       ntop = ntop + 1
*
*  now list channels if requested
      if (bastst) then
        write (6,1255)
        write (9,1255)
1255    format(/' ASYMPTOTIC ATOMIC ENERGY LEVELS'
     +    /'   N   S   J      EINT(CM-1)')
        do 1265 i = 1, nlevel
          write (6,2260) i, ishold(i), jhold(i), ehold(i) * econv
          write (9,2260) i, ishold(i), jhold(i), ehold(i) * econv
2260      format (3i4, f13.3)
1265    continue
        write (6, 255)
        write (9, 255)
255     format(/' CHANNEL BASIS FUNCTIONS'
     +      /'   N   S   J   J12   L      EINT(CM-1)')
        do 265  i = 1, n
          write (6, 260) i, is(i), j(i), j12(i), l(i), eint(i) * econv
          write (9, 260) i, is(i), j(i), j12(i), l(i), eint(i) * econv
260       format (3i4, 2i5, f13.3)
265     continue
      end if
*
*  now calculate coupling matrix elements
      if (bastst) then
        write (6, 280)
        write (9, 280)
280     format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
      end if
*
*  i counts v2 elements
*  inum counts v2 elements for given lambda
*  ilam counts number of v2 matrices
*  ij is address of given v2 element in present v2 matrix
      i = 0
      ilam=0
      do 320 lb = 1, nlam
        ilam = ilam + 1
        inum = 0
        do 310  icol= 1, n
          do 300  irow = icol, n
            ij = ntop * (icol - 1) +irow
            call vlm3p2s (j(irow), j12(irow), l(irow),
     :        j(icol), j12(icol), l(icol), jtot, lb, vee)
            if (vee .eq. 0) goto 300
              i = i + 1
              inum = inum + 1
              if (i .gt. nv2max) goto 300
                v2(i) = vee
                iv2(i) = ij
                if (bastst) then
                  write (6, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
                  write (9, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
290               format (i4, 2i7, 2i6, i6, g17.8)
                endif
300       continue
310     continue
410   if(ilam.gt.nlammx) then
        write(6,311) ilam
311     format(/' ILAM.GT.NLAMMX IN BA3P2S')
        call exit
      end if
      lamnum(ilam) = inum
      if (bastst) then
        write (6, 315) ilam, lamnum(ilam)
        write (9, 315) ilam, lamnum(ilam)
315     format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
      end if
320   continue
      if (i .gt. nv2max) then
        write (6, 350) i, nv2max
        write (9, 350) i, nv2max
350     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (bastst) then
        write (6, 360) i
        write (9, 360) i
360     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i6)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine vlm3p2s (j1, j12_1, l1, j2, j12_2, l2, jtot, lb, vee)
* --------------------------------------------------------------------
*  subroutine to evaluate the angular coupling matrix element for
*  inelastic collisions of a 3P atom with a 2S atom

*  author:  paul dagdigian
*  current revision date: 4-sep-2018
* --------------------------------------------------------------------
*  variables in call list:
*  j1,j12_1,l1:    initial electronic and orbital angular momenta
*  j2,j12_2,l2:    final electronic and orbital angular momenta
*  jtot:     total angular momentum
*  lb:       value of expansion index:
*            lb=1:  X2Pi potential
*            lb=2:  4Pi potential
*            lb=3:  2Sigma- potential
*            lb=4:  4Sigma- potential
*            lb=5:  apar spin-orbit matrix element
*            lb=6:  aperp spin-orbit matrix element
*
*  the matrix element is evaluated using Eq. 31 in the report "The case (e) 
*  scattering basis for O(3P)-H(2S) collisions," along with the matrix elements
*  of the electrostatic and spin-orbit interactions in the fully coupled basis
*  given by Eqs. 187-190 and Eqs.. 181-186.
*
*  vee:      on return:  contains desired coupling matrix element
*  subroutines called:
*  xf3j:     evaluates 3j symbol
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      vee = 0.d0
      xl1 = l1
      xl2 = l2
      xj1 = j1
      xj2 = j2
      xj12_1 = j12_1 + 0.5d0
      xj12_2 = j12_2 + 0.5d0
      xjtot = jtot + 0.5d0
      xl12 = sqrt((2.d0*xl1 + 1.d0) * (2.d0*xl2 + 1.d0))
      xfj12 = xf3j(xl1, xj12_1, xjtot, 0.d0,-0.5d0, 0.5d0)
     +  * xf3j(xl2, xj12_2, xjtot, 0.d0,-0.5d0, 0.5d0)
      xfj32 = xf3j(xl1, xj12_1, xjtot, 0.d0, -1.5d0, 1.5d0)
     +  * xf3j(xl2, xj12_2, xjtot, 0.d0,-1.5d0, 1.5d0)
      xfj52 = xf3j(xl1, xj12_1, xjtot, 0.d0, -2.5d0, 2.5d0)
     +  * xf3j(xl2, xj12_2, xjtot, 0.d0,-2.5d0, 2.5d0)
*
*  spin-orbit parameters assumed to be indpendent of R
*  so matrix elements of apar and aperp not required
*  since Hso computed as Hso(R) - Hso(R=inf)
      iacnst = 1
*
      if (j1 .eq. 2 .and. xj12_1 .eq. 2.5d0) ist1 = 1
      if (j1 .eq. 2 .and. xj12_1 .eq. 1.5d0) ist1 = 2
      if (j1 .eq. 1 .and. xj12_1 .eq. 1.5d0) ist1 = 3
      if (j1 .eq. 1 .and. xj12_1 .eq. 0.5d0) ist1 = 4
      if (j1 .eq. 0 .and. xj12_1 .eq. 0.5d0) ist1 = 5
*
      if (j2 .eq. 2 .and. xj12_2 .eq. 2.5d0) ist2 = 1
      if (j2 .eq. 2 .and. xj12_2 .eq. 1.5d0) ist2 = 2
      if (j2 .eq. 1 .and. xj12_2 .eq. 1.5d0) ist2 = 3
      if (j2 .eq. 1 .and. xj12_2 .eq. 0.5d0) ist2 = 4
      if (j2 .eq. 0 .and. xj12_2 .eq. 0.5d0) ist2 = 5
*
      vee12 = 0.d0
      vee32 = 0.d0
      vee52 = 0.d0
      goto (10,20,30,40,50,60),lb
      goto 1500
*
*  Matrix elements below are in j_, j, Omega basis
*
*  X2Pi state
10    if (ist1 .eq. 2 .and. ist2 .eq. 2) vee32 = 25.d0 / 30.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 2)) vee32 = - sqrt(5.d0) * 5.d0 / 30.d0
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee32 = 1.d0 / 6.d0
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee12 = 25.d0 / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 2)) then
          vee12 = - sqrt(5.d0) * 5.d0 / 90.d0
      end if
      if ((ist1 .eq. 2 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 2)) then
          vee12 = sqrt(10.d0) * 10.d0 / 90.d0
      end if
      if ((ist1 .eq. 2 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 2)) then
          vee12 = - sqrt(5.d0) * 5.d0 / 45.d0
      end if
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee12 = 1.d0 / 18.d0
      if ((ist1 .eq. 3 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 3)) then
          vee12 = - sqrt(2.d0) * 2.d0 / 18.d0
      end if
      if ((ist1 .eq. 3 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 3)) vee12 = 1.d0 / 9.d0
      if (ist1 .eq. 4 .and. ist2 .eq. 4) vee12 = 4.d0 / 9.d0
      if ((ist1 .eq. 4 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 4)) vee12 = - sqrt(2.d0) * 2.d0 / 9.d0
      if (ist1 .eq. 5 .and. ist2 .eq. 5) vee12 = 2.d0 / 9.d0
      goto 1000
*
*  4Pi state
20    if (ist1 .eq. 1 .and. ist2 .eq. 1) vee52 = 1.d0
      if (ist1 .eq. 1 .and. ist2 .eq. 1) vee32 = 3.d0 / 5.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee32 = 1.d0 / 5.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq.1)) vee32 = 1.d0 / sqrt(5.d0)
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee32 = 2.d0 / 30.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 2)) vee32 = sqrt(5.d0) * 2.d0 / 30.d0
      if (ist1 .eq. 3 .and. ist2 .eq. 2) vee32 = 2.d0 / 6.d0
      if (ist1 .eq. 1 .and. ist2. eq. 1) vee12 = 2.d0 / 5.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.   
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee12 = sqrt(6.d0) / 30.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 3) .or.   
     +  (ist1 .eq. 3 .and. ist2 .eq. 1)) vee12 = 1.d0 / sqrt(30.d0)
      if ((ist1 .eq. 1 .and. ist2 .eq. 4) .or.   
     +  (ist1 .eq. 4 .and. ist2 .eq. 1)) vee12 = 1.d0 / sqrt(15.d0)
      if ((ist1 .eq. 1 .and. ist2 .eq. 5) .or.   
     +  (ist1 .eq. 5 .and. ist2 .eq. 1)) vee12 = sqrt(2.d0 / 15.d0)
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee12 = 14.d0 / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 3) .or.   
     +  (ist1 .eq. 3 .and. ist2 .eq. 2)) vee12 = sqrt(5.d0) * 14.d0 / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 4) .or.   
     +  (ist1 .eq. 4 .and. ist2 .eq. 2)) vee12 = - sqrt(10.d0) / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 5) .or.   
     +  (ist1 .eq. 5 .and. ist2 .eq. 2)) vee12 = - sqrt(5.d0) / 45.d0
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee12 = 14.d0 / 18.d0 
      if ((ist1 .eq. 3 .and. ist2 .eq. 4) .or.   
     +  (ist1 .eq. 4 .and. ist2 .eq. 3)) vee12 = - sqrt(2.d0) / 18.d0
      if ((ist1 .eq. 3 .and. ist2 .eq. 5) .or.   
     +  (ist1 .eq. 5 .and. ist2 .eq. 3)) vee12 = - 1.d0 / 9.d0
      if (ist1 .eq. 4 .and. ist2 .eq. 4) vee12 = 2.d0 / 9.d0 
      if ((ist1 .eq. 4 .and. ist2 .eq. 5) .or.   
     +  (ist1 .eq. 5 .and. ist2 .eq. 4)) vee12 = sqrt(2.d0) 
     +  * 2.d0 / 9.d0      
      if (ist1 .eq. 5 .and. ist2 .eq. 5) vee12 = 4.d0 / 9.d0 
      goto 1000
*
*  2Sigma- state
30    if (ist1 .eq. 2 .and. ist2 .eq. 2) vee12 = 50.d0 / 90.d0  
      if ((ist1 .eq. 2 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 2)) vee12 = - sqrt(5.d0)  
     +  * 10.d0 / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 2)) vee12 = - sqrt(10.d0) 
     +  * 10.d0 / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 2)) vee12 = sqrt(5.d0) * 5.d0 / 45.d0
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee12 = 2.d0 / 18.d0  
      if ((ist1 .eq. 3 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 3)) vee12 = sqrt(2.d0) * 2.d0 / 18.d0
      if ((ist1 .eq. 3 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 3)) vee12 = - 1.d0 / 9.d0
      if (ist1 .eq. 4 .and. ist2 .eq. 4) vee12 = 2.d0 / 9.d0  
      if ((ist1 .eq. 4 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 4)) vee12 = - sqrt(2.d0) / 9.d0
      if (ist1 .eq. 5 .and. ist2 .eq. 5) vee12 = 1.d0 / 9.d0  
      goto 1000
*
*  4Sigma- state
40    if (ist1 .eq. 1 .and. ist2 .eq. 1) vee32 = 2.d0 / 5.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee32 = - 1.d0 / 5.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 1)) vee32 = - 1.d0 / sqrt(5.d0)
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee32 = 3.d0 / 30.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 2)) then
          vee32 = sqrt(5.d0) * 3.d0 / 30.d0
      end if
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee32 = 3.d0 / 6.d0
      if (ist1 .eq. 1 .and. ist2 .eq. 1) vee12 = 3.d0 / 5.d0  
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee12 = sqrt(6.d0) / 30.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 1)) vee12 = - 1.d0 / sqrt(30.d0)
      if ((ist1 .eq. 1 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 1)) vee12 = - 1.d0 / sqrt(15.d0)
      if ((ist1 .eq. 1 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 1)) vee12 = - sqrt(2.d0 / 15.d0)
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee12 = 1.d0 / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 3) .or.
     +  (ist1 .eq. 3 .and. ist2 .eq. 2)) vee12 = sqrt(5.d0) / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 2)) vee12 = sqrt(10.d0) / 90.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 2)) vee12 = sqrt(5.d0) / 45.d0
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee12 = 1.d0 / 18.d0
      if ((ist1 .eq. 3 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 3)) vee12 = sqrt(2.d0) / 18.d0
      if ((ist1 .eq. 3 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 3)) vee12 = 1.d0 / 9.d0
      if (ist1 .eq. 4 .and. ist2 .eq. 4) vee12 = 1.d0 / 9.d0
      if ((ist1 .eq. 4 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 4)) vee12 = sqrt(2.d0) / 9.d0
      if (ist1 .eq. 5 .and. ist2 .eq. 5) vee12 = 2.d0 / 9.d0
      goto 1000
*
*  apar spin-orbit matrix element
50    if (ist1 .eq. 1 .and. ist2 .eq. 1) vee52 = - 0.5d0
      if (ist1 .eq. 1 .and. ist2 .eq. 1) vee32 = - 1.d0 / 10.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee32 = - 1.d0 / 5.d0
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee32 = - 4.d0 / 10.d0
      if (ist1 .eq. 3 .and. ist2 .eq. 2) vee32 = 0.5d0
      if (ist1 .eq. 1 .and. ist2 .eq. 1) vee12 = 1.d0 / 10.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee12 = - sqrt(6.d0) / 30.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 1)) vee12 = 1.d0 / sqrt(30.d0)
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee12 = - 2.d0 / 30.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 2)) vee12 = - 1.0 / sqrt(45.d0)
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee12 = 2.d0 / 6.d0
      if ((ist1 .eq. 3 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 3)) vee12 = - 1.0 / sqrt(18.d0)
      if (ist1 .eq. 4 .and. ist2 .eq. 4) vee12 = 1.d0 / 6.d0
      if (ist1 .eq. 5 .and. ist2 .eq. 5) vee12 = 1.d0 / 3.d0
      goto 1000
*
*  aperp spin-orbit matrix element
60    if (ist1 .eq. 1 .and. ist2 .eq. 1) vee32 = 4.d0 / 10.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee32 = 1.d0 / 5.d0
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee32 = - 1.d0 / 10.d0
      if (ist1 .eq. 1 .and. ist2 .eq. 1) vee12 = - 6.d0 / 10.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 2) .or.
     +  (ist1 .eq. 2 .and. ist2 .eq. 1)) vee12 = sqrt(6.d0) / 30.d0
      if ((ist1 .eq. 1 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 1)) vee12 = - 1.d0 / sqrt(30.d0)
      if (ist1 .eq. 2 .and. ist2 .eq. 2) vee12 = 17.d0 / 30.d0
      if ((ist1 .eq. 2 .and. ist2 .eq. 5) .or.
     +  (ist1 .eq. 5 .and. ist2 .eq. 2)) vee12 = 1.0 / sqrt(45.d0)
      if (ist1 .eq. 3 .and. ist2 .eq. 3) vee12 = 1.d0 / 6.d0
      if ((ist1 .eq. 3 .and. ist2 .eq. 4) .or.
     +  (ist1 .eq. 4 .and. ist2 .eq. 3)) vee12 = 1.0 / sqrt(18.d0)
      if (ist1 .eq. 4 .and. ist2 .eq. 4) vee12 = 2.d0 / 6.d0
      if (ist1 .eq. 5 .and. ist2 .eq. 5) vee12 = 2.d0 / 3.d0
      goto 1000
c
1000  vee = 2.d0 * xl12 * (vee12 * xfj12 + vee32 * xfj32
     +  + vee52 * xfj52)
*
*  provision for constant spin-orbit parameters
      if (iacnst .eq. 1) then
        if (lb .eq. 5 .or. lb .eq.6) vee = 0.d0
      end if
1500  return
      end
* ------------------------------eof-----------------------------------
