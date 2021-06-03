* --------------------------------------------------------------------
      subroutine ba1d3p (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine coupling potential (electrostatic + spin-orbit)
*  for collision involving the 1D and 3P states of atom in a p^2 or p^4
*  electronic configuration with a structureless atom
*  author:  paul dagdigian
*  (based on the ba13p basis routine)
*
*  the spin-orbit part of the W matrix is computed as W(R) - W(R=inf)
*
*  current revision date:  23-sep-2014 by p.dagdigian
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains electronic angular momentum quantum number
*              for each channel (j = 0, 1, 2 for the 3P spin-orbit levels
*              and j = 10 for the 1D2 level)
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
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
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
*    isrcod:   number of real system dependent variables
*    en1d:     asymptotic energy of the 1D state (cm-1) relative to the center
*              gravity of the 3P sate
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different types of electronic coupling terms
*              this should be 1 here
*    nstate:   number of electronic states included
*              nstate=0:   just 1D state
*              nstate=1:   just 3P state
*              nstate=2:   both 1D and 3P states
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of case(a) interaction potentials actually used
*               this is :  nlam = nlam0 + nlam1
*               if only 1D channels, nlam = 3 (1Sigma+, 1Pi, 1Delta PE curves)
*               if only 3P channels, nlam = 6 (3Sigma- and 3Pi PE curves
*               and two spin-orbit matrix elements Axy and Axz, and the
*               R=inf values of Axy and Azy)
*               if both 1D and 3P channels, nlam= 19 (5 PE curves mentioned above,
*               2 spin-orbit 3P matrix elements, and 5 spin-orbit matrix
*               elements describing 1D-3P coupling, and the R-inf values of
*               7 spin-orbit matrix elements)
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
* ------------------------------------------------------------
      use mod_coeig, only: c0, c1, c2
      use mod_cov2, only: nv2max, junkv => ndummy, v2
      use mod_coiv2, only: iv2
      use mod_cocent, only: cent
      use mod_coeint, only: eint
      use
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"

      common /cosysi/ nscode, isicod, nterm, nstate
      common /cosysr/ isrcod, junkr, en1d
      common /conlam/ nlam, nlammx, lamnum(12)
      common /coered/ ered, rmu
      common /coskip/ nskip, iskip
      common /covvl/  vvl(19)
*   econv is conversion factor from cm-1 to hartrees
*   xmconv is converson factor from amu to atomic units
      common /coconv/ econv, xmconv
*   eigenvectors for the atomic Hamiltonian
      dimension j(1), l(1), jhold(1), ehold(1),
     :          ishold(1), is(1), ieig(0:2)
*  scratch arrays for computing asymmetric top energies and wave fns.
      dimension en0(4), en1(3), en2(2), vec(4,4), work(288)
      zero = 0.d0
      two = 2.d0
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (6, 5)
        write (9, 5)
5       format (' *** FLAGHF = .TRUE. FOR 1D/3P ATOM; ABORT ***' )
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (csflag) then
        write (6, 8)
        write (9, 8)
8      format
     :   ('  *** CSFLAG SET .FALSE. FOR 1D/3P ATOM CALCULATION ***')
        csflag=.false.
      end if
      nsum = 0
      if(nterm.ne.1) then
         write(6,9) nterm
         write(9,9) nterm
9        format(' *** NTERM = ',i3,' .NE. 1 FOR 1D/3P ATOM; ABORT')
         call exit
      end if
      if (nstate .eq. 0) then
        nsum = 3
      else if (nstate .eq. 1) then
        nsum = 6
      else if (nstate .eq. 2) then
        nsum = 19
      else
        write (6, 11) nstate
        write (9, 11) nstate
11      format(' *** NSTATE = ',i3,' .NE. 0, 1, OR 2; ABORT')
        call exit
      endif
      if (nlam .ne. nsum) then
        if (bastst) write (6, 14) nsum
        write (9, 14) nsum
14      format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
        nlam = nsum
      end if
      if (bastst) then
        if (flagsu) then
          write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
          write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16        format(/' **  1D/3P ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4,
     :      ' E=', f7.2, /, '         JTOT=', i5, 2x,' JLPAR=',i2)
        else
          write (6,20) npot,
     :        rmu * xmconv, ered * econv, jtot, jlpar
          write (9,20) npot,
     :        rmu * xmconv, ered * econv, jtot, jlpar
20        format(/,' **  CC 1D/3P ATOM ; NPOT =',i2,' ** RMU=', f9.4,
     :        ' E=', f10.2, /, '         JTOT=', i5, 2x,' JLPAR=',i2)
        end if
        if (.not. flagsu) write (9,30) rcut
30      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      endif
*
*  get 3P and 1D-3P spin-orbit parameters, for calculating asymptotic
*  electrostatic + spin-orbit atomic energies
      call pot (vv0, 8.5d0)
      axy = vvl(6) * econv
      byx = vvl(9) * econv
      bss = byx * sqrt(4.d0/3.d0)
      bxs = byx / sqrt(6.d0)
      bsy = byx / sqrt(2.d0)
      bxd = bsy
*
*  omega=2 atomic Hamiltonian.  basis:  3P2, 1D2 (in order of increasing energy)
      isize = 2
      c2(1,1) = - axy
      c2(1,2) = 2.d0 * bxd
      c2(2,1) = c2(1,2)
      c2(2,2) = en1d
cstart unix-darwin .or. unix-x86
      lwork = 144
      call dsyev('V','L',isize,c2,isize,en2,work,lwork,ierr)
cend
cstart .not. unix-darwin .and. .not. unix-x86
c;      write (6,1225)
c;1225    format('  *** DSYEV NOT ACCESSIBLE IN THIS VERSION OF CODE.')
c;        call exit
cend
*  make sure largest element of e.vectors are positive
      do ii = 1,isize
        imax = 1
        if (abs(c2(2,ii)) .gt. abs(c2(imax,ii))) imax = 2
        if (c2(imax,ii) .lt. 0.d0) then
          do jj = 1,isize
            c2(jj,ii) = -c2(jj,ii)
          enddo
        endif
      enddo
*
*  omega=1 atomic Hamiltonian.  basis:  3P2, 3P1, 1D2 (in order of
*  incrasing energy)
      isize = 3
      c1(1,1) = - axy
      c1(1,2) = 0.d0
      c1(2,1) = c1(1,2)
      c1(1,3) = byx / sqrt(2.d0) + bsy
      c1(3,1) = c1(1,3)
      c1(2,2) = axy
      c1(2,3) = byx / sqrt(2.d0) - bsy
      c1(3,2) = c1(2,3)
      c1(3,3) = en1d
      lwork = 144
      call dsyev('V','L',isize,c1,isize,en1,work,lwork,ierr)
*  make sure largest element of e.vectors are positive
      do ii = 1,isize
        imax = 1
        do ii1 = 2,isize
          if (abs(c1(ii1,ii)) .gt. abs(c1(imax,ii))) imax = ii1
        enddo
        if (c1(imax,ii) .lt. 0.d0) then
          do jj = 1,isize
            c1(jj,ii) = -c1(jj,ii)
          enddo
        endif
      enddo
*
*  omega=0 atomic Hamiltonian.  basis:  3P2, 3P1, 3P0, 1D2 (in order of
*  increasing energy)
      isize = 4

      c0(1,1) = - axy
      c0(1,2) = 0.d0
      c0(2,1) = c0(2,1)
      c0(1,3) = 0.d0
      c0(3,1) = c0(1,3)
      c0(1,4) = bxs * 2.d0 / sqrt(3.d0) + sqrt(2.d0/3.d0) * bss
      c0(4,1) = c0(1,4)
      c0(2,2) = axy
      c0(2,3) = 0.d0
      c0(3,2) = c0(2,3)
      c0(3,3) = 2.d0 * axy
      c0(2,4) = 0.d0
      c0(4,2) = c0(2,4)
      c0(3,4) = bxs * 2.d0 * sqrt(2.d0/3.d0) - bss / sqrt(3.d0)
      c0(4,3) = c0(3,4)
      c0(4,4) = en1d
      lwork = 144
      call dsyev('V','L',isize,c0,isize,en0,work,lwork,ierr)
*  make sure largest element of e.vectors are positive
      do ii = 1,isize
        imax = 1
        do ii1 = 2,isize
          if (abs(c0(ii1,ii)) .gt. abs(c0(imax,ii))) imax = ii1
        enddo
        if (c0(imax,ii) .lt. 0.d0) then
          do jj = 1,isize
            c0(jj,ii) = -c0(jj,ii)
          enddo
        endif
      enddo
*
*  assign quantum numbers and energies
*  zero of energy is energy of lowest (3P2) level
      n=0
      nlevel = 0
* here if 3P state is included
      if (nstate .ge. 1) then
        if (axy .gt. 0.d0) then
* state ordering for p^4 configuration
          j1 = 2
          j2 = 0
          idelj = -1
        else
* state ordering for p^2 configuration
          j1 = 0
          j2 = 2
          idelj = 1
        endif
        do 120 ji = j1, j2, idelj
          lmin=iabs(jtot-ji)
          lmax=jtot+ji
          do 110 li=lmin, lmax
            ix = (-1) ** (li - jtot)
            if (ix .eq. jlpar) then
*  here for correct orbital angular momentum
              n = n + 1
              if (n .gt. nmax) go to 130
              l(n) = li
              cent(n) = li * (li + 1.)
              is(n) = 3
              j(n) = ji
              if (axy .gt. 0.d0) then
*  here for p^4 configuration
                goto (102,104,106),ji + 1
102               eint(n) = (en0(3) - en0(1))/econv
                  goto 108
104               eint(n) = (en0(2) - en0(1))/econv
                  goto 108
106               eint(n) = 0.d0
108             continue
              else
*  here for p^2 configuration
                goto (202,204,206),ji + 1
202               eint(n) = 0.d0
                  goto 208
204               eint(n) = (en0(2) - en0(1))/econv
                  goto 208
206               eint(n) = (en0(3) - en0(1))/econv
208             continue
              endif
            endif
110       continue
          nlevel = nlevel + 1
          if (axy .gt. 0.d0) then
*  here for p^4 configuration
            goto (112,114,116),ji + 1
112           ehold(nlevel) = (en0(3) - en0(1))/econv
              goto 115
114           ehold(nlevel) = (en0(2) - en0(1))/econv
              goto 115
116           ehold(nlevel) = 0.d0
115         continue
          else
            goto (212,214,216),ji + 1
212           ehold(nlevel) = 0.d0
              goto 215
214           ehold(nlevel) = (en0(2) - en0(1))/econv
              goto 215
216           ehold(nlevel) = (en0(3) - en0(1))/econv
215         continue
          endif
*  here for p^2 configuration
          jhold(nlevel) = ji
          ishold(nlevel) = 3
120     continue
      endif
*  here if 1D (j=2) state is included
      if (nstate .ne. 1) then
        ji = 2
        lmin=iabs(jtot - ji)
        lmax=jtot + ji
        do  125 li=lmin,lmax
          ix = (-1) ** (li - jtot)
          if (ix .eq. jlpar) then
*  here for correct orbital angular momentum
            n = n + 1
            if (n .gt. nmax) go to 130
            l(n) = li
            cent(n) = li * (li + 1.)
            is(n) = 1
            j(n) = 2
            eint(n) = (en0(4) - en0(1))/econv
          end if
125     continue
        nlevel=nlevel+1
        ehold(nlevel) = (en0(4) - en0(1))/econv
        jhold(nlevel) = 2
        ishold(nlevel) = 1
      endif
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
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions or for rcut < 0
      if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
        emin = 1.e+7
        do 145  i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is
            if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this
*  condition is met
            end if
          end if
 145    continue
*  now eliminate all channels with eint .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn = 0
          do 150 i = 1, n
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              nn = nn + 1
              eint(nn) = eint(i)
              is(nn) = is(i)
              j(nn) = j(i)
              cent(nn) = cent(i)
              l(nn) = l(i)
            end if
150       continue
*  reset number of channels
          n = nn
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
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
185     format('*** NO OPEN LEVELS IN BA1D3P; ABOST')	
        if (bastst) return
        call exit
      endif
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
      ntop = max(n, nlevop)
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
*
*  now list channels if requested
      if (bastst) then
        write (6,1255)
        write (9,1255)
1255    format(/' ASYMPTOTIC ATOMIC ENERGY LEVELS'
     +    /'   N   S   J      EINT(CM-1)  ',
     +    'OMEGA   EFN: 3P2       3P1       3P0       1D2')
        do 1265 i = 1, nlevel
          write (6,2260) i, ishold(i), jhold(i), ehold(i) * econv
          write (9,2260) i, ishold(i), jhold(i), ehold(i) * econv
2260      format (3i4, f13.3)
          if (ishold(i) .eq. 1) then
*  eigenvalue numbers of 1D state vs. omega
            ieig(2) = 2
            ieig(1) = 3
            ieig(0) = 4
          else
            isj = jhold(i) + 1
            goto (2262,2264,2266), isj
2262          ieig(2) = 0
              ieig(1) = 0
              ieig(0) = 3
              goto 2270
2264          ieig(2) = 0
              ieig(1) = 2
              ieig(0) = 2
              goto 2270
2266          ieig(2) = 1
              ieig(1) = 1
              ieig(0) = 1
2270        continue
          endif
          write (6,2272) (c0(ii,ieig(0)), ii=1,4)
          write (9,2272) (c0(ii,ieig(0)), ii=1,4)
2272      format(32x,'0',6x,4f10.6)
          if (ieig(1) .ne. 0) then
          write (6,2274) (c1(ii,ieig(1)), ii=1,3)
          write (9,2274) (c1(ii,ieig(1)), ii=1,3)
2274      format(32x,'1',6x,2f10.6,10x,f10.6)
          endif
          if (ieig(2) .ne. 0) then
          write (6,2276) (c2(ii,ieig(2)), ii=1,2)
          write (9,2276) (c2(ii,ieig(2)), ii=1,2)
2276      format(32x,'2',6x,f10.6,20x,f10.6)
          endif
1265    continue
        write (6, 255)
        write (9, 255)
255     format(/' CHANNEL BASIS FUNCTIONS'
     +      /'   N   S   J    L      EINT(CM-1)')
        do 265  i = 1, n
          write (6, 260) i, is(i), j(i), l(i), eint(i) * econv
          write (9, 260) i, is(i), j(i), l(i), eint(i) * econv
260       format (3i4, i5, f13.3)
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
* i counts v2 elements
* inum counts v2 elements for given lambda
* ilam counts number of v2 matrices
* ij is address of given v2 element in present v2 matrix
      i = 0
      ilam=0
      do 320 il = 1, 19
        lb = il
        ilam=ilam+1
        inum = 0
        if (nstate.eq.0 .and. lb.gt.3) goto 410
        if (nstate.eq.1) then
          goto (410,410,410,405,405,
     +      405,405,410,410,410,410,410,
     +      405,405,410,410,410,410,410), lb
        endif
405     ij=0
        do 310  icol= 1, n
          do 300  irow = icol, n
            ij = ntop * (icol - 1) +irow
            call vlm1d3p (j(irow), l(irow), is(irow), j(icol),
     :        l(icol), is(icol), jtot, lb, vee)
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
410     if(ilam.gt.nlammx) then
        write(6,311) ilam
311     format(/' ILAM.GT.NLAMMX IN BA1D3P')
        call exit
      end if
      lamnum(ilam) = inum
      if (bastst) then
        write (6, 315) ilam, lamnum(ilam)
        write (9, 315) ilam, lamnum(ilam)
315     format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
      end if
320   continue
      nlam = ilam
      if ( i.gt. nv2max) then
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
      subroutine vlm1d3p (j1, l1, i1, j2, l2, i2, jtot, lb, vee)
* --------------------------------------------------------------------
*  subroutine to evaluate the angular coupling matrix element for
*  rotationally/electronically inelastic collisions of a 1D and/or 3P
*  atom (p^2 or p^4 electron configuration) with a structureless target

*  author:  paul dagdigian
*  current revision date: 23-apr-2013
* --------------------------------------------------------------------
*  variables in call list:
*  j1,l1,i1:    initial electronic orbital, orbital, and spin angular momenta
*  j2,l2,i1:    final electronic orbital, orbital, and spin angular momenta
*  jtot:     total angular momentum
*  lb:       value of expansion index:
*            lb=1:  singlet Sigma potential
*            lb=2:  singlet Pi potential
*            lb=3:  singlet Delta potential
*            lb=4:  triplet Sigma potential
*            lb=5:  triplet Pi potential
*            lb=6:  Axy spin-orbit matrix element (couples 3P states)
*            lb=7:  Azy spin-orbit matrix element (couples 3P states)
*            lb=8:  Bss spin-orbit matrix element (couples 3P and 1D states)
*            lb=9:  Byx spin-orbit matrix element (couples 3P and 1D states)
*           lb=10:  Bxs spin-orbit matrix element (couples 3P and 1D states)
*           lb=11:  Bsy spin-orbit matrix element (couples 3P and 1D states)
*           lb=12:  Bxd spin-orbit matrix element (couples 3P and 1D states)
*           lb=13:  Axy spin-orbit matrix element for large R
*           lb=14:  Azy spin-orbit matrix element for large R
*           lb=15:  Bss spin-orbit matrix element for large R
*           lb=16:  Byx spin-orbit matrix element for large R
*           lb=17:  Bxs spin-orbit matrix element for large R
*           lb=18:  Bsy spin-orbit matrix element for large R
*           lb=19:  Bxd spin-orbit matrix element for large R
*  vee:      on return:  contains desired coupling matrix element
*  subroutines called:
*  xf3j:     evaluates 3j symbol
*  xf6j:     evaluates 6j symbol
* --------------------------------------------------------------------
      use mod_coeig, only: c0, c1, c2
      implicit double precision (a-h,o-z)
      data onth, twth, frth, sqrt2, onsqt3 /0.333333333333333d0,
     +  0.666666666666667d0, 1.333333333333333d0,
     +  1.414213562373095d0, 0.577350269189626d0 /
      vee = 0.d0
      xl1 = l1
      xl2 = l2
      xj1 = j1
      xj2 = j2
      xjtot = jtot
      xl12 = sqrt((2.d0*xl1 + 1.d0) * (2.d0*xl2 + 1.d0))
      xf0 = xf3j(xl1, xj1, xjtot, 0.d0, 0.d0, 0.d0)
     +  * xf3j(xl2, xj2, xjtot, 0.d0, 0.d0, 0.d0)
      xf1 = xf3j(xl1, xj1, xjtot, 0.d0, 1.d0, -1.d0)
     +  * xf3j(xl2, xj2, xjtot, 0.d0, 1.d0, -1.d0)
      xfm1 = xf3j(xl1, xj1, xjtot, 0.d0, -1.d0, 1.d0)
     +  * xf3j(xl2, xj2, xjtot, 0.d0, -1.d0, 1.d0)
      xf2 = xf3j(xl1, xj1, xjtot, 0.d0, 2.d0, -2.d0)
     +  * xf3j(xl2, xj2, xjtot, 0.d0, 2.d0, -2.d0)
      if (i1 .eq. 1) then
*  eigenvalue numbers of 1D state vs. omega
        ieig2_1 = 2
        ieig1_1 = 3
        ieig0_1 = 4
      else
        ij1 = j1 + 1
        goto (1110,1120,1130), ij1
*  eigenvalue numbers of 3PJ states vs. omega
1110      ieig2_1 = 0
          ieig1_1 = 0
          ieig0_1 = 3
          goto 1140
1120      ieig2_1 = 0
          ieig1_1 = 2
          ieig0_1 = 2
          goto 1140
1130      ieig2_1 = 1
          ieig1_1 = 1
          ieig0_1 = 1
1140    continue
      endif
      if (i2 .eq. 1) then
*  eigenvalue numbers of 1D state vs. omega
        ieig2_2 = 2
        ieig1_2 = 3
        ieig0_2 = 4
      else
        ij2 = j2 + 1
        goto (1210,1220,1230), ij2
*  eigenvalue numbers of 3PJ states vs. omega
1210      ieig2_2 = 0
          ieig1_2 = 0
          ieig0_2 = 3
          goto 1240
1220      ieig2_2 = 0
          ieig1_2 = 2
          ieig0_2 = 2
          goto 1240
1230      ieig2_2 = 1
          ieig1_2 = 1
          ieig0_2 = 1
1240    continue
      endif
      sign = 1.d0
      goto (10,20,30,40,50,60,70,80,90,100,110,120,
     +  130,140,150,160,170,180,190),lb
*
*  1Sigma+ PE curve
10    om0 = c0(4,ieig0_1) * c0(4,ieig0_2)
      vee = xl12 * xf0 * om0
      goto 1000
*
*  1Pi PE curve
20    if (ieig1_1 .ne. 0 .and. ieig1_2 .ne. 0) then
        om1 = c1(3,ieig1_1) * c1(3,ieig1_2)
*  factor of 2 in next line and below takes account of the degeneracy of
*  the omega > 0 levels
        vee = 2.d0 * xl12 * xf1 * om1
        goto 1000
      endif
      goto 1000
*
*  1Delta PE curve
30    if (ieig2_1 .ne. 0 .and. ieig2_2 .ne. 0) then
        om2 = c2(2,ieig2_1) * c2(2,ieig2_2)
        vee = 2.d0 * xl12 * xf2 * om2
        goto 1000
      endif
      goto 1000
*
*  3Sigma- PE curve
*
*   3P0 <-> 3P0
40    om0 = c0(3,ieig0_1) * c0(3,ieig0_2)
      vee = onth * xl12 * xf0 * om0
*   3P1 <-> 3P1
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(2,ieig1_1) * c1(2,ieig1_2)
        vee = vee + 2.d0 * 0.5d0 * xl12 * xf1 * om1
      endif
*   3P2 <-> 3P2
      om0 = c0(1,ieig0_1) * c0(1,ieig0_2)
      vee = vee + twth * xl12 * xf0 * om0
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(1,ieig1_1) * c1(1,ieig1_2)
        vee = vee + 2.d0 * 0.5d0 * xl12 *xf1 * om1
      endif
*   3P2 <-> 3P1
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(1,ieig1_1) * c1(2,ieig1_2)
        vee = vee - 0.5d0 * xl12 * (xf1 - xfm1) * om1
      endif
*   3P1 <-> 3P2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(2,ieig1_1) * c1(1,ieig1_2)
        vee = vee - 0.5d0 * xl12 * (xf1 - xfm1) * om1
      endif
*   3P2 <-> 3P0
      om0 = c0(1,ieig0_1) * c0(3,ieig0_2)
      vee = vee - (sqrt2/3.d0) * xl12 * xf0 * om0
*   3P0 <-> 3P2
      om0 = c0(3,ieig0_1) * c0(1,ieig0_2)
      vee = vee - (sqrt2/3.d0) * xl12 * xf0 * om0
      goto 1000
*
*  3Pi PE curve
*
*   3P0 <-> 3P0
50    om0 = c0(3,ieig0_1) * c0(3,ieig0_2)
      vee = twth * xl12 * xf0 * om0
*   3P1 <-> 3P1
      om0 = c0(2,ieig1_1) * c0(2,ieig1_2)
      vee = vee + xl12 * xf0 * om0
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(2,ieig1_1) * c1(2,ieig1_2)
        vee = vee + 2.d0 * 0.5d0 * xl12 * xf1 * om1
      endif
*   3P2 <-> 3P2
      om0 = c0(1,ieig0_1) * c0(1,ieig0_2)
      vee = vee + onth * xl12 * xf0 * om0
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(1,ieig1_1) * c1(1,ieig1_2)
        vee = vee + 2.d0 * 0.5d0 * xl12 * xf1 * om1
      endif
      if (ieig2_1.ne.0 .and. ieig2_2.ne.0) then
        om2 = c2(1,ieig2_1) * c2(1,ieig2_2)
        vee = vee + 2.d0 * xl12 * xf2 * om2
      endif
*   3P2 <-> 3P1
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(1,ieig1_1) * c1(2,ieig1_2)
        vee = vee + 2.d0 * 0.5d0 * xl12 * xf1 * om1
      endif
*   3P1 <-> 3P2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(2,ieig1_1) * c1(1,ieig1_2)
        vee = vee + 2.d0 * 0.5d0 * xl12 * xf1 * om1
      endif
*   3P2 <-> 3P0
      om0 = c0(1,ieig0_1) * c0(3,ieig0_2)
      vee = vee + (sqrt2/3.d0) * xl12 * xf0 * om0
*   3P0 <-> 3P2
      om0 = c0(3,ieig0_1) * c0(1,ieig0_2)
      vee = vee + (sqrt2/3.d0) * xl12 * xf0 * om0
      goto 1000
*
*  Axy spin-orbit matrix element
60    continue
*   3P0 <-> 3P0
      om0 = c0(3,ieig0_1) * c0(3,ieig0_2)
      vee = twth * xl12 * xf0 * om0
*   3P1 <-> 3P1
      om0 = c0(2,ieig1_1) * c0(2,ieig1_2)
      vee = vee + xl12 * xf0 * om0
*   3P2 <-> 3P2
      om0 = c0(1,ieig0_1) * c0(1,ieig0_2)
      vee = vee + onth * xl12 * xf0 * om0
      if (ieig2_1.ne.0 .and. ieig2_2.ne.0) then
        om2 = c2(1,ieig2_1) * c2(1,ieig2_2)
        vee = vee - 2.d0 * xl12 *xf2 * om2
      endif
*   3P2 <-> 3P0
      om0 = c0(1,ieig0_1) * c0(3,ieig0_2)
      vee = vee + (sqrt2/3.d0) * xl12 * xf0 * om0
*   3P0 <-> 3P2
      om0 = c0(3,ieig0_1) * c0(1,ieig0_2)
      vee = vee + (sqrt2/3.d0) * xl12 * xf0 * om0
      vee = vee * sign
      goto 1000
*  Axy for large R
130   sign = -1.d0
      goto 60
*
*  Azy spin-orbit matrix element
70    continue
*   3P0 <-> 3P0
      om0 = c0(3,ieig0_1) * c0(3,ieig0_2)
      vee = frth * sqrt2 * xl12 * xf0 * om0
*   3P1 <-> 3P1
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(2,ieig1_1) * c1(2,ieig1_2)
        vee = vee + 2.d0 * sqrt2 * xl12 * xf1 * om1
      endif
*   3P2 <-> 3P2
      om0 = c0(1,ieig0_1) * c0(1,ieig0_2)
      vee = vee - frth * sqrt2 * xl12 * xf0 * om0
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(1,ieig1_1) * c1(1,ieig1_2)
        vee = vee - 2.d0 * sqrt2 * xl12 * xf1 * om1
      endif
*   3P2 <-> 3P0
      om0 = c0(1,ieig0_1) * c0(3,ieig0_2)
      vee = vee - twth * xl12 * xf0 * om0
*   3P0 <-> 3P2
      om0 = c0(3,ieig0_1) * c0(1,ieig0_2)
      vee = vee - twth * xl12 * xf0 * om0
      vee = vee * sign
      goto 1000
*  Azy for large R
140   sign = -1.d0
      goto 70
*
*  Bss spin-orbit matrix element
80    continue
*   3P0 <-> 1D2
      om0 = c0(3,ieig0_1) * c0(4,ieig0_2)
      vee = - onsqt3 * xl12 * xf0 * om0
*   1D2 <-> 3P0
      om0 = c0(4,ieig0_1) * c0(3,ieig0_2)
      vee = vee - onsqt3 * xl12 * xf0 * om0
*   3P2 <-> 1D2
      om0 = c0(1,ieig0_1) * c0(4,ieig0_2)
      vee = vee + sqrt2 * onsqt3 * xl12 * xf0 * om0
*   1D2 <-> 3P2
      om0 = c0(4,ieig0_1) * c0(1,ieig0_2)
      vee = vee + sqrt2 * onsqt3 * xl12 * xf0 * om0
      vee = vee * sign
      goto 1000
*  Bss for large R
150   sign = -1.d0
      goto 80
*
*  Byx spin-orbit matrix element
90    continue
*   3P1 <-> 1D2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(2,ieig1_1) * c1(3,ieig1_2)
        vee = (2.0d0 / sqrt2) * xl12 * xf1 * om1
      endif
*   1D2 <-> 3P1
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(3,ieig1_1) * c1(2,ieig1_2)
        vee = vee + (2.0d0 / sqrt2) * xl12 * xf1 * om1
      endif
*   3P2 <-> 1D2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(1,ieig1_1) * c1(3,ieig1_2)
        vee = vee + (2.0d0 / sqrt2) * xl12 * xf1 * om1
      endif
*   1D2 <-> 3P2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(3,ieig1_1) * c1(1,ieig1_2)
        vee = vee + (2.0d0 / sqrt2) * xl12 * xf1 * om1
      endif
      vee = vee * sign
      goto 1000
*  Byx for large R
160   sign = -1.d0
      goto 90
*
*  Bxs spin-orbit matrix element
100   continue
*   3P0 <-> 1D2
      om0 = c0(3,ieig0_1) * c0(4,ieig0_2)
      vee = 2.d0 * sqrt2 * onsqt3 * xl12 * xf0 *om0
*   1D2 <-> 3P0
      om0 = c0(4,ieig0_1) * c0(3,ieig0_2)
      vee = vee + 2.d0 * sqrt2 * onsqt3 * xl12 * xf0 *om0
*   3P2 <-> 1D2
      om0 = c0(1,ieig0_1) * c0(4,ieig0_2)
      vee = vee + 2.d0 * onsqt3 * xl12 * xf0 * om0
*   1D2 <-> 3P2
      om0 = c0(4,ieig0_1) * c0(1,ieig0_2)
      vee = vee + 2.d0 * onsqt3 * xl12 * xf0 * om0
      vee = vee * sign
      goto 1000
*  Bxs for large R
170   sign = -1.d0
      goto 100
*
*  Bsy spin-orbit matrix element
110   continue
*   3P1 <-> 1D2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(2,ieig1_1) * c1(3,ieig1_2)
        vee = - 2.0d0 * xl12 * xf1 * om1
      endif
*   1D2 <-> 3P1
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(3,ieig1_1) * c1(2,ieig1_2)
        vee = vee - 2.0d0 * xl12 * xf1 * om1
      endif
*   3P2 <-> 1D2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(1,ieig1_1) * c1(3,ieig1_2)
        vee = vee + 2.d0 *xl12 * xf1 * om1
      endif
*   1D2 <-> 3P2
      if (ieig1_1.ne.0 .and. ieig1_2.ne.0) then
        om1 = c1(3,ieig1_1) * c1(1,ieig1_2)
        vee = vee + 2.d0 *xl12 * xf1 * om1
      endif
      vee = vee * sign
      goto 1000
*  Bsy for large R
180   sign = -1.d0
      goto 110
*
*  Bxd spin-orbit matrix element
120   continue
*   1D2 <-> 3P2
      if (ieig2_1.ne.0 .and. ieig2_2.ne.0) then
        om2 = c2(2,ieig2_1) * c2(1,ieig2_2)
        vee = 2.d0 * 2.d0 * xl12 * xf2 * om2
      endif
*   3P2 <-> 1D2
      if (ieig2_1.ne.0 .and. ieig2_2.ne.0) then
        om2 = c2(1,ieig2_1) * c2(2,ieig2_2)
        vee = vee + 2.d0 * 2.d0 * xl12 * xf2 * om2
      endif
      vee = vee * sign
      goto 1000
*  Bxd for large R
190   sign = -1.d0
      goto 120
*
1000  return
      end
* ------------------------------eof-----------------------------------
