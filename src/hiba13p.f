* --------------------------------------------------
      subroutine ba13p (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential
*  for collision of a singlet and triplet atom with a structureless atom
*  authors:  brigitte pouilly and millard alexander
*  current revision date:  5-apr-2003 by mha
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains electronic angular momentum quantum number
*              for each channel
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
*    sc1,sc2:  scratch vectors of length at least nmax
*    sc3,sc4:  scratch vectors of length at least nmax
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
*    en(1):    asymptotic energy of j=0 fine-structure level of triplet
*    en(2):    asymptotic energy of j=1 fine-structure level of triplet
*    en(3):    asymptotic energy of j=2 fine-structure level of triplet
*    en(4):    asymptotic energy of singlet state (cm-1)
*    de(4):    morse De for 3Pi, 3Sig, 1Pi, 1Sig
*    re(4):    morse re for 3Pi, 3Sig, 1Pi, 1Sig
*    be(4):    morse beta for 3Pi, 3Sig, 1Pi, 1Sig
*    rl(4):    msv RL for 3Pi, 3Sig, 1Pi, 1Sig
*    cl(4):    msv C6 for 3Pi, 3Sig, 1Pi, 1Sig
*    cmix:     mixing coefficient of J=1 singlet and triplet levels
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different types of electronic coupling terms
*              this should be 1 here
*    nstate:   number of electronic states included
*              nstate=0:   just singlet state
*              nstate=1:   just triplet state
*              nstate=2:   both singlet and triplet state
*    ipol:     =0 if no polarization state preparation
*    ipol:     =1 if polaratization state preparation
*    npot:     potential designator (1-4 currently allowed)
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
*               if only singlet channels or only triplet channels, nlam= 2
*               if both singlet and triplet channels, nlam= 4
*    nlammx:    the maximum number of angular coupling terms
*  variables in common block /cosgpi/
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*    nlm(1):    the total number of terms in the expansion of the v-sig
*               potential (nlamsg)
*    nlm(2):    the total number of terms in the expansion of the v-pi
*               potential (nlampi)
*    nlm(3):    the total number of terms in the expansion of the v1
*               potential (nlam1)
*    nlm(4)     the total number of terms in the expansion of the v2
*               potential (nlam2)
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
      use mod_cov2, only: nv2max, junkv => ndummy, v2
      use mod_coiv2, only: iv2
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"

      common /cosysi/ nscode, isicod, nterm, nstate,ipol, npot
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coskip/ nskip, iskip
      common /coconv/ econv, xmconv
      dimension j(1), l(1), jhold(1), ehold(1), sc1(1), sc2(1), sc3(1),
     :          sc4(1), ishold(1), is(1)
*   econv is conversion factor from cm-1 to hartrees
*   xmconv is converson factor from amu to atomic units
      zero = 0.d0
      two = 2.d0
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (6, 5)
        write (9, 5)
5       format (' *** FLAGHF = .TRUE. FOR 1/3 ATOM; ABORT ***' )
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (jlpar .eq. -1) then
        if (ipol .ne. 1 .and. ipol .ne. 0) then
           write (6, 6) ipol
           write (9, 6) ipol
6          format(' ** IPOL =',i3, ' SET EQUAL TO ZERO IN BA13P')
           ipol=0
        endif
      else
        if (ipol .ne. 0) then
           write (6, 7) ipol
           write (9, 7) ipol
7          format(' ** IPOL =',i3, ' SET EQUAL TO ZERO IN BA13P')
           ipol=0
        endif
      endif
      if (csflag) then
        write (6, 8)
        write (9, 8)
8      format
     :   ('  *** CSFLAG SET .FALSE. FOR 1/3 ATOM CALCULATION ***')
        csflag=.false.
      end if
      nsum = 0
      if(nterm.ne.1) then
         write(6,9) nterm
         write(9,9) nterm
9        format(' *** NTERM = ',i3,' .NE. 1 FOR 1/3 ATOM; ABORT')
         call exit
      end if
      if (nstate .eq. 0) then
        nsum=2
      else if (nstate .eq. 1) then
        nsum=2
      else if (nstate .eq. 2) then
        nsum=4
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
      if (clist) then
        if (flagsu) then
          write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
          write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16        format(/' **  1/3 ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4,
     :      '             E=', f7.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)
        else
            write (6,20) npot,
     :          rmu * xmconv, ered * econv, jtot, jlpar
            write (9,20) npot,
     :          rmu * xmconv, ered * econv, jtot, jlpar
20          format(/,' **  CC 1/3 ATOM ; NPOT =',i2,' ** RMU=', f9.4,
     :           '       E=',f7.2,'   JTOT=', i5, 2x,' JLPAR=',i2)
        end if
        if (.not. flagsu) write (9,30) rcut
30      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      endif
*  assign quantum numbers and energies
      n=0
      nlevel = 0
* here if triplet state is included
      if (nstate .ge. 1) then
        jmin=0
        jmax=2
        do 120 ji=jmin, jmax
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
            is(n) = 1
            j(n) = ji
            eint(n) = en(ji+1)/econv
          end if
110       continue
          nlevel = nlevel + 1
          ehold(nlevel) = en(ji+1)/econv
          jhold(nlevel) = ji
          ishold(nlevel) = 1
120     continue
      endif
*  here if singlet state is included
      if (nstate .ne. 1) then
        lmin=iabs(jtot-1)
        lmax=jtot+1
        do  125 li=lmin,lmax
          ix = (-1) ** (li - jtot)
          if (ix .eq. jlpar) then
*  here for correct orbital angular momentum
            n = n + 1
            if (n .gt. nmax) go to 130
            l(n) = li
            cent(n) = li * (li + 1.)
            if (ipol .eq. 0) then
              is(n) = 0
            else if (ipol.eq.1) then
              if (li .lt. jtot) is(n)=-2
              if (li .gt. jtot) is(n)=2
            endif
            j(n) = 1
            eint(n) = en(4)/econv
          end if
125     continue
        if (ipol .eq. 0) then
          nlevel=nlevel+1
          ehold(nlevel) = en(4)/econv
          jhold(nlevel) = 1
          ishold(nlevel) = 0
        else if (ipol .eq. 1) then
          do 126 ii = -2,2,4
            nlevel=nlevel+1
            ehold(nlevel) = en(4)/econv
            jhold(nlevel) = 1
            ishold(nlevel) = ii
126       continue
        endif
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
      nlevop = nlevel
      ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1

*  now list channels if requested
      if (clist) then
        write (6, 255)
        write (9, 255)
255     format(/'   N   S   J    L      EINT(CM-1)')
        do 265  i = 1, n
          write (6, 260) i, is(i), j(i), l(i), eint(i) * econv
          write (9, 260) i, is(i), j(i), l(i), eint(i) * econv
260       format (3i4, i5, f13.3)
265     continue
      end if
*  now calculate coupling matrix elements
      if (bastst) then
        write (6, 280)
        write (9, 280)
280     format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
      end if
* i counts v2 elements
* inum counts v2 elements for given lambda
* ilam counts number of v2 matrices
* ij is address of given v2 element in present v2 matrix
      i = 0
      ilam=0
      do 320 il = 0, 6, 2
        lb = il
        ilam=ilam+1
        inum = 0
        ij=0
        do 310  icol= 1, n
          do 300  irow = icol, n
            ij = ntop * (icol - 1) +irow
              call vlm13p (j(irow), l(irow), is(irow), j(icol),
     :                     l(icol), is(icol), jtot, lb, cmix, vee)
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
      if(ilam.gt.nlammx) then
        write(6,311) ilam
311     format(/' ILAM.GT.NLAMMX IN BA13P')
        call exit
      end if
      lamnum(ilam) = inum
      if (bastst) then
        write (6, 315) ilam, lamnum(ilam)
        write (9, 315) ilam, lamnum(ilam)
315     format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
      end if
320   continue
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
      if (clist) then
        write (6, 360) i
        write (9, 360) i
360     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i6)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine vlm13p (j1, l1, i1, j2, l2, i2, jtot, lb, cmix, vee)
* --------------------------------------------------------------------
*  subroutine to evaluate the angular coupling matrix element for rotationally
*  inelastic collisions of a singlet and/or triplet P atom with a
*  structureless target

*  authors:  brigitte pouilly and millard alexander
*  current revision date: 23-jul-90
* --------------------------------------------------------------------
*  variables in call list:
*  j1,l1,i1:    initial electronic orbital, orbital, and spin angular momenta
*  j2,l2,i1:    final electronic orbital, orbital, and spin angular momenta
*  jtot:     total angular momentum
*  lb:       value of expansion index:
*            lb=0:  triplet Pi potential
*            lb=2:  triplet Sigma potential
*            lb=4:  singlet Pi potential
*            lb=6:  singlet Sigma potential
*  cmix:     cosine of mixing angle of singlet and triplet states
*  v:        on return:  contains desired coupling matrix element
*  subroutines called:
*  xf3j:     evaluates 3j symbol
*  xf6j:     evaluates 6j symbol
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension vterm(4),cs(4),sn(4)
      data  zero, one, two, three /0.d0, 1.d0, 2.d0, 3.d0/
      if (i1 .eq. 1) then
        ind1=j1+1
      else
        ind1=4
      endif
      if (i2 .eq. 1) then
        ind2=j2+1
      else
        ind2=4
      endif
      do  10 i=0,6,2
        is=1
        inu=i
        if (i .gt. 2) then
          is=0
          inu=i-4
        endif
        xnu=inu
        xjtot=jtot
        xl1=l1
        xj1=j1
        xs=is
        xl2=l2
        xj2=j2

        phase=(-1)**(jtot - is)
        xnorm=three*sqrt((two*l1+one)*(two*l2+one)*
     :                   (two*j1+one)*(two*j2+one))
        term1=xf3j(one,one,xnu,zero,zero,zero)
        term2=xf3j(xl1,xl2,xnu,zero,zero,zero)
        term3=xf6j(xj1,xj2,xnu,xl2,xl1,xjtot)
        term4=xf6j(xj1,xj2,xnu,one,one,xs)
        vterm(i/2+1)=phase*xnorm*term1*term2*term3*term4
10    continue
      smix=sqrt(one-cmix*cmix)
      do 20 i=1,4
        if (i .eq. 1 .or. i .eq. 3) then
          cs(i)=zero
          sn(i)=one
        else if (i .eq. 2) then
          cs(i)=-smix
          sn(i)=cmix
        else if (i .eq. 4) then
          cs(i)=cmix
          sn(i)=smix
        endif

20    continue
* cs is the singlet amplitude in each channel
* sn is the triplet amplitude in each channel
* note that only the j=1 channels are mixed
      if (lb .le. 2) then
* here for 3Pi or 3Sigma coupling
        vee=vterm(lb/2+1)*sn(ind1)*sn(ind2)
      else
        vee=vterm(lb/2+1)*cs(ind1)*cs(ind2)
      endif

      return
      end
* -----------------------------------------
      subroutine tcasea(j,jlpar)
* -----------------------------------------
*   matrix to transform from case a to case e, solely for 13P scattering
*   this matrix, tatoe, for total-j=j is returned
*   author:  thierry duhoo and millard alexander
*   latest revision date:  30-dec-1995
* -----------------------------------------
      implicit double precision (a-h,o-z)
      common /cotrans/ ttrans
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix
      dimension tatoe(6,6), cmat(6,6), ttrans(6,6)
      data zero, one,two ,three/0.d0, 1.d0, 2.d0, 3.d0/
      if (j .lt. 2) then
* error if jtot is less than 2
        write (6, 5) jtot
5       format (' *** JTOT = ',i2,' .LT. 2 IN TCASEA; ABORT **')
        stop
      endif
      if (cmix .le. 0.8d0) then
        write (6, 8) cmix
8       format (' *** WARNING; CMIX =',f6.3,'; POSSIBLY TOO SMALL **')
      endif
      smix=sqrt(one-cmix*cmix)
* initialization of the matrix tatoe
      call dset(36,zero,tatoe,1)
*
* initialize matrix which transforms from case e to spin-orbit diagonalized
* basis
      call dscal(36,zero,cmat,1)
      xj=j
      if(jlpar.gt.0) then
* here for f-labelled states
*       singlet state (assumed state 6) is not coupled
*
        tatoe(1,1)=-sqrt(1.d0/3.d0)
        tatoe(1,3)=sqrt(2.d0/3.d0)
        tatoe(2,2)=-sqrt(1.d0/2.d0)
        tatoe(2,4)=sqrt(1.d0/2.d0)
        denom=(two*xj+one)*(two*xj-one)
        tatoe(3,1)=sqrt(xj*(xj-one)/denom)
        tatoe(3,2)=sqrt((xj+one)*(xj-one)/denom)
        tatoe(3,3)=sqrt(xj*(xj-one)/(two*denom))
        tatoe(3,4)=sqrt((xj+one)*(xj-one)/denom)
        tatoe(3,5)=sqrt((xj+one)*(xj+two)/(two*denom))
        denom=(two*xj+three)*(two*xj-one)
        tatoe(4,1)=-sqrt((two*xj*(xj+one))/(three*denom))
        tatoe(4,2)=-sqrt(three/(two*denom))
        tatoe(4,3)=-sqrt((xj*(xj+one))/(three*denom))
        tatoe(4,4)=-sqrt(three/(two*denom))
        tatoe(4,5)=sqrt((three*(xj+two)*(xj-one))/(denom))
        denom=(two*xj+three)*(two*xj+one)
        tatoe(5,1)=sqrt(((xj+two)*(xj+one))/(denom))
        tatoe(5,2)=-sqrt((xj*(xj+two))/(denom))
        tatoe(5,3)=sqrt(((xj+two)*(xj+one))/(two*denom))
        tatoe(5,4)=-sqrt((xj*(xj+two))/(denom))
        tatoe(5,5)=sqrt((xj*(xj-one))/(two*denom))
        tatoe(6,6)=one
* now matrix which transforms from case e to spin-orbit diagonalized
* basis
        cmat(1,1)=one
        cmat(2,2)=cmix
        cmat(2,6)=-smix
        cmat(6,2)=smix
        cmat(3,3)=one
        cmat(4,4)=one
        cmat(5,5)=one
        cmat(6,6)=cmix
* here for e-labelled states
      else
        denom=two*(two*xj+one)
        tatoe(1,1)=-sqrt((xj+one)/denom)
        tatoe(1,2)=sqrt((two*xj)/denom)
        tatoe(1,3)=sqrt((xj+one)/denom)
        tatoe(2,1)=-sqrt((xj)/denom)
        tatoe(2,2)=-sqrt((two*(xj+one))/denom)
        tatoe(2,3)=sqrt((xj)/denom)
        tatoe(3,1)=sqrt((xj-one)/denom)
        tatoe(3,3)=sqrt((xj-one)/denom)
        tatoe(3,4)=sqrt(two*(xj+two)/denom)
        tatoe(4,1)=-sqrt((xj+two)/denom)
        tatoe(4,3)=-sqrt((xj+two)/denom)
        tatoe(4,4)=sqrt((two*(xj-one))/denom)
        denom=two*xj+one
        tatoe(5,5)=sqrt(xj/denom)
        tatoe(5,6)=sqrt((xj+one)/denom)
        tatoe(6,5)=-sqrt((xj+one)/denom)
        tatoe(6,6)=sqrt(xj/denom)
* now matrix which transforms from case e to spin-orbit diagonalized
* basis
        cmat(1,1)=cmix
        cmat(1,5)=-smix
        cmat(5,1)=smix
        cmat(2,2)=cmix
        cmat(2,6)=-smix
        cmat(6,2)=smix
        cmat(3,3)=one
        cmat(4,4)=one
        cmat(5,5)=cmix
        cmat(6,6)=cmix
      endif
      call mxma(cmat,1,6,tatoe,1,6,ttrans,1,6,6,6,6)
*      call mxoutd(94, cmat, 6,6,0,.false.)
*      call mxoutd(94, tatoe, 6,6,0,.false.)
*      call mxoutd(94, ttrans, 6,6,0,.false.)
*      call exit
      return
      end
