* --------------------------------------------------
      subroutine bah3p (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  isc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential
*  for collision of a triplet atom in a P state and a homonuclear molecule
*  authors:  millard alexander
*
*  moved j12 (renamed isc1 array) to common block /coj12/ to pass array
*  to other subrs (p.dagdigian)
*
*  revised routine so that CC scattering calculations work.  it appears that
*  this routine had previously been used only for bound-state calculations
*
*  current revision date:  27-mar-2013 by p. dagdigian
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains diatom angular momentum quantum number
*              for each channel (jmol)
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains atomic angular momentum (ja) of each channel
*              (0, 1, or 2)
*    jhold:    on return contains diatom angular momentum quantum number
*              for each level
*    ehold:    on return contains energy in hartrees of each level
*    ishold:   on return contains atomic angular momentum of each level
*    nlevel:   on return contains number of energetically distinct
*              levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              levels used in channel basis which are open
*              asymptotically
*    isc1,sc2: scratch vectors of length at least nmax
*    sc3,sc4:  scratch vectors of length at least nmax
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*    flaghf:   if .true., then system with half-integer spin (this is not the
*              case here)
*              if .false., then system with integer spin
*    flagsu:   if .true., then molecule-surface collisons (not applicable)
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then Dubernet and Hutson case 1C
*              if .false., then Dubernet and Hutson case 1A
*    nu:       coupled-states projection index (called P by dubernet and
*              hutson
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(l)=jlpar
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
*    brot:     rotational constant of molecule
*    aso1:     spin-orbit constant of 3P atom j=1, relative to j=0 (cm-1)
*    aso2:     spin-orbit constant of 3P atom j=2, relative to j=0 (cm-1)
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different types of electronic coupling terms
*              this should be 1 here
*    iop:      ortho (iop=1) or para (iop=0)
*    jmax:     maximum rotational quantum number of diatomic
*  variable in common block /cocent/
*    cent:     array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:     array containing channel energies (in hartree)
*  variable in common block /coj12/
*    j12:      array containing vector sum of ja + j (similar situation
*              with molecule-molecule scattering, see hibastpln basis
*              routine
*  variables in common block /coered/
*    ered:     collision energy in atomic units (hartrees)
*    rmu:      collision reduced mass in atomic units
*              (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:     the number of case(a) interaction potentials actually used
*              this is 5 here
*    nlammx:   the maximum number of angular coupling terms
*  variable in common block /cov2/
*    nv2max:   maximum core memory allocated for the v2 matrix
*    v2:       lower triangle of nonzero angular coupling matrix elements
*              only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:       matrix address of v2 matrix for each non-zero element
*  variable in common block /coconv/
*     econv:   conversion factor from cm-1 to hartrees
*     xmconv:  converson factor from amu to atomic units
*  subroutines called:
*   vlmh3p:    returns angular coupling coefficient for particular
*              choice of channel index
* ------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysr/ isrcod, junkr, brot,aso1, aso2
      common /cosysi/ nscode, isicod, nterm, iop,jmax
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /coipar/ iiipar(9), iprint
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(5)
      common /coeint/ eint(5)
      common /coj12/  j12(1)
      common /coered/ ered, rmu
      common /coskip/ nskip, iskip
      common /coconv/ econv, xmconv
      common /cojtot/ jjtot,jjlpar
      dimension j(9), l(9), jhold(9), ehold(9), isc1(9), sc2(9), sc3(9),
     :          sc4(9), ishold(9), is(9)
      dimension lamr(9),lama(9),lam12(9), mu(9)
      data lamr /2,0,2,2,2,4,4,4,4/
      data lama /2,2,0,2,2,0,2,2,2/
      data lam12 /0,2,2,2,4,4,2,4,6/
      data mu / 0,0,0,2,1,0,0,2,1/
*   econv is conversion factor from cm-1 to hartrees
*   xmconv is converson factor from amu to atomic units
      zero = 0.d0
      half=0.5d0
      two = 2.d0
*  check for consistency in the values of flaghf and csflag
      if (flaghf .and. (.not.csflag .or. csflag .and. ihomo))  then
        write (6, 5)
        write (9, 5)
5       format
     : (' *** FLAGHF = .TRUE. FOR CC OR CS1C MOL+ 3P ATOM; ABORT ***' )
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (flaghf) then
        xnu=nu+half
      else
        xnu=nu
      endif
      if(nterm.ne.1) then
        write(6,9) nterm
        write(9,9) nterm
9       format(' *** NTERM = ',i3,' .NE. 1 FOR MOL + 3P ATOM; ABORT')
        call exit
      end if
      nsum=5
      if (bastst) write (6, 14) nsum
      write (9, 14) nsum
14    format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
      nlam = nsum
      if (bastst) then
        if (flagsu) then
          write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
          write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16        format(/
     :  ' **  MOL + 3P ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4,
     :  '             E=', f9.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)
        else
          if(.not.csflag) then
            write (6,20)
     :          rmu * xmconv, ered * econv, jtot, jlpar
            write (9,20)
     :          rmu * xmconv, ered * econv, jtot, jlpar
20          format(/,' **  CC MOL + 3P ATOM',' ** RMU=', f9.4,
     :           '       E=',f9.2,'   JTOT=', i5, 2x,' JLPAR=',i2)
          else if (csflag .and. ihomo) then
            write (6,22)
     :          rmu * xmconv, ered * econv, jtot, xnu
            write (9,22)
     :          rmu * xmconv, ered * econv, jtot, xnu
22          format(/,' **  CD 1C MOL + 3P ATOM',' ** RMU=', f9.4,
     :           '       E=',f9.2,'   JTOT=', i5, 2x,' P=',f4.1)
          else if (csflag .and. .not.ihomo) then
            write (6,23)
     :          rmu * xmconv, ered * econv, jtot, xnu
            write (9,23)
     :          rmu * xmconv, ered * econv, jtot, xnu
23          format(/,' **  CD 1A MOL + 3P ATOM',' ** RMU=', f9.4,
     :           '       E=',f9.2,'   JTOT=', i5, 2x,' P=',f7.2)
          endif
        end if
        if (.not. flagsu) write (9,30) rcut
30      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      endif
* check for consistency in projection quantum number
      if (nu .gt. jtot .and. csflag) then
        write (6, 35) nu, jtot
        write (9, 35) nu, jtot
35      format (' NU =',i3,' .GT. JTOT =',i3,'; JTOT RESET TO NU')
        jtot=nu
      endif

*  save jtot and jlpar for use later in transformatin
      jjtot=jtot
      jjlpar=jlpar
      if (iop .eq. 0) jmolmin=0
      if (iop .eq. 1) jmolmin=1
      if (.not. csflag) then
*  assign quantum numbers and energies for CC calculations
        n=0
        nlevel = 0
* sum over molecular levels
        do 105 jmol=jmolmin,jmax,2
        erot = brot*jmol*(jmol + 1)
* sum over atomic levels
        do 100 jai=0,2
          nlevel = nlevel+1
          if (jai .eq. 0) ehold(nlevel)=zero
          if (jai .eq. 1) ehold(nlevel)=+aso1/econv
          if (jai .eq. 2) ehold(nlevel)=aso2/econv
          ehold(nlevel) = ehold(nlevel) + erot/econv
          ishold(nlevel)=jai
          jhold(nlevel)=jmol
          xjai=jai
          j12min=abs(xjai-jmol)
          j12max=jai+jmol
          do 90 j12i=j12min,j12max
            lmin=iabs(jtot-j12i)
            lmax=jtot+j12i
            do 80 li=lmin,lmax
              if ((-1)**li .eq. jlpar*(-1)**(jmolmin-1)) then
                n = n + 1
                if (n .gt. nmax) go to 180
                l(n) = li
                is(n)=jai
                j(n)=jmol
                j12(n)=j12i
* centrifugal contribution to hamiltonian
                cent(n) = li*(li+1)
* constant channel energy
                sc2(n) = ehold(nlevel)
              endif
80          continue
90        continue
100     continue
105     continue
*  sort energies in order ot increasing energy
        if (n .gt. 1) then
          do 520 i1 = 1, n-1
            esave = sc2(i1)
            do 530 i2 = i1+1, n
              if (sc2(i2) .lt. esave) then
*  level i2 has lower energy than level i1, switch them
                esave = sc2(i2)
                sc2(i2) = sc2(i1)
                sc2(i1) = esave
                jsave = j(i2)
                j(i2) = j(i1)
                j(i1) = jsave
                jsave = is(i2)
                is(i2) = is(i1)
                is(i1) = jsave
                jsave = l(i2)
                l(i2) = l(i1)
                l(i1) = jsave
                jsave = j12(i2)
                j12(i2) = j12(i1)
                j12(i1) = jsave
              endif
530         continue
520       continue
        endif
*
* here for CS calculations
      elseif (csflag) then
*  assign quantum numbers and energies for CS calculations
*  here for Dubernet and Hutson case 1C
        if (ihomo) then
          n=0
          nlevel = 0
          xnu=nu
* sum over molecular levels
          do 125 jmol=jmolmin,jmax,2
          erot = brot*jmol*(jmol + 1)
* sum over atomic levels
          do 120 jai=0,2
            nlevel = nlevel+1
            if (jai .eq. 0) ehold(nlevel)=-zero/econv
            if (jai .eq. 1) ehold(nlevel)=+aso1/econv
            if (jai .eq. 2) ehold(nlevel)=aso2/econv
            ehold(nlevel) = ehold(nlevel) + erot/econv
            ishold(nlevel)=jai
            jhold(nlevel)=jmol
            j12min=iabs(jai-jmol)
            j12min=max(j12min,nu)
            j12max=jai+jmol
            if (j12max .ge. j12min) then
              do 110 j12i=j12min,j12max
                li=jtot
                xjtot=jtot
                n = n + 1
                if (n .gt. nmax) go to 180
                l(n) = li
                is(n)=jai
                j(n)=jmol
                j12(n)=j12i
                xj12=j12i
* centrifugal contribution to hamiltonian
                cent(n) = xjtot*(xjtot+1)+xj12*(xj12+1)-2*xnu**2
* constant channel energy
                 sc2(n) = ehold(nlevel)
110           continue
            endif
120       continue
125       continue
        else if (.not.ihomo) then
*
*  here for Dubernet and Hutson case 1A
*          if (iop .eq. 0) jmol=0
*          if (iop .eq. 1) jmol=1
*  including spin-orbit coupling
          if (.not.flaghf) then
            n=0
            nlevel = 0
            do 155 jmol=jmolmin,jmax,2
              erot = brot*jmol*(jmol + 1)
* sum over atomic levels
* if aso1 > 10000, assume aso = 0 and ja=0.5
              do 150 jai=0,2
                if (aso1 .lt. 10000d0) then
                  nlevel = nlevel+1
                  if (jai .eq. 0) ehold(nlevel)=-zero/econv
                  if (jai .eq. 1) ehold(nlevel)=+aso1/econv
                  if (jai .eq. 2) ehold(nlevel)=aso2/econv
                else
                  if (jai. eq. 1) goto 150
                  nlevel=nlevel+1
                  ehold(nlevel)=zero
                endif
                ehold(nlevel) = ehold(nlevel) + erot/econv
                ishold(nlevel)=jai
                jhold(nlevel)=jmol
                xjai=jai
                mink=-jmol
                maxk=-mink
                do 140 ik=mink,maxk
                  xk=ik
                  xomega=xnu-xk
                  if (abs(xomega) .le. xjai) then
                    li=jtot
                    n = n + 1
                    if (n .gt. nmax) go to 180
                    l(n) = li
                    xjtot=jtot
                    is(n)=jai
                    j(n)=jmol
                    j12(n)=ik
* centrifugal contribution to hamiltonian
                    cent(n)=xjtot*(xjtot+1)+xjai*(xjai+1)+
     :                  jmol*(jmol+1)-2*xnu**2+2*ik*xomega
* constant channel energy
                    sc2(n)=brot*jmol*(jmol+1)
                    if (aso1 .lt. 10000d0) then
                      if (jai .eq. 0) sc2(n)=sc2(n)-zero
                      if (jai .eq. 1) sc2(n)=sc2(n)+aso1
                      if (jai .eq. 2) sc2(n)=sc2(n)+aso2
                    else
                      sc2(n)=sc2(n)
                    endif
                    sc2(n)=sc2(n)/econv
                  endif
140             continue
150           continue
155         continue
          else
* here for spin-free basis
            n=0
            nlevel = 0
            do 165 jmol=jmolmin,jmax,2
* sum over molecular projection levels
              if (iop .eq. 1) kmin=1
              if (iop .eq. 0) kmin=0
              do 160 ik=-kmin,kmin
                nlevel = nlevel+1
                ehold(nlevel)=zero
                ishold(nlevel)=1
                iomega=nu-ik
                if (iabs(iomega) .le. 1) then
                  jhold(nlevel)=ik
                  li=jtot
                  n = n + 1
                  if (n .gt. nmax) go to 180
                  j12(n)=jmol
                  l(n) = li
                  is(n)=1
                  j(n)=ik
* centrifugal contribution to hamiltonian
                  cent(n) = jtot*(jtot+1)+2+jmol*(jmol+1)
     :                     -2*nu**2+2*ik*iomega
* constant channel energy
                  sc2(n)=brot*jmol*(jmol+1)
                  sc2(n)=sc2(n)/econv
                endif
160           continue
165         continue
         endif
        endif
      endif
*
180   if (n .gt. nmax) then
        write (9, 185) n, nmax
        write (6, 185) n, nmax
185     format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF',
     :         i4,'; ABORT')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions or for rcut < 0
*  or for bound state calculations!
      if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
        emin = 1.e+7
        do 190  i = 1, n
          if (sc2(i) .le. ered) then
*  here if channel is
            if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut)
     :          .gt. (ered - sc2(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (sc2(i) .lt. emin) emin = sc2(i)
*  emin now contains the lowest channel energy for which this
*  condition is met
            end if
          end if
 190    continue
*  now eliminate all channels with sc2(eint) .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn = 0
          do 195 i = 1, n
            if (sc2(i) .lt. emin) then
*  here if this channel is to be included
              nn = nn + 1
              sc2(nn) = sc2(i)
              is(nn) = is(i)
              j(nn) = j(i)
              cent(nn) = cent(i)
              l(nn) = l(i)
              j12(nn) = j12(i)
            end if
195       continue
*  reset number of channels
          n = nn
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
*  form list of all energetically distinct rotational levels included in the
*  channel basis and their energies and
*  sort this list to put closed levels at end
*  also determine number of levels which are open
*
*  find lowest energy
      emin = 1.e+7
      do 188 i = 1, nlevel
        if (ehold(i) .lt. emin) emin = ehold(i)
188   continue
*  shift energies so that lowest level has zero energy
      do 189 i = 1, n
        eint(i) = sc2(i) - emin
189   continue
      do 191 i = 1, nlevel
        ehold(i) = ehold(i) - emin
191   continue
*
      nlevop = 0
      do 580  i = 1, nlevel - 1
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
              go to 580
            end if
75        continue
        end if
580    continue
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
*
      ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      if (iop .eq. 0) jmol=0
      if (iop .eq. 1) jmol=1

*      emin=0
*      emin=sc2(idmin(n,sc2,1))*econv
*      erot=brot*(jmol+1)*jmol
*      do 205 i=1,n
*        eint(i)=sc2(i)-(erot+emin)/econv
*205   continue

*  now list channels if requested
      if (bastst) then
        if (ihomo) then
          if (.not.csflag) then
            write (6, 210)
            write (9, 210)
210         format(/'   N   JA(IS) JMOL(J)  J12    L      EINT(CM-1)')
            do 220  i = 1, n
              write (6, 215) i,is(i),j(i),j12(i),l(i),
     :             eint(i)*econv
              write (9, 215) i,is(i),j(i),j12(i),l(i),
     :             eint(i)*econv
215           format (i4, i7,i6, i8, i6, f14.3)
220         continue
          else
            write (6, 225)
            write (9, 225)
225         format
     :       (/'   N   JA(IS) JMOL(J)  J12   L  CENT    EINT(CM-1)')
            do 235  i = 1, n
              icent=cent(i)
              write (6, 230) i,is(i),j(i),j12(i),l(i),icent,
     :             eint(i)*econv
              write (9, 230) i,is(i),j(i),j12(i),l(i),icent,
     :             eint(i)*econv
230           format (i4, i7,i6,i8, i5,i6, f14.3)
235         continue
          endif
        elseif (.not.ihomo) then
          if (.not.flaghf) then
            write (6, 255)
            write (9, 255)
255         format(/'   N   JA(IS) JMOL(J)   K   CENT    EINT(CM-1)')
            do 265  i = 1, n
              icent=cent(i)
              write (6, 260) i,is(i),j(i),j12(i),icent,eint(i)*econv
              write (9, 260) i,is(i),j(i),j12(i),icent,eint(i)*econv
260           format (i4, i7,i6,i8,i6, f14.3)
265         continue
          else
            write (6, 270)
            write (9, 270)
270         format(/'   N   LA(IS)  JMOL  K(J)   CENT    EINT(CM-1)')
            do 275  i = 1, n
              write (6, 272) i,is(i),isc1(i),j(i),idnint(cent(i)),
     :                       eint(i)*econv
              write (9, 272) i,is(i),isc1(i),j(i),idnint(cent(i)),
     :                       eint(i)*econv
272           format (i4, i8,3i6, f14.3)
275         continue
          end if

        endif
      end if
*
*  now calculate coupling matrix elements
*  ordering of terms is as follows:
*  CC calculations (csflag=.false.) and CS case 1C calculations
*    lam=1  V220(R)
*    lam=2  V022(R)
*    lam=3  V202(R)
*    lam=4  V222(R)
*    lam=5  V224(R)
*    lam=6  V404(R)
*    lam=7  V422(R)
*    lam=8  V424(R)
*    lam=9  V426(R)
*   CS case 1A calculations
*    lam=1  v220(R)
*    lam=2  v020(R)
*    lam=3  v200(R)
*    lam=4  v222(R)
*    lam=5  v221(R)
*    lam=6  v400(R)
*    lam=7  v420(R)
*    lam=8  v422(R)
*    lam=9  v421(R)
      if (bastst .and. iprint .gt. 1) then
        if (.not.csflag .or. (csflag .and. ihomo)) then
          write (6, 280)
          write (9, 280)
280       format (/' ILAM  LR LA L12  ICOL  IROW   I    IV2    VEE')
        else
          write (6, 281)
          write (9, 281)
281       format (/' ILAM  LR LA MU   ICOL  IROW   I    IV2    VEE')
        endif
      end if
* i counts v2 elements
* inum counts v2 elements for given lambda
* ilam counts number of v2 matrices
* ij is address of given v2 element in present v2 matrix
      i = 0
      ilam=0
      do 320 il = 1,lammax(1)
*      do 320 il = 4,4
        ilam=ilam+1
        inum = 0
        ij=0
        ilamr=lamr(il)
        ilama=lama(il)
        ilam12=lam12(il)
        imu=mu(il)
        do 310  icol= 1, n
          do 300  irow = icol, n
*        do 310  icol= 2, 3
*          do 300  irow = icol, 3
            ij = ntop * (icol - 1) +irow
            j12row=j12(irow)
            j12col=j12(icol)
            if (.not. csflag .or. (csflag .and. ihomo)) then
              call vlmh3p (irow,icol,jtot,jlpar,j(irow),j(icol),
     :        is(irow),
     :        is(icol), j12row, j12col, l(irow), l(icol), ilamr,
     :        ilama, ilam12, nu, csflag, vee)
            else if (csflag .and. .not.ihomo) then
* NB j12row and j12col are krow and kcol here
              call vlmh3pc(irow,icol,jtot,jlpar,j(irow),j(icol),
     :        is(irow),
     :        is(icol), j12row, j12col, ilamr,
     :        ilama, imu, nu, jmol, flaghf, vee)
            endif
*           write (6,291) irow,icol,jtot,jlpar,j(irow),j(icol),is(irow),
*     :         is(icol),j12row,j12col,nu,ilamr,ilama,ilam12,vee
291    format(14i3,g17.8)
            if (vee .eq. 0) goto 300
              i = i + 1
              inum = inum + 1
              if (i .gt. nv2max) goto 300
                v2(i) = vee
                iv2(i) = ij
                if (bastst .and. iprint .gt. 1) then
                  if (.not. csflag .or. (csflag .and. ihomo)) then
                    write (6, 290) ilam, ilamr,ilama,ilam12, icol, irow,
     :                           i, iv2(i), vee
                    write (9, 290) ilam, ilamr,ilama,ilam12, icol, irow,
     :                           i, iv2(i), vee
290                 format (i4, i5,2i3,i5, 2i6, i6, g17.8)
                  elseif (csflag .and. .not.ihomo) then
                    write (6, 290) ilam, ilamr,ilama,imu, icol, irow,
     :                           i, iv2(i), vee
                    write (9, 290) ilam, ilamr,ilama,imu, icol, irow,
     :                           i, iv2(i), vee
                  endif
                endif
300       continue
310     continue
      if(ilam.gt.nlammx) then
        write(6,311) ilam
311     format(/' ILAM.GT.NLAMMX IN BAH3P')
        call exit
      end if
      lamnum(ilam) = inum
      if (bastst .and. iprint .gt. 1) then
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
      subroutine vlmh3p (irow, icol, jtot, jlpar, j, jp, ja, jap,
     :            j12, j12p, l, lp, lamr, lama, lam12, nu, csflag, vee)
* --------------------------------------------------------------------
*  subroutine to evaluate the angular coupling matrix element for rotationally
*  inelastic collisions of a homonuclear molecule (j=1) and a 3P atom
*  in full closecoupling (csflag = .false.) or dubernet-hutson case 1C
*  (csflag = .true.)

*  authors:  millard alexander
*  current revision date: 6-jul-1995
* --------------------------------------------------------------------
*  variables in call list:
*  irow, icol:  row and column of fully coupled states
*  jtot:      total angular momentum
*  jlpar:     parity (-1 for e, +1 for f)
*  j, jp      bra and ket values of molecular angular momentum
*  ja, jap:   bra and ket values of atomic angular momentum
*  j12, j12p: bra and ket values of total internal angular momentum
*  l, lp:     bra and ket values of orbital angular momentum
*  lamr, lama, lam12:      value of SF expansion indicies
*  nu         projection index (called P by Dubernet and Hutson)
*  vee:        on return:  contains desired coupling matrix element
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag
      data  zero, one, two, three /0.d0, 1.d0, 2.d0, 3.d0/
      data half /0.5d0/
      vee = zero
      if (j.eq. 0 .and. jp.eq. 0 .and. lamr .ne. 0) return
      iphase=ja+jp+j12-jtot-1
      if (csflag) iphase=ja+j+j12p+lam12+nu-1
      iphase=(-1)**iphase
      xja=ja
      xjap=jap
      xjtot=jtot
      xj12=j12
      xj12p=j12p
      xnu=nu
      xj=j
      xjp=jp
      xlama=lama
      xlamr=lamr
      xlam12=lam12
      xl=l
      xlp=lp
      t1=0.d0
      t2=0.d0
      t3=0.d0
      t4=0.d0
      t5=0.d0
      tt=0.d0

      xnorm=three*sqrt((two*j+1)*(two*jp+1)*(2*xlam12+1)*(2*xja+1)*
     :     (2*xjap+1)*(2*xj12+1)*(2*xj12p+1))
      if (.not.csflag) xnorm=xnorm*sqrt((2*xl+1)*(2*xlp+1))

      t1=xf3j(xj,xlamr,xjp,zero,zero,zero)
      t2=xf3j(one,xlama,one,zero,zero,zero)
      term=t1*t2
      if (term .ne. zero) then
        t3=xf6j(xja,xlama,xjap,one,one,one)
      else
        return
      endif
      if (t3 .ne. zero) then
        term=term*t3
        if (.not.csflag)
     :      t4=xf3j(xl,xlam12,xlp,zero,zero,zero)
        if (csflag)
     :      t4=xf3j(xj12,xlam12,xj12p,-xnu,zero,xnu)
      else
        return
      endif
      if (t4 .ne. zero) then
        term=term*t4
        if (.not.csflag)
     :     t5=xf6j(xj12,xl,xjtot,xlp,xj12p,xlam12)
        if (csflag) t5=one
      else
        return
      endif
      if (t5 .ne. zero) then
        term=term*t5
        tt=
     :    xf9j(xj12,xlam12,xj12p,xja,xlama,xjap,xj,xlamr,xjp)
        term=term*tt
      else
        return
      endif
      vee=term*xnorm*iphase
      return
      end
* --------------------------------------------------------------------
      subroutine vlmh3pc (irow, icol, jtot, jlpar, j, jp, ja, jap,
     :            k, kp, lamr, lama, mu, nu, jmol, flaghf, vee)
* --------------------------------------------------------------------
*  subroutine to evaluate the angular coupling matrix element for rotationally
*  inelastic collisions of a homonuclear molecule (j=1) and a 3P atom
*  dubernet and hutson case 1A
*  if flaghf = .false., spin included
*  if flaghf = .true., spinless atom (i'm not sure this is working)

*  authors:  millard alexander
*  current revision date: 13-jun-1997
* --------------------------------------------------------------------
*  variables in call list:
*  irow, icol:  row and column of fully coupled states
*  jtot:      total angular momentum
*  j, jp      bra and ket values of molecular angular momentum
*  ja, jap:   bra and ket values of atomic angular momentum
*  k, kp:      bra and ket values of projection of jmolecular
*  lamr, lama, mu:      value of BF expansion indicies
*  nu         projection index (called P by Dubernet and Hutson)
*  vee:        on return:  contains desired coupling matrix element
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical flaghf
      data  zero, one, two /0.d0, 1.d0, 2.d0/
      data half /0.5d0/
      vee = zero
*      if (jmol .eq. 0 .and. lamr .ne. 0) return
      if (j .eq. 0 .and. jp .eq. 0 .and. lamr .ne. 0) return
      iomeg=nu-k
      iomegp=nu-kp
      iphase=(-1)**(k+iomeg-1)
      if (.not.flaghf) then
        xomeg=iomeg
        xomegp=iomegp
        xja=ja
        xjap=jap
        xjtot=jtot
        xnu=nu
      else
        xomeg=iomeg
        xomegp=iomegp
        xja=ja
        xjap=jap
        xjtot=jtot
        xnu=nu
      endif
      xj=j
      xjp=jp
      xlama=lama
      xlamr=lamr
      xk=k
      xkp=kp
      t1=0.d0
      t2=0.d0
      t3=0.d0
      t4=0.d0
      t5=0.d0
      if (.not.flaghf) then
        xnorm=3*sqrt((two*j+1)*(2*jp+1)*(2*xja+1)*(2*xjap+1))
        if (mu .eq. 0) then
          xmu=zero
          t1=xf3j(xj,xlamr,xjp,-xk,xmu,xkp)
        else
          iden=-k+mu+kp
          if (iden .eq. 0) then
            xmu=mu
            t1=xf3j(xj,xlamr,xjp,-xk,xmu,xkp)
          else
            idenp=-k-mu+kp
            if (idenp .eq. 0) then
              xmu=-mu
              t1=xf3j(xj,xlamr,xjp,-xk,xmu,xkp)
            else
              t1=zero
            endif
          endif
        endif
        t2=xf3j(xj,xlamr,xjp,zero,zero,zero)
        t3=xf3j(one,xlama,one,zero,zero,zero)
        term=t1*t2*t3
        t4=zero
        t5=zero
        if (term .ne. zero) then
          t4=xf3j(xjap,xlama,xja,xomegp,-xmu,-xomeg)
        endif
        term=term*t4
        if (t4 .ne. zero) then
          t5=xf6j(xja,xlama,xjap,one,one,one)
        endif
        vee=t5*term*xnorm*iphase
      else
        xnorm=9
        if (mu .eq. 0) then
          xmu=zero
          t1=xf3j(xj,xlamr,xjp,-xk,xmu,xkp)
        else
          iden=-k+mu+kp
          if (iden .eq. 0) then
            xmu=mu
            t1=xf3j(xj,xlamr,xjp,-xk,xmu,xkp)
          else
            idenp=-k-mu+kp
            if (idenp .eq. 0) then
              xmu=-mu
              t1=xf3j(xj,xlamr,xjp,-xk,xmu,xkp)
            else
              t1=zero
            endif
          endif
        endif
        t2=xf3j(xj,xlamr,xjp,zero,zero,zero)
        t3=xf3j(one,xlama,one,zero,zero,zero)
        t4=xf3j(one,xlama,one,-xomeg,-xmu,xomegp)
        vee=t1*t2*t3*t4*xnorm*iphase
      endif
      return
      end
