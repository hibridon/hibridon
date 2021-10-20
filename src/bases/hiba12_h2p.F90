! syh2p (savh2p/ptrh2p) defines, saves variables and reads               *
!                  potential for homonuclear+2P atom scattering          *
! --------------------------------------------------
#include "assert.h"
module mod_hiba12_h2p
contains
subroutine bah2p (j, l, is, jhold, ehold, ishold, nlevel, nlevop, &
                  sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential
!  for collision of a doublet atom in a P state and a homonuclear molecule
!  authors:  millard alexander
!
!  revised routine so that CC scattering calculations work.  it appears that
!  this routine had previously been used only for bound-state calculations
!  (4-mar-2013 by p.dagdigian)
!
!  edited list of SF lambda values
!
!  current revision date:  3-sep-2020 by p.dagdigian
! --------------------------------------------------------------------
!  variables in call list:
!    j:        for CC calculation:  returns diatom angular momentum (j) for ea
!              for CD-case 1A and 1C calculation:  on return contains atom+dia
!              (j12)
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains atomic angular momentum (ja) of each channel
!              (0 or 1)
!    jhold:    for CC calculations:  on return contains diatom angular momentu
!              for each level
!              for CD calculations:  on return contains atom+diatom angular mo
!              for each level
!    ehold:    on return contains energy in hartrees of each level
!    ishold:   on return contains electronic angular momentum of each level
!    nlevel:   on return contains number of energetically distinct
!              levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              levels used in channel basis which are open
!              asymptotically
!    sc1,sc2: scratch vectors of length at least nmax
!    sc3,sc4:  scratch vectors of length at least nmax
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!              note that this is ignored in molecule-surface collisions!!!
!    jtot:     total angular momentum
!    flaghf:   if .true., then system with half-integer spin (this is the
!              case here)
!              if .false., then system with integer spin
!    flagsu:   if .true., then molecule-surface collisons (not applicable)
!    csflag:   if .true., then coupled-states calculation
!              if .false., the close-coupled calculation
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true. , then Dubernet and Hutson case 1C
!              if .false., then Dubernet and Hutson case 1A
!    nu:       coupled-states projection index (called P by dubernet and
!              hutson
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              the parity is [Dubernet and Hutson] (-1)**(j+l+L), where
!              L (=1 here) is the electronic orbital angular momentum.
!              hence only those channels are included for which
!                  (-1)**(j+l+1)=jlpar
!              in cs calculation jlpar is set equal to 1 in calling program
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to subroutines
!              propag and soutpt.  ntop is set in basis only if nu = numin
!              otherwise it is unchanged from the value supplied by the
!              calling program
!    note!!!   if flaghf = .true., then the true values of the rotational
!    quantum numbers, the total angular momentum, and the coupled-states
!    projection index are equal to the values stored in j, jtot, and nu
!    plus 1/2
!  variables in common block /cosysr/
!    isrcod:   number of real system dependent variables
!    brot:     rotational constant of molecule
!    aso:      spin-orbit constant of 2P atom (cm-1)
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different types of electronic coupling terms
!              this should be 1 here
!    iop:      ortho (iop=1) or para (iop=0)
!    jmax:     maximum rotational quantum number of diatomic
!  variables in common block /coered/
!    ered:     collision energy in atomic units (hartrees)
!    rmu:      collision reduced mass in atomic units
!              (mass of electron = 1)
!  variable in module mod_conlam
!    nlam:     the number of case(a) interaction potentials actually used
!              if this has not been set in the pot subroutine, it it
!              here set to 5
!  variable in common block /coconv/
!     econv:   conversion factor from cm-1 to hartrees
!     xmconv:  converson factor from amu to atomic units
!  subroutines called:
!   vlmh2p:    returns angular coupling coefficient for particular
!              choice of channel index
! ------------------------------------------------------------
use mod_cov2, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coj12, only: j12
use mod_conlam, only: nlam
implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable :: v2
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst
#include "common/parbas.F90"
#include "common/parbasl.F90"
common /cosysr/ isrcod, junkr, brot,aso
common /cosysi/ nscode, isicod, nterm, iop,jmax
common /coipar/ iiipar(9), iprint
common /coered/ ered, rmu
common /coskip/ nskip, iskip
common /coconv/ econv, xmconv
common /cojtot/ jjtot,jjlpar
dimension j(9), l(9), jhold(9), ehold(9), sc2(9), sc3(9), &
          sc4(9), ishold(9), is(9), isc8(9)
dimension lamr(13),lama(13),lam12(13), mu(13)
dimension lamrold(9),lamaold(9),lam12old(9), muold(9)
!  edited the list below (3-sep-2020)
data lamr  /2,4,6,0,2,4,6,2,4,6,2,4,6/
data lama  /0,0,0,2,2,2,2,2,2,2,2,2,2/
data lam12 /2,4,6,2,0,2,4,2,4,6,4,6,8/
data mu    /0,0,0,0,0,0,0,2,2,2,1,1,1/
!*  original list
!      data lamr  /2,4,6,8,0,2,4,6,8,2,4,6,8,2,4,6,8/
!      data lama  /0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2/
!      data lam12 /2,4,6,8,2,0,2,4,6,2,4,6,8,4,6,8,10/
!      data mu    /0,0,0,0,0,0,0,0,0,2,2,2,2,1,1,1,1/
!  here is information on lambda terms for compatability with older
!  pot subroutines
data lamrold /2,0,2,2,2,4,4,4,4/
data lamaold /2,2,0,2,2,0,2,2,2/
data lam12old /0,2,2,2,4,4,2,4,6/
data muold / 0,0,0,2,1,0,0,2,1/

!   econv is conversion factor from cm-1 to hartrees
!   xmconv is converson factor from amu to atomic units
zero = 0.d0
half=0.5d0
two = 2.d0
!  check for consistency in the values of flaghf and csflag
if (.not.flaghf .and. (.not.csflag .or. csflag .and. ihomo))  then
  write (6, 5)
  write (9, 5)
5   format &
 (' *** FLAGHF = .FALSE. FOR CC OR CS1C MOL+ 2P ATOM; ABORT ***' )
  if (bastst) then
    return
  else
    stop
  end if
end if
if (flaghf) then
  xnu=nu+half
else
  xnu=nu
endif
!  check that iop equals 0 or 1 only
if (iop.ne.0 .and. iop.ne.1) then
  write (6,6) iop
  write (6,9) iop
6   format (' *** IOP =',i2,'.  THIS PARAMETER MUST EQUAL 0 OR 1;', &
    ' ABORT ***')
  stop
end if
if(nterm.ne.1) then
  write(6,9) nterm
  write(9,9) nterm
9   format(' *** NTERM = ',i3,' .NE. 1 FOR MOL + 2P ATOM;', &
    ' ABORT ***')
  stop
end if
if (nlam.eq.0) nlam=5
nsum=nlam
if (bastst) write (6, 14) nsum
write (9, 14) nsum
14 format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
if (bastst) then
  if (flagsu) then
    write (6,16) rmu * xmconv, ered * econv, jtot+half, jlpar
    write (9,16) rmu * xmconv, ered * econv, jtot+half, jlpar
16     format(/ &
  ' **  MOL + 2P ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4, &
  '             E=', f9.2, /, ' JTOT=', f5.1, 2x,' JLPAR=',i2)
  else
    if(.not.csflag) then
      write (6,20) &
          rmu * xmconv, ered * econv, jtot+half, jlpar
      write (9,20) &
          rmu * xmconv, ered * econv, jtot+half, jlpar
20       format(/,' **  CC MOL + 2P ATOM',' ** RMU=', f9.4, &
           '       E=',f9.2,'   JTOT=', f5.1, 2x,' JLPAR=',i2)
    else if (csflag .and. ihomo) then
      write (6,22) &
          rmu * xmconv, ered * econv, jtot+half, xnu
      write (9,22) &
          rmu * xmconv, ered * econv, jtot+half, xnu
22       format(/,' **  CD 1C MOL + 2P ATOM',' ** RMU=', f9.4, &
           '       E=',f9.2,'   JTOT=', f5.1, 2x,' P=',f4.1)
    else if (csflag .and. .not.ihomo) then
      write (6,23) &
          rmu * xmconv, ered * econv, jtot+half, xnu
      write (9,23) &
          rmu * xmconv, ered * econv, jtot+half, xnu
23       format(/,' **  CD 1A MOL + 2P ATOM',' ** RMU=', f9.4, &
           '       E=',f9.2,'   JTOT=', f5.1, 2x,' P=',f7.2)
    endif
  end if
  if (.not. flagsu) write (9,30) rcut
30   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f6.2)
endif
! check for consistency in projection quantum number
if (nu .gt. jtot .and. csflag) then
  write (6, 35) nu, jtot
  write (9, 35) nu, jtot
35   format (' NU =',i3,' .GT. JTOT =',i3,'; JTOT RESET TO NU')
  jtot=nu
endif
!
!  save jtot and jlpar for use later in transformation
jjtot=jtot
jjlpar=jlpar
jmolmin=0
if (iop .eq. 1) jmolmin=1
if (.not. csflag) then
!  assign quantum numbers and energies for CC calculations
  n=0
  nlevel = 0
! sum over molecular levels
  do 105 jmol=jmolmin,jmax,2
! sum over atomic levels
  do 100 jai=0,1
    nlevel = nlevel+1
! rotational energy
    erot = brot*jmol*(jmol+1)
    ehold(nlevel)=(erot-aso)/econv
    if (jai .eq. 1) ehold(nlevel)=(erot+aso/two)/econv
    jhold(nlevel)=jmol
    ishold(nlevel)=jai
    xjai=jai+half
    j12min=abs(xjai-jmol)-half
    j12max=jai+jmol
    do 90 j12i=j12min,j12max
      lmin=iabs(jtot-j12i)
      lmax=jtot+j12i+1
      do 80 li=lmin,lmax
! parity of space-frame states is (-1)**(L-atom+jmol+L-orbital)
        if ((-1)**li .eq. jlpar*(-1)**(jmolmin-1)) then
          n = n + 1
          if (n .gt. nmax) go to 180
          l(n) = li
          is(n)=jai
          j(n)=jmol
          j12(n)=j12i
          isc8(n)=j12i
! centrifugal contribution to hamiltonian
          cent(n) = li*(li+1)
          sc2(n)=ehold(nlevel)
        endif
80       continue
90     continue
100   continue
105   continue
!
! here for CS calculations
elseif (csflag) then
!  assign quantum numbers and energies for CS calculations
!  here for Dubernet and Hutson case 1C
  if (ihomo) then
    n=0
    nlevel = 0
    xnu=nu+half
! sum over molecular levels
    do 125 jmol=jmolmin,jmax,2
    erot = brot*jmol*(jmol + 1)
! sum over atomic levels
    do 120 jai=0,1
      xjai=jai+half
      j12min=abs(xjai-jmol)-half
      j12min=max(j12min,nu)
      j12max=jai+jmol
      if (j12max .ge. j12min) then
        do 110 j12i=j12min,j12max
          nlevel = nlevel+1
          if (jai .eq. 0) ehold(nlevel)=-aso/econv
          if (jai .eq. 1) ehold(nlevel)=aso/(econv*two)
          ehold(nlevel) = ehold(nlevel) + erot/econv
          ishold(nlevel)=jai
          jhold(nlevel)=j12i
          li=jtot
          xjtot=jtot+half
          n = n + 1
          if (n .gt. nmax) go to 180
          l(n) = li
          is(n)=jai
          j(n)=j12i
          j12(n)=jmol
          isc8(n)=jmol
          xj12=j12i+half
! centrifugal contribution to hamiltonian
          cent(n) = xjtot*(xjtot+1)+xj12*(xj12+1)-2*xnu**2
! constant channel energy
          sc2(n)=brot*jmol*(jmol+1)
          if (jai .eq. 0) sc2(n)=sc2(n)-aso
          if (jai .eq. 1) sc2(n)=sc2(n)+aso/two
          sc2(n)=sc2(n)/econv
110         continue
      endif
120     continue
125     continue
  else if (.not.ihomo) then
!  here for Dubernet and Hutson case 1A
!          if (iop .eq. 0) jmol=0
!          if (iop .eq. 1) jmol=1
!  including spin-orbit coupling
    if (flaghf) then
      n=0
      nlevel = 0
      do 155 jmol=jmolmin,jmax,2
      erot = brot*jmol*(jmol + 1)
! sum over atomic levels
! if aso > 10000, assume aso = 0 and ja=0.5
        do 150 jai=0,1
          if (aso .lt. 10000d0) then
            nlevel = nlevel+1
            if (jai .eq. 0) ehold(nlevel)=-aso/econv
            if (jai .eq. 1) ehold(nlevel)=aso/(econv*two)
          else
            if (jai .eq. 1) goto 150
            nlevel=nlevel+1
            ehold(nlevel)=zero
          endif
          ehold(nlevel) = ehold(nlevel) + erot/econv
          ishold(nlevel)=jai
          xjai=jai+half
          mink=-jmol
          maxk=-mink
          do 140 ik=mink,maxk
            xk=ik
            xomega=xnu-xk
            if (abs(xomega) .le. xjai) then
              jhold(nlevel)=ik
              li=jtot
              n = n + 1
              if (n .gt. nmax) go to 180
              j12(n)=jmol
              isc8(n)=jmol
              l(n) = li
              xjtot=jtot+half
              is(n)=jai
              j(n)=ik
! centrifugal contribution to hamiltonian
              cent(n)=xjtot*(xjtot+1)+xjai*(xjai+1)+ &
                  jmol*(jmol+1)-2*xnu**2+2*ik*xomega
! constant channel energy
              sc2(n)=brot*jmol*(jmol+1)
              if (aso .lt. 10000d0) then
                if (jai .eq. 0) sc2(n)=sc2(n)-aso
                if (jai .eq. 1) sc2(n)=sc2(n)+aso/two
              else
                sc2(n)=sc2(n)
              endif
              sc2(n)=sc2(n)/econv
            endif
140           continue
150         continue
155       continue
    else
! here for spin-free basis
      n=0
      nlevel = 0
      do 165 jmol=jmolmin,jmax,2
! sum over molecular projection levels
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
            isc8(n)=jmol
            l(n) = li
            is(n)=1
            j(n)=ik
! centrifugal contribution to hamiltonian
            cent(n) = jtot*(jtot+1)+2+jmol*(jmol+1) &
                     -2*nu**2+2*ik*iomega
! constant channel energy
            sc2(n)=brot*jmol*(jmol+1)
            sc2(n)=sc2(n)/econv
          endif
160         continue
165       continue
   endif
  endif
endif
!
180 if (n .gt. nmax) then
  write (9, 185) n, nmax
  write (6, 185) n, nmax
185   format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF', &
         i4,'; ABORT')
  if (bastst) then
    return
  else
    stop
  end if
end if
!
!  find lowest energy
emin = 1.e+7
do 188 i = 1, n
  if (sc2(i) .lt. emin) emin = sc2(i)
188 continue
!  shift energies so that lowest level has zero energy
do 189 i = 1, n
  eint(i) = sc2(i) - emin
189 continue
do 191 i = 1, nlevel
  ehold(i) = ehold(i) - emin
191 continue
!  now check to see if any of the open channels are closed at r=rcut
!  this is not done for molecule-surface collisions or for rcut < 0
!  or for bound state calculations!
if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
  emin = 1.e+7
  do 190  i = 1, n
    if (eint(i) .le. ered) then
!  here if channel is
      if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut) &
          .gt. (ered - eint(i)) ) then
!  here if channel is open asymptotically but closed at r = rcut
        if (eint(i) .lt. emin) emin = eint(i)
!  emin now contains the lowest channel energy for which this
!  condition is met
      end if
    end if
190   continue
!  now eliminate all channels with eint .ge. emin if any of the channels
!  are open asymptotically but closed at r = rcut
  if (emin.lt.ered) then
    nn = 0
    do 195 i = 1, n
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        eint(nn) = eint(i)
        is(nn) = is(i)
        j(nn) = j(i)
        cent(nn) = cent(i)
        l(nn) = l(i)
      end if
195     continue
!  reset number of channels
    n = nn
  end if
end if
!  return if no channels
if (n .eq. 0) return
!  form list of all energetically distinct rotational levels included in the
!  channel basis and their energies and
!  sort this list to put closed levels at end
!  also determine number of levels which are open
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
75     continue
  end if
580  continue
if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
!
ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1
!  now list channels if requested
if (bastst) then
  if (ihomo) then
    if (.not.csflag) then
      write (6, 210)
      write (9, 210)
210       format(/'   N   JA(IS) JMOL  J12(J)    L    EINT(CM-1)')
      do 220  i = 1, n
        write (6, 215) i,is(i)+half,j(i),j12(i)+half,l(i), &
             eint(i)*econv
        write (9, 215) i,is(i)+half,j(i),j12(i)+half,l(i), &
             eint(i)*econv
215         format (i4, f8.1,i5, f7.1, i7, f12.3)
220       continue
    else
      write (6, 225)
      write (9, 225)
225       format(/'   N   JA(IS) JMOL  J12(J)   CENT   EINT(CM-1)', &
            ';  L = JTOT')
      do 235  i = 1, n
        xs=is(i)+half
        xj=j(i)+half
        icent=cent(i)
        write (6, 230) i,xs,j12(i),xj,icent,eint(i)*econv
        write (9, 230) i,xs,j12(i),xj,icent,eint(i)*econv
230         format (i4, f8.1,i6,f7.1, i8, f12.3)
235       continue
    endif
  elseif (.not.ihomo) then
    if (flaghf) then
      write (6, 255)
      write (9, 255)
255       format( &
       /'   N   JA(IS)  JMOL  K(J)    L    CENT    EINT(CM-1)')
      do 265  i = 1, n
        xs=is(i)+half
        icent=cent(i)
        write (6, 260) i,xs,j12(i),j(i),l(i),icent,eint(i)*econv
        write (9, 260) i,xs,j12(i),j(i),l(i),icent,eint(i)*econv
260         format (i4, f8.1,4i6, f14.3)
265       continue
    else
      write (6, 270)
      write (9, 270)
270       format(/'   N   LA(IS)  JMOL  K(J)   CENT    EINT(CM-1)')
      do 275  i = 1, n
        write (6, 272) i,is(i),j12(i),j(i),idnint(cent(i)), &
                       eint(i)*econv
        write (9, 272) i,is(i),j12(i),j(i),idnint(cent(i)), &
                       eint(i)*econv
272         format (i4, i8,3i6, f14.3)
275       continue
    end if
  endif
end if
!
!  now calculate coupling matrix elements
!  ordering of terms is as follows:
!  CS calculations case 1A (csflag=.true. and ihomo=.false)
!  body frame expansion coefficients
!    lam=1  v200(R)
!    lam=2  v400(R)
!    lam=3  v600(R)
!    lam=4  v800(R)
!    lam=5  v020(R)
!    lam=6  v220(R)
!    lam=7  v420(R)
!    lam=8  v620(R)
!    lam=9  v820(R)
!    lam=10  v222(R)
!    lam=11  v422(R)
!    lam=12  v622(R)
!    lam=13  v822(R)
!    lam=14  v221(R)
!    lam=15  v421(R)
!    lam=16  v621(R)
!    lam=17  v821(R)
!  CC calculations (csflag=.false.)
!    lam=1  V202(R)
!    lam=2  V402(R)
!    lam=3  V602(R)
!    lam=4  V802(R)
!    lam=5  V020(R)
!    lam=6  V220(R)
!    lam=7  V422(R)
!    lam=8  V624(R)
!    lam=9  V826(R)
!    lam=10  V222(R)
!    lam=11  V424(R)
!    lam=12  V626(R)
!    lam=13  V828(R)
!    lam=14  V224(R)
!    lam=15  V426(R)
!    lam=16  V628(R)
!    lam=17  V8210(R)
!
!  print out matrix elements in BASTST = T and CHLIST = T
if (bastst .and. clist) then
  if (.not.csflag .or. (csflag .and. ihomo)) then
    write (6, 280)
    write (9, 280)
280     format (/' ILAM  LR LA L12  ICOL  IROW   I    IV2    VEE')
  else
    write (6, 281)
    write (9, 281)
281     format (/' ILAM  LR LA MU   ICOL  IROW   I    IV2    VEE')
  endif
end if
! i counts v2 elements
! inum counts v2 elements for given lambda
! ilam counts number of v2 matrices
! ij is address of given v2 element in present v2 matrix
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=lammax(1), num_channels=ntop)
ilam=0
do 320 il = 1,lammax(1)
!      do 320 il = 4,4
  ilam=ilam+1
  ancouma => v2%get_angular_coupling_matrix(il)
  inum = 0
  ij=0
  if (lammax(1) .eq. 13) then
! here for new (14-feb-1998) generation of pot subroutines
    ilamr=lamr(il)
    ilama=lama(il)
    ilam12=lam12(il)
    imu=mu(il)
  else
! here for older pot subroutines
    ilamr=lamrold(il)
    ilama=lamaold(il)
    ilam12=lam12old(il)
    imu=muold(il)
  endif
  do 310  icol= 1, n
    do 300  irow = icol, n
!        do 310  icol= 2, 3
!          do 300  irow = icol, 3
      ij = ntop * (icol - 1) +irow
!6/13/2002:  reverse definition of j and j12 for CD calculations
      if (csflag) then
         jrow=j12(irow)
         jcol=j12(icol)
         j12row=j(irow)
         j12col=j(icol)
      else
         jrow=j(irow)
         jcol=j(icol)
         j12row=j12(irow)
         j12col=j12(icol)
      endif
      if (.not. csflag .or. (csflag .and. ihomo)) then
        call vlmh2p (irow, icol, jtot, jlpar, jrow, jcol,is(irow), &
        is(icol), j12row, j12col, l(irow), l(icol), ilamr, &
        ilama, ilam12, nu, csflag, vee)
      else if (csflag .and. .not.ihomo) then
        call vlmh2pc(irow, icol, jtot, jlpar, jrow, jcol,is(irow), &
        is(icol), j(irow), j(icol), ilamr, &
        ilama, imu, nu, jmol, flaghf, vee)
      endif
!           write (6,291) irow,icol,jtot,jlpar,jrow,jcol,is(irow),is(icol),
!    : j(irow),j(icol),nu,ilamr,ilama,ilam12,vee
291  format(14i3,g17.8)
      if (vee .eq. 0) goto 300
        inum = inum + 1
          call ancouma%set_element(irow, icol, vee)
          if (bastst .and. clist) then
            if (.not. csflag .or. (csflag .and. ihomo)) then
              write (6, 290) ilam, ilamr,ilama,ilam12, icol, irow, &
                           i, ij, vee
              write (9, 290) ilam, ilamr,ilama,ilam12, icol, irow, &
                           i, ij, vee
290               format (i4, i5,2i3,i5, 2i6, i6, g17.8)
            elseif (csflag .and. .not.ihomo) then
              write (6, 290) ilam, ilamr,ilama,imu, icol, irow, &
                           i, ij, vee
              write (9, 290) ilam, ilamr,ilama,imu, icol, irow, &
                           i, ij, vee
            endif
          endif
300     continue
310   continue

if (bastst .and. iprint .gt. 1) then
  write (6, 315) ilam, ancouma%get_num_nonzero_elements()
  write (9, 315) ilam, ancouma%get_num_nonzero_elements()
315   format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
end if
320 continue
if (clist) then
  write (6, 360) v2%get_num_nonzero_elements()
  write (9, 360) v2%get_num_nonzero_elements()
360   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', &
           i6)
end if
return
end
! --------------------------------------------------------------------
subroutine vlmh2p (irow, icol, jtot, jlpar, j, jp, ja, jap, &
            j12, j12p, l, lp, lamr, lama, lam12, nu, csflag, vee)
! --------------------------------------------------------------------
!  subroutine to evaluate the angular coupling matrix element for rotationally
!  inelastic collisions of a homonuclear molecule (j=1) and a 2P atom
!  in full closecoupling (csflag = .false.) or dubernet-hutson case 1C
!  (csflag = .true.)

!  authors:  millard alexander
!  current revision date: 6-jul-1995
! --------------------------------------------------------------------
!  variables in call list:
!  irow, icol:  row and column of fully coupled states
!  jtot:      total angular momentum
!  jlpar:     parity (-1 for e, +1 for f)
!  j, jp      bra and ket values of molecular angular momentum
!  ja, jap:   bra and ket values of atomic angular momentum
!  j12, j12p: bra and ket values of total internal angular momentum
!  l, lp:     bra and ket values of orbital angular momentum
!  lamr, lama, lam12:      value of SF expansion indicies
!  nu         projection index (called P by Dubernet and Hutson)
!  vee:        on return:  contains desired coupling matrix element
! --------------------------------------------------------------------
implicit double precision (a-h,o-z)
logical csflag
data  zero, one, two, three /0.d0, 1.d0, 2.d0, 3.d0/
data half /0.5d0/
vee = zero
if (j.eq. 0 .and. jp.eq. 0 .and. lamr .ne. 0) return
iphase=ja+jp+j12-jtot
if (csflag) iphase=ja+j+j12p+lam12+nu+1
iphase=(-1)**iphase
xja=ja+half
xjap=jap+half
xjtot=jtot+half
xj12=j12+half
xj12p=j12p+half
xnu=nu+half
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

xnorm=three*sqrt((two*j+1)*(two*jp+1)*(2*xlam12+1)*(2*xja+1)* &
     (2*xjap+1)*(2*xj12+1)*(2*xj12p+1))
if (.not.csflag) xnorm=xnorm*sqrt((2*xl+1)*(2*xlp+1))

t1=xf3j(xj,xlamr,xjp,zero,zero,zero)
t2=xf3j(one,xlama,one,zero,zero,zero)
term=t1*t2
if (term .ne. zero) then
  t3=xf6j(xja,xlama,xjap,one,half,one)
else
  return
endif
if (t3 .ne. zero) then
  term=term*t3
  if (.not.csflag) &
      t4=xf3j(xl,xlam12,xlp,zero,zero,zero)
  if (csflag) &
      t4=xf3j(xj12,xlam12,xj12p,-xnu,zero,xnu)
else
  return
endif
if (t4 .ne. zero) then
  term=term*t4
  if (.not.csflag) &
     t5=xf6j(xj12,xl,xjtot,xlp,xj12p,xlam12)
  if (csflag) t5=one
else
  return
endif
if (t5 .ne. zero) then
  term=term*t5
  tt= &
    xf9j(xj12,xlam12,xj12p,xja,xlama,xjap,xj,xlamr,xjp)
  term=term*tt
else
  return
endif
vee=term*xnorm*iphase
return
end
! --------------------------------------------------------------------
subroutine vlmh2pc (irow, icol, jtot, jlpar, j, jp, ja, jap, &
            k, kp, lamr, lama, mu, nu, jmol, flaghf, vee)
! --------------------------------------------------------------------
!  subroutine to evaluate the angular coupling matrix element for rotationally
!  inelastic collisions of a homonuclear molecule (j=1) and a 2P atom
!  dubernet and hutson case 1A
!  if flaghf = .true., spin included
!  if flaghf = .false., spinless atom

!  authors:  millard alexander
!  current revision date: 14-apr-1997
! --------------------------------------------------------------------
!  variables in call list:
!  irow, icol:  row and column of fully coupled states
!  jtot:      total angular momentum
!  j, jp      bra and ket values of molecular angular momentum
!  ja, jap:   bra and ket values of atomic angular momentum
!  k, kp:      bra and ket values of projection of jmolecular
!  lamr, lama, mu:      value of BF expansion indicies
!  nu         projection index (called P by Dubernet and Hutson)
!  vee:        on return:  contains desired coupling matrix element
! --------------------------------------------------------------------
implicit double precision (a-h,o-z)
logical flaghf
data  zero, one, two /0.d0, 1.d0, 2.d0/
data half /0.5d0/
vee = zero
!      if (jmol .eq. 0 .and. lamr .ne. 0) return
if (j .eq. 0 .and. jp .eq. 0 .and. lamr .ne. 0) return
iomeg=nu-k
iomegp=nu-kp
iphase=(-1)**(k+iomeg)
if (flaghf) then
  xomeg=iomeg+half
  xomegp=iomegp+half
  xja=ja+half
  xjap=jap+half
  xjtot=jtot+half
  xnu=nu+half
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
if (flaghf) then
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
    t5=xf6j(xja,xlama,xjap,one,half,one)
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
!  -----------------------------------------------------------------------
subroutine syh2p (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for homonuclear
!   + 2P atom scattering
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  current revision date: 29-nov-1995 by mha
! NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
!  -----------------------------------------------------------------------
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     rotational constant of molecule
!    aso:      spin-orbit constant of atom
!  variables in common block /cosysi/
!    nscod:    total number of variable names which are passed to HINPUT
!              nscod must equal isrcod + isicod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
!              para states will be included if iop=1 and only ortho states if
!              iop=-1
!    jmax:     the maximum rotational angular momenta for the diatomic
!  variable in common /cosys/
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters.  Note that the ordering
!             of the variable names in scod must correspond to the ordering
!             of the variable names in cosysi followed by the ordering of
!             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
logical readpt, existf
double precision brot, aso
character*1 dot
character*8 scod
character*(*) fname
character*60 line, filnam, potfil, filnm1
parameter (icod=3, ircod=2)
parameter (lencod = icod + ircod + 3)
#include "common/parbas.F90"
common /cosys/ scod(lencod)
common /cosysi/ nscode, isicod, nterm, iop, jmax
common /cosysr/ isrcod, junkr, brot, aso
save potfil
!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
!  NTERM must be the first variable
!  followed by the system dependent real variables
!  in the same order as in the common block /cosysr/
!  then the three variable names LAMMIN, LAMMAX, MPROJ, in that order
#include "common/comdot.F90"
scod(1)='NTERM'
scod(2)='IOP'
scod(3)='JMAX'
scod(4)='BROT'
scod(5)='ASO'
scod(6)='LAMMIN'
scod(7)='LAMMAX'
scod(8)='MPROJ'
nscode = lencod
isicod = icod
isrcod = ircod
irpot = 1
!  set default values for homonuclear+2P atom
  nterm = 1
  mproj(1) = 0
if (iread .eq. 0) then
  lammin(1)= 1
  lammax(1) = 5
  jmax = 0
  iop = 1
endif
potfil=' '
if (iread .eq. 0) return
!  line 18
read (8, *, err=80) iop, jmax
!  line 21
read (8, *, err=80) brot, aso
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  return
endif
read (8, 60, end=100) line
60 format (a)
goto 100
! here if read error occurs
80 write(6,90)
90 format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
return
!
entry ptrh2p (fname,readpt)
line = fname
readpt = .true.
100 if (readpt) then
  l=1
  call parse(line,l,filnam,lc)
  if(lc.eq.0) then
    write(6,102)
102     format(' FILENAME MISSING FOR POTENTIAL INPUT')
  end if
  j=index(filnam(1:lc),dot)
  if(j.eq.0) then
     call gennam(potfil,filnam,0,'BIN',lc)
     filnam = potfil
  end if
  potfil=filnam
  filnm1 = 'potdata/'//filnam
  inquire(file=filnm1,exist=existf)
  if(.not.existf) then
    write(6,105) filnam(1:lc)
105    format(' FILE ',(a),' NOT FOUND')
   return
  end if
! now call loapot(iunit,filnam) routine to read potential parameters
  call loapot(1,filnam)
endif
close (8)
return
!
entry savh2p (readpt)
!  save input parameters for symmetric top + atom scattering
!  the order of the write statements should be identical to the read statement
!  above. for consistency with the data file written by gendat, format
!  statements should reserve the first 30 spaces for data, spaces 31-33 should
!  be left blank, and the names of the variables should be printed in spaces
!  34-80
!  line 18:
write (8, 220) iop, jmax
220 format (4i4, 14x,'   iop, jmax')
!  line 21
write (8, 250) brot, aso
250 format(f12.4,f14.4, 6x, '   brot, aso')
write (8, 60) potfil
return
end
end module mod_hiba12_h2p