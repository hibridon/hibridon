#include "assert.h"
#include "unused.h"
module mod_hiba15_diat2p
  use mod_assert, only: fassert
contains
! sydiat2p (savdiat2p/ptrdiat2p) defines, saves variables and read       *
!                  potential for heteronuclear + 2P atom scattering      *
! --------------------------------------------------
subroutine badiat2p (bqs, jhold, ehold, ishold, nlevel, nlevop, &
                  isc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential
!  for collision of a doublet atom in a P state and a homonuclear molecule
!  authors:  millard alexander
!  current revision date:  7-apr-2003 by mha
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains atom+diatom angular momentum quantum number
!              for each channel (j12)
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains atomic angular momentum (ja) of each channel
!              (0 or 1)
!    jhold:    on return contains j12 angular momentum quantum number
!              for each level
!    ehold:    on return contains energy in hartrees of each level
!    ishold:   on return contains electronic angular momentum of each level
!    nlevel:   on return contains number of energetically distinct
!              levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              levels used in channel basis which are open
!              asymptotically
!    isc1,sc2: scratch vectors of length at least nmax
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
!
!  variable in module mod_conlam
!    nlam:      the number of case(a) interaction potentials actually used
!               if this has not been set in the pot subroutine, it it
!               here set to 5
!  variable in common block /coconv/
!     econv:    conversion factor from cm-1 to hartrees
!     xmconv:   converson factor from amu to atomic units
!  subroutines called:
!   vlmh2p:    returns angular coupling coefficient for particular
!              choice of channel index
! ------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam
use mod_hiba12_h2p, only: vlmh2p, vlmh2pc
use mod_cosysi, only: ispar
use mod_cosysr, only: rspar
use constants, only: econv, xmconv
use mod_par, only: iprint
use mod_parbas, only: lammax
use mod_par, only: boundc
use mod_ered, only: ered, rmu
use mod_jtot, only: jjtot, jjlpar
use mod_hitypes, only: bqs_type
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable, target :: v2
type(bqs_type), intent(out) :: bqs
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst
dimension  jhold(40), ehold(40), isc1(40), &
          sc2(40), sc3(40), &
          sc4(40), ishold(40)
dimension lamr(48),lama(48),lam12(48), mu(48)
dimension lamrold(40),lamaold(40),lam12old(40), muold(40)
data lamr  /1,2,3,4,5,6,7,8,9,10,11,12, &
            0,1,2,3,4,5,6,7,8,9,10,11,12, &
            1,2,3,4,5,6,7,8,9,10,11,12, &
             2,3,4,5,6,7,8,9,10,11,12/
data lama  /0,0,0,0,0,0,0,0,0,0,0,0, &
            2,2,2,2,2,2,2,2,2,2,2,2,2, &
            2,2,2,2,2,2,2,2,2,2,2,2, &
            2,2,2,2,2,2,2,2,2,2,2/
data lam12 /1,2,3,4,5,6,7,8,9,10,11,12, &
            2,1,3,0,2,4,1,3,5,2,4, &
            6,3,5,7,4,6,8,5,7,9,6,8,10, &
            7,9,11,8,10,12,9,11,13,10,12,14/
data mu    /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
            0,0,0,0,0,0,0,0,0,0, &
            1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2/
!  here is information on lambda terms for compatability with older
!  pot subroutines
data lamrold /1,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,10, &
              1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10/
data lamaold /0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
data lam12old /1,2,3,4,5,6,7,8,9,10,2,1,3,0,2,4,1,3,5,2,4, &
  6,3,5,7,4,6,8,5,7,9,6,8,10,7,9,11,8,10,12/
data muold /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
            1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2/

integer, pointer :: nterm, iop, jmax
real(8), pointer :: brot, aso
nterm=>ispar(1); iop=>ispar(2); jmax=>ispar(3)
brot=>rspar(1); aso=>rspar(2)
UNUSED_DUMMY(numin)
UNUSED_DUMMY(sc3)
UNUSED_DUMMY(sc4)


!   econv is conversion factor from cm-1 to hartrees
!   xmconv is converson factor from amu to atomic units
zero = 0.d0
half=0.5d0
two = 2.d0

if (.not.flaghf .and. (.not.csflag .or. csflag .and. ihomo))  then
  write (6, 5)
  write (9, 5)
5   format &
 (' *** FLAGHF = .FALSE. FOR CC OR CS1C MOL+ 2P ATOM; ABORT ***' )
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
9   format(' *** NTERM = ',i3,' .NE. 1 FOR MOL + 2P ATOM; ABORT')
  call exit
end if
if (nlam.eq.0) nlam=5
nsum=nlam
if (bastst) write (6, 14) nsum
write (9, 14) nsum
14 format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
if (clist) then
  if (flagsu) then
    write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
    write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16     format(/ &
  ' **  MOL + 2P ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4, &
  '             E=', f9.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)
  else
    if(.not.csflag) then
      write (6,20) &
          rmu * xmconv, ered * econv, jtot, jlpar
      write (9,20) &
          rmu * xmconv, ered * econv, jtot, jlpar
20       format(/,' **  CC MOL + 2P ATOM',' ** RMU=', f9.4, &
           '       E=',f9.2,'   JTOT=', i5, 2x,' JLPAR=',i2)
    else if (csflag .and. ihomo) then
      write (6,22) &
          rmu * xmconv, ered * econv, jtot, xnu
      write (9,22) &
          rmu * xmconv, ered * econv, jtot, xnu
22       format(/,' **  CD 1C MOL + 2P ATOM',' ** RMU=', f9.4, &
           '       E=',f9.2,'   JTOT=', i5, 2x,' P=',f4.1)
    else if (csflag .and. .not.ihomo) then
      write (6,23) &
          rmu * xmconv, ered * econv, jtot, xnu
      write (9,23) &
          rmu * xmconv, ered * econv, jtot, xnu
23       format(/,' **  CD 1A MOL + 2P ATOM',' ** RMU=', f9.4, &
           '       E=',f9.2,'   JTOT=', i5, 2x,' P=',f7.2)
    endif
  end if
  if (.not. flagsu) write (9,30) rcut
30   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
endif
! check for consistency in projection quantum number
if (nu .gt. jtot .and. csflag) then
  write (6, 35) nu, jtot
  write (9, 35) nu, jtot
35   format (' NU =',i3,' .GT. JTOT =',i3,'; JTOT RESET TO NU')
  jtot=nu
endif

!  save jtot and jlpar fr use later in transformatin
jjtot=jtot
jjlpar=jlpar
jmolmin=0
call bqs%init(nmax)
if (.not. csflag) then
!  assign quantum numbers and energies for CC calculations
  n=0
  nlevel = 0
! sum over molecular levels
  do 105 jmol=jmolmin,jmax
! sum over atomic levels
  do 100 jai=0,1
    nlevel = nlevel+1
    if (jai .eq. 0) ehold(nlevel)=-aso/econv
    if (jai .eq. 1) ehold(nlevel)=aso/(econv*two)
    ishold(nlevel)=jai
    xjai=jai+half
    j12min=abs(xjai-jmol)-half
    j12max=jai+jmol
    do 90 j12i=j12min,j12max
!repl 5/5/98 jhold(nlevel)=j12i
      jhold(nlevel)=jmol
      lmin=iabs(jtot-j12i)
      lmax=jtot+j12i+1
      do 80 li=lmin,lmax
! parity of space-frame states bqs%inq(-1)**(L-atom+jmol+L-orbital)
        if ((-1)**li .eq. jlpar*(-1)**(jmolmin-1)) then
          n = n + 1
          if (n .gt. nmax) go to 180
          bqs%lq(n) = li
          bqs%inq(n)=jai
!repl 5/5/98     bqs%jq(n)=j12i
          bqs%jq(n)=jmol
          bqs%length = n
!repl 5/5/98     isc1(n)=jmol
          isc1(n)=j12i
! centrifugal contribution to hamiltonian
          cent(n) = li*(li+1)
! constant channel energy
          sc2(n)=brot*jmol*(jmol+1)
          if (jai .eq. 0) sc2(n)=sc2(n)-aso
          if (jai .eq. 1) sc2(n)=sc2(n)+aso/two
          sc2(n)=sc2(n)/econv
        endif
        if (n>0) ehold(nlevel)=sc2(n)
80       continue
90     continue
100   continue
105   continue
! here for CS calculations
elseif (csflag) then
!  assign quantum numbers and energies for CS calculations
!  here for Dubernet and Hutson case 1C
  if (ihomo) then
    n=0
    nlevel = 0
    xnu=nu+half
! sum over molecular levels
    do 125 jmol=jmolmin,jmax
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
          ishold(nlevel)=jai
          jhold(nlevel)=j12i
          li=jtot
          xjtot=jtot+half
          n = n + 1
          if (n .gt. nmax) go to 180
          bqs%lq(n) = li
          bqs%inq(n)=jai
          bqs%jq(n)=j12i
          bqs%length = n
          isc1(n)=jmol
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
      do 155 jmol=jmolmin,jmax
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
              isc1(n)=jmol
              bqs%lq(n) = li
              xjtot=jtot+half
              bqs%inq(n)=jai
              bqs%jq(n)=ik
              bqs%length = n
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
      do 165 jmol=jmolmin,jmax
! sum over molecular projection levels
        kmin=0
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
            isc1(n)=jmol
            bqs%lq(n) = li
            bqs%inq(n)=1
            bqs%jq(n)=ik
            bqs%length = n
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

180 if (n .gt. nmax) then
  write (9, 185) n, nmax
  write (6, 185) n, nmax
185   format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF', &
         i4,'; ABORT')
  if (bastst) then
    return
  else
    call exit
  end if
end if
!  now check to see if any of the open channels are closed at r=rcut
!  this is not done for molecule-surface collisions or for rcut < 0
!  or for bound state calculations!
if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
  emin = 1.e+7
  do 190  i = 1, n
    if (sc2(i) .le. ered) then
!  here if channel is
      if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut) &
          .gt. (ered - sc2(i)) ) then
!  here if channel is open asymptotically but closed at r = rcut
        if (sc2(i) .lt. emin) emin = sc2(i)
!  emin now contains the lowest channel energy for which this
!  condition is met
      end if
    end if
190   continue
!  now eliminate all channels with sc2(eint) .ge. emin if any of the channels
!  are open asymptotically but closed at r = rcut
  if (emin.lt.ered) then
    nn = 0
    do 195 i = 1, n
      if (sc2(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        sc2(nn) = sc2(i)
        bqs%inq(nn) = bqs%inq(i)
        bqs%jq(nn) = bqs%jq(i)
        cent(nn) = cent(i)
        bqs%lq(nn) = bqs%lq(i)
      end if
195     continue
!  reset number of channels
    n = nn
    bqs%length = n
  end if
end if
!  return if no channels
if (n .eq. 0) return
!  eliminate closed channels from level list
! sort level list in terms of increasing energy
ilow = 1
if (nlevel.gt.1) then
   do i=1,nlevel-1
      elow=ehold(i)
      do ji=i+1,nlevel
         ex=ehold(ji)
         if (ex.lt.elow) then
           ilow=ji
           elow=ex
         endif
      enddo
! ilow is now index of lowest energy
! swap levels
      if (ilow.ne.i) then
         eswap=ehold(i)
         ehold(i)=ehold(ilow)
         ehold(ilow)=eswap
         jswap=jhold(i)
         jhold(i)=jhold(ilow)
         jhold(ilow)=jswap
         iswap=ishold(i)
         ishold(i)=ishold(ilow)
         ishold(ilow)=iswap
      endif
   enddo
endif
nlevop=0
do i=1,nlevel
   if (ered-ehold(i).gt.0d0) then
      nlevop=nlevop+1
   endif
enddo

ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1
!  shift so that zero of energy is average spin-orbit state
!  e.g. j=3/2 are shifted by aso/2 and j=1/2 by -aso
 jmol=0
emin=0
erot=brot*(jmol+1)*jmol
do 205 i=1,n
  eint(i)=sc2(i)-(erot+emin)/econv
205 continue
!  now list channels if requested
if (clist) then
  if (ihomo) then
    if (.not.csflag) then
      write (6, 210)
      write (9, 210)
210       format(/'   N   JA(IS) JMOL  J12(J)    L      EINT(CM-1)')
      do 220  i = 1, n
        write (6, 215) i,bqs%inq(i)+half,bqs%jq(i),isc1(i)+half,bqs%lq(i), &
             eint(i)*econv
        write (9, 215) i,bqs%inq(i)+half,bqs%jq(i),isc1(i)+half,bqs%lq(i), &
             eint(i)*econv
215         format (i4, f8.1,i5, f7.1, i5, f14.3)
220       continue
    else
      write (6, 225)
      write (9, 225)
225       format(/'   N   JA(IS) JMOL  J12(J)   CENT   EINT(CM-1)', &
            ';  L = JTOT')
      do 235  i = 1, n
        xs=bqs%inq(i)+half
        xj=bqs%jq(i)+half
        icent=cent(i)
        write (6, 230) i,xs,isc1(i),xj,icent,eint(i)*econv
        write (9, 230) i,xs,isc1(i),xj,icent,eint(i)*econv
230         format (i4, f8.1,i6,f7.1, i6, f14.3)
235       continue
    endif
  elseif (.not.ihomo) then
    if (flaghf) then
      write (6, 255)
      write (9, 255)
255       format( &
       /'   N   JA(IS)  JMOL  K(J)    L    CENT    EINT(CM-1)')
      do 265  i = 1, n
        xs=bqs%inq(i)+half
        icent=cent(i)
        write (6, 260) i,xs,isc1(i),bqs%jq(i),bqs%lq(i),icent,eint(i)*econv
        write (9, 260) i,xs,isc1(i),bqs%jq(i),bqs%lq(i),icent,eint(i)*econv
260         format (i4, f8.1,4i6, f14.3)
265       continue
    else
      write (6, 270)
      write (9, 270)
270       format(/'   N   LA(IS)  JMOL  K(J)   CENT    EINT(CM-1)')
      do 275  i = 1, n
        write (6, 272) i,bqs%inq(i),isc1(i),bqs%jq(i),idnint(cent(i)), &
                       eint(i)*econv
        write (9, 272) i,bqs%inq(i),isc1(i),bqs%jq(i),idnint(cent(i)), &
                       eint(i)*econv
272         format (i4, i8,3i6, f14.3)
275       continue
    end if

  endif
end if
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
if (bastst .and. iprint .gt. 1) then
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
i = 0
ilam=0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 320 il = 1,lammax(1)
!      do 320 il = 4,4
  ilam=ilam+1
  ancouma => v2%get_angular_coupling_matrix(ilam)
  inum = 0
  ij=0
  if (lammax(1) .eq. 48) then
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
  do icol= 1, n
    do irow = icol, n
      ij = ntop * (icol - 1) +irow
      !repl 5/5/98jrow=isc1(irow)
      !repl 5/5/98jcol=isc1(icol)
      ! 6/13/2002:  reverse definition of j and j12 for CD calculations
      if (csflag) then
         jrow=isc1(irow)
         jcol=isc1(icol)
         j12row=bqs%jq(irow)
         j12col=bqs%jq(icol)
      else
         jrow=bqs%jq(irow)
         jcol=bqs%jq(icol)
         j12row=isc1(irow)
         j12col=isc1(icol)
      endif
      if (.not. csflag .or. (csflag .and. ihomo)) then
        call vlmh2p (jtot, jrow, jcol,bqs%inq(irow), &
        bqs%inq(icol), j12row, j12col, bqs%lq(irow), bqs%lq(icol), ilamr, &
        ilama, ilam12, nu, csflag, vee)
      else if (csflag .and. .not.ihomo) then
        call vlmh2pc(jtot, jrow, jcol,bqs%inq(irow), &
        bqs%inq(icol), bqs%jq(irow), bqs%jq(icol), ilamr, &
        ilama, imu, nu, flaghf, vee)
      endif
      if (vee .ne. 0) then
        i = i + 1
        inum = inum + 1
        call ancouma%set_element(irow=irow, icol=icol, vee=vee)
        if (bastst .and. iprint .gt. 1) then
          if (.not. csflag .or. (csflag .and. ihomo)) then
            write (6, 290) ilam, ilamr,ilama,ilam12, icol, irow, i, vee
            write (9, 290) ilam, ilamr,ilama,ilam12, icol, irow, i, vee
290               format (i4, i5,2i3,i5, 2i6, i6, g17.8)
          elseif (csflag .and. .not.ihomo) then
            write (6, 290) ilam, ilamr,ilama,imu, icol, irow, i, vee
            write (9, 290) ilam, ilamr,ilama,imu, icol, irow, i, vee
          endif
        endif
      end if
    end do
  end do
if (bastst .and. iprint .gt. 1) then
  write (6, 315) ilam, ancouma%get_num_nonzero_elements()
  write (9, 315) ilam, ancouma%get_num_nonzero_elements()
315   format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
end if
320 continue
if (clist) then
  write (6, 360) i
  write (9, 360) i
360   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', &
           i6)
end if
return
end
!  -----------------------------------------------------------------------
subroutine sydiat2p (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for heteronuclear
!   + 2P atom scattering
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  current revision date: 22-Sep-2005 by Jacek Klos
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
!    jmax:     the maximum rotational angular momenta for the diatomic
!  variable in common bloc /cosys/
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters.  Note that the ordering
!             of the variable names in scod must correspond to the ordering
!             of the variable names in cosysi followed by the ordering of
!             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
!  -----------------------------------------------------------------------
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, rspar
use funit, only: FUNIT_INP
use mod_parbas, only: lammin, lammax, mproj
use mod_hiutil, only: gennam, get_token
use mod_hipot, only: loapot
implicit none
integer, intent(inout) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: icod, ircod, j, l, lc
logical existf
character*1 dot
character*(*) fname
character*60 line, filnam, potfil
character*68 filnm1
parameter (icod=3, ircod=2)
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

integer, pointer, save :: nterm, iop, jmax
real(8), pointer, save :: brot, aso
nterm=>ispar(1); iop=>ispar(2); jmax=>ispar(3)
brot=>rspar(1); aso=>rspar(2)


scod(1)='NTERM'
scod(2)='IOP'
scod(3)='JMAX'
scod(4)='BROT'
scod(5)='ASO'
scod(6)='LAMMIN'
scod(7)='LAMMAX'
scod(8)='MPROJ'
nscode = icod + ircod + 3
isicod = icod
isrcod = ircod
irpot = 1
!  set default values for heteronuclear+2P atom
  nterm = 1
  mproj(1) = 0
if (iread .eq. 0) then
  lammin(1)= 1
  lammax(1) = 40
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
entry ptrdiat2p (fname,readpt)
line = fname
readpt = .true.
100 if (readpt) then
  l=1
  call get_token(line,l,filnam,lc)
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
entry savdiat2p ()
!  save input parameters for heteronuclear diatom + 2P atom scattering
!  the order of the write statements should be identical to the read statement
!  above. for consistency with the data file written by gendat, format
!  statements should reserve the first 30 spaces for data, spaces 31-33 should
!  be left blank, and the names of the variables should be printed in spaces
!  34-80
!  line 18:
write (FUNIT_INP, 220) iop, jmax
220 format (4i4, 14x,' iop, jmax')
!  line 21
write (FUNIT_INP, 250) brot, aso
250 format(f12.4,f14.4, 6x, '   brot, aso')
write (FUNIT_INP, 60) potfil
return
end
end module mod_hiba15_diat2p
