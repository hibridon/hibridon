#include "assert.h"
module mod_hiba10_22p
contains
! sy22p (sav22p/ptr22p) defines, save variables and reads                *
!                  potential for 2S / 2P atom scattering                 *
! --------------------------------------------------
subroutine ba22p (j, l, is, jhold, ehold, ishold, nlevel, nlevop, &
                  isc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential
!  for collision of a doublet atom in an S state and a doublet atom in a
!  P state
!  authors:  brigitte pouilly and millard alexander
!  current revision date:  13-may-1997 by mha
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains electronic angular momentum quantum number
!              for each channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains electronic spin of each channel (0 or 1)
!    jhold:    on return contains electronic angular momentum quantum number
!              for each level
!    ehold:    on return contains energy in hartrees of each level
!    ishold:   on return contains spin multiplicity of each level
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
!              case here since we are labelling asymptotic state)
!              if .false., then system with integer spin
!    flagsu:   if .true., then molecule-surface collisons
!    csflag:   if .true., then coupled-states calculation
!              if .false., the close-coupled calculation
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true. , then homonuclear molecule
!              if .false., then heteronuclear molecule
!              if the molecule is homonuclear (ihomo = .true.), the
!              rotational levels included go from jmin to jmax in steps
!              of 2 and only even lambda terms in the anisotropy are
!              included
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              only those channels are included for which
!                  (-1)**(l-jtot)=jlpar
!              n.b. jlpar=+1 corresponds to f levels, jlpar=-1, to e levels
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
!    aso:      spin-orbit constant of 2P atom (cm-1)
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different types of electronic coupling terms
!              this should be 1 here
!               in case (a) basis
!
!  variables in common block /coered/
!    ered:      collision energy in atomic units (hartrees)
!    rmu:       collision reduced mass in atomic units
!               (mass of electron = 1)
!  variable in module mod_conlam
!    nlam:      the number of case(a) interaction potentials actually used
!               this is 14 here
!  variable in common block /coconv/
!     econv:    conversion factor from cm-1 to hartrees
!     xmconv:   converson factor from amu to atomic units
!  subroutines called:
!   vlm22p:    returns angular coupling coefficient for particular
!              choice of channel index
! ------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use constants, only: econv, xmconv
#include "common/parbasl.F90"
implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst
#include "common/parbas.F90"

common /coered/ ered, rmu
common /coskip/ nskip, iskip
integer :: nskip, iskip
common /cojtot/ jjtot,jjlpar
dimension j(9), l(9), jhold(9), ehold(9), isc1(9), sc2(9), sc3(9), &
          sc4(9), ishold(9), is(9)
integer, pointer :: nterm, nphoto
real(8), pointer :: aso
nterm=>ispar(1); nphoto=>ispar(2)
aso=>rspar(1)

zero = 0.d0
half=0.5d0
two = 2.d0
!  check for consistency in the values of flaghf and csflag
if (.not.flaghf) then
  write (6, 5)
  write (9, 5)
5   format (' *** FLAGHF = .FALSE. FOR 2/2 ATOM; ABORT ***' )
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (csflag) then
  write (6, 8)
  write (9, 8)
8  format &
   ('  *** CSFLAG SET .FALSE. FOR 2/2 ATOM CALCULATION ***')
  csflag=.false.
end if
if(nterm.ne.1) then
   write(6,9) nterm
   write(9,9) nterm
9    format(' *** NTERM = ',i3,' .NE. 1 FOR 2/2 ATOM; ABORT')
   call exit
end if
nsum=14
if (bastst) write (6, 14) nsum
write (9, 14) nsum
14 format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
nlam = nsum
if (clist) then
  if (flagsu) then
    write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
    write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16     format(/' **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4, &
      '             E=', f9.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)
  else
      write (6,20) &
          rmu * xmconv, ered * econv, jtot, jlpar
      write (9,20) &
          rmu * xmconv, ered * econv, jtot, jlpar
20       format(/,' **  CC 2/2 ATOM',' ** RMU=', f9.4, &
           '       E=',f9.2,'   JTOT=', i5, 2x,' JLPAR=',i2)
  end if
  if (.not. flagsu) write (9,30) rcut
30   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
endif
!  assign quantum numbers and energies in case (e) basis
n=0
nlevel = 0
jmin=0
jmax=1
!  save jtot and jlpar fr use later in transformatin
jjtot=jtot
jjlpar=jlpar
do 120 ji=jmin, jmax
  nlevel = nlevel+1
  if (ji .eq. 0) ehold(nlevel)=aso/econv
  if (ji .eq. 1) ehold(nlevel)=-aso/(econv*two)
  jhold(nlevel)=ji
  ishold(nlevel)=1
  if (jlpar .eq. -1) then
! here for e levels
    j12max=ji+1
    j12min=j12max-ji
    lmin=iabs(jtot-1)
    lmax=jtot+1
    do 90 j12=j12min,j12max
      do 80 li=lmin, lmax, 2
        n = n + 1
        if (n .gt. nmax) go to 180
        l(n) = li
! centrifugal contribution to hamiltonian in case (a) basis
        cent(n) = jtot*(jtot+1)
! add on constant terms in various cases
! here for 3Pi0 and 3Pi1
        if (n .eq. 2 .or. n .eq. 3) cent(n)=cent(n)+1
! here for 3Pi2
        if (n .eq. 4) cent(n)=cent(n)-3
! here for 1Pi1
        if (n .eq. 6) cent(n)=cent(n)-1
        is(n)=1
        isc1(n) = j12
        j(n) = ji
! constant channel energy is zero in case a
        eint(n)=zero
        if (ji .eq. 0) sc2(n)=aso/econv
        if (ji .eq. 1) sc2(n)=-aso/(econv*two)
80       continue
90     continue
! here for f levels
  else
    j12max=ji+1
    j12min=j12max-1
    do 110 j12=j12min,j12max
      if (j12 .lt. 2) then
        lmin=jtot
        lmax=jtot
      else
        lmin=iabs(jtot-2)
        lmax=jtot+2
      endif
      do 100 li=lmin, lmax,2
        n = n + 1
        if (n .gt. nmax) go to 180
        l(n) = li
! centrifugal contribution to hamiltonian in case (a) basis
        cent(n) = jtot*(jtot+1)
! add on constant terms in various cases
! here for 3Pi0 and 3Pi1
        if (n .eq. 2 .or. n .eq. 3) cent(n)=cent(n)+1
! here for 3Pi2
        if (n .eq. 4) cent(n)=cent(n)-3
! here for 1Pi1
        if (n .eq. 6) cent(n)=cent(n)-1
        is(n)=1
        isc1(n) = j12
        j(n) = ji
! constant channel energy is zero in case a
        eint(n)=zero
        if (ji .eq. 0) sc2(n)=aso/econv
        if (ji .eq. 1) sc2(n)=-aso/(econv*two)
100        continue
110     continue
  endif
120 continue
180 if (n .gt. nmax) then
  write (9, 185) n, nmax
  write (6, 185) n, nmax
185   format(/' *** NCHANNELS=', i3,' .GT. MAX DIMENSION OF', &
         i3,'; ABORT')
  if (bastst) then
    return
  else
    call exit
  end if
end if
!  now check to see if any of the open channels are closed at r=rcut
!  this is not done for molecule-surface collisions or for rcut < 0
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
nlevop = nlevel
ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1

!  now list channels if requested
if (clist) then
  write (6, 255)
  write (9, 255)
255   format(/'   N    J   J12    L      EINT(CM-1)')
  do 265  i = 1, n
    write (6, 260) i, j(i)+half, isc1(i), l(i), sc2(i) * econv
    write (9, 260) i, j(i)+half, isc1(i), l(i), sc2(i) * econv
260     format (i4, f5.1,i6, i5, f14.3)
265   continue
end if
!  now calculate coupling matrix elements
!  this assumes case (a) basis which is, for jlpar = -1 (e levels)
!  3Sig1, 3Pi0, 3Pi1, 3Pi2, 1Sig, 1Pi1
!  and, for jlpar = +1 (f levels)
!  3Sig1, 3Pi0, 3Pi1, 3Pi2, 3Sig0, 1Pi1
!  ordering of terms is as follows:
!  lam = 1; W(1Sig)
!  lam = 2; W(1Pi)
!  lam = 3; W(3Sig)
!  lam = 4; W(3Pi)
!  lam = 5; BLperp^2 for 1sig
!  lam = 6; BLperp^2 for 1Pi
!  lam = 7; BLperp^2 for 3sig
!  lam = 8; BLperp^2 for 3Pi
!  lam = 9 ; <1Sig|l+|1Pi> term
!  lam = 10; <3Sig|l+|3Pi> term
!  lam = 11; constant spin-orbit term
!  lam = 12; <1Sig|l+|1Pi> spin-orbit term
!  lam = 13; <3Sig|l+|3Pi> spin-orbit term
!  lam = 14: term proportional just to B
if (bastst) then
  write (6, 280)
  write (9, 280)
280   format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
end if
! i counts v2 elements
! inum counts v2 elements for given lambda
! ilam counts number of v2 matrices
i = 0
ilam=0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 320 il = 1,lammax(1)
  lb = il
  ilam=ilam+1
  ancouma => v2%get_angular_coupling_matrix(ilam)
  inum = 0
  do icol= 1, n
    do irow = icol, n
      call vlm22p (irow, icol, jtot, jlpar, lb, vee)
      if (vee .ne. 0) then
        i = i + 1
        inum = inum + 1
          call ancouma%set_element(irow=irow, icol=icol, vee=vee)
          if (bastst) then
            write (6, 290) ilam, lb, icol, irow, i, vee
            write (9, 290) ilam, lb, icol, irow, i, vee
290             format (i4, 2i7, 2i6, g17.8)
          end if
        end if
      end do
    end do
if (bastst) then
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
! --------------------------------------------------------------------
subroutine vlm22p (irow, icol, jtot, jlpar, lam, vee)
! --------------------------------------------------------------------
!  subroutine to evaluate the angular coupling matrix element for rotationally
!  inelastic collisions of a doublet S and double P atom
!  structureless target

!  authors:  brigitte pouilly and millard alexander
!  current revision date: 29-sep-92
! --------------------------------------------------------------------
!  variables in call list3:
!  irow, icol:  row and column of case a state
!  jlpar:  parity (-1 for e, +1 for f)
!  jtot:     total angular momentum
!  lam:      value of expansion index lambda:
!  this assumes case (a) basis which is, for jlpar = -1 (e levels)
!  3Sig1, 3Pi0, 3Pi1, 3Pi2, 1Sig, 1Pi1
!  and, for jlpar = +1 (f levels)
!  3Sig1, 3Pi0, 3Pi1, 3Pi2, 3Sig0, 1Pi1
!  ordering of terms is as follows:
!  lam = 1; W(1Sig)
!  lam = 2; W(1Pi)
!  lam = 3; W(3Sig)
!  lam = 4; W(3Pi)
!  lam = 5; BLperp^2 for 1sig
!  lam = 6; BLperp^2 for 1Pi
!  lam = 7; BLperp^2 for 3sig
!  lam = 8; BLperp^2 for 3Pi
!  lam = 9 ; <1Sig|l+|1Pi> term
!  lam = 10; <3Sig|l+|3Pi> term
!  lam = 11; constant spin-orbit term
!  lam = 12; <1Sig|l+|1Pi> spin-orbit term
!  lam = 13; <3Sig|l+|3Pi> spin-orbit term
!  lam = 14:  term proportional just to B
!  vee:        on return:  contains desired coupling matrix element
! --------------------------------------------------------------------
implicit double precision (a-h,o-z)
data  zero, one, two, three /0.d0, 1.d0, 2.d0, 3.d0/
data half /0.5d0/
sq2=sqrt(two)
xjjp1=jtot*(jtot+1)
sq=sqrt(xjjp1)
sqm2=sqrt(xjjp1-2)
vee = zero
!  here for f states
if (jlpar .gt. 0) then
  if (lam .eq. 2) then
!  here for diagonal terms for singlet pi
    if (irow .eq. 6 .and. icol .eq. 6) vee=one
!  here for diagonal terms for triplet sigma
  else if (lam .eq. 3) then
!  3Sig1 - 3Sig1
    if (irow .eq. 1 .and. icol .eq. 1) vee=one
!  3Sig0 - 3Sig0
    if (irow .eq. 5 .and. icol .eq. 5) vee=one
!  here for diagonal terms for triplet pi
  else if (lam .eq. 4) then
    if (irow .eq. icol .and. (irow .eq. 2 .or. irow .eq. 3 &
        .or. irow .eq. 4)) vee=one
!  here for diagonal terms of l^2 in singlet pi
  else if (lam .eq. 6) then
    if (irow .eq. icol .and. irow .eq. 6) vee=one
!  here for diagonal terms of l^2 in triplet sigma
  else if (lam .eq. 7) then
!  3Sig1 - 3Sig1
    if (irow .eq. 1 .and. icol .eq. 1) vee=one
!  3Sig0 - 3Sig0
    if (irow .eq. 5 .and. icol .eq. 5) vee=one
!  here for diagonal terms of l^2 in triplet pi
  else if (lam .eq. 8) then
    if (irow .eq. icol .and. (irow .eq. 2 .or. irow .eq. 3 &
        .or. irow .eq. 4)) vee=one
!  here for triplet sigma / triplet pi rotational coupling
  else if (lam .eq. 10) then
!  3Sig1 - 3Pi0
    if (icol .eq. 1 .and. irow .eq. 2) vee=two*sq
!  3Sig1 - 3Pi1
    if (icol .eq. 1 .and. irow .eq. 3) vee=-two*sq2
!  3Sig1 - 3Pi2
    if (icol .eq. 1 .and. irow .eq. 4) vee=two*sqm2
!  3Pi0 - 3Sig0
    if (icol .eq. 2 .and. irow .eq. 5) vee=-4.d0
!  3P11 - 3Sig0
    if (icol .eq. 3 .and. irow .eq. 5) vee=two*sq2*sq
!  here for constant spin-orbit term
  else if (lam .eq. 11) then
!  3Pi1 - 1Pi1
    if (icol .eq. 3 .and. irow .eq. 6) vee=-half
!  3Pi0 - 3Pi0
    if (icol .eq. 2 .and. irow .eq. 2) vee=half
!  3Pi2 - 3Pi2
    if (icol .eq. 4 .and. irow .eq. 4) vee=-half
!  here for triplet spin-orbit term
  else if (lam .eq. 13) then
!  3Sig1 - 3Pi1
    if (icol .eq. 1 .and. irow .eq. 3) vee=one/sq2
!  3Sig1 - 1Pi1
    if (icol .eq. 1 .and. irow .eq. 6) vee=-one/sq2
!  3Pi0 - 3Sig0
    if (icol .eq. 2 .and. irow .eq. 5) vee=one
  else if (lam .eq. 14) then
!  3Sig1 - 3Sig0
    if (icol .eq. 1 .and. irow .eq. 5) vee=-two*sq
!  3Pi0 - 3Pi1
    if (icol .eq. 2 .and. irow .eq. 3) vee=-sq2*sq
!  3Pi1 - 3Pi2
    if (icol .eq. 3 .and. irow .eq. 4) vee=-sq2*sqm2
endif
! here for e-labelled levels
else if (jlpar .lt. 0) then
!  here for diagonal terms for singlet sigma
  if (lam .eq. 1) then
!  1Sig - 1Sig electrostatic potential
    if (icol .eq. 5 .and. irow .eq. 5) vee=one
  else if (lam .eq. 2) then
!  here for diagonal terms for singlet pi
    if (irow .eq. 6 .and. icol .eq. 6) vee=one
!  here for diagonal terms for triplet sigma
  else if (lam .eq. 3) then
    if (irow .eq. 1 .and. icol .eq. 1) vee=one
!  here for diagonal terms for triplet pi
  else if (lam .eq. 4) then
    if (irow .eq. icol .and. (irow .eq. 2 .or. irow .eq. 3 &
        .or. irow .eq. 4)) vee=one
!  here for diagonal terms of l^2 in singlet sigma
  else if (lam .eq. 5) then
    if (irow .eq. icol .and. irow .eq. 5) vee=one
!  here for diagonal terms of l^2 in singlet pi
  else if (lam .eq. 6) then
    if (irow .eq. icol .and. irow .eq. 6) vee=one
!  here for diagonal terms of l^2 in triplet sigma
  else if (lam .eq. 7) then
    if (irow .eq. icol .and. irow .eq. 1) vee=one
!  here for diagonal terms of l^2 in triplet pi
  else if (lam .eq. 8) then
    if (irow .eq. icol .and. (irow .eq. 2 .or. irow .eq. 3 &
        .or. irow .eq. 4)) vee=one
!  here for singlet sigma / singlet pi rotational coupling
  else if (lam .eq. 9) then
    if (icol .eq. 5 .and. irow .eq. 6) vee=two*sq2*sq
!  here for triplet sigma / triplet pi rotational coupling
  else if (lam .eq. 10) then
!  3Sig1 - 3Pi0
    if (icol .eq. 1 .and. irow .eq. 2) vee=-two*sq
!  3Sig1 - 3Pi1
    if (icol .eq. 1 .and. irow .eq. 3) vee=-two*sq2
!  3Sig1 - 3Pi2
    if (icol .eq. 1 .and. irow .eq. 4) vee=two*sqm2
!  here for constant spin-orbit term
  else if (lam .eq. 11) then
!  3Pi1 - 1Pi1
    if (icol .eq. 3 .and. irow .eq. 6) vee=-half
!  3Pi0 - 3Pi0
    if (icol .eq. 2 .and. irow .eq. 2) vee=half
!  3Pi2 - 3Pi2
    if (icol .eq. 4 .and. irow .eq. 4) vee=-half
!  here for singlet spin-orbit term
  else if (lam .eq. 12) then
!  3Pi0 - 1Sig0
    if (icol .eq. 2 .and. irow .eq. 5) vee=one
!  here for triplet spin-orbit term
  else if (lam .eq. 13) then
!  3Sig1 - 3Pi1
    if (icol .eq. 1 .and. irow .eq. 3) vee=one/sq2
!  3Sig1 - 1Pi1
    if (icol .eq. 1 .and. irow .eq. 6) vee=-one/sq2
  else if (lam .eq. 14) then
!  3Pi0 - 3Pi1
    if (icol .eq. 2 .and. irow .eq. 3) vee=-sq2*sq
!  3Pi1 - 3Pi2
    if (icol .eq. 3 .and. irow .eq. 4) vee=-sq2*sqm2
  end if
endif
return
end
! -----------------------------------------
subroutine tcasea22(j,jlpar)
! -----------------------------------------
!   matrix to transform from case a to case e, solely for 22P scattering
!   this matrix, tatoe, for total-j=j is returned
!   author:  brigitte pouilly and millard alexander
!   latest revision date:  30-dec-1995
! -----------------------------------------
implicit double precision (a-h,o-z)
common /cotrans/ t(6,6)
data zero, one,two ,three,six/0.d0, 1.d0, 2.d0, 3.d0, 6.d0/
if (j .lt. 2) then
! error if jtot is less than 2
  write (6, 5) jtot
5   format (' *** JTOT = ',i2,' .LT. 2 IN TCASEA; ABORT **')
  stop
endif
! initialization of the matrix tatoe
call dset(36,zero,t,1)
!
xj=j
xjp1=j+1
xjp2=j+2
xjm1=j-1
if(jlpar.lt.0) then
! here for e-labelled states
! transformation from asymptotic states (rows) into case (a) states (columns)
! with ordering (of columns)
!  3Sig1, 3Pi0, 3Pi1, 3Pi2, 1Sig, 1Pi1
  t(1,1)=-sqrt(xjp1/three)
  t(1,2)= sqrt(two*xj/three)
  t(1,3)= sqrt(xjp1/three)
  t(1,5)=-sqrt(xj/three)
  t(1,6)=-sqrt(xjp1/three)
  t(2,1)=t(1,5)
  t(2,2)=-sqrt(two*xjp1/three)
  t(2,3)= sqrt(xj/three)
  t(2,5)= sqrt(xjp1/three)
  t(2,6)=-sqrt(xj/three)
  t(3,5)= sqrt(two*xj/three)
  t(3,6)= sqrt(two*xjp1/three)
  t(3,1)=-sqrt(xjp1/six)
  t(3,2)= sqrt(xj/three)
  t(3,3)= sqrt(xjp1/six)
  t(4,5)=-sqrt(2*xjp1/three)
  t(4,6)= sqrt(2*xj/three)
  t(4,1)=-sqrt(xj/six)
  t(4,2)=-sqrt(xjp1/three)
  t(4,3)= sqrt(xj/six)
  t(5,1)= sqrt(xjm1/two)
  t(5,3)=t(5,1)
  t(5,4)= sqrt(xj+two)
  t(6,1)=-sqrt((xj+two)/two)
  t(6,3)=t(6,1)
  t(6,4)= sqrt(xjm1)
  denom=one/sqrt(two*xj+one)
  call dscal(36,denom,t,1)
! here for f-labelled states
!  3Sig1, 3Pi0, 3Pi1, 3Pi2, 3Sig0, 1Pi1

else
  t(1,5)=-sqrt(one/three)
  t(1,2)= sqrt(two/three)
  t(2,1)=t(1,5)
  t(2,3)=-t(1,5)
  t(2,6)=t(1,5)
  t(3,1)=-sqrt(one/six)
  t(3,3)= sqrt(one/six)
  t(3,6)=sqrt(two/three)
  t(4,5)=sqrt(xj*xjm1)
  t(4,1)= sqrt(xjp1*xjm1)
  t(4,2)=sqrt(xj*xjm1/two)
  t(4,3)= sqrt(xjp1*xjm1)
  t(4,4)=sqrt(xjp1*xjp2/two)
  denom=one/sqrt((two*xj+one)*(two*xj-one))
  call dscal(6,denom,t(4,1),6)
  t(5,5)=-sqrt(two*xj*xjp1/three)
  t(5,1)=-sqrt(three/two)
  t(5,2)=-sqrt(xjp1*xj/three)
  t(5,3)=-sqrt(three/two)
  t(5,4)=sqrt(three*xjp2*xjm1)
  denom=one/sqrt((two*xj+three)*(two*xj-one))
  call dscal(6,denom,t(5,1),6)
  t(6,5)=sqrt(xjp1*xjp2)
  t(6,1)=-sqrt(xj*xjp2)
  t(6,2)= sqrt(xjp1*xjp2/two)
  t(6,3)=-sqrt(xjp2*xj)
  t(6,4)=sqrt(xj*xjm1/two)
  denom=one/sqrt((two*xj+one)*(two*xj+three))
  call dscal(6,denom,t(6,1),6)
endif
return
end
! -----------------------------------------
subroutine trans22(w,n,nmax)
!  --------------------------------------------
!  to transform w-matrix from case (a) to case (e)
!  this is returned in w
!  b:  scratch matrices
!  also restore correct case (e) channel energies
!  this is -aso for channels 1 and 2 and +0.5 aso for channels (3-6)

!  author:  millard alexander
!  latest revision date:  4-oct-1992
!  --------------------------------------------
use mod_coeint, only: eint
use mod_cosysr, only: isrcod, junkr, rspar
use constants, only: econv, xmconv
implicit double precision (a-h,o-z)
#if defined(HIB_UNIX_IBM)
character*1 forma, formb
#endif
common /cotrans/ t(6,6)
common /cojtot/ j,jlpar
dimension w(42), sc(49)
data one,half /1.d0,0.5d0/
data nnmax /6/

real(8), pointer :: aso
aso=>rspar(1)

call tcasea22(j,jlpar)
if (n .ne. 6) then
  write (6, 10) n
  write (9, 10) n
10   format ('*** N =',i3,' .NE. 6 IN TRANS22; ABORT ***')
  stop
endif
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
 call mxma (t,1,nnmax,w,1,nmax,sc,1,nmax,n,n,n)
 call mxma (sc,1,nmax,t,nnmax,1,w,1,nmax,n,n,n)
#endif
#if defined(HIB_UNIX_IBM)
forma='N'
formb='T'
call dgemul (t,nnmax,forma,w,nmax,forma,sc,nmax,n,n,n)
call dgemul (sc,nmax,forma,t,nnmax,formb,w,nmax,n,n,n)
#endif
return
!  w now contains desired product
! restore internal energies of each channel
entry energ22
eint(1)=one
eint(2)=one
eint(3)=-half
eint(4)=-half
eint(5)=-half
eint(6)=-half
fact=aso/econv
call dscal(6,fact,eint,1)

return
end
! -----------------------------------------------------------------------
subroutine sy22p (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for collisions of
! atom in doublet S state with atom in doublet P state
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  current revision date: 14-mar-1997 by mha
!  -----------------------------------------------------------------------
!  variable in common cosysi
!    nscode:  total number of system dependent parameters
!             nscode = isicod + isrcod +3
!    isicod:  number of integer system dependent parameters
!    nterm:    number of different electronic coupling terms
!              this should be 1 here
!    nphoto:  this should be 1 here
!    nvib:    vibrational quantum number of initial state (0-6)
!  variable in common /cosys
!    scod:    character*8 array of dimension lcode, which contains names
!             of all system dependent parameters
!  line 16:
!    filnam:  name of file containing potential parameters
!
!  subroutines called: loapot(iunit,filnam)
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_conlam, only: nlam
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, iscod=>ispar
use mod_cosysr, only: isrcod, junkr, rcod=>rspar
implicit none
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
real(8) :: aso
integer :: ibran, j, jlpar, jtot1, jtot2, jtotd, l, lc, nphoto, nterm, nvib
logical existf
character*1 dot
character*(*) fname
character*60 filnam, line, potfil, filnm1
common /coipar/ jtot1,jtot2,jtotd,jlpar
#include "common/parbas.F90"
common /coskip/ nskip,iskip
integer :: nskip, iskip
#include "common/comdot.F90"
save potfil
!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
isicod = 4
scod(1) = 'NTERM'
scod(2) = 'NPHOTO'
scod(3) = 'NVIB'
scod(4) = 'IBRAN'
!  then all the system dependent real variables
!  in the same order as in the common block /cosysr/
isrcod = 1
scod(5) = 'A-SO'
nscode = isicod + isrcod
!  set default values for 2S + 2P atom scattering
nterm = 1
iscod(1)=1
iscod(2)=1
iscod(3)=0
iscod(4)=0
lammin(1)= 1
lammax(1) = 13
rcod(1)=881.
mproj(1)=0
niout=2
indout(1)=0
indout(2)=1
potfil = ' '
if(iread.eq.0) return
read (8, *, err=888) nphoto, nvib, ibran, aso
iscod(2)=nphoto
iscod(3)=nvib
iscod(4)=ibran
if (nvib .gt. 6) then
  write (6, 280) nvib
280   format ('*** NVIB = ',i2, '; SHOULD BE 0-6')
  go to 888
endif
if (iabs(ibran) .gt. 1) then
  write (6, 281) ibran
281   format ('*** IBRAN = ',i2, '; SHOULD BE -1, 0, or 1')
  go to 888
endif
if (jlpar .eq. -1) then
  if (ibran .eq. 0) then
    write (6, 282) jlpar, ibran
282     format ('*** JLPAR = ',i2, '; IBRAN = ', i2, &
            'Q-BRANCH NOT ALLOWED')
    goto 888
  endif
else
  if (ibran .ne. 0) then
    write (6, 283) jlpar, ibran
283     format ('*** JLPAR = ',i2, '; IBRAN = ', i2, &
            'P,R-BRANCH NOT ALLOWED')
    goto 888
  endif
endif
print *, ibran, jlpar
rcod(1)=aso
line=' '
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  return
endif
read (8, 285, end=286) line
285 format (a)
goto 286
! here if read error occurs
888 write(6,1000)
1000 format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
return
! --------------------------------------------------------------
entry ptr22p (fname,readpt)
line = fname
readpt = .true.
286 if (readpt) then
  l=1
  call parse(line,l,filnam,lc)
  if(lc.eq.0) then
    write(6,1020)
1020     format(' FILENAME MISSING FOR POTENTIAL INPUT')
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
    write(6,1025) filnam(1:lc)
1025     format(' FILE ',(a),' NOT FOUND')
    return
  end if
! now call loapot(iunit,filnam) routine to read potential parameters
  call loapot(1,filnam)
end if
close (8)
irpot=1
return
!
entry sav22p (readpt)
!  save input parameters for singlet-sigma + atom scattering
nphoto=iscod(2)
nvib=iscod(3)
ibran=iscod(4)
aso=rcod(1)
write (8, 300) nphoto, nvib, ibran, aso
300 format(i4,2(2x,i4),1pg12.5,3x,'  nphoto, nvib, ibran, a-so')
write (8, 285) potfil
return
end
end module mod_hiba10_22p
