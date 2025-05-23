#include "assert.h"
#include "unused.h"
! sy1sg1sg (sav1sg1sg/ptr1sg1sg) defines, saves variables and reads      *
!                  potentials for 1sigma - 1sigma (different molecules)  *
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of two unlike clsoed-shell diatomic molecules.  the second molecule
!  can be homonuclear.
!
!  this subroutine is an extension of ba2mol, which only treats collisions
!  of identical molecules
!
!  pot matrix elements are computed assuming the angular dependence defined
!  in the appendix of green, jcp 62, 2271 (1975)
!
!  author:  paul dagdigian
!  current revision date:  19-may-2022 by pjd
! --------------------------------------------------------------------
!     This module contains (explictly) the number of terms and their
!     indices in the expansion of the PES.  Its contents should be
!     filled in the pot routine.
!
!     This module replaces lammin, lammax, mproj in hibridon.
!
module mod_1sg1sg
  use mod_assert, only: fassert
implicit none
!
type lm_type
integer :: l1, l2, ltot
end type lm_type
!
type(lm_type), dimension(:), allocatable :: lms
end module mod_1sg1sg
! --------------------------------------------------------------------
module mod_hiba25_1sg1sg
  contains
  
subroutine ba1sg1sg (bqs, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  variables in call list:
!    bqs%jq:   on return contains combined rotational quantum numbers for each
!              channel.  in terms of the rotational quantum numbers of each
!              molecule we have:  j = 10 * j1 + j2
!    bqs%lq:   on return contains orbital angular momentum for each
!              channel
!    bqs%inq:  additional quantum index of each channel:  on return this is
!              set equal to j12 (the vector resultant of j1 and j2) times nsym
!    jhold:    on return contains rotational quantum numbers for each
!              rotational level, irrespective of projection degeneracy and the
!              degeneracy associated with different values of j12
!              note that jhold = 10 * j1hold + j2hold
!    ehold:    on return contains energy in hartrees of each rotational
!              level
!    ishold:   on return this is set equal to nsym
!    nlevel:   on return contains number of energetically distinct
!              rotational levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              rotational levels used in channel basis which are open
!              asymptotically
!    sc1:      scratch vector (not used here)
!    sc2:      scratch vector (not used here)
!    sc3:      scratch vector (not used here)
!    sc4:      scratch vector (not used here)
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!    jtot:     total angular momentum
!              in cc calculation jtot+1/2 is the total angular momentum
!              in cs calculation jtot is the l-bar quantum number
!    flaghf:   if .true., then system has half-integer spin
!              if .false., then system has integer spin (this is the case
!              here)
!    flagsu:   if .true., then molecule-surface collisons (this variable
!              should always be false in the present case)
!    csflag:   if .true., then coupled-states calculation
!              if .false., the close-coupled calculation
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true. , then homonuclear molecule
!              if .false., then heteronuclear molecule
!              the homonuclear option is not specifically implemented here
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              parity=(-1)**jlpar (by definition parity=(-1)**(j1+j2+l)
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to subroutines
!              propag and soutpt.  ntop is set in basis only if nu = numin
!              otherwise it is unchanged from the value supplied by the
!              calling program
!  variables in common block /cosysr/
!    isrcod:   number of real system dependent variables
!    b1rot:    rotational constant in cm-1 of molecuole 1
!    d1rot:    centrifugal distortion constant in cm-1 of molecule 1
!    b2rot:    rotational constant in cm-1 of molecuole 2
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   numbr of integer system dependent variables
!    j1max:    maximum rotational angular momentum of molecule 1
!    j2min:    minimum rotational angular momentum of molecule 2
!    j2max:    maximum rotational angular momentum of molecule 2
!    ipotsy:   symmetry of potential.  set to 2 for homonuclear
!              molecule 2, set to 1 for heteronuclear molecule 2
!  subroutines called:
!   vlmlml:    returns molecule-molecule angular coupling coefficient for
!              particular choice of channel index
! --------------------------------------------------------------------
use mod_1sg1sg
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam
use mod_cosysi, only: ispar
use mod_cosysr, only: rspar
use constants, only: econv, xmconv
use mod_par, only: iprint
use mod_par, only: boundc
use mod_ered, only: ered, rmu
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
type(bqs_type), intent(out) :: bqs
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst
dimension jhold(1), ehold(1), &
    sc1(1), sc2(1), sc3(1), sc4(1), ishold(1)
data zero, ione, itwo, ithree / 0.d0, 1, 2, 3 /
integer, pointer :: j1max, j2min, j2max, ipotsy2
real(8), pointer :: b1rot, d1rot, b2rot
j1max=>ispar(1); j2min=>ispar(2); j2max=>ispar(3); ipotsy2=>ispar(4)
b1rot=>rspar(1); d1rot=>rspar(2); b2rot=>rspar(3)

UNUSED_DUMMY(sc1)
UNUSED_DUMMY(sc2)
UNUSED_DUMMY(sc3)
UNUSED_DUMMY(sc4)
UNUSED_DUMMY(ihomo)
UNUSED_DUMMY(nu)

!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (9, 5)
5   format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***' )
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (flagsu) then
  write (6, 10)
  write (9, 10)
10   format ('  *** FLAGSU = .TRUE. FOR', &
         ' MOLECULE-MOLECULE COLLISION; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (csflag) then
  write (6, 25)
  write (9, 25)
25   format (' *** CSFLAG=.TRUE.; NOT ALLOWED IN 1SIG-1SIG BASIS;', &
          ' ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
!  check that ipotsy2 is either 1 or 2
if (ipotsy2.lt.1 .or. ipotsy2.gt.2) then
  write (6, 27) ipotsy2
  write (9, 27) ipotsy2
27   format (' *** IPOTSY2 .NE. 1 .OR.2;', &
          ' ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (bastst) then
  write (6,  90) rmu * xmconv, ered * econv, jtot, jlpar
  write (9,  90) rmu * xmconv, ered * econv, jtot, jlpar
90   format(/2x, ' **  CC  1SIGMA-1SIGMA ** RMU=',f10.6, &
     '  E=', f9. 3, 2x, 'JTOT=', i4, &
     '  JLPAR=', i2)
end if
write (9, 95) rcut
95   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
!  assign quantum numbers and energies for rotational levels
nlevel = 0
do j1 = 0, j1max
  do j2 = j2min, j2max, ipotsy2
    nlevel = nlevel + 1
    jsub = 10 * j1 + j2
    jhold(nlevel) = jsub
    ishold(nlevel) = 0
    jj1 = j1 * (j1 + 1.d0)
!
!  subtract energy of lowest H2 level (13-jan-2022)
    ezero = b2rot * j2min * (j2min + 1d0)
    ehold(nlevel) = (b1rot * jj1 - d1rot * jj1**2 &
      + b2rot * j2 * (j2 + 1.d0) - ezero) / econv
  end do
end do
!
!  form list of all energetically open rotational levels included in the
!  calculations and their energies
nlevop = 0
do 80 i = 1, nlevel - 1
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
        goto 80
      end if
75     continue
  end if
80 continue
if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
!
!  set up CC scattering channel basis
call bqs%init(nmax)
n = 0
do i = 1, nlevel
  j1 = jhold(i)/10
  j2 = mod(jhold(i),10)
  j12min = abs(j1 - j2)
  j12max = j1 + j2
  do j12i = j12min, j12max
    lmin = iabs(jtot - j12i)
    lmax = jtot + j12i
    do li = lmin, lmax
      ix = (- 1) ** (j1 + j2 - jtot + li)
      if (ix .eq. jlpar) then
        n = n + 1
        if (n .le. nmax) then
          bqs%jq(n) = 10 * j1 + j2
          bqs%inq(n) = 0
          bqs%lq(n) = li
          bqs%j12(n) = j12i
          bqs%length = n
          eint(n) = ehold(i)
          cent(n) = li * (li + 1.d0)
        else
          write (9, 230) n, nmax
          write (6, 230) n, nmax
230           format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF', &
              i4,'; ABORT')
          if (bastst) then
            return
          else
            call exit
          end if
        end if
      end if
    end do
  end do
end do
!
!  order channel basis functions in order of increasing energy
!  sort algorithm from Numerical Recipes
if (n .gt. 1) then
  do 212 jj = 2, n
    ekeep = eint(jj)
    jkeep = bqs%jq(jj)
    ikeep = bqs%inq(jj)
    lkeep = bqs%lq(jj)
    j12kp = bqs%j12(jj)
    ckeep = cent(jj)
    do 211 ii = jj-1, 1, -1
      if (eint(ii).le.ekeep) goto 210
      eint(ii+1) = eint(ii)
      bqs%jq(ii+1) = bqs%jq(ii)
      bqs%inq(ii+1) = bqs%inq(ii)
      bqs%lq(ii+1) = bqs%lq(ii)
      bqs%j12(ii+1) = bqs%j12(ii)
      cent(ii+1) = cent(ii)
211     continue
    ii = 0
210     eint(ii+1) = ekeep
    bqs%jq(ii+1) = jkeep
    bqs%inq(ii+1) = ikeep
    bqs%lq(ii+1) = lkeep
    bqs%j12(ii+1) = j12kp
    cent(ii+1) = ckeep
212   continue
end if
!
!  now check to see if any of the open channels are closed at r = rcut
nn = 0
if (rcut .gt. 0.d0 .and. .not.boundc) then
  emin = 1.e+7
  do 120 i = 1, n
    if (eint(i) .le. ered) then
!  here if channel is open asymptotically
      if (jtot * (jtot + 1) / (2.d0 * rmu * rcut ** 2) &
          .gt. (ered - eint(i)) ) then
!  here is channel is open asymptotically but closed at r  = rcut
        if (eint(i) .lt. emin) emin = eint(i)
!  emin now contains the lowest rotational energy for which this condition is met
      end if
    end if
120   continue
!  now eliminate all channels with eint .gt. emin if any of the channels
!  are open asymptoptically but closed at r  = rcut
  if (emin .lt. ered) then
    do 130 i = 1, n
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        eint(nn) = eint(i)
        bqs%inq(nn) = bqs%inq(i)
        bqs%lq(nn) = bqs%lq(i)
        bqs%j12(nn) = bqs%j12(i)
        cent(nn) = cent(i)
      end if
130     continue
!  reset number of channels
    n = nn
    bqs%length = n
  end if
end if
!  return if no channels
if (n .eq. 0) return
if (nn .eq. numin) then
  ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure that this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) then
    ntop = ntop + 1
  else
    if (n .gt. ntop) then
      write(6,160) nn, n, ntop
      write(9,160) nn, n, ntop
160       format (' *** NCH = ',i3,' AT NU = ',i2,' .GT. NTOP = ',i3, &
          '; ABORT ***',/, &
      '      CHECK RCUT')
      call exit
    end if
  end if
end if
!
!  now list channels if requested
if (clist) then
  if (bastst) write (6, 200)
  write(9, 200)
200   format(/'   N     J1     J2     J12     L   EINT(CM-1)')
  do i = 1, n
    ecm = eint(i) * econv
    j1 = bqs%jq(i)/10
    j2 = mod(bqs%jq(i),10)
    write (6, 220) i, j1, j2, bqs%j12(i), bqs%lq(i), ecm
    write (9, 220) i, j1, j2, bqs%j12(i), bqs%lq(i), ecm
220     format (i4, 4i7, f10.3)
  end do
end if
!
!     calculate coupling matrix elements
if (bastst .and. iprint.eq.2) then
  write (6, 285)
  write (9, 285)
285   format (/' ILAM  L1   L2  LTOT   ICOL IROW      I', &
    '      IV2         VEE')
end if
i = 0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 400 ilam = 1, nlam
!     ilam denotes a particular L1,L2,L term
  ancouma => v2%get_angular_coupling_matrix(ilam)
  ll1 = lms(ilam)%l1
  ll2 = lms(ilam)%l2
  lltot = lms(ilam)%ltot
  inum = 0
  do 355 icol = 1, n
    j1c = bqs%jq(icol)/10
    j2c = mod(bqs%jq(icol),10)
    j12c = bqs%j12(icol)
    lc = bqs%lq(icol)
    do 350 irow = icol, n
      j1r = bqs%jq(irow)/10
      j2r = mod(bqs%jq(irow),10)
      j12r = bqs%j12(irow)
      lr = bqs%lq(irow)
!     initialize potential to zero
      vee = zero
      call v1sgsg(j1r,j2r,j12r,lr,j1c,j2c,j12c, &
          lc,jtot,ll1,ll2,lltot,vee)
      if (vee .ne. zero) then
        i = i + 1
        inum = inum + 1
        call ancouma%set_element(irow=irow, icol=icol, vee=vee)
        if (bastst .and. iprint.ge.2) then
          write (6, 290) ilam, ll1, ll2, lltot, &
              icol, irow, i, vee
          write (9, 290) ilam, ll1, ll2, lltot, &
              icol, irow, i, vee
290             format (i4, 3i5, 2x, 2i5, i8, e20.7)
        end if
      end if
350     continue
355   continue
  if (bastst) then
    write (6, 370) ilam, ancouma%get_num_nonzero_elements()
    write (9, 370) ilam, ancouma%get_num_nonzero_elements()
370     format ('ILAM=',i4,' LAMNUM(ILAM) = ',i7)
  end if
400 continue
if (clist .and. bastst) then
  write (6, 420) v2%get_num_nonzero_elements()
  write (9, 420) v2%get_num_nonzero_elements()
420   format (' *** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', &
      i6)
end if
return
end
! --------------------------------------------------------------------
subroutine v1sgsg(j1r,j2r,j12r,lr,j1c,j2c,j12c, &
          lc,jtot,l1,l2,ltot,vee)
!  subroutine to calculate v-lambda matrices for close-coupled
!  treatment of collisions of a molecule in a 1sigma
!  electronic state with a (different) molecule in a 1sigma state
!
!  matrix elements originally derived by s. green, jcp 62, 2271 (1973);
!  see also lanza et al., jcp 140, 064316 (2014) for explicit formula
!
!  subroutines called:
!     xf3j, xf6j, xf9j
!
!  written  by p. dagdigian
!  current revision date:  2-jun-2017 (p.dagdigian)
! --------------------------------------------------------------------
use mod_1sg1sg
use mod_hiutil, only: xf3j, xf6j, xf9j
implicit double precision (a-h,o-z)
data zero, one, two/ 0.d0, 1.d0, 2.d0/, &
  sq4pi3 / 44.546623974653656d0 /
vee = zero
iph = jtot + j1c + j2c + j12r
phase = 1.d0
if (iph .ne. 2*(iph/2)) phase = -1.d0
xj1r = j1r
xj2r = j2r
xj12r = j12r
xlr = lr
xj1c = j1c
xj2c = j2c
xj12c = j12c
xlc = lc
xjtot = jtot
xl1 = l1
xl2 = l2
xltot = ltot
facj = (two*xltot+one) * sqrt((two*xl1+one) * (two*xl2+one) &
  * (two*xj1r+one) * (two*xj2r + one) * (two*xj12r + one) &
  * (two*xlr+one) * (two*xj1c+one) * (two*xj2c+one) &
  * (two*xj12c+one) * (two*xlc+one) )
cg1 = xf3j (xltot, xlr, xlc, zero, zero, zero)
if (cg1 .eq. zero) return
cg2 = xf3j (xl1, xj1r, xj1c, zero, zero, zero)
if (cg2 .eq. zero) return
cg3 = xf3j (xl2, xj2r, xj2c, zero, zero, zero)
if (cg3 .eq. zero) return
f6j = xf6j (xlr, xlc, xltot, xj12c, xj12r, xjtot)
if (f6j .eq. zero) return
f9j = xf9j (xj12r, xj2r, xj1r, xj12c, xj2c, xj1c, &
  xltot, xl2, xl1)
if (f9j .eq. zero) return
vee = phase * facj * cg1 * cg2 * cg3 * f6j * f9j / sq4pi3
return
end
! -----------------------------------------------------------------------
subroutine sy1sg1sg (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for collisions of
!  two different 1sigma linear molecules
!  current revision date: 23-may-2017 by p.dagdigian
!  -----------------------------------------------------------------------
!  variable in common cosysi
!    nscode:   total number of system dependent parameters
!              nscode = isicod + isrcod + 3
!    isicod:   number of integer system dependent parameters
!    nterm:    number of different electronic coupling terms
!              this should be 1 here
!    j1max:    maximum rotational angular momentum for molecule 1
!    j2min:    minimum rotational angular momentum for molecule 2
!    j2max:    maximum rotational angular momentum for molecule 2
!    ipotsy2:  symmetry of potential.  set to 2 for homonuclear
!              molecule 2, set to 1 for heteronuclear molecule 2
!    iop:      ortho/para label for levele of molecule 2. If ihomo=.true.
!              then only para states will be included if iop=1 and
!              only ortho states if iop=-1
!  variable in common  /cosys/
!    scod:    character*8 array of dimension lcode, which contains names
!             of all system dependent parameters
!  line 16:
!    filnam:  name of file containing potential parameters
!
!  subroutines called: loapot(iunit,filnam)
!  -----------------------------------------------------------------------
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, rspar
use funit, only: FUNIT_INP
use mod_hiutil, only: gennam, get_token
use mod_hipot, only: loapot
implicit none
integer, intent(inout) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: j, l, lc
logical existf
character*1 dot
character*(*) fname
character*60 filnam, line, potfil
character*68 filnm1
#include "common/comdot.F90"
save potfil

integer, pointer, save :: j1max, j2min, j2max, ipotsy2, iop
real(8), pointer, save :: b1rot, d1rot, b2rot
j1max=>ispar(1); j2min=>ispar(2); j2max=>ispar(3); ipotsy2=>ispar(4); iop=>ispar(5)
b1rot=>rspar(1); d1rot=>rspar(2); b2rot=>rspar(3)

!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
isicod = 4
scod(1) = 'J1MAX'
scod(2) = 'J2MIN'
scod(3) = 'J2MAX'
scod(4) = 'POTSY2'
!  then all the system dependent real variables
!  in the same order as in the common block /cosysr/
isrcod = 3
scod(5) = 'B1ROT'
scod(6) = 'D1ROT'
scod(7) = 'B2ROT'
nscode = isicod + isrcod + 3
if(iread.eq.0) return
!  read the last few lines of the input file
read (8, *, err=888) j1max, j2min, j2max, ipotsy2
read (8, *, err=888) b1rot, d1rot, b2rot
read (8, *, err=888) potfil
call loapot(10, potfil)
close(8)
return
! here if read error occurs
888 write(6,1000)
1000 format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
return
! --------------------------------------------------------------
entry ptr1sg1sg (fname,readpt)
line = fname
readpt = .true.
if (readpt) then
  l=1
  call get_token(line,l,filnam,lc)
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
! --------------------------------------------------------------
entry sav1sg1sg ()
!  save input parameters for two unlike 1sigma molecule scattering
write (FUNIT_INP, 310) j1max, j2min, j2max, ipotsy2
310 format(5i4,15x,'j1max, j2min, j2min, ipotsy2')
write (FUNIT_INP, 320) b1rot, d1rot, b2rot
320 format(f10.7, e12.5, f10.7, '   b1rot, d1rot, b2rot')
write (FUNIT_INP, 285) potfil
285 format (a)
return
end

! --------------------------------eof---------------------------------

end module mod_hiba25_1sg1sg
