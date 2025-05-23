#include "assert.h"
#include "unused.h"
! sy1sg (sav1sg/ptr1sg) defines, saves variables and reads              *
!                  potential for singlet sigma scattering               *
!************************************************************************
!                                                                       *
!                            basis  library                             *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!  1. basis      dispatcher routine to select basis                     *
!                see in header for basis subroutine for list of basis   *
!                routines currently available                           *
!  2. ba1sg      basis for singlet sigma scattering                     *
!  3. vlm1sg     calculates v-lamda matrices for above                  *
!  4. is_twomol check if a basis is for mol-mol collision (j=10j1+j2)   *
!  5. is_j12    check if j12 is used in a basis                         *
!                                                                       *
!     current revision:  24-jul-2019 (p.dagdigian)                       *
!                                                                       *
!************************************************************************
module mod_hiba01_1sg
  use mod_assert, only: fassert
contains
! --------------------------------------------------------------------
subroutine ba1sg (bqs, jhold, ehold, ishold, nlevel, nlevop, &
                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential
!  for collision of singlet-sigma molecule with a structureless atom
!  authors:  millard alexander and hans-joachim werner
!  current revision date:  18-jun-2006 by mha
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains rotational quantum numbers for each
!              channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains symmetry index for each channel
!    jhold:    on return contains rotational quantum numbers for each
!              rotational level
!    ehold:    on return contains energy in hartrees of each rotational
!              level
!    ishold:   on return contains symmetry index of each rotational level
!    nlevel:   on return contains number of energetically distinct
!              rotational levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              rotational levels used in channel basis which are open
!              asymptotically
!    sc1,sc2:  scratch vectors of length at least nmax
!    sc3,sc4:  scratch vectors of length at least nmax
!              these scratch vectors are not used here
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!              note that this is ignored in molecule-surface collisions!!!
!    jtot:     total angular momentum
!    flaghf:   if .true., then system with half-integer spin
!              if .false., then system with integer spin (this is the case
!              here)
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
!                  (-1)**(j+l-jtot)=jlpar
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
!    brot:     rotational constant in cm-1
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   numbr of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    jmin:     the minimum rotational angular momenta
!    jmax:     the maximum rotational angular momenta
!  subroutines called:
!   vlm1sg:    returns singlet-sigma angular coupling coefficient for
!              particular choice of initial and final channel quantum numbers
!
! --------------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam
use mod_cosysi, only: ispar, convert_ispar_to_mat
use mod_cosysr, only: convert_rspar_to_mat
use constants, only: econv, xmconv
use mod_par, only: iprint
use mod_parbas, only: maxvib, ntv, ivcol, ivrow, lammin, lammax
use mod_par, only: boundc
use mod_ered, only: ered, rmu
use mod_skip, only: nskip
use mod_vib, only: nvib=>nvibs, ivib=>ivibs
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
type(bqs_type), intent(out) :: bqs
integer, intent(out), dimension(:) :: jhold
real(8), intent(out), dimension(:) :: ehold
integer, intent(out), dimension(:) :: ishold
integer, intent(out) :: nlevel
integer, intent(out) :: nlevop
real(8), intent(out), dimension(:) :: sc1
real(8), intent(out), dimension(:) :: sc2
real(8), intent(out), dimension(:) :: sc3
real(8), intent(out), dimension(:) :: sc4
real(8), intent(in) :: rcut
integer, intent(in) :: jtot
logical, intent(in) :: flaghf
logical, intent(in) :: flagsu
logical, intent(in) :: csflag
logical, intent(in) :: clist
logical, intent(in) :: bastst
logical, intent(in) :: ihomo
integer, intent(in) :: nu
integer, intent(in) :: numin
integer, intent(in) :: jlpar
integer, intent(out) :: n
integer, intent(in) :: nmax
integer, intent(out) :: ntop
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
!   econv is conversion factor from cm-1 to hartrees
!   xmconv is converson factor from amu to atomic units
real(8), dimension(4, maxvib) :: rpar
integer, dimension(2, maxvib) :: iscod
integer, pointer :: nterm, nvmin, nvmax
nterm=>ispar(1); nvmin=>ispar(2); nvmax=>ispar(3)

UNUSED_DUMMY(sc1)
sc1(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(sc2)
sc2(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(sc3)
sc3(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(sc4)
sc4(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.

call convert_rspar_to_mat(4,maxvib, rpar)
call convert_ispar_to_mat(2, maxvib, 4, iscod)
zero = 0.d0
two = 2.d0
!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (6, 1)
  write (9, 1)
1   format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***' )
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (flagsu .and. .not. csflag) then
  write (6, 8)
  write (9, 8)
8  format &
   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
! check that jmin .ge. nu for bound state calculations (this shouldn't be a bu
! but is?)
nsum = 0
iva=0
ive=0
!  check that requested vib levels have been defined
do i=1,nvib
  if(ivib(i).ne.ivib(1)+i-1) then
    write(6,5) (ivib(k),k=1,nvib)
5     format(/' INPUT ERROR: NON-SEQUENTIAL VIBRATIONAL', &
            ' STATES: ',10i5)
    call exit
  end if
  if(ivib(i).eq.nvmin) iva=i
  if(ivib(i).eq.nvmax) ive=i
end do
if(iva.eq.0.or.ive.eq.0) then
  write(6,15) nvmin,nvmax,(ivib(i),i=1,nvib)
15     format(/' PARAMETERS UNDEFINED FOR VIBRATIONAL STATE'/ &
  1x,' REQUESTED STATES:',i3,'-',i2,' DEFINED STATES:',10i5)
  call exit
end if

if (clist) then
  if (flagsu) then
    write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
    write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16     format(/,' **  1-SIGMA ON UNCORRUGATED SURFACE ** RMU=', f9.4, &
      '             E=', f7.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)
  else
    if (.not. csflag) then
      write (6,20) rmu * xmconv, ered * econv, jtot, jlpar
      write (9,20) rmu * xmconv, ered * econv, jtot, jlpar
20       format(/,' **  CC SINGLET SIGMA DIATOMIC ** RMU=', f9.4, &
           '       E=',f7.2,'   JTOT=', i5, 2x,' JLPAR=',i2)
    else
      write (6,25) rmu * xmconv, ered * econv, jtot, nu
      write (9,25) rmu * xmconv, ered * econv, jtot, nu
25       format(/,' **  CS SINGLET SIGMA DIATOMIC ** RMU=', f9.4, &
           '       E=',f7.2,'   JTOT=', i5, 2x,' NU=',i2)
    end if
  end if
  if (.not. flagsu) write (9,30) rcut
30   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
!  assign quantum numbers and energies for rotational levels
  write(6,31) ' State    B(v)',' D(v)','H(v)','E(v)'
  write(9,31) ' State    B(v)',' D(v)','H(v)','E(v)'
31   format(/2(a,8x),a,10x,a)
  do iv=iva,ive
    write(6,35) ivib(iv),(rpar(jj,iv),jj=1,4)
    write(9,35) ivib(iv),(rpar(jj,iv),jj=1,4)
  end do
35   format(1x,i3,2x,3g12.5,f15.8)
end if
call bqs%init(nmax)
n=0
nskip = 1
if (ihomo) nskip = 2
do 120 iv=iva,ive
  jmin=iscod(1,iv)
  if (boundc.and.csflag) then
      if (jmin.lt.nu) then
         write(6, 7) jmin, nu
         write(9, 7) jmin, nu
7        format(/ &
      ' JMIN = ',i2,', .LT. NU = ',i2,' IN BASIS; JMIN RESET')
         jmin=nu
      endif
  endif
  jmax=iscod(2,iv)
  brot=rpar(1,iv)/econv
  drot=rpar(2,iv)/econv
  hrot=rpar(3,iv)/econv
  evib=rpar(4,iv)/econv
  do 120  ji = jmin, jmax, nskip
    jj1=ji*(ji+1)
    ee=brot*jj1 - drot*jj1**2 + hrot*jj1**3 + evib
    if (.not. csflag) then
    !  here for cc calculations
      lmax = jtot + ji
      lmin = iabs (jtot - ji)
      do li = lmin, lmax
        ix = (-1) ** (ji + li - jtot)
        if (ix .eq. jlpar) then
    !  here for correct combination of j and l
          n = n + 1
          if (n .gt. nmax) go to 130
          bqs%lq(n) = li
          cent(n) = li * (li + 1.)
          bqs%inq(n) = ivib(iv)
          bqs%jq(n) = ji
          bqs%length = n
          eint(n) = ee
        end if
      end do
    else
  !  here for cs calculations
      if (ji .lt. nu) go to 120
      n = n + 1
      if (n .gt. nmax) go to 130
      bqs%lq(n) = jtot
      if (.not.boundc) then
        cent(n) = jtot * (jtot + 1)
      else
        xjtot=jtot
        xj=bqs%jq(n)
        xnu=nu
        cent(n)=xjtot*(xjtot+1)+xj*(xj+1)-2*xnu*xnu
      endif
      bqs%inq(n) = ivib(iv)
      bqs%jq(n) = ji
      bqs%length = n
      eint(n) = ee
    end if
120 continue
130 if (n .gt. nmax) then
  write (9, 140) n, nmax
  write (6, 140) n, nmax
140   format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF', &
         i4,'; ABORT')
  if (bastst) then
    return
  else
    call exit
  end if
end if
!  now check to see if any of the open channels are closed at r=rcut
!  this is not done for molecule-surface collisions or for rcut < 0
!  and for bound state calculations
if (.not.flagsu .and. rcut .gt. 0.d0 .and. .not.boundc) then
  emin = 1.e+7
  do i = 1, n
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
  end do
!  now eliminate all channels with eint .ge. emin if any of the channels
!  are open asymptotically but closed at r = rcut
  if (emin.lt.ered) then
    nn=n
    n = 0
    do i = 1, nn
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        n = n + 1
        eint(n) = eint(i)
        bqs%inq(n) = bqs%inq(i)
        bqs%jq(n) = bqs%jq(i)
        cent(n) = cent(i)
        bqs%lq(n) = bqs%lq(i)
        bqs%length = n
      end if
    end do
!  reset number of channels
  end if
end if
!  return if no channels
if (n .eq. 0) return
if (nu .eq. numin) then
  ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1
else
  if (n.gt.ntop) then
    write (6, 160) nu, n, ntop
    write (9, 160) nu, n, ntop
160     format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3, &
            '; ABORT **',/, &
    '     CHECK RCUT')
    call exit
  endif
end if
nlevel = 0
nlevop = 0
!  form list of all energetically open rotational levels included in the
!  calculations and their energies
!  if homonuclear diatomic, skip space is two
do iv=iva,ive
  jmin=iscod(1,iv)
  jmax=iscod(2,iv)
  brot=rpar(1,iv)/econv
  drot=rpar(2,iv)/econv
  hrot=rpar(3,iv)/econv
  evib=rpar(4,iv)/econv
  do ji = jmin, jmax, nskip
    jj1=ji*(ji+1)
    ee=brot*jj1 - drot*jj1**2 + hrot*jj1**3 + evib
    nlevel = nlevel + 1
    ehold(nlevel) = ee
    jhold(nlevel) = ji
    ishold(nlevel) = ivib(iv)
  end do
end do
!  now sort this list to put closed levels at end
!  also determine number of levels which are open
nlevop = 0
if (nlevel .gt. 1) then
  do 80  i = 1, nlevel - 1
    if (ehold(i) .le. ered) then
       nlevop = nlevop + 1
    else
      do ii = i + 1, nlevel
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
          go to 80
        end if
      end do
    end if
80   continue
else
! here for only one level
  if (ehold(1) .le. ered) then
    nlevop=1
  else
    nlevop=0
  endif
endif
if (nlevop .le. 0) then
  write (9, 85)
  write (6, 85)
85   format('*** NO OPEN LEVELS IN BA1SG; ABORT')
  if (bastst) return
  call exit
endif
if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
if (nu .eq. numin) then
  ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1
end if
!  now list channels if requested
if (clist) then
  write (6, 255)
  write (9, 255)
255   format(/'   N   V   J    L      EINT(CM-1)',/)
  do 265  i = 1, n
    write (6, 260) i, bqs%inq(i), bqs%jq(i), bqs%lq(i), eint(i) * econv
    write (9, 260) i, bqs%inq(i), bqs%jq(i), bqs%lq(i), eint(i) * econv
260     format (3i4, i5, f13.3)
265   continue
  write (6, 256)
256   format(/' OPEN LEVELS:'// &
          '   N   V   J      EINT(CM-1)',/)
  do 266  i = 1, nlevop
    write (6, 261) i, ishold(i), jhold(i),  ehold(i) * econv
261     format (3i4, f13.3)
266   continue
end if
!  now calculate coupling matrix elements
if (bastst .and. iprint.gt.1) then
  write (6, 280)
  write (9, 280)
280   format (/'  VR VC  LAMBDA ILAM    I   ICOL  IROW', &
           '  IV2       VEE')
end if
!
!  the following only for safety. Checks that pot routine
!  contains all required vib. levels in correct order
!
do 295 irow=1,n
ivr=bqs%inq(irow)
do 290 icol=1,irow
ivc=bqs%inq(icol)
do 290 ivb=1,ntv(1)
if(ivrow(ivb,1).eq.ivr.and. &
   ivcol(ivb,1).eq.ivc) goto 295
290 continue
write(6,291) ivr,ivc
write(9,291) ivr,ivc
291 format(/' VIBRATIONAL STATE NOT DEFINED IN POT',2i5)
call exit
295 continue
! i counts v2 elements
! inum counts v2 elements for given lambda
! ilam counts numver of v2 matrices
! ij is address of given v2 element in present v2 matrix
i = 0
ilam=0
ASSERT(.not. allocated(v2))
num_lam = 0
do il = lammin(1), lammax(1), nskip
  num_lam = num_lam + ntv(1)
end do
v2 = ancou_type(nlam=num_lam, num_channels=ntop)
do 320 iv=1,ntv(1)
    ivr=ivrow(iv,1)
    ivc=ivcol(iv,1)
do 320 il = lammin(1), lammax(1), nskip
  lb=il
  ilam=ilam+1
  inum = 0
  ancouma => v2%get_angular_coupling_matrix(ilam)
  do 310  icol= 1, n
    do 300  irow = icol, n
      if(bqs%inq(irow).ne.ivr.or.bqs%inq(icol).ne.ivc) goto 300
!  here for coupling between molecular rotational levels
      call vlm1sg (bqs%jq(irow), bqs%lq(irow), bqs%jq(icol), bqs%lq(icol), jtot, &
                   nu, lb, vee, csflag)
      if (vee .eq. 0) goto 300
        i = i + 1
        inum = inum + 1
         if (bastst .and. iprint.gt.1) then
           write (6, 340) ivr,ivc,il,ilam,i,irow,icol, vee
           write (9, 340) ivr,ivc,il,ilam,i,irow,icol, vee
340            format (1x,2i3,5i6, g17.8)
         end if
         call ancouma%set_element(irow=irow, icol=icol, vee=vee)
300     continue
310   continue
  if (bastst .and. iprint.ge.1) then
    write (6, 420) ivr,ivc,il,inum
    write (9, 420) ivr,ivc,il,inum
420     format(' IVR=',i2,'  IVC=',i2,'  LAMBDA=',i2, &
           ' NUMBER OF NONZERO MATRIX ELEMENTS',i5)
  end if
  if(inum.ne.0) nlam=ilam
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
subroutine vlm1sg (j1, l1, j2, l2, jtot, nu, lb, v, csflag)
! --------------------------------------------------------------------
!  subroutine to evaluate the angular coupling matrix element for rotationally
!  inelastic collisions of a 1sigma molecule with a structureless target
!  current revision date: 5-sep-88
!  variables in call list:
!  j1,l1:    initial rotational and orbital angular momenta
!  j2,l2:    final rotational and orbital angular momenta
!  jtot:     total angular momentum
!  nu:       coupled states projection index
!  lb:       value of legendre expansion index
!  v:        on return:  contains desired coupling matrix element
!  csflag:   if .true., then coupled states calculation
!  subroutines called:
!  xf3j:     evaluates 3j symbol
!  xf6j:     evaluates 6j symbol
! --------------------------------------------------------------------
use mod_hiutil, only: xf3j, xf6j
implicit double precision (a-h,o-z)
!      real v, xj1, xj2, xjtot, xl1, xl2, xlb, xnu, zero
integer ij, il, j1, j2, jtot, l1, l2, lb, nu
logical csflag
zero = 0.
v = 0.
xj1 = j1
xl1 = l1
xj2 = j2
xl2 = l2
xjtot = jtot
xlb = lb
xnu = nu
if (.not. csflag) then
!  here for close-coupling, for the matrix element see w.lester, meth. comp.
!  phys. 10, 211 (1971)
  ij = (-1) ** (j1 + j2 + lb)
  il = (-1) ** (l1 + l2 + lb)
  if (ij .ne. -1  .and. il .ne. -1) then
    v = xf3j (xj1, xlb, xj2, zero, zero, zero) * &
        xf3j (xl1, xlb, xl2, zero, zero, zero) * &
        xf6j (xj1, xj2, xlb, xl2, xl1, xjtot) * &
        sqrt ( (2 * j1 + 1.) * (2 * j2 + 1.) * (2 * l1 + 1.) &
             * (2 * l2 + 1.) ) * (-1) ** ( j1 + j2 - jtot )
  end if
else
!  here for coupled-states, see p. mcguire and d.j. kouri, j. chem. phys.
!  60, 2488 (1974)
  v = xf3j (xj1, xlb, xj2, zero, zero, zero) * &
      xf3j (xj1, xlb, xj2, -xnu, zero, xnu) * &
      sqrt ( (2 * j1 + 1.) * (2 * j2 + 1.) ) * (-1) ** nu
end if
return
end
!  -----------------------------------------------------------------------
subroutine sy1sg (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for singlet-sigma
!   + atom scattering using werner-follmeg potential form
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  current revision date: 10-mar-1994 by mha
!  -----------------------------------------------------------------------
!  variable in common /cosys
!    scod:    character*8 array of dimension lcode, which contains names
!             of all system dependent parameters
!  line 13:
!    nlam:     total number of anistropic terms in potential
!    jmin:     minimum molecular rotational angular momenta included in
!              the channel basis
!    jmax      maximum molecular rotational angular momenta inncluded in
!              the channel basis
!  line 14:
!    brot:    rotational constant of the molecule in cm-1
!  line 16:
!    filnam:  name of file containing potential parameters
!
!  subroutines called: loapot(iunit,filnam)
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, iscod=>ispar
use mod_cosysr, only: isrcod, rcod => rspar
use mod_par, only: ihomo
use funit, only: FUNIT_INP
use mod_parbas, only: maxvib, lammin, lammax, mproj
use mod_skip, only: nskip
use mod_vib, only: nvib => nvibs, ivib => ivibs
use mod_hiutil, only: gennam, get_token
use mod_hipot, only: loapot
implicit none
integer, intent(inout) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: i, iofi, iofr
integer :: j, l, lc
logical existf
character*1 dot
character*4 char
character*(*) fname
character*60 filnam, line, potfil
character*68 filnm1

save potfil
!equivalence(iscod(1),nterm),(iscod(2),nvibmn),(iscod(3),nvibmx)
#include "common/comdot.F90"

!     number and names of system dependent parameters
!  set default values for singlet-sigma scattering

iscod(1) = 1
if (iread .eq. 0) then
  iscod(2) = 0
  iscod(3) = 0
  lammin(1)=0
  lammax(1)=-1
  mproj(1)=0
  niout=1
  indout(1)=0
  nvib = 1
endif
if (.not.readpt)irpot=1
if (iread .eq. 1) irpot=1
if (ihomo) nskip = 2
potfil = ' '
!  read number of vib states
if(iread.ne.0) read (8, *, err=88) nvib, iscod(2),iscod(3) ! nvib, vmin, vmax
if(nvib.gt.iscod(3)-iscod(2)+1) then
  write (6,40) nvib, iscod(2), iscod(3)
40   format(' ** PROBABLE VIBRATIONAL LEVEL NUMBERING ERROR:',/, &
         '   VMIN =',i2,', VMAX=',i3,', BUT NVIB =',i3)
  return
endif
if(nvib.gt.maxvib) stop 'nvib'
if(nvib.le.0) nvib=1
!  read data for each vib state
!  brot, drot, hrot are bv, dv, hv (see huber herzberg, page x)
scod(1)='NTERM'
scod(2)='VMIN'
scod(3)='VMAX'
isrcod=0
isicod=3
iofr=2*nvib+4-1
do i = 1,nvib
  if(iread.ne.0) then
    read (8, *, err=99) ivib(i),(iscod(isicod+j),j=1,2) ! iv, jmin, jmax
    read (8, *, err=99) (rcod(isrcod+j),j=1,3) ! brot, drot, hrot
    read (8, *, err=99) rcod(isrcod+4) ! evib
  end if
  char=' '
  if(nvib.gt.1.or.ivib(i).ne.0) then
    if(ivib(i).le.9) write(char,'(''('',i1,'')'')') ivib(i)
    if(ivib(i).gt.9) write(char,'(''('',i2,'')'')') ivib(i)
  end if
  scod(isicod+1)='JMIN'//char
  scod(isicod+2)='JMAX'//char
  scod(iofr+1)='BROT'//char
  scod(iofr+2)='DROT'//char
  scod(iofr+3)='HROT'//char
  scod(iofr+4)='EVIB'//char
  iofr=iofr+4
  isicod=isicod+2
  isrcod=isrcod+4
end do
if(isicod+isrcod+3.gt.size(scod,1)) stop 'lencod'
nscode=isicod+isrcod
line=' '
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  close (8)
  return
endif
read (8, 85, end=186) line
85 format (a)
goto 186
! here if read error occurs
88 write(6,89)
89 format(/' *** ERROR DURING READ FROM INPUT FILE ***')
return
99 write (6, 100) i
100 format(/' ** ERROR DURING READ:', &
  ' PROBABLY NOT ENOUGH VIBRATIONAL LEVELS SUPPLIED')
return
! --------------------------------------------------------------
entry ptr1sg (fname,readpt)
line = fname
readpt = .true.
186 if (readpt) then
  l=1
  call get_token(line,l,filnam,lc)
  if(lc.eq.0) then
    write(6,190)
190     format(' FILENAME MISSING FOR POTENTIAL INPUT')
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
    write(6,195) filnam(1:lc)
195     format(' FILE ',(a),' NOT FOUND')
    return
  end if
! now call loapot(iunit,filnam) routine to read potential parameters
  call loapot(1,filnam)
end if
!      close (8)
irpot=1
return
! --------------------------------------------------------------
entry sav1sg ()
!  save input parameters for singlet-sigma + atom scattering
if (iscod(3) .lt. iscod(2)) then
  write (6, 210) iscod(3), iscod(2)
210   format ('**  VMAX =',i3,' .LT. VMIN =',i3,' SET VMAX = VMIN')
iscod(3)=iscod(2)
endif
write (FUNIT_INP, 220) nvib, iscod(2),iscod(3)
220 format(3i4, t34,'nvib, vmin,vmax')
iofi=3
iofr=0
do 301 i=1,nvib
write (FUNIT_INP, 310) ivib(i),(iscod(iofi+j),j=1,2)
310 format (3i4, t50,'iv,jmin,jmax')
write (FUNIT_INP, 320) (rcod(iofr+j),j=1,4)
iofi=iofi+2
iofr=iofr+4
320 format(3g14.6,t50,'brot,drot,hrot'/f15.8,t50,'evib')
301 continue
write (FUNIT_INP, 85) potfil
return
end

end module mod_hiba01_1sg

module mod_basis
  use mod_assert, only: fassert
implicit none
contains

  logical function basis_exists(ibasty)
!     ------------------------------------------------------------------
  !
  !  returns true if the given base type exists
  !
  !  ibasty:    one of the supported basis types
  !
  !     ------------------------------------------------------------------
  use mod_comxbs, only: maxbas
  implicit none
  integer, intent(in) :: ibasty

  if ((ibasty >= 1) .and. (ibasty <= maxbas)) then
    basis_exists = .true.
  else
    if ((ibasty == 99) .or. (ibasty == 100)) then
      !  99 is for user defined base not for molecule-molecule collision
      ! 100 is for user defined base for molecule-molecule collision
      basis_exists = .true.
    else
      basis_exists = .false.
    end if
  end if
  end function

  integer function basis_get_isa(ibasty, ispar)
  !     ------------------------------------------------------------------
  !
  !     returns the value of the isa parameter if the base has an isa parameter. returns zero otherwise
  !
  !  ibasty:    one of the supported basis types
  ! 
  !  isparam:   the array containing the values of the integer typed base 
  !              specific parameters
  !
  !     ------------------------------------------------------------------
  implicit none
  integer, intent(in) :: ibasty
  integer, intent(in), dimension(:) :: ispar
  integer :: isa
  ASSERT( basis_exists(ibasty) )

  select case (ibasty)
  case (2)
    isa = ispar(6)
  case (3)
    isa = ispar(4)
  case (4)
    isa = ispar(5)
  case (5)
    isa = ispar(4)
  case (11)
    isa = ispar(4)
  case (14)
    isa = ispar(4)
  case (19)
    isa = ispar(3)
  case default
    isa = 0
  end select
  basis_get_isa = isa

  end function
  logical function basis_is_singlet(ibasty)
  !     ------------------------------------------------------------------
  !
  !     checks if a basis is for a singlet system
  !
  !     ------------------------------------------------------------------
  implicit none
  integer, intent(in) :: ibasty
  integer, parameter :: singlet_base_types(*) = [1, 6, 8, 9, 11, 16, 17, 18, 21, 24, 25, 27, 29, 30]
  integer :: i
  do i = 1, size(singlet_base_types)
    if ( singlet_base_types(i) == ibasty ) then
      basis_is_singlet = .true.
      return
    end if
  end do
  basis_is_singlet = .false.

  end function basis_is_singlet

end module mod_basis

