#include "assert.h"
#include "unused.h"
module mod_hiba06_stp
contains
! systp (savstp/ptrstp) defines, save variables and reads                *
!                  potential for symmetric top/atom scattering           *
!                  this routine for sym top with inversion doubling      *
! ----------------------------------------------------------------------
subroutine bastp (bqs, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of a symmetric top molecule possessing inversion doubling (e.g. NH3)
!  with a structureless atom or with an uncorrugated surface
!  (basis type = 6)
!
!  authors:  millard alexander
!  revision:  11-aug-2011 by p.j.dagdigian
!    rotational channel basis now chosen by j < jmax and ehold < emax
!    as in bastp1 basis routine.  also changed printout of basis in bastst
!    also, added inversion splitting, labeled as delta
!  revision:  4-jun-2013 by q. ma (fix a bug in counting anisotropic
!     terms)
!  revision:  22-oct-2013 by p.j.dagdigian (capable of setting up 3
!     nuclear spin modifications or nuclear spin = 1)
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains rotational quantum numbers for each
!              channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains symmetry of inversion vibration for
!              each channel times k
!  Note that we have adopted the following convention for the symmetry
!  index "is" so that on return the symmetric top molecular levels can be
!  uniquely identified by the two labels "j" and "is":
!           is = +/- k
!       where k is the projection quantum number
!       the plus sign refers to the "plus" inversion levels and the minus
!       sign to the "minus" inversion levels
!       for k=0, where only eps=+1 states exist, the index "is" is equal to zero
!    jhold:    on return contains rotational quantum numbers for each
!              rotational level
!    ehold:    on return contains energy in hartrees of each rotational
!              level
!    ishold:   on return contains symmetry index of each energetically
!              distinct level
!    nlevel:   on return contains number of energetically distinct
!              rotational levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              rotational levels used in channel basis which are open
!              asymptotically
!    ietmp:    scratch array used to create channel list
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!              note that this is ignored in molecule-surface collisions!!!
!    jtot:     total angular momentum
!              in cc calculation jtot is the total angular momentum
!              in cs calculation jtot is the l-bar quantum number
!    flaghf:   if .true., then system has half-integer spin
!              if .false., then system has integer spin
!    flagsu:   if .true., then molecule-surface collisons
!    csflag:   if .true., then coupled-states calculation
!              if .false., the close-coupled calculation
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true., then the molecule posseses interchange symmetry
!              (e.g. NH3), so that only the ortho or para levels will be
!              included depending on the value of the parameter iop in common
!              /cosysi/ (see below)
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              only those channels are included for which
!                  (-1)**(parity+l-jtot)=jlpar
!              where parity designates the parity of the molecular state
!              [by definition parity=isym (-1)^k - see Eq. (A21) of
!              S. Green, J. Chem. Phys. 73, 2740 (1980) ]
!              in cs calculation jlpar is set equal to 1 in calling program
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to subroutines
!              propag and soutpt.  ntop is set in basis only if nu = numin
!              otherwise it is unchanged from the value supplied by the
!              calling program
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     A=B rotational constant for prolate top
!    crot:     C rotational constant for prolate top
!    delta:    inversion splitting
!    emax:     the maximum rotational energy (in cm-1) for a channel to be
!              included in the basis
!  variables in common block /cosysi/
!    nscod:    total number of variable names which are passed to HINPUT
!              nscod must equal isrcod + isicod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    numpot:   the number of the potential used, this variable is passed
!              to the pot subroutine
!    ipotsy:   cylindrical symmetry of potential.  If ihomo = .true., only
!              terms with mu equal to an integral multiple of ipotsy
!              can be included in the potential.  Example:  for NH3, ipotsy = 3.
!    iop:      nuclear spin modification label for molecular levels. If ihomo=.true.
!              then for molecules with equivalent nuclei of spin = 1/2 (e.e. NH3) only
!              para levels with k a multiple of 3, and only symmectric inverson levels
!              with k = 0, will be included if iop=1; only ortho levels
!              with k not a multiple of 3 will be included if iop=-1.
!                For a molecule with nuclear spins = 1 (e.g. ND3), set iop as follows:
!              iop =  2 for rotational levels with A1 nuclear spin symmetry
!              iop = -1 for rotational levels with A2 nuclear spin symmetry
!              iop =  1 for rotational levels with E nuclear spin symmetry
!    jmax:     the maximum rotational angular momentum for the symmetric top
!               the zero of energy is assumed to be the j=0, k=0 level
!  variable in module mod_conlam




!  variable in common block /coconv/
!   econv:        conversion factor from cm-1 to hartrees
!   xmconv:       converson factor from amu to atomic units
!  subroutines called:
!   vlmstp:    returns angular coupling coefficient for particular
!              choice of channel index
! --------------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coamat, only: ietmp ! ietmp(1)
use mod_conlam, only: nlam
use mod_cosysi, only: ispar
use mod_cosysr, only: rspar
use mod_hibasutil, only: vlmstp
use constants, only: econv, xmconv
use mod_par, only: iprint
use mod_parbas, only: lammin, lammax, mproj
use mod_par, only: boundc
use mod_ered, only: ered, rmu
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
real(8), intent(out), dimension(:) :: sc1
real(8), intent(out), dimension(:) :: sc2
real(8), intent(out), dimension(:) :: sc3
real(8), intent(out), dimension(:) :: sc4
type(ancou_type), intent(out), allocatable, target :: v2
type(bqs_type), intent(out) :: bqs
type(ancouma_type), pointer :: ancouma
logical flaghf, csflag, clist, flagsu, ihomo, bastst
dimension jhold(1), ehold(1), &
          ishold(1)
integer, pointer :: nterm, numpot, ipotsy, iop, jmax
real(8), pointer :: brot, crot, delta, emax
integer, dimension(nmax) :: k ! on return contains projection quantum number for each channel
integer, dimension(nmax) :: ieps ! on return contains epsilon label for each channel
integer, dimension(nmax) :: jtemp ! scratch array used to create channel list
integer, dimension(nmax) :: ktemp ! scratch array used to create channel list

nterm=>ispar(1); numpot=>ispar(2); ipotsy=>ispar(3); iop=>ispar(4); jmax=>ispar(5)
brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); emax=>rspar(4)


UNUSED_DUMMY(sc1)
sc1(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(sc2)
sc2(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(sc3)
sc3(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(sc4)
sc4(1) = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.

zero = 0.d0
two = 2.d0
!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (6, 5)
  write (9, 5)
5  format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
  stop
end if
if (flagsu .and. .not. csflag) then
  write (6, 6)
  write (9, 6)
6 format &
   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
  stop
end if
do 25  i = 1, nterm
  if (ihomo .and. mod(mproj(i), ipotsy) .ne. 0 ) then
    write (20, 20) i, mproj(i), ipotsy
    write (6, 20) i, mproj(i), ipotsy
20     format(/' *** MPROJ(',i2,')=',i2, &
     ' .NE. INTEGER MULTIPLE OF', i2, ' FOR IHOMO=T; ABORT ***')
    stop
  end if
25 continue
xjtot = jtot
nsum = 0
do 35  i = 1, nterm
  if (mproj(i) .gt. lammin(i) ) then
    write (6, 30) i, mproj(i), lammin(i)
    write (9, 30) i, mproj(i), lammin(i)
30     format (' *** MPROJ=',i2,' > LAMMIN=',i2, &
            ' FOR TERM',i2,'; ABORT ***')
    stop
  end if
  nsum = nsum + lammax(i) - lammin(i) + 1
35 continue
if (nsum .ne. nlam) then
  write (6, 45) nsum, nlam
  write (9, 45) nsum, nlam
45   format (' *** NLAM IN BASIS=', i3,' .NE. NLAM FROM SYSDAT=', &
           i3, '; ABORT ***')
  stop
end if
if (bastst) write (6, 46) nsum
write (9, 46) nsum
46 format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS IN POTENTIAL =', &
        i3)
nlam = nsum
if (clist) then
  if (flagsu) then
    if (ihomo) then
      if (bastst) &
      write (6,55) rmu * xmconv, brot, crot, delta, ipotsy, iop, &
                   ered * econv, jtot, nu
      write (9,55) rmu * xmconv, brot, crot, delta, ipotsy, iop, &
                   ered * econv, jtot, nu
55     format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **', &
      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=', f7.3, &
           '  DELTA=',f7.3,/,'     POT-SYM=', i2,'  O/P=',i2, &
           '  E=', f7.2, '       LBAR=', i5, '  NU=', i3)
    else
      write (6,60) rmu * xmconv, brot, crot, delta, ipotsy, &
                   ered * econv, jtot, nu
      write (9,60) rmu * xmconv, brot, crot, delta, ipotsy, &
                   ered * econv, jtot, nu
60     format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **', &
      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=', f7.3, &
           '  DELTA=',f7.3,/,'     POT-SYM=', i2, &
           '     E=', f7.2, '       LBAR=', i5, '  NU=', i3)
    end if
  else
    if (csflag) then
      if (ihomo) then
        if (bastst) &
        write (6,65) rmu * xmconv, brot, crot, delta, ipotsy, iop, &
                   ered * econv, jtot, nu
        write (9,65) rmu * xmconv, brot, crot, delta, ipotsy, iop, &
                   ered * econv, jtot, nu
65       format(/,' **  CS SYMMETRIC TOP **', &
      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3, &
           '  DELTA=',f7.3,/,'     POT-SYM=', i2,'  O/P=',i2, &
           '  E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
      else
        if (bastst) &
        write (6,75) rmu * xmconv, brot, crot, delta, ipotsy, &
                   ered * econv, jtot, nu
        write (9,75) rmu * xmconv, brot, crot, delta, ipotsy, &
                   ered * econv, jtot, nu
75       format(/,' **  CS SYMMETRIC TOP **', &
      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3, &
           '  DELTA=',f7.3/,'     POT-SYM=', i2, &
           '     E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
      end if
    else
      if (ihomo) then
        if (bastst) &
        write (6,80) rmu * xmconv, brot, crot, delta, ipotsy, iop, &
                     ered * econv, jtot, jlpar
        write (9,80) rmu * xmconv, brot, crot, delta, ipotsy, iop, &
                     ered * econv, jtot, jlpar
80       format(/,' **  CC SYMMETRIC TOP **', &
      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3, &
           '  DELTA=',f7.3/,'     POT-SYM=', i2,'  O/P=',i2, &
           '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
      else
        if (bastst) &
        write (6,85) rmu * xmconv, brot, crot, delta, &
                     ered * econv, jtot, jlpar
        write (9,85) rmu * xmconv, brot, crot, delta, &
                     ered * econv, jtot, jlpar
85       format(/,' **  CC SYMMETRIC TOP **', &
      /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3, &
           '  DELTA=',f7.3/, &
           '     E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
      end if
    end if
  end if
  if (.not. flagsu) write (9,90) rcut
90   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
end if
!
!  set up rotational basis
!
!  first set up list of j,k levels within the requested
!  nuclear spin modification (iop value)
nlist = 0
do 105  ki = 0, jmax
  do 105  ji = ki, jmax
    numeps = 2
    if (ki.eq.0) numeps = 1
    do 105 iep = 1, numeps
!  check that level has correct nuclear permutation symmetry
      if (ihomo) then
!  levels with iop = -1 (called ortho levels for NH3) and
!  iop =2 (called para levels for NH3) have ki a multiple of 3
!  all other levels have iop = 1 (called para levels for NH3)
        if (iop.eq.-1 .and. mod(ki,ipotsy).ne.0) goto 105
        if (iop.eq.2 .and. mod(ki,ipotsy).ne.0) goto 105
        if (iop.eq.1 .and. mod(ki,ipotsy).eq.0) goto 105
      end if
!  check to make sure rotational energy is less than emax
      roteng = brot * ji * (ji + 1) + (crot - brot) * ki * ki
      if (roteng .gt. emax) go to 105
      nlist = nlist + 1
      jtemp(nlist) = ji
      ktemp(nlist) = ki
      ietmp(nlist) = 1 - 2*(iep - 1)
!  add inversion splitting to - inversion levels
      isym = -ietmp(nlist) * (-1) ** jtemp(nlist)
      if (iop.eq.2) isym = -isym
      if (isym .eq. -1) roteng = roteng + delta
      ehold(nlist)= roteng / econv
105 continue
!  now sort this list in terms of increasing energy
if (nlist .gt. 1) then
  do 120 i1 = 1, nlist - 1
    esave = ehold(i1)
    do 115 i2 = i1 + 1, nlist
      if (ehold(i2) .lt. esave) then
!  state i2 has a lower energy than state i1, switch them
        esave = ehold(i2)
        ehold(i2) = ehold(i1)
        ehold(i1) = esave
        jsave = jtemp(i2)
        jtemp(i2) = jtemp(i1)
        jtemp(i1) = jsave
        ksave = ktemp(i2)
        ktemp(i2) = ktemp(i1)
        ktemp(i1) = ksave
        isave = ietmp(i2)
        ietmp(i2) = ietmp(i1)
        ietmp(i1) = isave
      end if
115     continue
120   continue
end if
!  print this list if bastst = .true.
if (bastst) then
  write (6, 130)
130   format (/,10x, &
    'SORTED LEVEL LIST',/,'   N   J   K  EPS INV   EINT(CM-1)')
  do  140 i = 1, nlist
    isym = -ietmp(i) * (-1) ** jtemp(i)
    if (iop.eq.2) isym = -isym
    ecm = ehold(i) * econv
    write (6, 135) i, jtemp(i), ktemp(i), ietmp(i), isym, ecm
135     format (5i4,f10.3)
140   continue
end if
!  now set up channel and level list for scattering calculation
call bqs%init(nmax)
n = 0
nlevel = 0
do 170  njk = 1, nlist
  ji = jtemp(njk)
  ki = ktemp(njk)
  nlevel = nlevel + 1
  jhold(nlevel) = ji
!  ishold set equal to k times inversion symmetry
  ishold(nlevel) = -ietmp(njk) * (-1.d0)**ji * ki
!  originally, ishold set equal to k times epsilon
!        ishold(nlevel) = ietmp(njk) * ki
!
!  here for cs calculations; include levels only if j is at least
!  equal to coupled states projection index
  if (csflag) then
!  here for cs calculations; include state only if j at least equal to coupled
!  states projection index
    if (ji .ge. nu) then
      n = n + 1
      if (n .gt. nmax) then
        write (6, 150) n, nmax
        write (9, 150) n, nmax
150         format(/' *** NCHANNELS=', i5, &
            ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
        stop
      end if
      bqs%jq(n) = ji
      k(n) = ki
      ieps(n) = ietmp(njk)
!  is set equal to k times inversion symmetry
      bqs%inq(n) = -ietmp(njk) * (-1.d0)**ji * ki
!  originally, ishold set equal to k times epsilon
!            bqs%inq(n) = ietmp(njk) * ki
      eint(n) = ehold(njk)
      bqs%lq(n) = jtot
      bqs%length = n
      cent(n) = jtot * (jtot + 1)
    end if
  else if (.not. csflag) then
!
!  here for cc calculations.  first calculate range of orbital angular
!  momentum quantum numbers allowed for this state
!
!  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem. Phys.
!  73, 2740 (1980) ]
      isym = -ietmp(njk) * (-1) ** ji
      if (iop.eq.2) isym = -isym
      ipar = isym * (-1) ** ki
      lmax = jtot + ji
      lmin = iabs (jtot - ji)
      do 155  li = lmin, lmax
        lpar = (-1) ** (li - jtot)
!  check to see if this channel has the desired parity, if so, include it
        if (ipar * lpar .eq. jlpar) then
          n = n + 1
          if (n .gt. nmax) then
            write (6, 150) n, nmax
            write (9, 150) n, nmax
            stop
          end if
          bqs%jq(n) = ji
          k(n) = ki
          ieps(n) = ietmp(njk)
!  is set equal to k times inversion symmetry
          bqs%inq(n) = -ietmp(njk) * (-1.d0)**ji * ki
!  originally, ishold set equal to k times epsilon
!                bqs%inq(n) = ishold(nlevel)
          eint(n) = ehold(nlevel)
          bqs%lq(n) = li
          bqs%length = n
          cent(n) = li * (li + 1)
        end if
155       continue
    end if
170 continue
!  also determine number of levels which are open
nlevop = 0
do 250  i = 1, nlevel
  if (ehold(i) .le. ered) then
    nlevop = nlevop + 1
  end if
250 continue
!  now check to see if any of the open channels are closed at r=rcut
!  this is not done for molecule-surface collisions or for rcut < 0
if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
  emin = 1.e+7
  do 290  i = 1, n
    if (eint(i) .le. ered) then
!  here if channel is open asymptotically
      if ( jtot * (jtot + 1) / (two * rmu * rcut * rcut) &
          .gt. (ered - eint(i)) ) then
!  here if channel is open asymptotically but closed at r = rcut
        if (eint(i) .lt. emin) emin = eint(i)
!  emin now contains the lowest channel energy for which this condition is met
      end if
    end if
290   continue
!  now eliminate all channels with eint .ge. emin if any of the channels
!  are open asymptotically but closed at r = rcut
  if (emin.lt.ered) then
    nn = 0
    do 300 i = 1, n
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        eint(nn) = eint(i)
        bqs%jq(nn) = bqs%jq(i)
        ieps(nn) = ieps(i)
        bqs%inq(nn) = bqs%inq(i)
        cent(nn) = cent(i)
        k(nn) = k(i)
        bqs%lq(nn) = bqs%lq(i)
      end if
300     continue
!  reset number of channels
    n = nn
    bqs%length = n
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
    write (6, 303) nu, n, ntop
    write (9, 303) nu, n, ntop
303     format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3, &
            '; ABORT **',/, &
    '     CHECK RCUT')
    call exit
  end if
end if
!  now list channels if requested
if (clist) then
  if (.not.csflag) then
    if (bastst) write (6, 305)
    write (9,305)
305     format ( &
     /'   N   J   K  EPS INV  L    EINT(CM-1)')
  else
    if (bastst) write (6, 310) nu
    write (9,310) nu
310     format ( &
     /'   N   J   K  EPS INV  L    EINT(CM-1) ** NU = ',i2)
  end if
  do 330  i = 1, n
    if (k(i) .ne. 0) isym = bqs%inq(i) / k(i)
    ecm = eint(i) * econv
    if (bastst) &
      write (6, 320) i, bqs%jq(i), k(i), ieps(i), isym, bqs%lq(i), ecm
      write (9, 320) i, bqs%jq(i), k(i), ieps(i), isym, bqs%lq(i), ecm
320       format (6i4, f10.3)
330   continue
end if
!  now calculate coupling matrix elements
!  i counts v2 elements
!  inum counts v2 elements for given lambda
!  ilam counts numver of v2 matrices
i = 0
if (bastst.and. iprint.ge. 2) then
  write (6, 340)
  write (9, 340)
340   format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
end if
ilam = 0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 400 iterm = 1, nterm
  lbmin = lammin(iterm)
!  if bastst = .true., then get the matrix elements of the lb=0 term
!  in the potential
  if (bastst .and. iterm .eq. 1) lbmin = 0
  do 390 lb = lbmin, lammax(iterm)
!  ilam is the index for the next term in the potential matrix
!  lb is the actual value of lambda
    ilam = ilam + 1
    ancouma => v2%get_angular_coupling_matrix(ilam)
    mu = mproj(iterm)
    inum = 0
    do icol = 1, n
      do irow = icol, n
        lrow = bqs%lq(irow)
        if (csflag) lrow = nu
        call vlmstp (bqs%jq(irow), lrow, bqs%jq(icol), bqs%lq(icol), jtot, &
                     k(irow), k(icol), lb, mu, ieps(irow), &
                     ieps(icol), vee, csflag)
        if (vee .ne. zero) then
          i = i + 1
          inum = inum + 1
          call ancouma%set_element(irow=irow, icol=icol, vee=vee)
          if (bastst .and. iprint.ge.2) then
            write (6, 345) ilam, lb, icol, irow, i, vee
            write (9, 345) ilam, lb, icol, irow, i, vee
345               format (i4, 2i7, 2i6, i6, g17.8)
          end if
        end if
      end do
    end do
    if (bastst) then
      
      write (6, 370) ilam, lb, mu, ancouma%get_num_nonzero_elements()
      write (9, 370) ilam, lb, mu, ancouma%get_num_nonzero_elements()
370       format ('ILAM=',i3,' LAM=',i3,' MU=',i3, &
         ' LAMNUM(ILAM) = ',i6)
    end if
390   continue
400 continue
if (clist .and. bastst) then
  write (6, 460) v2%get_num_nonzero_elements()
  write (9, 460) v2%get_num_nonzero_elements()
460   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS ', &
           i6)
end if
return
end
!  -----------------------------------------------------------------------
subroutine systp (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for symmetric top
!   + atom scattering
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
! NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
!  current revision date: 22-jun-2011 by p.j.dagdigian
!    list of input parameters changed
!  -----------------------------------------------------------------------
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     A=B rotational constant for prolate top
!    crot:     C rotational constant for prolate top
!    delta:    inversion splitting
!    emax:     the maximum rotational energy (in cm-1) for a channel to be
!              included in the basis
!  variables in common block /cosysi/
!    nscode:   total number of variable names which are passed to HINPUT
!              nscode must equal isrcod + isicod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    numpot:   the number of the potential used, this variable is passed
!              to the pot subroutine
!    ipotsy:   cylindrical symmetry of potential.  Only terms with mu equal to
!              an integral multiple of ipotsy can be included in the potential
!              Example:  for NH3, ipotsy = 3
!    iop:      ortho/para label for molecular states. If ihomo=.true. then only
!              para states will be included if iop=1 and only ortho states if
!              iop=-1
!    jmax:     the maximum rotational quantum number for the asymmetric top
!  variable in common cosys/
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters.  Note that the ordering
!             of the variable names in scod must correspond to the ordering
!             of the variable names in cosysi followed by the ordering of
!             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, rspar
use funit, only: FUNIT_INP
use mod_parbas, only: lammin, lammax, mproj
use mod_hiutil, only: gennam, get_token
implicit none
integer, intent(inout) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: j, l, lc
logical existf
integer icod, ircod
character*1 dot
character*(*) fname
character*60 line, filnam, potfil
character*68 filnm1
parameter (icod = 5, ircod = 4)
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
integer, pointer, save :: nterm, numpot, ipotsy, iop, jmax
real(8), pointer, save :: brot, crot, delta, emax
nterm=>ispar(1); numpot=>ispar(2); ipotsy=>ispar(3); iop=>ispar(4); jmax=>ispar(5)
brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); emax=>rspar(4)

scod(1)='NTERM'
scod(2)='NUMPOT'
scod(3)='IPOTSY'
scod(4)='IOP'
scod(5)='JMAX'
scod(6)='BROT'
scod(7)='CROT'
scod(8)='DELTA'
scod(9)='EMAX'
scod(10)='LAMMIN'
scod(11)='LAMMAX'
scod(12)='MPROJ'
nscode = icod + ircod + 3
isicod = icod
isrcod = ircod
irpot = 1
!  set default values for symmetric top scattering
!  (symmetric AB3 molecules)
nterm = 4
if (iread .eq. 0) then
  mproj(1) = 0
  mproj(2) = 3
  mproj(3) = 6
  mproj(4) = 9
  lammin(1) = 1
  lammin(2) = 3
  lammin(3) = 6
  lammin(4) = 9
  lammax(1) = 8
  lammax(2) = 8
  lammax(3) = 8
  lammax(4) = 9
  ipotsy = 3
  iop = 1
  jmax = 5
  niout=11
!
!  INDOUT VALUES TO BE SET ***
  niout=11
  indout(1)=0
  indout(2)=1
  indout(3)=-1
  indout(4)=2
  indout(5)=-2
  indout(6)=3
  indout(7)=-3
  indout(8)=4
  indout(9)=-4
  indout(10)=5
  indout(11)=-5
endif
potfil=' '
if (iread .eq. 0) return
!  line 18
read (8, *, err=80) ipotsy, iop
!  line 19
read (8, *, err=80) jmax, emax
!  line 20
read (8, *, err=80) brot, crot, delta
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
entry ptrstp (fname,readpt)
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
entry savstp (readpt)
!  save input parameters for symmetric top + atom scattering
!  the order of the write statements should be identical to the read statement
!  above. for consistency with the data file written by gendat, format
!  statements should reserve the first 30 spaces for data, spaces 31-33 should
!  be left blank, and the names of the variables should be printed in spaces
!  34-80
!  line 18:
write (FUNIT_INP, 220) ipotsy, iop
220 format (2i4, 22x,'   ipotsy, iop')
!  line 20
write (FUNIT_INP, 230) jmax, emax
230 format (i4, f12.5,14x, '   jmax, emax')
!  line 21
write (FUNIT_INP, 250) brot, crot, delta
250 format(3f9.4, 6x, 'brot, crot, delta')
write (FUNIT_INP, 60) potfil
return
end
end module mod_hiba06_stp
