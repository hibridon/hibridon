#include "assert.h"
#include "unused.h"
module mod_hiba14_2del
  use mod_assert, only: fassert
contains
! sy2del (sav2de/ptr2de) defines, saves variables and reads              *
!                  potential for doublet-delta scattering                *
! --------------------------------------------------------------------
subroutine ba2del (bqs, jhold, ehold, ishold, nlevel, &
                  nlevop, ifi, c32, c52, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of a 2delta molecule in an intermediate coupling basis with a
!  structureless atom or with an uncorrugated surface
!  revised from 'ba2pi' basis routine (authors millard alexander
!  and hans-joachim werner
!  authors:  paul dagdigian and millard alexander
!  current revision date:  7-apr-2003 by pjd
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains rotational quantum numbers for each
!              channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains symmetry index of each channel
!  Note that we have adopted the following convention for the symmetry
!  index "is" so that on return the doublet delta molecular levels can be
!  uniquely identified by the two labels "j" and "is":

!           for F1 levels is = eps
!           for F2 levels is = 2*eps
!       where  eps = +1 or -1 is the "true" case (a) symmetry index
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
!    ifi:      on return contains Fi label for each channel
!    c32:      on return contains omega=3/2 component in mixed states
!    c52:      on return contains omega=5/2 component in mixed states
!    sc4:      scratch array
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!              note that this is ignored in molecule-surface collisions!!!
!    jtot:     total angular momentum
!              in cc calculation jtot+1/2 is the total angular momentum
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
!    ihomo:    if .true., then homonuclear molecule
!              only the s or a levels will be included depending on the
!              value of the parameter isa in common cosysi/ (see below)
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              only those channels are included for which
!                  (-1)**(parity+l-jtot)=jlpar
!              where parity designates the parity of the molecular state
!              (by definition parity=eps*(-1)**(j-1/2) )
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
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     rotational constant in cm-1
!    aso:      spin-orbit constant in cm-1
!    p, q:     lambda-doubling constants cm-1
!  variables in common block /cosysi/
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    jmax:     the maximum rotational angular momenta for each channel
!              in each spin-orbit manifold with convention
!              omega .le. j .le. jmax+0.5
!    igu:      permutation inversion symmetry of electronic state
!              igu=1 for gerade states, igu=-1 for ungerade states
!              for heteronuclear molecules igu should be +1
!    isa:      s/a label for molecular states. if ihomo=.true. then only
!              s states will be included if isa=1 and only a states if
!              isa=-1
!    npar:     number of symmetry doublets included (npar=2 will ensure
!              both lambda doublets; npar=1, just eps=1 levels, npar=-1,
!              just eps=-1 levels
!  subroutines called:
!   vlm2del:    returns angular coupling coefficient for particular
!              choice of channel index
! --------------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam
use mod_cosysi, only: ispar
use mod_cosysr, only: rspar
use constants, only: econv, xmconv
use mod_par, only: iprint
use mod_parbas, only: lammin, lammax, mproj
use mod_par, only: boundc
use mod_ered, only: ered, rmu
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable, target :: v2
type(bqs_type), intent(out) :: bqs
type(ancouma_type), pointer :: ancouma
logical flaghf, csflag, clist, flagsu, ihomo, bastst
dimension jhold(1), ehold(1), ifi(1), &
          c32(1), c52(1), ieps(2), ishold(1), sc4(1)
!   econv is conversion factor from cm-1 to hartrees
!   xmconv is converson factor from amu to atomic units
data ieps / -1, 1 / , ione, itwo, min10 / 1, 2, -10 /
integer, pointer :: nterm, jmax, igu, isa, npar
real(8), pointer :: brot, aso, p, q
nterm=>ispar(1); jmax=>ispar(2); igu=>ispar(3); isa=>ispar(4); npar=>ispar(5)
brot=>rspar(1); aso=>rspar(2); p=>rspar(3); q=>rspar(4)

UNUSED_DUMMY(sc4)

pi2 = 1.570796327d0
zero = 0.d0
one = 1.d0
two = 2.d0
four = 4.d0
half = 0.5d0
quart = 0.25d0
xjtot = jtot + half
xnu = nu + half
!  check for consistency in the values of flaghf and csflag
if (.not.flaghf) then
  write (6, 5)
  write (9, 5)
5  format (' *** FLAGHF = .FALSE. FOR DOUBLET SYSTEM; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (flagsu .and. .not. csflag) then
  write (6, 6)
  write (9, 6)
6 format &
   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
!  check for consistency in values of nterm, lammin, lammax, mproj
if (nterm .lt. 2) then
  write (6, 7) nterm
  write (9, 7) nterm
7   format (' *** NTERM=',i2, &
   ' .LT. 2 FOR 2-PI BASIS; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (nterm .gt. 2) then
  write (6, 8) nterm
  write (9, 8) nterm
8   format (' *** NTERM=',i2, &
   ' .GT. 2 FOR 2-PI BASIS; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
do 10  i = 1, nterm
  if (ihomo .and. lammin(i) .lt. 2) then
    write (9, 9) i, lammin(i)
    write (6, 9) i, lammin(i)
9     format(/' *** LAMMIN(',i2,')=',i2, &
           ' < 2 FOR IHOMO=T; ABORT ***')
    if (bastst) then
      return
    else
      call exit
    end if
  end if
10 continue
nsum = 0
do 13  i = 1, nterm
  if (mproj(i) .gt. lammin(i) ) then
    write (6, 11) i, mproj(i), lammin(i)
    write (9, 11) i, mproj(i), lammin(i)
11     format (' *** MPROJ=',i2,' > LAMMIN=',i2, &
            ' FOR TERM',i2,'; ABORT ***')
    if (bastst) then
      return
    else
      call exit
    end if
  end if
  if (ihomo) then
    if (mod(lammax(i)-lammin(i),2) .ne. 0) then
      write (6, 12) i, lammin(i), lammax(i)
      write (9, 12) i, lammin(i), lammax(i)
12       format &
         (' *** IHOMO=T BUT ODD NO. OF TERMS FOR I=',i2, &
        /,'     LAMMIN=',i2,' LAMMAX=',i2,'; ABORT ***')
      if (bastst) then
        return
      else
        call exit
      end if
    end if
  end if
  if (ihomo) then
      nsum = nsum + (lammax(i) - lammin(i) )/2 + 1
  else
      nsum = nsum +  lammax(i) - lammin(i) + 1
  end if
13 continue
if (nlam .ne. nsum) then
  write (6, 14) nlam, nsum
  write (9, 14) nlam, nsum
14   format (' ** NLAM=',i2, ' RESET TO NLAM=',i2, &
          '    IN BASIS')
  nlam = nsum
end if
if (clist) then
  if (flagsu) then
    if (ihomo) then
      if (bastst) &
      write (6,20) rmu * xmconv, brot, aso, igu, isa, &
                   ered * econv, jtot, xnu
      write (9,20) rmu * xmconv, brot, aso, igu, isa, &
                   ered * econv, jtot, xnu
20     format(/,' **  2DEL INT. COUPLING UNCORRUGATED SURFACE **', &
      '     RMU=', f9.4,'  BROT=', f7.3, '  A-SO=', f7.2,/, &
            '     g/u=', i2,'  s/a=',i2, &
           '  E=', f7.2, '       LBAR=', i5, '  NU=', f5.1)
    else
      write (6,21) rmu * xmconv, brot, aso, &
                   ered * econv, jtot, xnu
      write (9,21) rmu * xmconv, brot, aso, &
                   ered * econv, jtot, xnu
21     format(/,' **  2DEL INT. COUPLING UNCORRUGATED SURFACE **', &
      '     RMU=', f9.4,'  BROT=', f7.3, '  A-SO=', f7.2,/, &
           '     E=', f7.2, '       LBAR=', i5, '  NU=', f5.1)
    end if
  else
    if (csflag) then
      if (ihomo) then
        if (bastst) &
        write (6,25) rmu * xmconv, brot, aso, igu, isa, &
                   ered * econv, jtot, xnu
        write (9,25) rmu * xmconv, brot, aso, igu, isa, &
                   ered * econv, jtot, xnu
25       format(/,' **  CS 2DEL INT. COUPLING ** RMU=', f9.4, &
             '  BROT=', f7.3, '  A-SO=', f7.2,/, &
            '     g/u=', i2,'  s/a=',i2, &
             '  E=', f7.2,'  LBAR=', i5, 2x,' NU=', f5.1)
      else
        if (bastst) &
        write (6,26) rmu * xmconv, brot, aso, &
                   ered * econv, jtot, xnu
        write (9,26) rmu * xmconv, brot, aso, &
                   ered * econv, jtot, xnu
26       format(/,' **  CS 2DEL INT. COUPLING ** RMU=', f9.4, &
             '  BROT=', f7.3, '  A-SO=', f7.2,/ &
             '     E=', f7.2,'  LBAR=', i5, 2x,' NU=', f5.1)
      end if
    else
      if (ihomo) then
        if (bastst) &
        write (6,30) rmu * xmconv, brot, aso, igu, isa, &
                     ered * econv, xjtot, jlpar
        write (9,30) rmu * xmconv, brot, aso, igu, isa, &
                     ered * econv, xjtot, jlpar
30       format(/,' **  CC 2DEL INT. COUPLING ** RMU=', f9.4, &
             '  BROT=', f7.3, '  A-SO=', f7.2, /, &
            '     g/u=', i2,'  s/a=',i2, &
             '  E=', f7.2, '  JTOT=', f5.1, 2x,' JLPAR=', i2)
      else
        if (bastst) &
        write (6,31) rmu * xmconv, brot, aso, &
                     ered * econv, xjtot, jlpar
        write (9,31) rmu * xmconv, brot, aso, &
                     ered * econv, xjtot, jlpar
31       format(/,' **  CC 2DEL INT. COUPLING ** RMU=', f9.4, &
             '  BROT=', f7.3, '  A-SO=', f7.2, /, &
             '     E=', f7.2, '  JTOT=', f5.1, 2x,' JLPAR=', i2)
      end if
    end if
  end if
  if (.not. flagsu) write (9,35) rcut
35   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
end if
!  first set up list of all case (a) levels for omega=3/2
!  for homonuclear molecules in gerade electronic states the e levels
!  are s for j=1/2, a for j=3/2, etc. while the f levels are a for
!  j=1/2, s for j=3/2, etc.  this reverses for ungerade states.
!  see table i of alexander and corey, j. chem. phys. 84, 100 (1986)
call bqs%init(nmax)
n = 0
do 45 ji = 1, jmax
  npbot=1
  nptop=2
  if (npar .eq. 1) npbot=2
  if (npar .eq. -1) nptop=1
  do 45 ip = npbot,nptop
!  ipar = 1 for e levels, ipar=-1 for f levels
    ipar = ieps(ip) * (-1) ** ji
    if (.not. ihomo .or. isa.eq.0 .or. &
       (ihomo .and. (ipar*igu .eq. isa) ) ) then
      n = n + 1
      if (n .gt. nmax) then
        write (6, 40) n, nmax
        write (9, 40) n, nmax
40         format(/' *** NCHANNELS=', i4, &
           ' .GT. MAX DIMENSION OF',i4,' ABORT ***')
        stop
      end if
      bqs%inq(n) = ieps(ip)
      bqs%jq(n) = ji
      bqs%length = n
    end if
45 continue
!  now assign omega values and energies for case (a) levels
do 50 i = 1, n
!  initialize the array ifi to be 1 (corresponding to F1)
  ifi(i) = 1
!  now set up arrays of internal energies for the case (a) levels
!    eint contains the energies of the omega = 3/2 levels
!    cent contains the energies of the omega = 5/2 levels
!    the matrix elements are given by j. m. brown, a. s.-c. cheung,
!    and a. j. merer, j. mol. spectr. 124, 464 (1987)
!    z bqs%inq(j + 1/2)^2
  x = bqs%jq(i) + one
  z = x * x
  eint(i) = - aso + brot * (z + two) - half * bqs%inq(i) &
             * (p + four * q) * x * (z - one)
  cent(i) = aso + brot * (z - two)
!  now calculate the mixing angle due to the j.s term in the hamiltonian
!  store this angle temporarily in the array c32 and store the 3/2 - 5/2
!  coupling matrix element temporarily in the array c52
  if (bqs%jq(i) .eq. 1) then
!  no mixing if j = 3/2
    c52(i) = 0.
    c32(i) = 0.
  else if (bqs%jq(i) .gt. 1) then
      c52(i) = - sqrt(z - four) * (brot - bqs%inq(i) &
              * half * q * x * (z - one))
    adum1=two*c52(i)
    adum2=eint(i)-cent(i)
    c32(i) = half * atan2 (adum1,adum2) + pi2
  end if
50 continue
!  n now contains the number of omega = 3/2 channels levels
! check if j=3/2 is to be assigned to F1 or F2 manifold
y = aso/brot
if (y .le. 2) then
!  reorder omega=3/2 channels so that j=3/2 comes at the end, appropriate
!  to F2
  nmn = 1
  if (bqs%jq(2) .eq. 1) nmn=2
! nmn is the number of j=1 levels

  eh1=eint(1)
  eh2=eint(2)
  centh1=cent(1)
  centh2=cent(2)
  ish1=bqs%inq(1)
  ish2=bqs%inq(2)
  do 53 i = nmn+1, n+1
    bqs%jq(i-nmn)=bqs%jq(i)
    cent(i-nmn)=cent(i)
    bqs%inq(i-nmn)=bqs%inq(i)
    c32(i-nmn)=c32(i)
    c52(i-nmn)=c52(i)
    eint(i-nmn)=eint(i)
53   continue
  if (nmn .eq. 2) then
    eint(n-1)=eh1
    eint(n)=eh2
    bqs%jq(n-1)=1
    bqs%jq(n)=1
    ifi(n-1)=2
    ifi(n)=2
    bqs%inq(n-1)=ish1
    bqs%inq(n)=ish2
  else
    eint(n)=eh1
    bqs%jq(n)=1
    ifi(n)=2
    bqs%inq(n)=ish1
  endif
endif
!  now add the omega = 5/2 levels
nn = n
do 60 i = 1, n
if ( bqs%jq(i) .eq. 1) then
!  here for j = 3/2, which can exist only for omega = 3/2 and hence
!  is always a pure case (a) state
  c32(i) = 1.
  c52(i) = 0.
  cent(i) = bqs%lq(i) * ( bqs%lq(i) + one )
else if ( bqs%jq(i) .ge. 2) then
!  here for j .ge. 5/2
!  the index nn points to the nominal omega = 5/2 channels
  nn = nn + 1
  if (nn .gt. nmax) then
    write (6, 40) nn, nmax
    write (9, 40) nn, nmax
    stop
  end if
  bqs%jq(nn) = bqs%jq(i)
  bqs%lq(nn) = bqs%lq(i)
! assign F2 label to each level
  ifi(nn) = 2
  bqs%inq(nn) = bqs%inq(i)
!  now determine array of mixing coefficients:
!    the c32 array indicates the amount of omega = 3/2 component in each
!    mixed level and the c52 array indicates the amount of omega = 5/2
!    component
!  also save the j.s matrix element temporarily in the variable h12
!  and save the case(a) energies temporarily in the variables e32 and e52
  h12 = c52(i)
  e32 = eint(i)
  e52 = cent(i)
  sn = sin (c32(i))
  cs = cos (c32(i))
  c52(i) = sn
  c32(i) = cs
  c52(nn) = cs
  c32(nn) = - sn
!  now calculate the energies of the mixed states and store these in the
!  array eint
  eint(i) = e32 * c32(i)**2 + e52 * c52(i)**2 &
            + two *  c32(i) * c52(i) * h12
  eint(nn) = e32 * c32(nn)**2 + e52 * c52(nn)**2 &
            + two *  c32(nn) * c52(nn) * h12
end if
60 continue
!  set n to equal total number of levels
n = nn
!  find lowest energy
emin = 1.e+7
do 65  i = 1, n
  if (eint(i) .lt. emin) emin = eint(i)
65 continue
!  form list of all energetically distinct rotational levels included in the
!  channel basis and their energies (with zero of energy set at lowest level)
nlevel = 0
do 70  i = 1, n
  nlevel = nlevel + 1
  ehold(nlevel) = (eint(i) - emin) / econv
  jhold(nlevel) = bqs%jq(i)
  ishold(nlevel) = bqs%inq(i) * ifi(i)
70 continue
!  now sort this list to put closed levels at end
!  also determine number of levels which are open
nlevop = 0
do 80  i = 1, nlevel - 1
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
        go to 80
      end if
75     continue
  end if
80 continue
if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
!  set up coupled-states channel basis (if desired)
if (csflag) then
  nn = 0
  do 85  i = 1, n
!  include only channels with j at least equal to coupled states
!  projection index
    if (bqs%jq(i) .ge. nu) then
      nn = nn + 1
      if (nn .gt. nmax) then
        write (6, 40) nn, nmax
        write (9, 40) nn, nmax
        stop
      end if
      eint(nn) = (eint(i) - emin) / econv
      bqs%jq(nn) = bqs%jq(i)
      c32(nn) = c32(i)
      c52(nn) = c52(i)
      bqs%inq(nn) = bqs%inq(i)
      ifi(nn) = ifi(i)
      bqs%lq(nn) = jtot
      cent(nn) = jtot * (jtot + 1)
    end if
85   continue
!  set number of coupled states channels
  n = nn
  bqs%length = n
!  set up close-coupled channel basis (if desired)
else if (.not. csflag) then
!  first move all indices of rotational levels to the top of the vectors
!  e.g. move bqs%jq(n) to bqs%jq(nmax), bqs%jq(n-1) to bqs%jq(nmax-1),  etc
  do 90  i = 1, n
!  move (n-i+1)th element to (nmax-i+1)th element
    inew = nmax - i + 1
    iold = n - i + 1
    eint(inew) = eint(iold)
    bqs%jq(inew) = bqs%jq(iold)
    c32(inew) = c32(iold)
    c52(inew) = c52(iold)
    bqs%inq(inew) = bqs%inq(iold)
    ifi(inew) = ifi(iold)
90   continue
  nn = 0
  do 100  i = 1, n
!  now take (nmax-n+i)th element, duplicate it as many times as is
!  required for rotational degeneray, and store the new elements in
!  nn, nn+1, ....
    ipoint = nmax - n + i
    ji = bqs%jq(ipoint)
    lmax = jtot + ji + 1
    lmin = iabs (jtot - ji)
    do 95  li = lmin, lmax
      ix = (-1) ** (ji + li - jtot) * bqs%inq(ipoint)
      if (ix .eq. jlpar) then
        nn = nn + 1
        if (nn .gt. nmax) then
         write (6, 40) nn, nmax
         write (9, 40) nn, nmax
         stop
        end if
        eint(nn) = (eint(ipoint) - emin) / econv
        bqs%jq(nn) = ji
        c32(nn) = c32(ipoint)
        c52(nn) = c52(ipoint)
        bqs%inq(nn) = bqs%inq(ipoint)
        ifi(nn) = ifi(ipoint)
        bqs%lq(nn) = li
        cent(nn) = li * ( li + 1)
      end if
95     continue
100   continue
!  set number of close coupled channels
  n = nn
  bqs%length = n
end if
!  now check to see if any of the open channels are closed at r=rcut
!  this is not done for molecule-surface collisions or for rcut < 0
if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
  emin = 1.e+7
  do 120  i = 1, n
    if (eint(i) .le. ered) then
!  here if channel is open asymptotically
      if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut) &
          .gt. (ered - eint(i)) ) then
!  here if channel is open asymptotically but closed at r = rcut
        if (eint(i) .lt. emin) emin = eint(i)
!  emin now contains the lowest channel energy for which this condition is met
      end if
    end if
120   continue
!  now eliminate all channels with eint .ge. emin if any of the channels
!  are open asymptotically but closed at r = rcut
  if (emin.lt.ered) then
    nn = 0
    do 130 i = 1, n
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        eint(nn) = eint(i)
        bqs%jq(nn) = bqs%jq(i)
        c32(nn) = c32(i)
        c52(nn) = c52(i)
        bqs%inq(nn) = bqs%inq(i)
        cent(nn) = cent(i)
        ifi(nn) = ifi(i)
        bqs%lq(nn) = bqs%lq(i)
      end if
130     continue
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
    write (6, 160) nu, n, ntop
    write (9, 160) nu, n, ntop
160     format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3, &
            '; ABORT **',/, &
    '     CHECK RCUT')
    call exit
  endif
end if
!  now list channels if requested
if (clist) then
  if (.not.csflag) then
    if (bastst) write (6, 200)
    write (9,200)
200     format ( &
     /'   N   J  EPS OMEGA   L    EINT(CM-1)  C-3/2     C-5/2')
  else
    if (bastst) write (6,210) xnu
    write (9,210) xnu
210     format ( &
     /'   N   J  EPS OMEGA   L    EINT(CM-1)  C-3/2     C-5/2', &
       '    **  NU=', f4.1/)
  end if
  do 230  i = 1, n
    fj = bqs%jq(i) + half
    ecm = eint(i) * econv
    if (bastst) &
    write (6, 220) i, fj, bqs%inq(i), ifi(i), bqs%lq(i), ecm, &
                   c32(i), c52(i)
    write (9, 220) i, fj, bqs%inq(i), ifi(i), bqs%lq(i), ecm, &
                   c32(i), c52(i)
220     format (i4, f5.1, i4, i4, i6, 3f10.3)
230   continue
endif
!  now calculate coupling matrix elements
!  i counts v2 elements
!  inum counts v2 elements for given lambda
!  ilam counts number of v2 matrices
!  ij is address of given v2 element in present v2 matrix
i = 0
istep = 1
if (ihomo) istep = 2
nlam2 = (lammax(2) - lammin(2))/istep + 1
nlam0 = (lammax(1) - lammin(1))/istep + 1
if (bastst .and. iprint .gt. 1) then
  write (6, 285)
  write (9, 285)
285   format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
end if
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 400 ilam = 1, nlam
  ancouma => v2%get_angular_coupling_matrix(ilam)
  ! ilam is the angular expansion label
  ! here for l=0 term in electrostatic potential
  if (ilam .le. nlam0) then
    lb = lammin(1) + (ilam - 1) * istep
  else
    lb = lammin(2) + (ilam - nlam0 - 1) * istep
  end if
!  lb is the actual value of lambda
  ij=0
  inum=0
  do icol = 1, n
    do  irow = icol, n
      ij = ntop * (icol - 1) +irow
      vee=0
      lrow = bqs%lq(irow)
      if (csflag) lrow = nu
      ! here for l=0 terms in potential (average potential)
      if (ilam .le. nlam0) then
        call vlm2del (bqs%jq(irow), lrow, bqs%jq(icol), bqs%lq(icol), jtot, &
                     ione, ione, lb, bqs%inq(irow), bqs%inq(icol), &
                     v3232, csflag)
        call vlm2del (bqs%jq(irow), lrow, bqs%jq(icol), bqs%lq(icol), jtot, &
                     itwo, itwo, lb, bqs%inq(irow), bqs%inq(icol), &
                     v5252,csflag)
        vee = c32(irow) * c32(icol) * v3232 &
            + c52(irow) * c52(icol) * v5252
      else
        ! here for l=4 terms in potential (difference potential)
        call vlm2del (bqs%jq(irow), lrow, bqs%jq(icol), bqs%lq(icol), jtot, &
                     ione, itwo, lb, bqs%inq(irow), bqs%inq(icol), &
                     v3252, csflag)
        call vlm2del (bqs%jq(irow), lrow, bqs%jq(icol), bqs%lq(icol), jtot, &
                     itwo, ione, lb, bqs%inq(irow), bqs%inq(icol), &
                     v5232, csflag)
        vee = c32(irow) * c52(icol) * v3252 &
            + c52(irow) * c32(icol) * v5232
      end if
      if(vee .ne. zero) then
        i=i+1
        inum=inum+1
        call ancouma%set_element(irow=irow, icol=icol, vee=vee)
        if(bastst .and. iprint .gt. 1) then
          write (6, 290) ilam, lb, icol, irow, i, vee
          write (9, 290) ilam, lb, icol, irow, i, vee
290             format (i4, 2i7, 2i6, g17.8)
        end if
      end if
    end do
  end do
  if (bastst) then
    write (6, 360) ilam, ancouma%get_num_nonzero_elements()
    write (9, 360) ilam, ancouma%get_num_nonzero_elements()
360     format ('ILAM=', i3, ' LAMNUM(ILAM) =', i6)
  end if
400 continue
if (clist .and. bastst) then
  write (6, 430) v2%get_num_nonzero_elements()
  write (9, 430) v2%get_num_nonzero_elements()
430   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', i4)
end if
! now multiply is array by ifi array, so that index will equal +/- 1 for
! F1 levels and +/- 2 for F2 levels
if (bastst) then
  write (6, 440)
  write (9, 440)
440      format &
   ('  **  FINAL CHANNEL LIST AFTER BASIS', &
    ' WITH Fi LABEL FOR INDEX')
  if (.not.csflag) then
    write (6, 200)
    write (9,200)
  else
    write (6,210) xnu
    write (9,210) xnu
  end if
end if
do 450 i = 1, nn
    bqs%inq(i) = ifi(i) * bqs%inq(i)
    fj = bqs%jq(i) + half
    ecm = eint(i) * econv
    if (bastst) then
      write (6, 220) i, fj, bqs%inq(i), ifi(i), bqs%lq(i), ecm, &
                   c32(i), c52(i)
      write (9, 220) i, fj, bqs%inq(i), ifi(i), bqs%lq(i), ecm, &
                   c32(i), c52(i)
    end if
450 continue
return
end
! ----------------------------------------------------------------------
subroutine vlm2del (jp, lp, j, l, jtot, iomegp, iomeg, lambda, &
                   iepsp, ieps, v, csflag)
!  subroutine to calculate v-lambda matrices for close-coupled and
!  coupled-states treatments of collisions of a molecule in a 2delta
!  electronic state
!  cc and cs matrix elements derived by p.j. dagdigian, in an extension
!  of the treatment by m.h. alexander, chem. phys. 92, 337 (1985)
!
!  note that for cc collisions of a 2pi molecule with a flat surface, the
!  coupling matrix elements [m.h. alexander, j. chem. phys. 80, 3485 (1984)]
!  are identical to the cs matrix elements here
!  -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element (bra)
!              minus 1/2
!    lp:       orbital angular momentum of left side of matrix element (bra)
!    j:        rotational quantum number of right side of matrix element (ket)
!              minus 1/2
!    l:        orbital angular momentum of right side of matrix element (ket)
!    jtot:     total angular momentum
!    iomegp:   omega quantum number of bra
!    iomeg:    omega quantum number of ket
!    lambda:   order of legendre term in expansion of potential
!    iepsp:    symmetry index of bra
!    ieps:     symmetry index of ket
!    v:        on return, contains matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index) minus 1/2
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum minus 1/2
!    for collisions of a 2delta molecule with a surface, nu is equivalent
!    to m (the projection of j along the surface normal) minus 1/2
!  subroutines called:
!     xf3j, xf6j
!  -----------------------------------------------------------------------
use mod_hiutil, only: xf3j, xf6j
implicit double precision (a-h,o-z)
integer iomegp, iomeg, jp, j, jtot, lp, l, lambda, &
        iepsp, ieps, nu, iphase
!      real half, v, xinu, xomegp, xomeg, xjp, xj, xjtot, xlp,
!     :     xl, xlamda, xnorm, zero, two
logical csflag
data half, ione, two, four, zero/ 0.5d0, 1, 2.d0, 4.d0, 0.0d0/
v = zero
xjp = jp + half
xj = j + half
xomegp = iomegp + half
xomeg = iomeg + half
xjtot = jtot + half
if (csflag) then
  nu = lp
  xnu = nu + half
end if
xlp = lp
xl = l
xlamda = lambda
iphase = ieps * iepsp * ((-1) ** (jp + j + lambda + 1))
if (iphase .eq. 1) return
if (csflag) then
  iphase = nu - iomeg
  xnorm = (2. * xjp + 1.) * (2. * xj + 1.)
  xnu = nu + half
  x = xf3j (xjp, xlamda, xj, -xnu, zero, xnu)
else
  iphase = jp + j + ione + jtot - iomeg
  xnorm = (2. * xjp + 1.) * (2. * xj + 1.) * (2. * lp + 1.) &
        * (2. * l + 1.)
  x = xf3j (xlp, xlamda, xl, zero, zero, zero)
  if  (x .eq. zero) return
  x = x * xf6j (xjp, xlp, xjtot, xl, xj, xlamda)
end if
if (x .eq. zero) return
if (iomegp .eq. iomeg) then
  x = x * xf3j (xjp, xlamda, xj, -xomegp, zero, xomegp)
else
  x = - ieps * x * xf3j (xjp, xlamda, xj, -xomegp, four, -xomeg)
end if
iphase = (-1) ** iabs(iphase)
v = iphase * x * sqrt(xnorm)
return
end
! ----------------------------------------------------------------
subroutine sy2del (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for doublet-delta
!   + atom scattering using werner-follmeg potential form
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  current revision date: 13-oct-1999 by mha
!  NB, this subroutine is indentical to sy2pi except for name changes
!  -----------------------------------------------------------------------
!    nlam:      the total number of angular coupling terms
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     rotational constant in cm-1
!    aso:      spin-orbit constant in cm-1
!    p, q:     lambda-doubling constants cm-1
!  variables in common block /cosysi/
!    nscode:  total number of system dependent parameters
!             nscode = isicod + isrcod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:   number of different angular terms in potential
!              nb for singlet sigma molecule, nterm should be 1
!    jmax:     the maximum rotational angular momenta for each channel
!              in each spin-orbit manifold with convention
!              omega .le. j .le. jmax+0.5
!    igu:      permutation inversion symmetry of electronic state
!              igu=1 for gerade states, igu=-1 for ungerade states
!              for heteronuclear molecules igu should be +1
!    isa:      s/a label for molecular states. if ihomo=.true. then only
!              s states will be included if isa=1 and only a states if
!              isa=-1
!    npar:     number of symmetry doublets included (npar=2 will ensure
!              both lambda doublets)
!  variable in common bloc /cosys/
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters
!  subroutines called: loapot(iunit,filnam)
use mod_coiout, only: niout, indout
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
integer :: j, l, lc
logical existf
character*1 dot
character*(*) fname
character*60 filnam, line, potfil
character*68 filnm1
save potfil
#include "common/comdot.F90"

integer, pointer, save :: nterm, jmax, igu, isa, npar
real(8), pointer, save :: brot, aso, p, q
nterm=>ispar(1); jmax=>ispar(2); igu=>ispar(3); isa=>ispar(4); npar=>ispar(5)
brot=>rspar(1); aso=>rspar(2); p=>rspar(3); q=>rspar(4)

isicod = 5
isrcod = 4
nscode = isicod+isrcod
scod(1) = 'NTERM'
scod(2) = 'JMAX'
scod(3) = 'IGU'
scod(4) = 'ISA'
scod(5) = 'NPAR'
scod(6) = 'BROT'
scod(7) = 'ASO'
scod(8) = 'P'
scod(9) = 'Q'
potfil = ' '
!  set default values for doublet delta scattering
nterm = 2
mproj(2) = 2
if (iread .eq. 0) then
  mproj(1) = 0
  lammin(1)= 1
  lammax(1) = 1
  lammin(2) = 2
  lammax(2) = 2
  niout=4
  indout(1)=-1
  indout(2)=1
  indout(3)=-2
  indout(4)=2
  npar=2
  igu=1
  isa=0
  jmax=3
endif
if(iread.eq.0) return
!  line 13
read (8, *, err=888) jmax, igu, isa, npar
!  line 14
read (8, *, err=888) brot, aso, p, q
!  line 15 name of file containing potential parameters
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
!
entry ptr2del (fname,readpt)
line = fname
readpt = .true.
286 if (readpt) then
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
!
entry sav2del ()
!  save input parameters for doublet-delta + atom scattering
!  line 13:
write (FUNIT_INP, 315) jmax, igu, isa, npar
315 format(4i4,18x,'jmax, igu, isa, npar')
!  line 14
write (FUNIT_INP, 320) brot, aso, p, q
320 format(f10.5,f10.4,2(1pg11.3),' brot, aso, p, q')
!  line 15
write (FUNIT_INP, 285) potfil
return
end
end module mod_hiba14_2del
