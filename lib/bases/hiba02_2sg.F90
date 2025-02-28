#include "assert.h"
#include "unused.h"
module mod_hiba02_2sg
contains
! sy2sg (sav2sg/ptr2sg) defines, saves variables and reads          *
!                  potential for doublet sigma scattering                *
! --------------------------------------------------------------------
subroutine ba2sg (bqs, jhold, ehold, ishold, nlevel, &
                  nlevop, nrot, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of a 2sigma molecule in a hund's case(a) basis with a
!  structureless atom or with an uncorrugated surface
!  author:  millard alexander
!  fixed is label to take care of +/- symmetry of the sigma state
!  current revision date: 9-dec-2011 by p.j.dagdigian
!  --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains rotational quantum numbers for each
!              channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains symmetry index eps of each channel
!              eps = +1 or -1
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
!    nrot:     n quantum number of each level
!    sc2-sc4:  scratch arrays (only sc2 is used here)
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!              note that this is ignored in molecule-surface collisions!!!
!    jtot:     total angular momentum
!              in cc calculation jtot+1/2 is the total angular momentum
!              in cs calculation jtot is the l-bar quantum number
!    flaghf:   if .true., then system has half-integer spin
!              if .false., then system has integer spin (this is the case
!              here)
!    flagsu:   if .true., then molecule-surface collisons
!    csflag:   if .true., then coupled-states calculation
!              if .false., the close-coupled calculation
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true., then homonuclear molecule
!              if .false., then heteronuclear molecule
!              if the molecule is homonuclear (ihomo = .true.), either
!              only the s or the a rotational levels are included
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              only those channels are included for which
!                  (-1)**(parity+l-jtot)=jlpar
!              where parity is +1 for e levels and -1 for f levels
!              with the standard definition that e levels are those for
!              which eps*(-1)**(j-1/2-s) = 1
!              and f levels are those for which eps(-1)**(j-1/2-s) = -1
!              here s=0 for sigma-plus levels and s=1 for sigma-minus
!              levels
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
!    brotsg:   rotational constant in cm-1
!    gsr:      2sigma state spin-rotation constant in cm-1
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables (integer
!              plus real)
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    nrmax:    the maximum case (b) rotational angular momenta
!    npar:     number of symmetry doublets included (npar=2 will ensure
!              both spin doublets)
!    isym:     if isym=+1, then the electronic symmetry is sigma-plus
!              if isym=-1, then the electronic symmetry is sigma-minus
!    igu:      if igu=+1, then the inversion symmetry is gerade
!              if igu=-1, then the inversion symmetry is ungerade
!    isa:      s/a symmetry index, if the molecule is homonuclear (ihomo=t)
!              then, if isa=+1 then only the s-levels are included in the
!              basis, if isa=-1, then only the a-levels are included
!               zero of energy is taken to be n=0
!  variable in common block /coconv/
!   econv:        conversion factor from cm-1 to hartrees
!   xmconv:       converson factor from amu to atomic units
!  subroutines called:
!   vlm2sg:    returns angular coupling coefficient for particular
!              choice of channel index
! --------------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam
use mod_cosysi, only: ispar
use mod_cosysr, only: rspar
use mod_hibasutil, only: vlm2sg
use constants, only: econv, xmconv
use mod_par, only: iprint
use funit, only: FUNIT_INP
use mod_parbas, only: lammin, lammax, mproj
use mod_par, only: boundc
use mod_ered, only: ered, rmu
use mod_hitypes, only: bqs_type

implicit double precision (a-h,o-z)
integer, intent(out), dimension(:) :: nrot
real(8), intent(out), dimension(:) :: sc2
real(8), intent(out), dimension(:) :: sc3
real(8), intent(out), dimension(:) :: sc4
type(ancou_type), intent(out), allocatable, target :: v2
type(bqs_type), intent(out) :: bqs
type(ancouma_type), pointer :: ancouma
logical clist, csflag, flaghf, flagsu, ihomo, bastst
dimension jhold(1), ehold(1), &
          ieps(2), ishold(1)
data ieps / -1, 1 /
integer, pointer :: nterm, nrmax, npar, isym, igu, isa
real(8), pointer :: brotsg, gsr, drotsg, hrotsg
nterm=>ispar(1); nrmax=>ispar(2); npar=>ispar(3); isym=>ispar(4); igu=>ispar(5); isa=>ispar(6)
brotsg=>rspar(1); gsr=>rspar(2); drotsg=>rspar(3); hrotsg=>rspar(4)

UNUSED_DUMMY(sc2)
UNUSED_DUMMY(sc3)
UNUSED_DUMMY(sc4)

half = 0.5d0
zero = 0.d0
xjtot = jtot + half
xnu = nu + half
!  check for consistency in the values of flaghf and csflag
if (.not. flaghf) then
  write (6, 7)
  write (9, 7)
7   format (' *** FLAGHF = .FALSE. FOR DOUBLET SYSTEM; ABORT ***')
  stop
end if
if (flagsu .and. .not. csflag) then
  write (6, 8)
  write (FUNIT_INP, 8)
8   format &
   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
  stop
end if
!  check that isym equals +1 or -1
if (abs(isym).ne.1) then
  write (6, 112)
  write (6, 112)
112   format (' *** ISYM MUST EQUAL +1 OR -1; ABORT ***')
  call exit
end if
!  check for consistency in values of nterm, lammin, lammax, mproj
if (nterm .gt. 1) then
  write (6, 9) nterm
  write (9, 9) nterm
9   format (' *** NTERM=',i2, &
   ' .GT. 1 FOR 2-SIGMA BASIS; ABORT ***')
  stop
end if
nsum = 0
do 13  i = 1, nterm
  if (mproj(i) .gt. lammin(i) ) then
    write (6, 11) i, mproj(i), lammin(i)
    write (9, 11) i, mproj(i), lammin(i)
11     format (' *** MPROJ=',i2,' > LAMMIN=',i2, &
            ' FOR TERM',i2,'; ABORT ***')
    stop
  end if
  if (ihomo) then
    if (mod(lammax(i)-lammin(i),2) .ne. 0) then
      write (6, 12) i, lammin(i), lammax(i)
      write (9, 12) i, lammin(i), lammax(i)
12       format &
         (' *** IHOMO=T BUT ODD NO. OF TERMS FOR I=',i2, &
        /,'     LAMMIN=',i2,' LAMMAX=',i2,'; ABORT ***')
      stop
    end if
      nsum = nsum + (lammax(i) - lammin(i) )/2+1
  else
    nsum = nsum + lammax(i) - lammin(i)+1
  end if
13 continue
if (bastst) write (6, 15) nsum
if (bastst) write (9, 15) nsum
if (nsum.ne.nlam) write (6, 16) nsum, nlam
15 format (' ** TOTAL NUMBER OF ANISTROPIC TERMS IN POTENTIAL =', &
        i3)
16 format (' ** TOTAL NUMBER OF ANISTROPIC TERMS IN POTENTIAL =', &
        i3,' BUT NLAM=',i3,'  NLAM WILL BE REDEFINED')
nlam = nsum
if (clist) then
  if (flagsu) then
    if (ihomo) then
      write (9,20) rmu * xmconv, brotsg, gsr, isym, igu, isa, &
                   ered * econv, jtot, xnu
      if (bastst) &
      write (6,20) rmu * xmconv, brotsg, gsr, isym, igu, isa, &
                   ered * econv, jtot, xnu
20       format(/,' **  2SIGMA UNCORRUGATED SURFACE **', &
      '     RMU=', f9.4,'  BROT=', f7.3,'  G-SR=',g10.3, &
      '  +/-=', i2,' g/u=',i2, &
    '  s/a=', i2,/,'     E=', f7.2, '       LBAR=', i5, &
    '  NU=', f5.1)
    else
      write (9,24) rmu * xmconv, brotsg, gsr, isym, ered * econv, &
                 jtot, xnu
      if (bastst) &
      write (6,24) rmu * xmconv, brotsg, gsr, isym, ered * econv, &
                 jtot, xnu
24       format(/,' **  2SIGMA UNCORRUGATED SURFACE **', &
      '     RMU=', f9.4,'  BROT=', f7.3,'  G-SR=',g10.3, &
      '  +/-=', i2, &
       /,'     E=', f7.2, '       LBAR=', i5, &
    '  NU=', f5.1)
    end if
  else
    if (ihomo) then
      if (csflag) then
        write (9,25) rmu * xmconv, brotsg, gsr, isym, igu, isa, &
                     ered * econv, jtot, xnu
        if (bastst) &
        write (6,25) rmu * xmconv, brotsg, gsr, isym, igu, isa, &
                     ered * econv, jtot, xnu
25         format(/,' ** CS 2SIGMA ** RMU=', f8.4, &
             '  BROT=', f7.3,'  G-SR=',g10.3, &
         ' +/-=', i2,' g/u=', i2, &
        ' s/a=', i2,/,'     E=', f7.2,'  LBAR=', i5, 2x, &
        '  NU=', f5.1)
      else
        if (bastst) &
        write (6,30) rmu * xmconv, brotsg, gsr, isym, igu, isa, &
                     ered * econv, xjtot, jlpar
        write (9,30) rmu * xmconv, brotsg, gsr, isym, igu, isa, &
                     ered * econv, xjtot, jlpar
30         format(/,' **  CC 2SIGMA ** RMU=', f8.4, &
             '  BROT=', f7.3,'  G-SR=',g10.3, &
            ' +/-=', i2,' g/u=', i2, &
     ' s/a=', i2,/,'     E=', f7.2,'  JTOT=', f5.1, &
       2x,' JLPAR=', i2)
      end if
    else
      if (csflag) then
        write (9,31) rmu * xmconv, brotsg, gsr, isym, &
                     ered * econv, jtot, xnu
        if (bastst) &
        write (6,31) rmu * xmconv, brotsg, gsr, isym, &
                     ered * econv, jtot, xnu
31         format(/,' **  CS 2SIGMA ** RMU=', f9.4, &
             '  BROT=', f7.3,'  G-SR=',g10.3, '  +/-=', i2, &
               /,'     E=', f7.2,'  LBAR=', i5, 2x, &
     ' NU=', f5.1)
      else
        if (bastst) &
        write (6,32) rmu * xmconv, brotsg, gsr, isym, &
                     ered * econv, xjtot, jlpar
        write (9,32) rmu * xmconv, brotsg, gsr, isym, &
                     ered * econv, xjtot, jlpar
32         format(/,' **  CC 2SIGMA ** RMU=', f9.4, &
             '  BROT=', f7.3,'  G-SR=',g10.3, '  +/-=', i2, &
        /,'     E=', f7.2,'  JTOT=', f5.1, &
          2x,' JLPAR=', i2)
      end if
    end if
  end if
  if (.not. flagsu) write (9,35) rcut
35   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
end if
!  first set up list of all case (a) levels
call bqs%init(nmax)
n = 0
!  the n=0, eps=+1, j=1/2 level is
!     s for sigma-g-plus or for sigma-u-minus
!     a for sigma-g-minus or for sigma-g-plus
isplus = igu * isym * isa
do 45 nnrot = 0, nrmax
  do 40 ip = 1, npar
!  include only the eps=+1 level for n=0
    if (nnrot .eq. 0 .and. ip .eq. 1) go to 40
!  now calculate j for each level
    if (ip .eq. 1) then
      ji = nnrot -1
    else
      ji = nnrot
    end if
!  actual half integer value of j is ji + 1/2
!  if homonuclear, include level only if allowed
!  see table i of alexander and corey, j. chem. phys. 84, 100 (1986)
    if (.not.ihomo .or. isa.eq.0 .or. &
       (ihomo .and. ieps(ip)*(-1)**ji.eq.isplus)) then
       n = n + 1
       bqs%inq(n) = ieps(ip)
       bqs%jq(n) = ji
       nrot(n) = nnrot
!  now assign energies for case (a) level and store in array eint
!    the matrix elements are given by a. j. kotlar, r. w. field,
!    and j. i. steinfeld, j. mol. spectr. 80, 86 (1980)
       x = float(bqs%jq(n)) + 1.
       nn1 = x * (x - bqs%inq(n))
       eint(n) = brotsg * nn1 - drotsg*nn1**2 + hrotsg*nn1**3 &
                       - half * (1 - bqs%inq(n) * x) * gsr
!  next statement to take care of +/- symmetry of the sigma state
       bqs%inq(n) = bqs%inq(n) * isym
       bqs%length = n
    end if
40   continue
45 continue
!  n now contains the number of levels
!  assume zero of energy is n=0
emin = zero
!      emin = 1.e+7
!      do 65  i = 1, n
!        if (eint(i) .lt. emin) emin = eint(i)
!65    continue
!  form list of all energetically distinct rotational levels included in the
!  channel basis and their energies (with zero of energy set at lowest level)
nlevel = 0
do 70  i = 1, n
  nlevel = nlevel + 1
  ehold(nlevel) = (eint(i) - emin) / econv
  jhold(nlevel) = bqs%jq(i)
  ishold(nlevel) = bqs%inq(i)
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
      eint(nn) = (eint(i) - emin) / econv
      bqs%jq(nn) = bqs%jq(i)
      bqs%inq(nn) = bqs%inq(i)
      nrot(nn) = nrot(i)
      bqs%lq(nn) = jtot
      if (.not. boundc) then
           cent(nn) = jtot * (jtot + 1)
      else
           xjtot=jtot+0.5d0
           xj=bqs%jq(nn)+0.5d0
           xnu = nu+0.5d0
           cent(nn)=xjtot*(xjtot+1)+xj*(xj+1)-2*xnu*xnu
      endif

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
    bqs%inq(inew) = bqs%inq(iold)
    nrot(inew) = nrot(iold)
90   continue
  nn = 0
  do 100  i = 1, n
!  now take (nmax-n+i)th element, duplicate it as many times as is
!  required for rotational degeneray, and store the new elements in
!  nn, nn+1, ....
    ipoint = nmax - n + i
    ji = bqs%jq(ipoint)
!  ipar is the symmetry of the molecular wavefunction
!  eps(-1)**(j-1/2) for sigma plus states
!  eps(-1)**(j-1/2-1) for sigma minus state
!  ipar=1 for e levels and -1 for f levels
    ipar = (-1) ** (ji + (isym - 1)/2 ) * bqs%inq(ipoint)
    lmax = jtot + ji + 1
    lmin = iabs (jtot - ji)
    do 95  li = lmin, lmax
      ix = ipar * (-1) ** (li - jtot)
      if (ix .eq. jlpar) then
        nn = nn + 1
        eint(nn) = (eint(ipoint) - emin) / econv
        bqs%jq(nn) = ji
        bqs%inq(nn) = bqs%inq(ipoint)
        nrot(nn) = nrot(ipoint)
        bqs%lq(nn) = li
        cent(nn) = li * ( li + 1)
      end if
95     continue
100   continue
!  set number of close coupled channels
  n = nn
  bqs%length = n
end if
if (n .gt. nmax) then
  write (6, 110) n, nmax
  write (9, 110) n, nmax
110   format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF', &
         i4,' ABORT ***')
  stop
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
        bqs%inq(nn) = bqs%inq(i)
        nrot(nn) = nrot(i)
        cent(nn) = cent(i)
        bqs%lq(nn) = bqs%lq(i)
      end if
130     continue
!  reset number of channels
    n = nn
    bqs%length = n
  end if
end if
!  if no channels, return
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
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1
!  now list channels if requested
if (clist) then
  if (.not.csflag) then
    if (bastst) write (6, 200)
    write (9,200)
200     format ( &
     /'   N   J  EPS   N    L   EINT(CM-1)')
  else
    if (bastst) write (6, 210) xnu
    write (9,210) xnu
210     format &
     (/'   N   J  EPS   N    L   EINT(CM-1)', &
       '    **  NU=', f4.1/)
  end if
  do 230  i = 1, n
    fj = bqs%jq(i) + half
    ecm = eint(i) * econv
    if (bastst) &
    write (6, 220) i, fj, bqs%inq(i), nrot(i), bqs%lq(i), ecm
    write (9, 220) i, fj, bqs%inq(i), nrot(i), bqs%lq(i), ecm
220     format (i4, f5.1, 2i4, i6, 1f10.3)
230   continue
end if
!  now calculate coupling matrix elements
!  i counts v2 elements
!  inum counts v2 elements for given lambda
!  ilam counts numver of v2 matrices
i = 0
if (bastst.and.iprint.ge.2) then
  write (6, 285)
  write (9, 285)
285   format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
end if
lamsum = 0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 400 ilam = 1, nlam
!  ilam is the angular expansion label
  lb = ilam
  inum = 0
  ancouma => v2%get_angular_coupling_matrix(ilam)
  do icol = 1, n
    do irow = icol, n
      lrow = bqs%lq(irow)
      if (csflag) lrow = nu
      !  always initialize potential to zero
      vee = zero
      call vlm2sg (bqs%jq(irow), lrow, bqs%jq(icol), bqs%lq(icol), jtot, &
                   lb, bqs%inq(irow), bqs%inq(icol), vee, csflag)
      if (vee .ne. zero) then
        i = i + 1
        inum = inum + 1
        call ancouma%set_element(irow=irow, icol=icol, vee=vee)
        if (bastst.and.iprint.ge.2) then
          write (6, 290) ilam, lb, icol, irow, i, &
                         vee
          write (9, 290) ilam, lb, icol, irow, i, &
                         vee
290             format (i4, 2i7, 2i6, g17.8)
        end if
      end if
    end do  ! row
  end do  ! col
  if (bastst) then
    write (6, 370) ilam, ancouma%get_num_nonzero_elements()
    write (9, 370) ilam, ancouma%get_num_nonzero_elements()
370   format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
  end if
  lamsum = lamsum + ancouma%get_num_nonzero_elements()
400 continue
if (clist .and. bastst) then
  write (6, 420) v2%get_num_nonzero_elements()
  write (9, 420) v2%get_num_nonzero_elements()
420   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', &
           i6)
end if
if(bastst) call v2%print_summary(unit=9)
return
end
!  -----------------------------------------------------------------------
!  -----------------------------------------------------------------------
subroutine sy2sg (irpot, readpt, iread)
use mod_coiout, only: niout, indout
use mod_conlam, only: nlam
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, rspar
use mod_par, only: ihomo
use funit, only: FUNIT_INP
use mod_parbas, only: lammin, lammax, mproj
use mod_skip, only: nskip
use mod_hiutil, only: gennam, get_token
implicit none
!  subroutine to read in system dependent parameters for doublet-sigma
!   + atom scattering
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  current revision date: 8-apr-1997 by mha
!  -----------------------------------------------------------------------
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:   rotational constant for 2sigma state in cm-1
!    gsr:      2sigma state spin-rotation constant in cm-1
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different diabatic potentials, this is 1 here
!    nrmax:    the maximum case (b) rotational angular momenta
!    npar:     number of symmetry doublets included (npar=2 will ensure
!              both spin doublets)
!    isym:     if isym=+1, then the electronic symmetry is sigma-plus
!              if isym=-1, then the electronic symmetry is sigma-minus
!    igu:      if igu=+1, then the inversion symmetry is gerade
!              if igu=-1, then the inversion symmetry is ungerade
!    isa:      s/a symmetry index, if the molecule is homonuclear (ihomo=t)
!              then, if isa=+1 then only the s-levels are included in the
!              basis, if isa=-1, then only the a-levels are included
!    interp:   parameter to control spline interpolation in subroutine pot;
!              this should always be set equal to 1 for the first calculation
!  variable in common /cosys
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters.  note that the ordering
!             of the variable names in scod must correspond to the ordering
!             of the variable names in cosysi followed by the ordering of
!             variable names in cosysr followed by lammin, lammax, and mproj
! ------------------------------------------------------------------------
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: j, l, lc
logical existf
character*(*) fname
character*60 :: line,filnam,potfil
character*68 :: filnm1
character*1 dot
save potfil
#include "common/comdot.F90"
integer, pointer, save :: nterm, nrmax, npar, isym, igu, isa
real(8), pointer, save :: brot, gsr, drot, hrot
brot=>rspar(1); gsr=>rspar(2); drot=>rspar(3); hrot=>rspar(4)
nterm=>ispar(1); nrmax=>ispar(2); npar=>ispar(3); isym=>ispar(4); igu=>ispar(5); isa=>ispar(6)

irpot = 1
!     number and names of system dependent parameters
isicod = 6
isrcod = 4
nscode = isicod+isrcod
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
scod(1) = 'NTERM'
scod(2) = 'NRMAX'
scod(3) = 'NPAR'
scod(4) = 'ISYM'
scod(5) = 'IGU'
scod(6) = 'ISA'
!  then all the system dependent real variables
!  in the same order as in the common block /cosysr/
scod(7) = 'BROT'
scod(8) = 'GSR'
scod(9) = 'DROT'
scod(10) = 'HROT'
!  set default values for doublet-sigma scattering
potfil=' '
nterm = 1
if (iread .eq. 0) then
  nrmax = 3
  mproj(1) = 0
  lammin(1)= 1
  lammax(1) = -1
  isa=0
  isym=1
  npar=2
  niout=2
  indout(1)=-1
  indout(2)=1
endif
if(iread.eq.0) return
!  line 13
read (8, *,err=888) nrmax, npar, isym, igu, isa
!  line 14
read (8, *,err=888) brot, gsr, drot, hrot
!  line 16 name of file containing potential parameters
line=' '
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  return
endif
read (8, 285, end=286) line
potfil=line
285 format (a)
goto 286
! here if read error occurs
888 write(6,1000)
1000 format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
return
!
entry ptr2sg (fname,readpt)
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
if(nskip.eq.2) ihomo=.true.
nlam=0
nlam=nlam+(lammax(1)-lammin(1))/nskip+1
irpot=1
return
!
entry sav2sg (readpt)

ASSERT(nrmax .eq. ispar(2))
ASSERT(npar .eq. ispar(3))
ASSERT(isym .eq. ispar(4))
ASSERT(igu .eq. ispar(5))
ASSERT(isa .eq. ispar(6))

ASSERT(brot .eq. rspar(1))
ASSERT(gsr .eq. rspar(2))

!  save input parameters for doublet-sigma + atom scattering
!  line 13:
write (FUNIT_INP, 220) nrmax, npar, isym, igu, isa
220 format (5i4, t50, 'nrmax, npar, isym, igu, isa')
write (FUNIT_INP, 320) brot, gsr, drot, hrot
320 format (f12.6,3g12.4,t50,'brot, gsr, drot, hrot')
write (FUNIT_INP,330) potfil
330 format(a)
return
end
end module mod_hiba02_2sg
