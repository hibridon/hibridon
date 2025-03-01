#include "assert.h"
module mod_hiba09_stpln
logical :: twomol ! if .true. collision between symmetric top and linear
                  ! molecule, if .false. collision symmetric top-atom.
contains
! systpln (savstpln/ptrstpln) defines, save variables and reads          *
!                  potential for symmetric top and singlet sigma molecule*
!                  scattering                                            *
subroutine bastpln(bqs, jhold, ehold, ishold, nlevel, &
                  nlevop, ktemp, jtemp, &
                  ieps, isc1, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, &
                  jlpar, n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of a rigid symmetric top with a rigid diatomic molecule (basis type = 9)
!
!  authors:  millard alexander and peter vohralik
!  revision:  24-mar-1992 by c.rist
!
!  provided to pjd 30-aug-2011 by claire rist
!
!  included inversion splitting in symmetric top energies and changed
!  method of symmetric top level selection.  rotational channel basis
!  now chosen by j < jmax and ehold < emax (p.dagdigian)
!
!  revision:  24-aug-2012 by p.dagdigian
!  revision:  4-jun-2013 by q. ma (fix a bug in counting anisotropic
!     terms)
! --------------------------------------------------------------------
!  variables in call list:
!    bqs%jq:        on return contains combined rotational quantum numbers for each
!              channel.  in terms of the rotational quantum numbers of each
!              molecule we have:  j = 10 * j1 + j2
!    bqs%lq:        on return contains orbital angular momentum for each
!              channel
!    bqs%inq:       on return contains symmetry of inversion vibration for
!              each channel
!  Note that we have adopted the following convention for the symmetry
!  index "is" so that on return the symmetric top molecular levels can
!  be uniquely identified by the two labels "j" and "is":
!           is = +/- k
!       where k is the projection quantum number
!       the plus sign refers to the "plus" inversion levels and the
!       minus sign to the "minus" inversion levels
!       for k=0, where only eps=+1 states exist, the index "is" is equal
!       to zero.
!    jhold:    on return contains rotational quantum numbers for each
!              rotational level, irrespective of projection degeneracy and the
!              degeneracy associated with different values of j12
!              note that jhold = 10 * j1hold + j2hold
!    ehold:    on return contains energy in hartrees of each rotational
!              level
!    ishold:   on return contains symmetry index of each energetically
!              distinct level
!    nlevel:   on return contains number of energetically distinct
!              rotational levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              rotational levels used in channel basis which are open
!              asymptotically
!    ktemp:    scratch array used to create channel list
!    ieps:     on return contains epsilon label for each channel
!    isc1:     scratch array, not used
!    jtemp:    scratch array used to create channel list.
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
!              if ihomo.true., then symmetric top posseses interchange
!              symmetry (e.g. NH3),
!              so that only the ortho or para levels will be
!              included depending on the value of parameters iop
!              in common cosysi/ (see below)
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              jlpar (by definition jlpar=(-1)**(j2+l-jtot)*parity
!              jtot is added to be consistent with bastp routine but is
!              not needed here)
!              where parity is for the symmetric-top molecule
!              parity=isym*(-1)*k=-ieps*(-1)**(k+j1)
!
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to subroutines
!              propag and soutpt.  ntop is set in basis only if nu = numin
!              otherwise it is unchanged from the value supplied by the
!              calling program
!  variables in common block /cosysr/
!    isrcod:   number of real system dependent variables
!    brot:     A=B rotational constant for prolate top
!    crot:     C rotational constant for prolate top
!    delta:    inversion splitting
!    emax:     the maximum nh3 rotational energy (in cm-1) for a channel
!              to be included in the basis
!    drot:     rotational constant for linear molecule
!  variables in common block /cosysi/
!    nscod:    total number of variable names which are passed to HINPUT
!              nscod must equal isrcod + isicod + 3
!    isicod:   numbr of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    numpot:   the number of the potential used, this variable is passed
!              to the pot subroutine
!    ipotsy:   cylindrical symmetry of potential.  If ihomo = .true.,
!              only
!              terms with mu1 equal to an integral multiple of ipotsy
!              can be included in the potential.  Example:  for NH3,
!              ipotsy = 3
!    iop:      ortho/para label for molecular states. If ihomo=.true.
!              then only para states will be included if iop=1
!              and only ortho states if iop=-1
!    ninv:     number of inversion doublets included
!              if ninv = +1, only + inversion levels included
!              if ninv = -1, only - inversion levels included
!              if ninv = 2, both inversion levels included
!    jmax:     the maximum rotational angular momentum for the symmetric top
!    ipotsy2:  symmetry of potential. if linear molecule is homonuclear
!              then ipotsy=2 and only terms with lambda2  even can be
!              included in the potential,else ipotsy=1.
!    j2max:    the maximum rotational angular momentum for linear
!              molecule
!    j2min:    the minimum rotational angular momentum for linear
!              molecule
!  subroutines called:
!   vlmstpln:  returns molecule-molecule angular coupling coefficient for
!              particular choice of channel index
! --------------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coamat, only: ietmp ! ietmp(1)
use mod_conlam, only: nlam
use mod_cosysi, only: isicod, ispar
use mod_cosysr, only: isrcod, rspar
use mod_hibasutil, only: vlmstp, vlmstpln
use constants, only: econv, xmconv
use mod_par, only: iprint
use mod_parbas, only: ivrow, lammin, lammax, mproj, lam2, m2proj
use mod_ered, only: ered, rmu
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
type(bqs_type), intent(out) :: bqs
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst
character*40 fname
dimension jhold(1), ehold(1),&
          ishold(1)
dimension ieps(1), ktemp(1), jtemp(1), isc1(1)
data ione, itwo, ithree / 1, 2, 3 /
data one, zero, two / 1.0d0, 0.0d0, 2.d0 /
integer, pointer :: nterm, numpot, ipotsy, iop, ninv, jmax, ipotsy2, j2max, j2min
real(8), pointer :: brot, crot, delta, emax, drot
nterm=>ispar(1); numpot=>ispar(2); ipotsy=>ispar(3); iop=> ispar(4); ninv=>ispar(5)
jmax=>ispar(6); ipotsy2=>ispar(7); j2max=>ispar(8); j2min=>ispar(9)
brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); emax=>rspar(4); drot=>rspar(5)

twomol=.true.
!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (9, 5)
5   format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***' )
  if (bastst) then
    return
  else
    stop
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
    stop
  end if
end if
if (csflag) then
  write (6, 15)
  write (9, 15)
15   format (' *** CSFLAG=.TRUE.; NOT IMPLEMENTED IN BASIS;', &
          ' ABORT ***')
  if (bastst) then
    return
  else
    stop
  end if
end if
do 25  i = 1, nterm
  if (ihomo .and. mod(mproj(i), ipotsy) .ne. 0 ) then
    write (9, 20) i, mproj(i), ipotsy
    write (6, 20) i, mproj(i), ipotsy
20     format(/' *** MPROJ(',i2,')=',i2, &
     ' .NE. INTEGER MULTIPLE OF', i2, ' FOR IHOMO=T; ABORT ***')
    stop
  end if
  if (twomol) then
    if (mod(lam2(i), ipotsy2) .ne. 0 ) then
      write (9, 21) i, lam2(i), ipotsy2
      write (6, 21) i, lam2(i), ipotsy2
21       format(/' *** LAM2(',i2,')=',i2, &
      ' .NE. INTEGER MULTIPLE OF', i2, &
         ' FOR HOMONUCLEAR MOLECULE ABORT ***')
      stop
    endif
  endif
25 continue

nsum = 0
do 35  i = 1, nterm
  if (mproj(i) .gt. lammin(i) ) then
    write (6, 30) i, mproj(i), lammin(i)
    write (9, 30) i, mproj(i), lammin(i)
30     format (' *** MPROJ=',i2,' > LAMMIN=',i2, &
            ' FOR TERM',i2,'; ABORT ***')
    stop
  end if
  if (twomol) then
    if (m2proj(i).gt.lammin(i)) then
      write (6, 31) i, m2proj(i), lammin(i)
      write (9, 31) i, m2proj(i), lammin(i)
31       format (' *** M2PROJ=',i2,' > LAMMIN=',i2, &
            ' FOR TERM',i2,'; ABORT ***')
      stop
    endif
    if (m2proj(i).gt.lam2(i) ) then
      write (6, 32) i, m2proj(i), lam2(i)
      write (9, 32) i, m2proj(i), lam2(i)
32       format (' *** M2PROJ=',i2,' > LAM2=',i2, &
            ' FOR TERM',i2,'; ABORT ***')
      stop
    end if
  end if
  if(lammax(i) .lt. 0) goto 35
  nsum = nsum + lammax(i) - lammin(i) + 1
35 continue
if (nsum .ne. nlam) then
  write (6, 45) nsum, nlam
  write (9, 45) nsum, nlam
45   format (' *** NLAM IN BASIS=', i3,' .NE. NLAM FROM SYSDAT=', &
           i3, '; ABORT ***')
  stop
end if
nlam = nsum
if (bastst) then
  write (9, 50) nlam
  write (6, 50) nlam
50   format (' *** TOTAL NUMBER OF ANISOTROPIC TERMS =',i3)
end if
if (clist) then
  if (flagsu) then
    write(6,51)
    write(9,51)
51     format(/,' **  SYMMETRIC TOP - UNCORRUGATED SURFACE **', &
             ' NON IMPLEMENTED HERE, ABORT**')
    stop
  elseif (csflag) then
    write(6,52)
    write(9,52)
52     format(/,' **  SYMMETRIC TOP - COUPLED STATES **', &
             ' NON IMPLEMENTED HERE, ABORT**')
    stop
  elseif (twomol) then
    if (ihomo) then
      if (bastst) &
      write (6,80) rmu * xmconv, brot, crot, delta, drot, ipotsy, &
                     iop, ipotsy2, ered * econv, jtot, jlpar
      write (9,80) rmu * xmconv, brot, crot, delta, drot, ipotsy, &
                     iop, ipotsy2, ered * econv, jtot, jlpar
80       format(/,' **  CC SYMMETRIC TOP-LINEAR MOLECULE **', &
          /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3, &
            '  DELTA =',f7.3,'  DROT=',f7.3,/, &
            '     POT-SYM=', i2,'  O/P=',i2, ' HOMOLIN=', i2, &
             '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
    else
      if (bastst) &
      write (6,85) rmu * xmconv, brot, crot, delta, drot, &
                    ipotsy2, ered * econv, jtot, jlpar
      write (9,85) rmu * xmconv, brot, crot, delta, drot, &
                    ipotsy2, ered * econv, jtot, jlpar
85       format(/,' **  CC SYMMETRIC TOP **', &
          /,'     RMU=', f9.4,'  A=BROT=', f7.3, '  CROT=',f7.3, &
             '  DELTA = ',f7.3,'  DROT=',f7.3,/,'  HOMOLIN=', i2, &
             '     E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
    end if
  else
     write(6,86)
     write(9,86)
86      format(/,' **  SYMMETRIC TOP - ATOM **', &
            ' NON IMPLEMENTED HERE, ABORT**')
     stop
  end if
end if
if (.not. flagsu) write (9,90) rcut
90   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
!
!  first set up list of all j,k j2 states included in basis
nlist = 0
do 105  ki = 0, jmax
  do 105  ji = ki, jmax
    numeps = 2
    if (ki.eq.0) numeps = 1
    do 105 iep = 1, numeps
      if (ihomo) then
!  consider nuclear permutation symmetry for symmetric top
!  molecule:
!  ortho levels (iop = -1) have ki a multiple of 3
!  all other levels are para (iop = 1)
        if (iop.eq.-1 .and. mod(ki,ipotsy).ne.0) goto 105
        if (iop.eq.1 .and. mod(ki,ipotsy).eq.0) goto 105
      end if
!  check to make sure symmetric top rotational energy is less than emax
      roteng = brot * ji * (ji + 1) + (crot - brot) * ki * ki
      if (roteng .gt. emax) go to 105
!  add inversion splitting to - inversion levels
      iepsil = 1 - 2*(iep - 1)
      isym = -iepsil * (-1) ** ji
      if (isym .eq. -1) roteng = roteng + delta
!  include level only if matches ninv specification
      if (ninv.eq.2 .or. (ninv.eq.1 .and. isym.eq.1) .or. &
          (ninv.eq.-1 .and. isym.eq.-1)) then
        do ji2 = j2min, j2max, ipotsy2
          nlist = nlist+1
          jtemp(nlist) = 10*ji + ji2
          ktemp(nlist) = ki
          ietmp(nlist) = iepsil
          ehold(nlist)= (roteng + drot * ji2 * (ji2 + 1)) &
            / econv
        end do
      end if
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
  write (6, 125)
125   format (/,10x, &
      'SORTED LEVEL LIST',/,'   N   J   K  EPS INV J2   ', &
      'EINT(CM-1)')
  do  127 i = 1, nlist
    ji2 = mod(jtemp(i),10)
    ji1 = jtemp(i)/10
    isym = -ietmp(i) * (-1) ** ji1
    ecm = ehold(i) * econv
    write (6, 126) i, ji1, ktemp(i), ietmp(i), isym, &
        ji2, ecm
126     format (6i4,f10.3)
127   continue
end if
!  now set up channel and level list for scattering calculation
n = 0
nlevel = 0
call bqs%init(nmax)
do 170  njk = 1, nlist
  ki = ktemp(njk)
  ji = jtemp(njk)/10
  ji2 = mod(jtemp(njk),10)
  iepsil = ietmp(njk)
  nlevel = nlevel + 1
  jhold(nlevel) = 10*ji + ji2
!  ishold set equal to k times inversion symmetry
  ishold(nlevel) = -ietmp(njk) * (-1.d0)**ji * ki
!
!  construct channel basis for cc calculations.  first calculate range of orbital
!  angular momentum quantum numbers allowed for this state
!
!  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem. Phys.
!  73, 2740 (1980) ]
  isym = -ietmp(njk) * (-1) ** ji
  ipar = isym * (-1) ** ki
  j12max = ji + ji2
  j12min = iabs(ji - ji2)
  do 156 ji12 = j12min, j12max
    lmax = jtot + ji12
    lmin = iabs(jtot - ji12)
    do 155  li = lmin, lmax
      lpar = (-1) ** (li - jtot)
      j2par = (-1)**ji2
!  check to see if this channel has the desired parity, if so, include it
      if (ipar * lpar * j2par .eq. jlpar) then
        n = n + 1
        if (n .gt. nmax) then
          write (6, 150) n, nmax
          write (9, 150) n, nmax
150           format(/' *** NCHANNELS=', i5, &
              ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
          stop
        end if
        bqs%jq(n) = 10*ji + ji2
        bqs%inq(n) = ishold(nlevel)
        ieps(n) = iepsil
        bqs%j12(n) = ji12
        eint(n) = ehold(nlevel)
        bqs%lq(n) = li
        cent(n) = li * (li + 1)
        bqs%length = n
      end if
155     continue
156   continue
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
if (.not. flagsu .and. rcut .gt. 0.d0) then
  emin = 1.e+7
  do 290  i = 1, n
    if (eint(i) .le. ered) then
! here if channel is open asymptotically
      if ( jtot * (jtot + 1) / (two * rmu * rcut * rcut) &
          .gt. (ered - eint(i)) ) then
!  here if channel is open asymptotically but closed at r = rcut
        if (eint(i) .lt. emin) emin = eint(i)
!  emin now contains the lowest channel energy for which this condition
!  is met
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
        bqs%j12(nn)=bqs%j12(i)
        cent(nn) = cent(i)
        bqs%lq(nn) = bqs%lq(i)
      end if
300     continue
!  reset number of channels
    n = nn
    bqs%length = nn
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
    stop
  endif
end if
!  now list channels if requested
if (clist) then
  if (.not.csflag) then
    if (bastst) write (6, 305)
    write (9,305)
305     format ( &
     /'   N   J  EPS INV  K  J2  J12  L    EINT(CM-1)')
  else
    if (bastst) write (6, 310) nu
    write (9,310) nu
310     format ( &
     /'   N   J  EPS INV  K   L    EINT(CM-1) ** NU = ',i2)
  end if
  do 330  i = 1, n
    ki = iabs(bqs%inq(i))
    if (bqs%inq(i) .ne. 0) inv = bqs%inq(i) / ki
    ecm = eint(i) * econv
    if (twomol) then
      ji = bqs%jq(i)/10
      ji2 = mod(bqs%jq(i),10)
    else
      ji = bqs%jq(i)
      ji2 = 0
    endif
    if (bastst) &
      write (6, 320) i, ji, ieps(i), inv, ki, ji2, bqs%j12(i), &
                     bqs%lq(i), ecm
      write (9, 320) i, ji, ieps(i), inv, ki, ji2, bqs%j12(i), &
                     bqs%lq(i), ecm
320       format (8i4, f10.3)
330   continue
end if

!  now calculate coupling matrix elements
!  i counts v2 elements
!  inum counts v2 elements for given lambda
!  ilam counts number of v2 matrices
i = 0
if (bastst.and. iprint.ge. 2) then
  write (6, 340)
!        write (76,340)
  write (9, 340)
340   format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
end if
lamsum = 0
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
    if (twomol) then
      lb2 = lam2(iterm)
      mu2 = m2proj( iterm)
    endif
    inum = 0
    do 355  icol = 1, n
      do 350  irow = icol, n
        lr = bqs%lq(irow)
        lc = bqs%lq(icol)
        kc = iabs (bqs%inq(icol))
        kr = iabs (bqs%inq(irow))

        if (csflag) lrow = nu
        if (twomol) then
          jc = bqs%jq(icol)/10
          jr = bqs%jq(irow)/10
          j2c = mod(bqs%jq(icol), 10)
          j2r = mod(bqs%jq(irow), 10)
          kc = iabs (bqs%inq(icol))
          kr = iabs (bqs%inq(irow))
          j12r = bqs%j12(irow)
          j12c = bqs%j12(icol)
! this call to vlmstp is a test for sym-atom limit
!                call vlmstp (jr, lr, jc, lc, jtot,
!     :                     kr, kc, lb, mu, ieps(irow),
!     :                     ieps(icol), vee, csflag)

          call vlmstpln (jr, lr, jc, lc, j2r, j2c, &
                     j12r, j12c, jtot, kr, kc, lb, mu, lb2, mu2, &
                     ieps(irow), ieps(icol), vee, csflag)
        else
          call vlmstp (bqs%jq(irow), lr, bqs%jq(icol), lc, jtot, &
                     kr, kc, lb, mu, ieps(irow), &
                     ieps(icol), vee, csflag)
        endif
        if (vee .ne. zero) then
          i = i + 1
          inum = inum + 1
          call ancouma%set_element(irow=irow, icol=icol, vee=vee)
          if (bastst.and. iprint.ge.2) then
            write (6, 345) ilam, lb, icol, irow, i, vee
            write (9, 345) ilam, lb, icol, irow, i, vee
345               format (i4, 2i7, 2i6, i6, f17.10)
          end if
        end if
350       continue
355     continue
    if (bastst) then
      write (6, 370) ilam, lb, mu, lb2, mu2, ancouma%get_num_nonzero_elements()
      write (9, 370) ilam, lb, mu, lb2, mu2, ancouma%get_num_nonzero_elements()
370       format ('ILAM=',i3,' LAM=',i3,' MU=',i3, &
        ' LAM2=',i3,' MU2=',i3,' LAMNUM(ILAM) = ',i6)
    end if
390   continue
400 continue
if (clist .and. bastst) then
  write (6, 460) v2%get_num_nonzero_elements()
!        write(76, 460) lamsum
  write (9, 460) v2%get_num_nonzero_elements()
460   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', &
           i8)
end if
return
end

!  -----------------------------------------------------------------------
subroutine systpln (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for symmetric top
!   + linear molecule scattering
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  revision: 9-mar-1992 by c.rist
!  add delta (inversion splitting) parameter and
!  replace kmax,jmax0, ..., with jmax, emax parameters (p.j.dagdigian)
!  revision:  17-nov-2011 by p.j.dagdigian
!  current revision:  19-apr-2012 by q.ma (correct savestpln)
!
! NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
!  -----------------------------------------------------------------------
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     A=B rotational constant for prolate top
!    crot:     C rotational constant for prolate top
!    delta:    inversion splitting
!    emax:     the maximum symmetric top rotational energy (in cm-1)
!              for a channel to be included in the basis
!    drot:     rotational constant for linear molecule
!  variables in common block /cosysi/
!    nscod:    total number of variable names which are passed to HINPUT
!              nscod must equal isrcod + isicod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    numpot:   the number of the potential used, this variable is passed
!              to the pot subroutine
!    ipotsy:   cylindrical symmetry of potential.  Only terms with mu equal to
!              an integral multiple of ipotsy can be included in the potential
!              Example:  for NH3, ipotsy = 3
!    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
!              para states will be included if iop=1 and only ortho states if
!              iop=-1
!    ninv:     number of inversion doublets included
!              if ninv = +1, only + inversion levels included
!              if ninv = -1, only - inversion levels included
!              if ninv = 2, both inversion levels included
!    jmax:     the maximum rotational angular momentum for the symmetric top
!    ipotsy2:  symmetry of potential. if linear molecule is homonuclear
!              then ipotsy=2 and only terms with lambda2  even can be
!              included in the potential,else ipotsy=1.
!    j2max:    the maximum rotational angular momentum for linear
!              molecule
!    j2min:    the minimum rotational angular momentum for linear
!              molecule
!  variable in common /cosys
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters.  Note that the ordering
!             of the variable names in scod must correspond to the ordering
!             of the variable names in cosysi followed by the ordering of
!             variable names in cosysr followed by LAMMIN, LAMMAX, MPROJ,
!             LAM2 and MPROJ2.
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, rspar
use funit, only: FUNIT_INP
use mod_parbas, only: ivrow, lammin, lammax, mproj, lam2, m2proj
use mod_hiutil, only: gennam, get_token
implicit none
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
logical existf
integer icod, ircod
integer i, j, k, l, lc
character*1 dot
character*(*) fname
character*60 line, filnam, potfil, filnm1
parameter (icod=9, ircod=5)
#include "common/comdot.F90"
save potfil

integer, pointer, save :: nterm, numpot, ipotsy, iop, ninv, jmax, ipotsy2, j2max, j2min
real(8), pointer, save :: brot, crot, delta, emax, drot
nterm=>ispar(1); numpot=>ispar(2); ipotsy=>ispar(3); iop=> ispar(4); ninv=>ispar(5)
jmax=>ispar(6); ipotsy2=>ispar(7); j2max=>ispar(8); j2min=>ispar(9)
brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); emax=>rspar(4); drot=>rspar(5)


twomol = .true.
!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
!  NTERM must be the first variable
!  followed by the system dependent real variables
!  in the same order as in the common block /cosysr/
!  then the three variable names LAMMIN, LAMMAX, MPROJ, in that order
scod(1)='NTERM'
scod(2)='NUMPOT'
scod(3)='IPOTSY'
scod(4)='IOP'
scod(5)='NINV'
scod(6)='JMAX'
scod(7)='IPOTSY2'
scod(8)='J2MAX'
scod(9)='J2MIN'
scod(10)='BROT'
scod(11)='CROT'
scod(12)='DELTA'
scod(13)='EMAX'
scod(14)='DROT'
scod(15)='LAMMIN'
scod(16)='LAMMAX'
scod(17)='MPROJ'
scod(18)='LAM2'
scod(19)='M2PROJ'
nscode = icod + ircod + 5
isicod = icod
isrcod = ircod
irpot = 1
!  set default values for symmetric top scattering
nterm=12
do j=1, 4
  do i=1, 3
    k=3*(j-1)+i
    mproj(k) = 3*(i-1)
    if(j.le.2) m2proj(k) = 0
    if(j.eq.3) m2proj(k) = 1
    if(j.eq.4) m2proj(k) = 2
    if(j.eq.1) lam2(k) = 0
    if(j.ge.2) lam2(k) = 2
    lammin(k) = max(mproj(k), m2proj(k))
    lammax(k) = 6
  enddo
enddo
lammin(1)=1
if (iread .eq. 0) then

  jmax = 5
  ipotsy = 3
  iop = 1
  ninv = 2
  jmax =5
  ipotsy2=2
  j2min=0
  j2max=2
  niout=5
  indout(1)=0
  indout(2)=-3
  indout(3)=+3
  indout(4)=-6
  indout(5)=6
endif
potfil=' '
if (iread .eq. 0) return
!  line 18
read (8, *, err=80) ipotsy, iop, ninv, ipotsy2
!  line 19
!  line 20
read (8, *, err=80) jmax
!  line 21
read (8, *, err=80) j2min, j2max
!  line 22
read (8, *, err=80) brot, crot, delta, emax
!  line 23
read(8, *, err=80) drot
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
!************************************************************
entry ptrstpln (fname,readpt)
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
entry savstpln (readpt)
!  save input parameters for symmetric top + linear molecule scattering
!  the order of the write statements should be identical to the read statement
!  above. for consistency with the data file written by gendat, format
!  statements should reserve the first 30 spaces for data, spaces 31-33 should
!  be left blank, and the names of the variables should be printed in spaces
!  34-80
!  line 18:
write (FUNIT_INP, 220) ipotsy, iop, ninv, ipotsy2
220 format (4i4, 14x,'   ipotsy, iop, ninv, ipotsy2')
!  line 20
write (FUNIT_INP, 230) jmax
230 format (i4, 26x, '   jmax')
!  line 21
write (FUNIT_INP,231) j2min, j2max
231 format (2i4, 22x,'   j2min, j2max')
write (FUNIT_INP, 250) brot, crot, delta, emax
250 format(3f8.4, f8.2, ' brot, crot, delta, emax' )
write (FUNIT_INP, 251) drot
251 format(f12.6, 18x,'   drot')
write (FUNIT_INP, 60) potfil
return
end
end module mod_hiba09_stpln
