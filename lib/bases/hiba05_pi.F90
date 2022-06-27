! sypi (savpi/ptrpi) defines, save variables and reads                   *
!                  potential for general pi scattering                   *
! ----------------------------------------------------------------------
subroutine bapi  (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, c0, c1, c2, cf, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop)
! ----------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of a 1pi in a case (a) basis, or of a 2pi or 3pi molecule in an
!  intermediate coupling basis, with a structureless atom or with
!  an uncorrugated surface, including possible stark mixing of
!  the lambda-doublets by an external field
!  author:  didier lemoine
!  current revision date:  25-feb-2004 by mha
! ----------------------------------------------------------------------
!  variables call list:
!    j:        on return contains rotational quantum numbers for each
!              channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains symmetry index of each channel
!  note that we have adopted the following convention for the symmetry
!  index "is" so that on return the multiplet pi molecular levels can
!  be uniquely identified by the two labels "j" and "is":
!           for Fi levels is = i*eps
!  where i = 1,2,3 and eps = +/- 1 is the "true" case (a) symmetry index
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
!    ci:       on return contains omega=i(+1/2 for 2pi) component in the
!              states mixed by the j.s term in the diatomic hamiltonian
!              where i = 1,2,3
!    cf:       on return contains f component in the states mixed
!              by an external stark field
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!              note that this is ignored in molecule-surface collisions!
!    jtot:     total angular momentum
!              in cc calculation jtot is the total angular momentum
!              in cs calculation jtot is the l-bar quantum number
!    flaghf:   if .true., then system has half-integer spin
!              if .false., then system has integer spin
!    flagsu:   if .true., then molecule-surface collisions
!    csflag:   if .true., then coupled-states calculation
!              if .false., then close-coupled calculation
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true., then homonuclear molecule
!              only the s or a levels will be included depending on the
!              value of the parameter isa in common cosysi/ (see below)
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0
!              by calling program
!    jlpar:    total parity of included channels in cc calculation
!              only those channels are included for which
!                  (-1)**(parity+l-jtot)=jlpar
!              where parity designates the parity of the molecular state
!              (by definition parity=eps*(-1)**(J-S))
!              in cs calculation jlpar is set = 1 in calling program
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to
!              subroutines propag and soutpt. ntop is set in basis
!              only if nu = numin. otherwise it is unchanged from
!              the value supplied by the calling program
!    note!!   if flaghf = .true., then the true values of the rotational
!    quantum numbers, the total angular momentum, and the coupled-states
!    projection index are equal to the values stored in j, jtot, and nu
!    plus 1/2
!  variables in common block /cosysi/
!    nscode:   total number of system dependent parameters
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    jmax:     the maximum rotational angular momenta (-1/2 for 2pi) for
!              each channel in each spin-orbit manifold with convention
!              omega .le. j .le. jmax(+1/2 for 2pi)
!    igu:      permutation inversion symmetry of electronic state
!              igu=1 for gerade states, igu=-1 for ungerade states
!              for heteronuclear molecules igu is not used
!    isa:      s/a label for molecular states. if ihomo=.true.,
!              then only s states will be included if isa=1 and
!              only a states if isa=-1
!    npar:     number of symmetry doublets included
!              (npar=2 will ensure both lambda doublets)
!    imult:    spin multiplicity of pi state (imult = 2*S+1)
!    nman:     number of spin-orbit manifolds included
!              for a 2pi state nman=1 selects the
!              omega=1/2 manifold in a case (a) basis
!              otherwise, nman is set = imult
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     rotational constant in cm-1
!    aso:      spin-orbit constant in cm-1
!    o, p, q:  lambda-doubling constants cm-1
!    dmom:     dipole moment of molecule in Debye
!    efield:   stark field in kV / cm
!    eint:     array containing channel energies (in hartree)
!  variables in common block /coered/
!    ered:     collision energy in atomic units (hartrees)
!    rmu:      collision reduced mass in atomic units
!              (mass of electron = 1)
!  variable in common block /coconv/
!   econv:     conversion factor from cm-1 to hartrees
!   xmconv:    conversion factor from amu to atomic units
!  subroutines called:
!   rs:        diagonalizes 3*3 case (a) hamiltonian matrix
!              for a 3pi state
!   vlmpi:     returns angular coupling coefficient for
!              particular choice of channel index
! ----------------------------------------------------------------------
use mod_cov2, only: nv2max, junkv => ndummy, v2
use mod_coiv2, only: iv2
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam, nlammx, lamnum
use constants, only: econv, xmconv, ang2c
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use mod_par, only: iprint, rendai=>scat_rendai
#include "common/parbasl.F90"
implicit double precision (a-h,o-z)
logical flaghf, csflag, clist, flagsu, ihomo, bastst
character*80 string
character*27 case
character*2 chf
#include "common/parbas.F90"
common /coered/ ered, rmu
dimension j(1), is(1), l(1), jhold(1), ishold(1), ieps(2)
dimension c0(1), c1(1), c2(1), cf(1), ehold(1)
dimension e(3,3), eig(3), sc1(3), sc2(3), vec(3,3), vii(0:2)
!  cvtown: conversion factor from coulomb.volts to cm-1 (wavenumbers)
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
! work vector for dsyev
dimension work(9)
#endif
data cvtown / 0.0167917d0/, ieps / -1, 1/
integer, pointer :: nterm, jmax, igu, isa, npar, imult, nman
real(8), pointer :: brot, aso, o, p, q, dmom, efield
nterm=>ispar(1); jmax=>ispar(2); igu=>ispar(3); isa=>ispar(4); npar=>ispar(5); imult=>ispar(6); nman=>ispar(7)
brot=>rspar(1); aso=>rspar(2); o=>rspar(3); p=>rspar(4); q=>rspar(5); dmom=>rspar(6); efield=>rspar(7)

zero = 0.d0
one = 1.d0
two = 2.d0
four = 4.d0
half = 0.5d0
quart = 0.25d0
!  xhf is the half-integer part of most quantum numbers
!  ispin is the integer part of the spin (= 1 for 3pi and 0 otherwise)
xhf = half * mod(imult-1,2)
ispin = (imult - 1) / 2
!  check for consistency in the values of flaghf, flagsu and csflag
if (.not. flaghf .and. xhf .eq. half &
    .or. flaghf .and. xhf .eq. zero) then
  write (6, 10)
  write (9, 10)
10   format (/' *** MISMATCH BETWEEN FLAGHF AND IMULT; ABORT ***')
  stop
end if
if (flagsu .and. .not. csflag) then
  write (6, 20)
  write (9, 20)
20   format (/ &
  ' *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
  stop
end if
!  check for consistency in the values of nterm, lammin, lammax, mproj
if (nterm .gt. 2) then
  write (6, 30) nterm
  write (9, 30) nterm
30   format (/' *** NTERM=', i2, &
          ' > 2 FOR PI STATE BASIS; ABORT ***')
  stop
end if
do 70  i = 1, nterm
  if (ihomo .and. lammin(i) .lt. 2) then
    write (9, 40) i, lammin(i)
    write (6, 40) i, lammin(i)
40     format (/ &
    ' *** LAMMIN(', i2, ')=', i2, ' < 2 FOR IHOMO=T; ABORT ***')
    stop
  end if
  if (mproj(i) .gt. lammin(i)) then
    write (6, 50) i, mproj(i), lammin(i)
    write (9, 50) i, mproj(i), lammin(i)
50     format (/' *** MPROJ(', i2, ')=', i2, &
            ' > LAMMIN=', i2, '; ABORT ***')
    stop
  end if
  if (ihomo) then
    if (mod(lammax(i)-lammin(i),2) .ne. 0) then
      write (6, 60) i, lammax(i)-lammin(i)
      write (9, 60) i, lammax(i)-lammin(i)
60       format (/' *** IHOMO=T BUT (LAMMAX-LAMMIN)(', i2, ')=', &
              i2, ' NOT EVEN; ABORT ***')
      stop
    end if
  end if
70 continue
!  determine number of anisotropic terms
istep = 1
if (ihomo) istep = 2
nlam0 = (lammax(1) - lammin(1)) / istep + 1
nlam = nlam0
if (.not. (imult .eq. 2 .and. nman .eq. 1)) nman = imult
if (nman .ne. imult) nterm = 1
if (nterm .eq. 2) &
    nlam = nlam0 + (lammax(2) - lammin(2)) / istep + 1
if (nlammx .lt. nlam) then
  write (6, 80) nlam, nlammx
  write (9, 80) nlam, nlammx
80   format (/' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', &
          i2, ' > NLAMMX=', i2,'; ABORT')
  stop
end if
if (clist) then
  if (bastst) write (6, 90) nlam
  write (9, 90) nlam
90   format (/' ** TOTAL NUMBER OF ANISTROPIC TERMS =', i3/)
  chf = '  '
  if (flaghf) chf = '.5'
  case = 'STATE INTERMEDIATE COUPLING'
  if (nman .ne. imult) case = 'OMEGA=1/2 MANIFOLD CASE (A)'
  i = 27
  if (imult .eq. 1) i = 5
  if (flagsu) then
    if (bastst) write (6, 100) imult, case(1:i), rmu*xmconv, brot, &
                               aso, ered*econv, jtot, nu, chf
    write (9, 100) imult, case(1:i), rmu*xmconv, brot, &
                   aso, ered*econv, jtot, nu, chf
100     format (' **', i2, 'PI ', a, ' + FLAT SURFACE **'/, &
            '     RMU=', f9.4, '  BROT=', f7.3, '  A-SO=', f7.2/, &
            '     E=', f8.2, '     LBAR=', i4, '     M=', i3, a)
  else if (csflag) then
    if (bastst) write (6, 110) imult, case(1:i), rmu*xmconv, brot, &
                               aso, ered*econv, jtot, nu, chf
    write (9, 110) imult, case(1:i), rmu*xmconv, brot, &
                   aso, ered*econv, jtot, nu, chf
110     format (' ** CS', i2, 'PI ', a,  ' **'/, &
            '     RMU=', f9.4, '  BROT=', f7.3, '  A-SO=', f7.2/, &
            '     E=', f8.2, '     LBAR=', i4, '     NU=', i3, a)
  else
    if (bastst) write (6, 120) imult, case(1:i), rmu*xmconv, brot, &
                               aso, ered*econv, jlpar, jtot, chf
    write (9, 120) imult, case(1:i), rmu*xmconv, brot, &
                   aso, ered*econv, jlpar, jtot, chf
120     format (' ** CC', i2, 'PI ', a,  ' **'/, &
            '     RMU=', f9.4, '  BROT=', f7.3, '  A-SO=', f7.2/, &
            '     E=', f8.2, '    JLPAR=', i3, 5x, 'JTOT=', i4, a)
  end if
  if (ihomo) then
    if (bastst) write (6, 130) igu, isa
    write (9, 130) igu, isa
130     format ('     g/u=', i2, '  s/a=', i2)
  end if
  if (.not. flagsu) write (9, 140) rcut
140   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
end if
!  calculation of the mixing of the rotational basis by a stark field is
!  only valid for a molecule-surface calculation because the presence of
!  the electric field breaks the isotropy of space. one therefore cannot
!  use a total-J representation
if (dmom .eq. zero .or. npar .eq. 1 .or. nu+xhf .eq. zero) &
    efield = zero
if (efield .ne. zero .and. .not. flagsu) then
  write (6, 150)
  write (9, 150)
150   format (/' *** FLAGSU = .FALSE. FOR STARK MIXING; ABORT ***')
  stop
end if
!  first set up list of all case (a) levels for first spin-orbit state
n = 0
iommin = 0
if (imult .eq. 1) iommin = 1
do 160 ji = iommin, jmax
do 160 ip = 1, npar
!  for 1pi and 2pi states eps = +1 for e levels, eps = -1 for f levels
!  while for 3pi states eps = -1 for e levels, eps = +1 for f levels
!  for homonuclear molecules in gerade electronic states the e levels
!  are s for even j, a for odd j, while the f levels are a for even j,
!  s for odd j, where j = J-1/2 for 2pi. this reverses for ungerade
!  states. see alexander and corey, j. chem. phys. 84, 100 (1986)
  if (ihomo .and. igu*ieps(ip)*(-1)**(ji-ispin) .ne. isa) &
      go to 160
  n = n + 1
  j(n) = ji
  is(n) = ieps(ip)
  if (imult .eq. 1) eint(n) = (brot + half * is(n) * q) &
                              * ji * (ji + 1)
  if (imult .eq. 2) then
    x = ji + one
    eint(n) = - half * (nman - 1) * aso + brot * x * x + &
       half * (one - is(n) * x) * (p + q * (one - is(n) * x))
  end if
160 continue
!  n now contains the number of levels for the first spin-orbit state
!  now add the levels from the second spin-orbit state if any,
!  determine arrays of mixing coefficients (c0, c1 and c2),
!  and calculate internal energies of the mixed states
!  iso = +/- 1 designates a regular/inverted 2pi or 3pi state
iso = isign(1,nint(aso))
iso01 = (1 - iso) / 2
nn = n
do 170 i = 1, n
  if (j(i) .eq. 0 .or. nman .eq. 1) then
!  here for a pure case (a) state
!  for 3pi (2pi) j = 0 (1/2) can exist only for omega = 0 (1/2)
!  and hence is always a pure case (a) state
    if (nman .eq. 2) is(i) = (1 + iso01) * is(i)
    if (imult .eq. 3) eint(i) = - aso + two * brot + &
                                is(i) * (o + p + q)
    c0(i) = one
    c1(i) = zero
    c2(i) = zero
  else
!  here for j .ge. 1
!  the index nn points to the nominal omega=1(+1/2 for 2pi) levels
    nn = nn + 1
    j(nn) = j(i)
    if (imult .eq. 2) then
!  here for 2pi, the matrix elements are given by a.j. kotlar,
!  r.w. field and j.i. steinfeld, j. mol. spectr. 80, 86 (1980)
!  where x = J+1/2
!  save the case (a) energies temporarily in the variables e12 and e32,
!  and evaluate 1/2-3/2 coupling matrix element (h1232) and mixing
!  angle (angle) due to the j.s term in the hamiltonian
      e12 = eint(i)
      x = j(i) + one
      e32 = half * aso + (brot + half * q) * (x * x - one) - brot
      h1232 = - (brot + quart * p + half * (one - is(i) * x) * q) &
              * sqrt(x*x-one)
      angle = half * atan(two*h1232/(e12-e32))
!  c0 and c1 indicate the amount of omega = 1/2 and 3/2
!  components, respectively, in each mixed level
      c0(i) = cos(angle)
      c1(i) = sin(angle)
      c2(i) = zero
      c0(nn) = - c1(i)
      c1(nn) = c0(i)
      c2(nn) = zero
!  now store the energies of the mixed states in the array eint
      eint(i) = e12 * c0(i)**2 + e32 * c1(i)**2 &
                + two * c0(i) * c1(i) * h1232
      eint(nn) = e12 * c0(nn)**2 + e32 * c1(nn)**2 &
                 + two * c0(nn) * c1(nn) * h1232
!  assign symmetry index: i*eps for Fi levels
      is(nn) = (2 - iso01) * is(i)
      is(i) = (1 + iso01) * is(i)
    else if (imult .eq. 3) then
      is(nn) = 2 * is(i)
!  here for 3pi, the matrix elements are given by
!  j.m. brown and a.j. merer, j. mol. spectr. 74, 488 (1979), or by c.r.
!  brazier, r.s. ram and p.f. bernath, j. mol. spectr. 120, 381 (1986)
!  save the case (a) energies temporarily in the array e
      x = float(j(i)*(j(i)+1))
      e(1,1) = - aso + brot * (x + two) + is(i) * (o + p + q)
      e(1,2) = - (brot + half * is(i) * (p + two * q)) &
               * sqrt(two*x)
      e(1,3) = half * is(i) * q * sqrt(x*(x-two))
      e(2,1) = e(1,2)
      e(2,2) = brot * (x + two) + half * is(i) * q * x
      e(2,3) = - brot * sqrt(two*(x-two))
      e(3,1) = e(1,3)
      e(3,2) = e(2,3)
      e(3,3) = aso + brot * (x - two)
!  c0, c1 and c2 indicate the amount of omega = 0, 1 and 2
!  components, respectively, in each mixed level
!  store the energies of the mixed states in the array eint
      if (j(i) .eq. 1) then
!  only the omega=0 and omega=1 states are mixed for j = 1
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
        lwork=9
        call dsyev('V','L',2,e,3,eig,work,lwork,ierr)
        do iv=1,2
           do jv=1,2
              vec(iv,jv)=e(iv,jv)
           enddo
        enddo
#endif
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
        call rs (3, 2, e, eig, 1, vec, sc1, sc2, ierr)
#endif
        eint(i) = eig(1+iso01)
        c0(i) = vec(1,1+iso01)
        c1(i) = vec(2,1+iso01)
        c2(i) = zero
        eint(nn) = eig(2-iso01)
        c0(nn) = vec(1,2-iso01)
        c1(nn) = vec(2,2-iso01)
        c2(nn) = zero
      else if (j(i) .ge. 2) then
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
        lwork=9
        call dsyev('V','L',3,e,3,eig,work,lwork,ierr)
        call dcopy(9,e,1,vec,1)
#endif
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
        call rs (3, 3, e, eig, 1, vec, sc1, sc2, ierr)
#endif
        eint(i) = eig(2-iso)
        c0(i) = vec(1,2-iso)
        c1(i) = vec(2,2-iso)
        c2(i) = vec(3,2-iso)
        eint(nn) = eig(2)
        c0(nn) = vec(1,2)
        c1(nn) = vec(2,2)
        c2(nn) = vec(3,2)
!  at this step put all indices of omega=2 levels
!  at the top of the vectors
        eint(nmax-n+i) = eig(2+iso)
        c0(nmax-n+i) = vec(1,2+iso)
        c1(nmax-n+i) = vec(2,2+iso)
        c2(nmax-n+i) = vec(3,2+iso)
      end if
    end if
  end if
170 continue
!  nn now contains the number of levels for up to two spin-orbit states
nnn = nn
if (imult .eq. 3) then
!  here for 3pi add the omega=2 levels
  do 180 i = 1, n
    if (j(i) .le. 1) then
!  here for j = 0 or 1, which can exist only for omega = 0 or 1
      is(i) = (2 - iso) * is(i)
    else
!  here for j .ge. 2
!  the index nnn points to the nominal omega=2 levels
      nnn = nnn + 1
      j(nnn) = j(i)
      is(nnn) = (2 + iso) * is(i)
      is(i) = (2 - iso) * is(i)
      eint(nnn) = eint(nmax-n+i)
      c0(nnn) = c0(nmax-n+i)
      c1(nnn) = c1(nmax-n+i)
      c2(nnn) = c2(nmax-n+i)
    end if
180   continue
end if
!  set n to equal total number of levels
n = nnn
if (n .gt. nmax) then
  write (6, 190) n, nmax
  write (9, 190) n, nmax
190   format (/' *** NLEVELS =', i4, ' > MAX DIMENSION OF', &
          i4, '; ABORT ***')
  stop
end if
!  find lowest energy
emin = 1.e+7
do 195  i = 1, n
  if (eint(i) .lt. emin) emin = eint(i)
195 continue
!  now shift energies so lowest energy is zero
do 200 i = 1, n
  eint(i)=eint(i)-emin
200 continue
!  now evaluate mixing of lambda-doublets by an external stark field
!  assume that the field is weak enough that the rotational states
!  are not mixed
!  (see g.c. corey and d. lemoine, chem. phys. letters 160, 3 (1989))
!  ce and cf are the coefficients of the e/f field-free basis functions
!  for 1pi and 2pi states i is the nominal f state (eps = -1) and i+1 is
!  the nominal e state (eps = +1) (if npar = 2); cf. do 160 loop over ip
!  this reverses for 3pi states
if (efield .eq. zero) then
  do 205 i = 1, n, npar
    cf(i) = one - ispin
    if (npar .eq. 2) cf(i+1) = float(ispin)
205   continue
else
  write (9, 210)
210   format (/' ** STARK MIXING OF THE LAMBDA-DOUBLETS', &
          ' BY AN EXTERNAL ELECTRIC FIELD **')
  xommin = iommin + xhf
  fstark = two * cvtown * dmom * efield * (nu + xhf)
  angmix = quart * sign(one,fstark) * acos(-one)
  do 220 i = 1, n, 2
    if (j(i) .eq. zero .and. imult .eq. 3) then
!  for 3pi no stark mixing for j = 0, which can exist only for omega = 0
      cf(i) = zero
      cf(i+1) = one
      go to 220
    end if
!  save the e/f field-free energies temporarily in the variables
!  ee and ef, and evaluate stark coupling matrix element (hef/2)
!  and mixing angle (angle)
    ie = i - ispin + 1
    if = i + ispin
    ee = eint(ie)
    ef = eint(if)
    xomega = xommin + (one + xhf) * c1(i)**2 + two * c2(i)**2
    hef = fstark * xomega / ((j(i) + xhf) * (j(i) + xhf + one))
    if (ee .ne. ef) then
      angle = half * atan(hef/(ee-ef))
    else
      angle = angmix
    end if
!  determine mixing coefficients and calculate energies of mixed states
    ceie = cos(angle)
    cf(ie) = sin(angle)
    ceif = - cf(ie)
    cf(if) = ceie
    eint(ie) = ceie**2 * ee + cf(ie)**2 * ef + ceie * cf(ie) * hef
    eint(if) = ceif**2 * ee + cf(if)**2 * ef + ceif * cf(if) * hef
220   continue
end if
!  form list of all energetically and/or symmetrically distinct
!  rotational levels included in the channel basis and their energies
nlevel = 0
do 230 i = 1, n
  nlevel = nlevel + 1
  jhold(nlevel) = j(i)
  ishold(nlevel) = is(i)
  eint(i) = eint(i) / econv
  ehold(nlevel) = eint(i)
230 continue
!  now sort this list to put closed levels at end
!  also determine number of levels which are open
nlevop = 0
do 250 i = 1, nlevel-1
  if (ehold(i) .le. ered) then
    nlevop = nlevop + 1
  else
    do 240 ii = i+1, nlevel
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
        go to 250
      end if
240     continue
  end if
250 continue
if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
if (csflag) then
!  set up coupled-states channel basis (if desired)
  nn = 0
  do 260 i = 1, n
!  include only channels with j at least equal to coupled states
!  projection index
    if (j(i) .ge. nu) then
      nn = nn + 1
      j(nn) = j(i)
      is(nn) = is(i)
      eint(nn) = eint(i)
      c0(nn) = c0(i)
      c1(nn) = c1(i)
      c2(nn) = c2(i)
      cf(nn) = cf(i)
      l(nn) = jtot
      cent(nn) = jtot * (jtot + 1)
    end if
260   continue
else if (.not. csflag) then
!  set up close-coupled channel basis (if desired)
!  first move all indices of rotational levels to the top of the vectors
!  e.g. move j(n) to j(nmax), j(n-1) to j(nmax-1),  etc
  do 270 i = 1, n
!  move (n-i+1)th element to (nmax-i+1)th element
    inew = nmax - i + 1
    iold = n - i + 1
    j(inew) = j(iold)
    is(inew) = is(iold)
    eint(inew) = eint(iold)
    c0(inew) = c0(iold)
    c1(inew) = c1(iold)
    c2(inew) = c2(iold)
270   continue
  nn = 0
  ihf = nint(xhf+xhf)
  do 290 i = 1, n
!  now take (nmax-n+i)th element, duplicate it as many times as is
!  required for rotational degeneracy, and store the new elements
!  in nn, nn+1, ....
    ipoint = nmax - n + i
    ji = j(ipoint)
    lmax = jtot + ji + ihf
    lmin = iabs(jtot-ji)
!  in the e/f field-free basis:  is = i*eps for Fi level
!  where i = 1,2,3 and eps = +/- 1 is the "true" case (a) symmetry index
    if (isign(1,is(ipoint))*(-1)**(ji-ispin+lmin-jtot) .ne. jlpar) &
        lmin = lmin + 1
    do 280 li = lmin, lmax, 2
      nn = nn + 1
      j(nn) = ji
      is(nn) = is(ipoint)
      eint(nn) = eint(ipoint)
      c0(nn) = c0(ipoint)
      c1(nn) = c1(ipoint)
      c2(nn) = c2(ipoint)
      l(nn) = li
      cent(nn) = li * (li + 1)
280     continue
290   continue
end if
!  set number of coupled channels
n = nn
if (n .gt. nmax) then
  write (6, 300) n, nmax
  write (9, 300) n, nmax
300   format (/' *** NCHANNELS =', i4, ' > MAX DIMENSION OF', &
          i4, '; ABORT ***')
  stop
end if
!  now check to see if any of the open channels are closed at r=rcut
!  this is NOT done for molecule-surface collisions or bound-state calcs!
if (.not. flagsu .and. .not.boundc) then
  ecut = float(jtot) * (jtot + one) / (two * rmu * rcut**2)
  emin = 1.e+7
  do 310  i = 1, n
    if (eint(i) .le. ered) then
!  here if channel is open asymptotically
      if (ecut .gt. ered-eint(i)) then
!  here if channel is open asymptotically but closed at r = rcut
        if (eint(i) .lt. emin) emin = eint(i)
      end if
    end if
310   continue
!  emin now contains the lowest channel energy for which this condition
!  is met. eliminate all channels with eint .ge. emin if any of the
!  channels are open asymptotically but closed at r = rcut
  if (emin .lt. ered .and. rcut .lt. rendai) then
    nn = 0
    do 320 i = 1, n
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        j(nn) = j(i)
        is(nn) = is(i)
        eint(nn) = eint(i)
        c0(nn) = c0(i)
        c1(nn) = c1(i)
        c2(nn) = c2(i)
        cent(nn) = cent(i)
        l(nn) = l(i)
      end if
320     continue
!  reset number of channels
    n = nn
  end if
end if
!  return if no channels
if (n .eq. 0) return
!      if (nu .eq. numin) then
  ntop = max(n,nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on scalar machines such as a vax or an hp
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) ntop = ntop + 1
!      end if
!  now list channels if requested
if (clist) then
  string = '   N   J   EPS  Fi   L   EINT(CM-1)'
  if (nman .eq. 2) string(36:49) = '  C-1/2  C-3/2'
  if (imult .eq. 3) string(36:56) = '   C-0    C-1    C-2 '
  if (efield .ne. zero) then
    rdtodg = 180.d0 / acos(-one)
    index = 35
    if (nman .ne. 1) index = 35 + nman * 7 + 3
    string(index+1:index+21) = '  C-e    C-f   ANGMSF'
  end if
  if (bastst) write (6, 330) string
  write (9, 330) string
330   format (/a)
  do 370 i = 1, n
!  in the e/f field-free basis:  is = i*eps for Fi level
!  where i = 1,2,3 and eps = +/- 1 is the "true" case (a) symmetry index
!  therefore eps is the sign of is and Fi is the magnitude of is
    if (nman .eq. 1) then
      if (bastst) write (6, 333) i, j(i), chf, isign(1,is(i)), &
                               iabs(is(i)), l(i), eint(i)*econv
      write (9, 333) i, j(i), chf, isign(1,is(i)), &
                   iabs(is(i)), l(i), eint(i)*econv
333       format (2i4, a, i4, i3, i5, f10.3)
    endif
    if (nman .eq. 2) then
      if (bastst) write (6, 335) i, j(i), chf, isign(1,is(i)), &
                               iabs(is(i)), l(i), eint(i)*econv, &
                               c0(i), c1(i)
      write (9, 335) i, j(i), chf, isign(1,is(i)), &
                   iabs(is(i)), l(i), eint(i)*econv, &
                               c0(i), c1(i)
335       format (2i4, a, i4, i3, i5, f10.3,3x,3f7.3)
    else if (imult .eq. 3) then
      if (bastst) write (6, 335) i, j(i), chf, isign(1,is(i)), &
                               iabs(is(i)), l(i), eint(i)*econv, &
                               c0(i), c1(i), c2(i)
      write (9, 335) i, j(i), chf, isign(1,is(i)), &
                   iabs(is(i)), l(i), eint(i)*econv, &
                               c0(i), c1(i), c2(i)
    end if
    if (efield .ne. zero) then
      cei = (-1)**(ispin+i) * cf(i-(-1)**i)
!  angmsf is the orientation angle between the
!  molecular axis and the space-fixed axis Z
      angmsf = acos(two*(xommin+(one+xhf)*c1(i)**2+two*c2(i)**2)* &
               cei*cf(i)*(nu+xhf)/((j(i)+xhf)*(j(i)+xhf+one)))
      if (bastst) write (6, 360) cei, cf(i), rdtodg*angmsf
      write (9, 360) cei, cf(i), rdtodg*angmsf
    end if
360     format (2x, 2f7.3, f8.2)
370   continue
end if
!  now calculate coupling matrix elements
!  i counts v2 elements
!  inum counts v2 elements for given lambda
!  ilam counts number of v2 matrices
if (bastst.and. iprint.ge. 2) then
  write (6, 380)
  write (9, 380)
380   format (/' ILAM  LAMBDA  ICOL  IROW    I   IV2      VLAM')
end if
i = 0
lamsum = 0
if (csflag) lcol = nu
vii(0) = zero
vii(1) = zero
vii(2) = zero
v11 = zero
do 450 ilam = 1, nlam
!  ilam is the angular expansion label
  inum = 0
  do 430 icol = 1, n
    if (efield .eq. zero) then
      iepsc = isign(1,is(icol))
    else
      iepsc = 1
      ceicol = (-1)**(ispin+icol) * cf(icol-(-1)**icol)
    end if
    if (.not. csflag) lcol = l(icol)
    do 420 irow = icol, n
      if (efield .eq. zero) then
        iepsr = isign(1,is(irow))
      else
        iepsr = 1
        ceirow = (-1)**(ispin+irow) * cf(irow-(-1)**irow)
      end if
390       if (ilam .le. nlam0) then
!  here for l=0 terms in potential (average potential)
!  diagonal in omega
        lb = lammin(1) + (ilam - 1) * istep
!  lb is the actual value of lambda
        if (iabs(j(irow)-j(icol)) .gt. lb) go to 420
        ii = iommin - 1
400         ii = ii + 1
        call vlmpi (j(irow), l(irow), j(icol), lcol, jtot, ii, &
                     ii, lb, iepsr, iepsc, vii(ii), xhf, csflag)
        if (ii .lt. iommin+nman-1) go to 400
        vlam = vii(1)
        if (imult .ne. 1) vlam = &
            vii(0) * c0(irow) * c0(icol) + vii(1) * c1(irow) &
                   * c1(icol) + vii(2) * c2(irow) * c2(icol)
        k = 1
      else
!  here for l=2 terms in potential (difference potential)
!  off-diagonal in omega for 2pi states
        lb = lammin(2) + (ilam - nlam0 - 1) * istep
!  lb is the actual value of lambda
        if (iabs(j(irow)-j(icol)) .gt. lb) go to 420
        if (imult .ne. 2) &
        call vlmpi (j(irow), l(irow), j(icol), lcol, jtot, &
                     1, 3, lb, iepsr, iepsc, v11, xhf, csflag)
        vlam = v11
        if (imult .ne. 1) then
          ii = imult - 1
          call vlmpi (j(irow), l(irow), j(icol), lcol, jtot, &
                       0, ii, lb, iepsr, iepsc, v0i, xhf, csflag)
          call vlmpi (j(irow), l(irow), j(icol), lcol, jtot, &
                       ii, 0, lb, iepsr, iepsc, vi0, xhf, csflag)
          if (imult .eq. 2) then
            vlam = v0i * c0(irow) * c1(icol) + &
                   vi0 * c1(irow) * c0(icol)
          else if (imult .eq. 3) then
            vlam = v0i * c0(irow) * c2(icol) + v11 * c1(irow) &
                   * c1(icol) + vi0 * c2(irow) * c0(icol)
          end if
        end if
        k = - 1
      end if
      if (efield .ne. zero) then
!  here for potential matrix elements in basis mixed by a stark field
        if (iepsr .eq. iepsc) then
!  terms diagonal in eps
          vepsd = vlam
          iepsr = - 1
          go to 390
        else
!  terms off-diagonal in eps
          vepso = vlam
        end if
!  for 1pi and 2pi states eps = +1 for e levels, eps = -1 for f levels
!  while for 3pi states eps = -1 for e levels, eps = +1 for f levels
        vlam = vepsd * (ceirow * ceicol + k * cf(irow) * cf(icol)) &
             + vepso * (ceirow * cf(icol) + k * cf(irow) * ceicol)
        if (imult .eq. 3 .and. k .eq. -1) vlam = - vlam
      end if
      if (vlam .ne. zero) then
        i = i + 1
        if (i .le. nv2max) then
          inum = inum + 1
          v2(i) = vlam
          iv2(i) = ntop * (icol - 1) + irow
          if (bastst.and. iprint .ge. 2) then
            write (6, 410) ilam, lb, icol, irow, i, iv2(i), vlam
            write (9, 410) ilam, lb, icol, irow, i, iv2(i), vlam
410             format (i4, 2i7, 3i6, g17.8)
          end if
        end if
      end if
420     continue
430   continue
  if (i .le. nv2max) lamnum(ilam) = inum
  if (bastst) then
    write (6, 440) ilam, lamnum(ilam)
    write (9, 440) ilam, lamnum(ilam)
440     format ('ILAM=', i3, ' LAMNUM(ILAM) =', i6)
  end if
  lamsum = lamsum + lamnum(ilam)
450 continue
if (i .gt. nv2max) then
  write (6, 460) i, nv2max
  write (9, 460) i, nv2max
460   format (' *** NUMBER OF NONZERO V2 ELEMENTS = ', i6, &
          ' > NV2MAX=', i6, '; ABORT ***')
  stop
end if
if (clist .and. bastst) then
  write (6, 470) lamsum
  write (9, 470) lamsum
470   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', i6)
end if
return
end
!  ---------------------------------------------------------------------
subroutine vlmpi (jp, lp, j, l, jtot, iomegp, iomeg, lambda, &
                   iepsp, ieps, v, xhf, csflag)
!  ---------------------------------------------------------------------
!  subroutine to calculate v-lambda matrices for close-coupled and
!  coupled-states treatments of collisions of a molecule in a 1pi,
!  2pi or 3pi electronic state in a Hund's case (a) basis
!  the cc matrix elements are given in eqs. (27), (29) (except that
!  the sign in front of the second term in square brackets must be
!  minus) and (30) of m.h. alexander, chem. phys. 92, 337 (1985)
!  the cs matrix elements are given in eq. (14) of t. orlikowski and
!  m.h. alexander, j. chem. phys. 79, 6006 (1983), for a 2pi state
!  note that for cc collisions with a flat surface, the coupling
!  matrix elements [see m.h. alexander, j. chem. phys. 80, 3485 (1984)
!  for a 2pi molecule] are identical to the cs matrix elements here
!  author:  didier lemoine
!  current revision date:  9-jun-90
!  ---------------------------------------------------------------------
!  variables in call list:
!    jp:      rotational quantum number of left side of v element (bra)
!    lp:      orbital angular momentum of left side of v element (bra)
!    j:       rotational quantum number of right side of v element (ket)
!    l:       orbital angular momentum of right side of v element (ket)
!    jtot:    total angular momentum
!    iomegp:  omega quantum number of bra
!    iomeg:   omega quantum number of ket
!             iomeg should be set to 3 for coupling within the omega=1
!             manifold of a 1pi or 3pi state by the difference potential
!    lambda:  order of legendre term in expansion of potential
!    iepsp:   symmetry index of bra
!    ieps:    symmetry index of ket
!    v:       on return, contains matrix element
!    xhf:     equals 0 for a 1pi or 3pi molecule calculation
!             equals 1/2 for a 2pi molecule calculation
!    csflag:  if .true., then cs calculation; in this case:
!               jtot in call list is cs lbar
!               l is nu (cs projection index)
!               lp is not used
!             if .false., then cc calculation; in this case:
!               jtot is total angular momentum
!    note!!   if xhf = 1/2, then the true values of the rotational
!    quantum numbers, the total angular momentum, the omega quantum
!    numbers and the coupled-states projection index are equal to the
!    values stored in j, jp, jtot, iomeg, iomegp and nu plus 1/2
!    for collisions with a surface, nu is equivalent to m
!    (the projection of j along the surface normal)
!  functions called:
!     xf3j, xf6j
!  ---------------------------------------------------------------------
!      implicit none
integer ieps, iepsp, iomeg, iomegp, iphase, &
        j, jp, jtot, l, lp, lambda, nu
double precision v, x, xf3j, xf6j, xhf, xhf2, xj, xjp, xjtot, &
     xl, xlb, xlp, xnorm, xnu, xomeg, xomegp, zero, one, two
logical csflag
zero = 0.d0
one = 1.d0
two = 2.d0
v = zero
if (ieps*iepsp*(-1)**(j+jp+lambda+1) .eq. 1) return
xlb = float(lambda)
xjp = jp + xhf
xj = j + xhf
xhf2 = xhf + xhf
xnorm = (xjp + xjp + one) * (xj + xj + one)
if (csflag) then
  nu = l
  xnu = nu + xhf
  x = xf3j (xjp, xlb, xj, -xnu, zero, xnu)
  iphase = nu - iomeg
else
  xlp = float(lp)
  xl = float(l)
  x = xf3j (xlp, xlb, xl, zero, zero, zero)
  if  (x .eq. zero) return
  xjtot = jtot + xhf
  x = x * xf6j (xjp, xlp, xjtot, xl, xj, xlb)
  iphase = jtot + j + jp - iomeg + nint(xhf2)
  xnorm = xnorm * (xlp + xlp + one) * (xl + xl + one)
end if
if (x .eq. zero) return
xomegp = iomegp + xhf
xomeg = iomeg + xhf
if (xomeg .eq. 3.) xomeg = one
if (iomegp .eq. iomeg) then
  x = x * xf3j (xjp, xlb, xj, -xomegp, zero, xomegp)
else
  x = (one - xhf2 - xhf2) * ieps * x * &
      xf3j (xjp, xlb, xj, -xomegp, two, -xomeg)
end if
if (x .ne. zero) v = (-1)**iabs(iphase) * x * sqrt(xnorm)
return
end
!  ---------------------------------------------------------------------
subroutine sypi (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for singlet,
!  doublet or triplet pi molecule + atom or flat surface scattering
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  author:  didier lemoine and millard alexander
!  current revision date:  4-mar-1996 by mha
!  ---------------------------------------------------------------------
!  variable in common /cosys
!    scod:    character*8 array of dimension nscode, which contains
!             names of all system dependent parameters. note that the
!             ordering of the variable names in scod must correspond
!             to the ordering of the variable names in cosysi, cosysr,
!             cosysl and cobspt respectively
!  variables in common block /cosysi/
!    nscode:   total number of system dependent parameters
!              nscode = isicod + isrcod + islcod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different angular terms in potential
!              NB for pi molecule, nterm should be 1 or 2
!    jmax:     the maximum rotational angular momenta (-1/2 for 2pi) for
!              each channel in each spin-orbit manifold with convention
!              omega .le. j .le. jmax(+1/2 for 2pi)
!    igu:      permutation inversion symmetry of electronic state
!              igu=1 for gerade states, igu=-1 for ungerade states
!              for heteronuclear molecules igu is not used
!    isa:      s/a label for molecular states. if ihomo=.true.,
!              then only s states will be included if isa=1 and
!              only a states if isa=-1
!    npar:     number of symmetry doublets included
!              (npar=2 will ensure both lambda doublets)
!    imult:    spin multiplicity of pi state (imult = 2*S+1)
!    nman:     number of spin-orbit manifolds included
!              for a 2pi state nman=1 selects the
!              omega=1/2 manifold in a case (a) basis
!              otherwise, nman is set = imult in subroutine bapi
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    brot:     rotational constant in cm-1
!    aso:      spin-orbit constant in cm-1
!    o, p, q:  lambda-doubling constants cm-1
!    dmom:     dipole moment of molecule in Debye
!    efield:   stark field in kV / cm
!  variables in common bloc /cosysl/
!    islcod:   total number of logical system dependent variables
! ----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_conlam, only: nlam
use mod_cosys, only: scod
use mod_cosysl, only: islcod, lspar
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
implicit none
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: j, l, lc
logical existf
character*(*) fname
character*60 filnam, line, potfil, filnm1
character*1 dot
#include "common/parbas.F90"
common /coskip/ nskip,iskip
integer :: nskip, iskip
save potfil
#include "common/comdot.F90"
integer, pointer :: nterm, jmax, igu, isa, npar, imult, nman
real(8), pointer :: brot, aso, o, p, q, dmom, efield
nterm=>ispar(1); jmax=>ispar(2); igu=>ispar(3); isa=>ispar(4); npar=>ispar(5); imult=>ispar(6); nman=>ispar(7)
brot=>rspar(1); aso=>rspar(2); o=>rspar(3); p=>rspar(4); q=>rspar(5); dmom=>rspar(6); efield=>rspar(7)

potfil = ' '
!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
isicod = 7
scod(1) = 'NTERM'
scod(2) = 'JMAX'
scod(3) = 'IGU'
scod(4) = 'ISA'
scod(5) = 'NPAR'
scod(6) = 'IMULT'
scod(7) = 'NMAN'
!  then all the system dependent real variables
!  in the same order as in the common block /cosysr/
isrcod = 7
scod(8) = 'BROT'
scod(9) = 'ASO'
scod(10) = 'O'
scod(11) = 'P'
scod(12) = 'Q'
scod(13) = 'DMOM'
scod(14) = 'EFIELD'
!  then all the system dependent logical variables
!  in the same order as in the common block /cosysl/
islcod = 0
nscode = 14
!  set default values for doublet pi scattering
nterm = 2
if (iread .eq. 0) then
  jmax = 3
  igu = 1
  isa = 0
  npar = 2
  nman = 0
  imult = 2
  mproj(1) = 0
  mproj(2) = 2
  lammin(1)= 1
  lammax(1) = 1
  lammin(2) = 2
  lammax(2) = 2
  niout=4
  indout(1)=-1
  indout(2)=1
  indout(3)=-2
  indout(4)=2
endif
if (iread .eq. 0) return
!  line 19
!  line 20
read (8, *, err=888) jmax, igu, isa, npar, imult, nman
!  line 21
read (8, *, err=888) brot, aso
!  line 22
read (8, *, err=888) o, p, q
!  line 23
read (8, *, err=888) dmom, efield
!  line 16 name of file containing potential parameters
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
entry ptrpi (fname,readpt)
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
entry savpi (readpt)
!  save input parameters for singlet, doublet or
!  triplet-pi molecule + atom scattering
!  the order of the write statements should be identical to the read
!  statements above. for consistency with the data file written by
!  gendat, format statements should reserve the first 30 spaces for
!  data, spaces 31-33 should be left blank, and the names of the
!  variables should be printed in spaces 34-80
!  line 18:
!  line 19
!  line 20
write (8, 130) jmax, igu, isa, npar, imult, nman
130 format (6i4, 6x, '   jmax, igu, isa, npar, imult, nman')
!  line 21
write (8, 140) brot, aso
140 format (f10.5,f10.4, 10x, '   brot, aso')
!  line 22
write (8, 150) o, p, q
150 format (3(1pg12.4), ' o, p, q')
!  line 23
write(8, 160) dmom, efield
160 format (2f10.5, 10x, '   dmom, efield')
!  line 16
write (8, 285) potfil
return
end
