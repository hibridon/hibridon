! syastp1 (savastp1/ptrastp1) defines, saves variables and reads         *
!                  potential for C2v asymmetric top-atom scattering      *
! ----------------------------------------------------------------------
subroutine baastp1 (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, etemp, fjtemp, fktemp, fistmp, &
                  rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, &
                  jlpar, n, nmax, ntop)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of an asymmetric top molecule of C2v symmetry with a structureless atom
!  or with an uncorrugated surface
!
!  modification of hibaastp basis routine
!
!  here, the body-frame z system is taken to lie along the C2 axis of the molecule,
!  following green's convention
!
!  author:  paul dagdigian
!  current revision:  21-sep-2017 by pjd
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains rotational quantum number for each
!              channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains symmetry index (ieps * kp) for each
!              channel
!  Note that we have adopted the following convention for the symmetry
!  index "is" so that on return is = 0/1, with the symmetric top basis
!  functions written as [|j,kp> + (-1)^is*|j,-kp)\, respectively.
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
!    etemp:    scratch array used to create channel list
!    fjtemp:   scratch array used to create channel list
!    fktemp:   scratch array used to create channel list
!    fistmp:   scratch array used to create channel list
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
!              (e.g. H2O, CH2, H2CO), so that only the ortho or para levels will be
!              included depending on the value of the parameter iop in common
!              /cosysi/ (see below)
!              subroutine currently to allow only BA2 or A2BC type molecules,
!              where the a inertial axis is an axis of C2 symmetry.
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              only those channels are included for which
!                  (-1)**(is+j+kp+l-jtot)=jlpar
!              where parity designates the parity of the molecular state
!              in cs calculation jlpar is set equal to 1 in calling program
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to subroutines
!              propag and soutpt.  ntop is set in basis only if nu = numin
!              otherwise it is unchanged from the value supplied by the
!              calling program
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    arot:     A rotational constant for the asymmetric top
!    brot:     B rotational constant for the asymmetric top
!    crot:     C rotational constant for the asymmetric top
!    emax:     the maximum rotational energy (in cm^-1) for a channel to be
!              included in the basis
!  variables in common block /cosysi/
!    nscode:   total number of variable names which are passed to HINPUT
!              nscod must equal isrcod + isicod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    numpot:   the number of the potential used, this variable is passed
!              to the pot subroutine
!    iop:      ortho/para label for molecular states. If ihomo=.true. then only
!              para states will be included if iop=1 and only ortho states if
!              iop=-1
!    jmax:     the maximum rotational angular momentum for the asymmetric top
!  variables in common block /coered/
!    ered:      collision energy in atomic units (hartrees)
!    rmu:       collision reduced mass in atomic units
!               (mass of electron = 1)
!  variables in common block /coconv/
!   econv:      conversion factor from cm-1 to hartrees
!   xmconv:     converson factor from amu to atomic units
!  subroutines called:
!   vlmatp:     returns angular coupling coefficient for particular
!               choice of channel index
!   prmatp:     computes primitive cc and cs v-lambda matrix elements
!               between symmetric top basis fns.
!   rotham:     computes matrix elements of asymmetric top hamiltonian
! --------------------------------------------------------------------
use mod_cov2, only: nv2max, junkv => ndummy, v2
use mod_coiv2, only: iv2
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coatpi, only: narray, isiz
use mod_coatpr, only: c
use mod_coatp1, only: ctemp
use mod_coatp2, only: chold
use mod_coatp3, only: isizh
use mod_conlam, only: nlam, nlammx, lamnum
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use constants, only: econv, xmconv
use mod_par, only: iprint
#include "common/parbasl.F90"
implicit double precision (a-h,o-z)
logical flaghf, csflag, clist, flagsu, ihomo, bastst
#include "common/parbas.F90"
common /coered/ ered, rmu
dimension j(1), l(1), is(1), jhold(1), ehold(1), &
          ishold(1), etemp(1), fjtemp(1), fktemp(1), &
          fistmp(1)
!  scratch arrays for computing asymmetric top energies and wave fns.
dimension e(narray,narray), eig(narray), vec(narray,narray), &
  sc1(narray), sc2(narray), work(1000)
!
integer, pointer :: nterm, numpot, iop, jmax
real(8), pointer :: arot, brot, crot, emax
nterm=>ispar(1); numpot=>ispar(2); iop=>ispar(3); jmax=>ispar(4)
arot=>rspar(1); brot=>rspar(2); crot=>rspar(3); emax=>rspar(4)

zero = 0.d0
two = 2.d0
!
!  this basis routine applies to symmetrical molecules.  so set ihomo =.true.
ihomo = .true.
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
if (flagsu) then
  write (6, 10)
  write (9, 10)
10   format ('  *** FLAGSU = .TRUE. FOR ASYMMETRIC TOP', &
     ' COLLISIONS; NOT IMPLEMENTED.  ABORT ***')
  call exit
end if
xjtot = jtot
nsum = 0
do 35  i = 1, nterm
  if (mproj(i) .gt. lammin(i) ) then
    write (6, 30) mproj(i), lammin(i), i
    write (9, 30) mproj(i), lammin(i), i
30     format (' *** MPROJ=',i2,' > LAMMIN=',i2, &
            ' FOR TERM',i2,'; ABORT ***')
    stop
  end if
  nsum = nsum + (lammax(i) - lammin(i)) + 1
35 continue
if (nlammx .lt. nsum) then
  write (6, 40) nsum, nlammx
  write (9, 40) nsum, nlammx
40   format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2, &
          ' .GT. NLAMMX=', i2,'; ABORT')
  stop
end if
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
  if (csflag) then
    if (ihomo) then
      if (bastst) then
        write (6,65) rmu * xmconv, arot, brot, crot, &
             iop, ered * econv, jtot, nu
      end if
      write (9,65) rmu * xmconv, arot, brot, crot, &
             iop,ered * econv, jtot, nu
65       format(/,' **  CS ASYMMETRIC TOP **', &
        /,'     RMU=', f9.4,'  AROT=', f7.3,'  BROT=', &
        f7.3,/,'  CROT=',f7.3,/, &
       '  O/P=',i2,'  E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
    else
      if (bastst) then
        write (6,75) rmu * xmconv, arot, brot, crot, &
              ered * econv, jtot, nu
      end if
      write (9,75) rmu * xmconv, arot, brot, crot, &
              ered * econv, jtot, nu
75       format(/,' **  CS ASYMMETRIC TOP **', &
        /,'     RMU=', f9.4,'  AROT=', f7.3,'  BROT=', &
        '  CROT=',f7.3,/, &
        '     E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
    end if
  else
    if (ihomo) then
      if (bastst) then
        write (6,80) rmu * xmconv, arot, brot, crot, &
             iop, ered * econv, jtot, jlpar
      end if
      write (9,80) rmu * xmconv, arot, brot, crot, &
              iop, ered * econv, jtot, jlpar
80       format(/,' **  CC ASYMMETRIC TOP **', &
        /,'     RMU=', f9.4,'  AROT=', f7.3, '  BROT=',f7.3, &
        '  CROT=',f7.3,/,'    O/P=',i2, &
        '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
    else
      if (bastst) then
        write (6,85) rmu * xmconv, arot, brot, crot, &
             ered * econv, jtot, jlpar
      end if
      write (9,85) rmu * xmconv, arot, brot, crot, &
             ered * econv, jtot, jlpar
85       format(/,' **  CC ASYMMETRIC TOP **', &
        /,'     RMU=', f9.4,'  AROT=', f7.3, '  BROT=',f7.3, &
       '  CROT=',f7.3,/,  '     E=', f7.2, '  JTOT=', i4, 2x, &
        ' JLPAR=', i2)
    end if
  end if
  if (.not. flagsu) write (9,90) rcut
90   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
end if
!  first set up list of all j(kp,ko) states included in basis
nlist = 0
do ji = 0, jmax
! set up and diagonalize rotational hamiltonian for each value of ji
! hamiltonian block diagonalize into 1-4 subblocks
!
! E+ submatrix
  isize = (ji + 2)/2
  do mm = 1, isize
    do nn = 1, isize
      e(mm,nn) = 0.d0
    end do
  end do
  e(1,1) = rotham1(ji,0,ji,0)
  if (isize .gt. 1) then
    do mm = 2, isize
      kk = 2*mm - 2
      e(mm,mm) = rotham1(ji, kk, ji, kk)
      e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
      e(mm-1,mm) = e(mm,mm-1)
    end do
    e(2,1) = sqrt(2.d0)*e(2,1)
    e(1,2) = e(2,1)
  end if
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
  lwork = 144
  call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
#endif
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
write (6,1225)
1225   format('  *** DSYEV NOT ACCESSIBLE IN THIS VERSION OF CODE.')
  call exit
#endif
  do mm = 1, isize
    nlist = nlist + 1
!  eigenvalues stored in order of increasing energy (lowest KP first)
    etemp(nlist) = eig(mm)
    fjtemp(nlist) = ji
    fktemp(nlist) = 2*mm - 2
    fistmp(nlist) = 1
!  eigenfunctions expressed over basis KP = 0, 2, 4, ...
    nbas = int(fktemp(nlist)/2) + 1
    do nn = 1, isize
      isub = (nlist - 1)*narray + nn
!  make sure coeff for prolate k is positive
      if (e(mm, nbas) .gt. zero) then
        ctemp(isub) = e(mm, nn)
      else
        ctemp(isub) = -e(mm, nn)
      end if
    end do
  end do
!
! E- submatrix
  isize = ji/2
  if (isize .gt. 0) then
    do mm = 1, isize
      do nn = 1, isize
        e(mm,nn) = 0.d0
      end do
    end do
    e(1,1) = rotham1(ji,2,ji,2)
    if (isize .gt. 1) then
    do mm = 2, isize
      kk = 2*mm
      e(mm,mm) = rotham1(ji, kk, ji, kk)
      e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
      e(mm-1,mm) = e(mm,mm-1)
    end do
    end if
    lwork = 144
    call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
    do mm = 1, isize
      nlist = nlist + 1
!  eigenvalues stored in order of increasing energy (lowest KP first)
      etemp(nlist) = eig(mm)
      fjtemp(nlist) = ji
      fktemp(nlist) =  2*mm
      fistmp(nlist) = -1
!  eigenfunctions expressed over basis KP = 0, 2, 4, ...
      nbas = int(fktemp(nlist)/2)
      do nn = 1, isize
        isub = (nlist - 1)*narray + nn
!  make sure coeff for prolate k is positive
        if (e(mm, nbas) .gt. zero) then
          ctemp(isub) = e(mm, nn)
        else
          ctemp(isub) = -e(mm, nn)
        end if
      end do
    end do
  end if
!
! O+ submatrix
  isize = (ji + 1)/2
  if (isize .gt. 0) then
    do mm = 1, isize
      do nn = 1, isize
        e(mm,nn) = 0.d0
      end do
    end do
    e(1,1) = rotham1(ji,1,ji,1) + rotham1(ji,1,ji,-1)
    if (isize .gt. 1) then
    do mm = 2, isize
      kk = 2*mm - 1
      e(mm,mm) = rotham1(ji, kk, ji, kk)
      e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
      e(mm-1,mm) = e(mm,mm-1)
    end do
    end if
    lwork = 144
    call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
    do mm = 1, isize
      nlist = nlist + 1
!  eigenvalues stored in order of increasing energy (lowest KP first)
      etemp(nlist) = eig(mm)
      fjtemp(nlist) = ji
      fktemp(nlist) = 2*mm - 1
      fistmp(nlist) = 1
!  eigenfunctions expressed over basis KP = 1, 3, 5, ...
      nbas = int((fktemp(nlist) - 1)/2) + 1
      do nn = 1, isize
        isub = (nlist - 1)*narray + nn
!  make sure coeff for prolate k is positive
        if (e(mm, nbas) .gt. zero) then
          ctemp(isub) = e(mm, nn)
        else
          ctemp(isub) = -e(mm, nn)
        end if
      end do
    end do
  end if
!
! O- submatrix
  isize = (ji + 1)/2
  if (isize .gt. 0) then
    do mm = 1, isize
      do nn = 1, isize
        e(mm,nn) = 0.d0
      end do
    end do
    jx = jj
    e(1,1) = rotham1(ji,1,ji,1) - rotham1(ji,1,ji,-1)
    if (isize .gt. 1) then
      do mm = 2, isize
        kk = 2*mm - 1
        e(mm,mm) = rotham1(ji, kk, ji, kk)
        e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
        e(mm-1,mm) = e(mm,mm-1)
      end do
    end if
    lwork = 144
    call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
    do mm = 1, isize
      nlist = nlist + 1
!  eigenvalues stored in order of increasing energy (lowest KP first)
      etemp(nlist) = eig(mm)
      fjtemp(nlist) = ji
      fktemp(nlist) = 2*mm - 1
      fistmp(nlist) = -1
!  eigenfunctions expressed over basis KP = 1, 3, 5, ...
      nbas = int((fktemp(nlist) - 1)/2) + 1
      do nn = 1, isize
        isub = (nlist - 1)*narray + nn
!  make sure coeff for prolate k is positive
        if (e(mm, nbas) .gt. zero) then
          ctemp(isub) = e(mm, nn)
        else
          ctemp(isub) = -e(mm, nn)
        end if
      end do
    end do
  end if
end do
!
!  now sort this list in terms of increasing energy
if (nlist .gt. 1) then
  do 120 i1 = 1, nlist - 1
    esave = etemp(i1)
    do 115 i2 = i1 + 1, nlist
      if (etemp(i2) .lt. esave) then
!  state i2 has a lower energy than state i1, switch them
        esave = etemp(i2)
        etemp(i2) = etemp(i1)
        etemp(i1) = esave
        fjsave = fjtemp(i2)
        fjtemp(i2) = fjtemp(i1)
        fjtemp(i1) = fjsave
        fksave = fktemp(i2)
        fktemp(i2) = fktemp(i1)
        fktemp(i1) = fksave
        fissav = fistmp(i2)
        fistmp(i2) = fistmp(i1)
        fistmp(i1) = fissav
!  also move e.fn coeffs (don't worry about size of e.fn)
        do mm = 1, narray
          isub1 = (i1 - 1)*narray + mm
          isub2 = (i2 - 1)*narray + mm
          sc1(mm) = ctemp(isub2)
          ctemp(isub2) = ctemp(isub1)
          ctemp(isub1) = sc1(mm)
        end do
      end if
115     continue
120   continue
end if
!
!  now set up channel and level list
!  print this list if bastst = .true. or if clist = .true.
if (bastst .or. clist) then
  write (6, 130)
130   format (/,2x, &
   'LEVEL LIST SORTED BY ENERGY',/,'   N   J  ', &
     'IS  KP  KO      EINT(CM-1)')
end if
n = 0
nlevel = 0
do 170  njk = 1, nlist
  ki = fktemp(njk)
  ji = fjtemp(njk)
  isi = fistmp(njk)
!
!  delete state if energy > emax
  if (etemp(njk) .gt. emax) go to 170
!
!  check to see if (ki/isi) corresponds to an allowed ortho/para level for
!  a symmetric molecule.  the statements below correspond to BA2 and A2BC
!  type molecules like H2O, CH2, and H2CO
!
!  first compute kp and ko projection quantum numbers
  if (ji .eq. 2*(ji/2)) then
!  even j
    if (isi .eq. 1) then
!  isi = +1
      kp = 2 * ceiling(float(ki) * 0.5d0)
      ko = ji - ki
      goto 9900
    else
!  isi = -1
      kp = 2 * ceiling(float(ki - 2) * 0.5d0) + 1
      ko = ji + 1 - ki
      goto 9900
    end if
  else
!  odd j
    if (isi .eq. 1) then
!  isi = +1
      kp = 2 * ceiling(float(ki - 1) * 0.5d0) + 1
      ko = ji - ki
      goto 9900
    else
!  isi = -1
      kp = 2 * ceiling(float(ki - 1) * 0.5d0)
      ko = ji + 1 - ki
      goto 9900
    end if
  end if
9900   continue
  if ((kp+ko) .eq. 2*((kp+ko)/2)) then
    if (iop .eq. -1) goto 170
  else
    if (iop .eq. 1) goto 170
  end if
!  here if this state is to be included
  nlevel = nlevel + 1
  ehold(nlevel) = etemp(njk) / econv
  jhold(nlevel) = ji
  ishold(nlevel) = isi * ki
!  also move e.fn coeffs
  do mm = 1, narray
    isub = (nlevel - 1)*narray + mm
    isub1 = (njk - 1)*narray + mm
    chold(isub) = ctemp(isub1)
  end do
!
!  print this level if bastst = .true.
  if (bastst .or. clist) then
    ecm = ehold(nlevel) * econv
    write (6, 135) nlevel, jhold(nlevel), ishold(nlevel), &
       kp, ko, ecm
135     format (5i4, 3x, f10.3)
  end if
!
!  here for cs calculations; include state only if j at least equal to coupled
!  states projection index
  if (csflag) then
    if (ji .ge. nu) then
      n = n + 1
      if (n .gt. nmax) then
        write (6, 150) n, nmax
        write (9, 150) n, nmax
150         format(/' *** NCHANNELS=', i5, &
             ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
        stop
      end if
      is(n) = ishold(nlevel)
      j(n) = ji
      eint(n) = ehold(nlevel)
      l(n) = jtot
      cent(n) = jtot * (jtot + 1)
!  move e.fn also
      isiz(n) = isizh(nlevel)
      do mm = 1, narray
        isub = (n - 1)*narray + mm
        isub1 = (nlevel - 1)*narray + mm
        c(isub) = chold(isub1)
      end do
    end if
  else if (.not. csflag) then
!
!  here for cc calculations.  first calculate range of orbital angular
!  momentum quantum numbers allowed for this state
!  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem. Phys.
!  73, 2740 (1980) and Townes and Schawlow, Microwave Spectroscopy, Eq. (3-27),
!  p. 64.]  See Eq. (A3) of Green.
    ji = jhold(nlevel)
    ki = abs(ishold(nlevel))
    iss = sign(1, ishold(nlevel))
    ipar = (-1) ** (ji + ki) * iss
    lmax = jtot + ji
    lmin = iabs (jtot - ji)
    do 155  li = lmin, lmax
      ix = ipar * (-1) ** (li - jtot)
!  original defn
!            lpar = (-1) ** (li - jtot)
!
!  check to see if this channel has the desired parity, if so, include it
!  original defn
!            if (ipar * lpar .eq. jlpar) then
!
      if (ix .eq. jlpar) then
        n = n + 1
        if (n .gt. nmax) then
          write (6, 150) n, nmax
          write (9, 150) n, nmax
          stop
        end if
        is(n) = ishold(nlevel)
        j(n) = ji
        eint(n) = ehold(nlevel)
        l(n) = li
        cent(n) = li * (li + 1)
!  move e.fn also
        isiz(n) = isizh(nlevel)
        do mm = 1, narray
          isub = (n - 1)*narray + mm
          isub1 = (nlevel - 1)*narray + mm
          c(isub) = chold(isub1)
        end do
      end if
155     continue
  end if
170 continue
!
!  also determine number of levels which are open
nlevop = 0
do 250  i = 1, nlevel
  if (ehold(i) .le. ered) then
    nlevop = nlevop + 1
  end if
250 continue
!  now check to see if any of the open channels are closed at r=rcut
!  this is not done for molecule-surface collisions or for rcut < 0
if (.not.flagsu .and. rcut .gt. 0.d0 .and. .not.boundc) then
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
        j(nn) = j(i)
        is(nn) = is(i)
        cent(nn) = cent(i)
        l(nn) = l(i)
        isiz(nn) = isiz(i)
        do mm = 1, narray
          isub = (nn - 1)*narray + mm
          isub1 = (i - 1)*narray + mm
          c(isub) = c(isub1)
        end do
      end if
300     continue
!  reset number of channels
    n = nn
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
  endif
end if
!  now list channels if requested
if (clist) then
  if (.not.csflag) then
    if (bastst .or. clist) write (6, 305)
    write (9,305)
305     format ( &
    /,2x,'CHANNEL LIST SORTED BY ENERGY',/, &
     '   N   J  IS   L   EINT(CM-1)')
  else
    if (bastst .or. clist) write (6, 310) nu
    write (9,310) nu
310     format ( &
     /,2x,'CHANNEL LIST SORTED BY ENERGY',/, &
     '   N   J  IS   L   EINT(CM-1) ** NU = ',i2)
  end if
  do 330  i = 1, n
    ecm = eint(i) * econv
    if (bastst .or. clist) then
      isize = isiz(i)
      isub = (i - 1)*narray
      if (isize .gt. 12) then
        write (6, 320) i, j(i), is(i), l(i), ecm
        write (9, 320) i, j(i), is(i), l(i), ecm
320         format (4i4, f10.3)
      else
        if (isize .gt. 6) then
          write (6, 3201) i, j(i), is(i), l(i), ecm
          write (9, 3201) i, j(i), is(i), l(i), ecm
3201            format (4i4, f10.3)
        else
          write (6, 3202) i, j(i), is(i), l(i), ecm
          write (9, 3202) i, j(i), is(i), l(i), ecm
3202           format (4i4, f10.3)
        end if
      end if
    end if
330   continue
end if
!
!  now calculate coupling matrix elements
!  i counts v2 elements
!  inum counts v2 elements for given lambda
!  ilam counts number of v2 matrices
!  ij is address of given v2 element in present v2 matrix
i = 0
if (bastst.and. iprint.ge. 2) then
  write (6, 340)
  write (9, 340)
340   format (/' ILAM  LAMBDA   MU    ICOL  IROW    I    IV2    VEE')
end if
lamsum = 0
ilam = 0
do 400 iterm = 1, nterm
  lbmin = lammin(iterm)
!  if bastst = .true., then get the matrix elements of the lb=0 term
!  in the potential
  if (bastst .and. iterm .eq. 1) lbmin = 0
  do 390 lb = lbmin, lammax(iterm)
!  ilam is the index for the next term in the potential matrix
!  lb is the actual value of lambda
    ilam = ilam + 1
    mu = mproj(iterm)
    inum = 0
    ij = 0
    do 355  icol = 1, n
      do 350  irow = icol, n
        ij = ntop * (icol - 1) + irow
        lrow = l(irow)
        if (csflag) lrow = nu
        call vlmatp1 (j(irow), lrow, j(icol), l(icol), jtot, &
                     is(irow), is(icol), lb, mu, vee, &
                     csflag, irow, icol)

!  change check for nonzero v2 matrix element
!              if (vee .ne. zero) then
        if (abs(vee) .gt. 1.d-15) then

          i = i + 1
          if (i .le. nv2max) then
            inum = inum + 1
            v2(i) = vee
            iv2(i) = ij
            if (bastst.and. iprint.ge.2) then
              write (6, 345) ilam, lb, mu, icol, irow, i, iv2(i), &
                             vee
              write (9, 345) ilam, lb, mu, icol, irow, i, iv2(i), &
                             vee
345               format (i4, 3i7, 2i6, i6, g17.8)
            end if
          end if
        end if
350       continue
355     continue
    if (i .le. nv2max) lamnum(ilam) = inum
    if (bastst) then
      write (6, 370) ilam, lamnum(ilam)
      write (9, 370) ilam, lamnum(ilam)
370       format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
    end if
    lamsum = lamsum + lamnum(ilam)
390   continue
400 continue
if ( i.gt. nv2max) then
   write (6, 450) i, nv2max
   write (9, 450) i, nv2max
450    format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6, &
           ' .GT. NV2MAX=',i6,'; ABORT ***')
   stop
end if
if (clist .and. bastst) then
  write (6, 460) lamsum
  write (9, 460) lamsum
460   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS ', &
           i5)
end if
return
end
! --------------------------------------------------------------------
subroutine vlmatp1 (jp, lp, j, l, jtot, isp, is, lambda, mu, &
                   v, csflag, indp, ind)
!  subroutine to calculate v-lambda matrices for close-coupled or coupled-states
!  treatment of collisions of an asymmetric top with an atom
!  the angular dependence of the (lambda,mu) term in the potential is given by
!          Y(lambda,mu) + (-1)^mu * Y(lambda,-mu)
!  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
!  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
!  the expressions for the full matrix elements for definite-parity symmetric
!  top wavefunctions are given in eqs. (46-48) and (50) of the same article.
!  note, that in this article the bra indices (left side of matrix elements)
!  are primed, while in the conventions of the present subroutine the bra
!  indices are primed and the ket indices (right-hand side of matrix elements),
!  unprimed.
!
!  author of vlmstp (vlmstp1 is a copy of vlmstp):  millard alexander
!  revised from vlmstp subr for asymmetric top levels by paul dagdigian
!  current revision date:  3-sep-2009 by pjd
!  -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element (bra)
!    lp:       orbital angular momentum of left side of matrix element (bra)
!    j:        rotational quantum number of right side of matrix element (ket)
!    l:        orbital angular momentum of right side of matrix element (ket)
!    jtot:     total angular momentum
!    isp:      symmetry index of bra
!    is:       symmetry index of ket
!    lambda:   order of legendre term in expansion of potential
!    mu:       absolute value of index of legendre term in expansion of
!              potential
!    v:        on return, contains desired matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index)
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum and l and lp correspond
!                to the orbital angular momenta
!    indp:     number of left-side level in channel list
!    ind:      number of right-side level in channel list
!  subroutines called:
!    xf3j, xf6j, prmatp
!  -----------------------------------------------------------------------
use mod_coatpi, only: narray
use mod_coatpr, only: c 
implicit double precision (a-h,o-z)
logical csflag
!
data one, zero, two / 1.0d0, 0.0d0, 2.0d0 /
v = zero
!  determine k_prolate and symmetry index of the asymmetric top wave functions
kp = abs(isp)
iepsp = sign(1, isp)
k = abs(is)
ieps = sign(1, is)
!  check phase relations first
iphase = ieps * iepsp * (-1) ** (jp + j + lambda + mu)
if (iphase .eq. -1) return
!  first determine size and symmetry block of wave function expansions
!  ket
if (k .eq. 2*(k/2)) then
  if (ieps .eq. 1) then
!  E+ block
    isize = (j + 2)/2
    kmin = 0
  else
!  E- block
    isize = j/2
    kmin = 2
  end if
else
!  O+, O- blocks
    isize = (j + 1)/2
    kmin = 1
end if
!  bra
if (kp .eq. 2*(kp/2)) then
  if (iepsp .eq. 1) then
!  E+ block
    isizp = (jp + 2)/2
    kminp = 0
  else
!  E- block
    isizp = jp/2
    kminp = 2
  end if
else
!  O+,O- blocks
    isizp = (jp + 1)/2
    kminp = 1
end if
!
v = 0.d0
do 100 mm = 1,isize
do 100 nn = 1,isizp
  isubp = (indp - 1)*narray + nn
  isub = (ind - 1)*narray + mm
  kbas = kmin + (mm - 1)*2
  kbasp = kminp + (nn - 1)*2
  kdif = kbas - kbasp
  if (iabs(kdif) .eq. mu) then
    omeg = one
    if (csflag) then
!  the signed value of mu in the cs matrix elements is given by the
!  3j symbol (j' lambda j / -k' mu k) so that mu = k' - k = - kdif
      musign = - kdif
!  the multiplicative factor is given by Eq. (52) of S. Green, j. chem. phys.
!  64, 3463 (1976)
      if (kdif .gt. 0) omeg =  (-1) ** mu
    else if (.not. csflag) then
!  the signed value of mu in the cc matrix elements is given by the
!  3j symbol (j' j lambda / k' -k mu) so that mu = k - k'= kdif
      musign = kdif
!  the multiplicative factor is given by Eq. (48) of S. Green, j. chem. phys.
!  64, 3463 (1976)
      if (kdif .lt. 0)  omeg =  (-1) ** mu
    end if
!  contribution from (jp, kp, lp / Y(lambda, mu) / j, k, l), that is
!  the first term in Eq. (46) of S. Green, j. chem. phys. 64, 3463 (1976)
    call prmatp1 (jp, lp, j, l, jtot, kbasp, kbas, lambda, musign, &
               vprm, csflag)
    v = v + omeg * vprm * c(isubp) * c(isub)
  end if
  if (kbas + kbasp .eq. mu) then
!  n.b. for kbas = 0 and/or kbasp = 0, we recompute the same primitive matrix
!  element (here we follow MOLSCAT, although this might be somewhat inefficient
!  this is the second term in Eq. (46) of S. Green, j. chem. phys. 64, 3463
    if (.not.csflag) then
!  cc contribution from (jp, -kp, lp / Y(lambda, mu) / j, k, l)
      call prmatp1 (jp, lp, j, l, jtot, -kbasp, kbas, lambda, mu, &
                 vprm, csflag)
      v = v + vprm * iepsp * c(isubp) * c(isub)
    else if (csflag) then
!  cs contribution from (jp, kp, lp / Y(lambda, mu) / j, -k, l)
      call prmatp1 (jp, lp, j, l, jtot, kbasp, -kbas, lambda, mu, &
                 vprm, csflag)
      v = v + vprm * ieps * c(isubp) * c(isub)
    end if
  end if
100 continue
return
end
! ----------------------------------------------------------------------
subroutine prmatp1 (jp, lp, j, l, jtot, kp, k, lambda, mu, &
                   vprm, csflag)
!  subroutine to calculate primitive v-lambda matrix elements for close-coupled
!  treatment of collisions of a symmetric top with an atom
!  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
!  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
!  note, that in this article the bra indices are unprimed and the ket indices
!  primed, while in the conventions of the present subroutine the bra indices
!  are primed and the ket indices, unprimed.
!
!  duplicate of prmstp
!  author:  millard alexander
!  current revision date:  27-mar-90
!  -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element (bra)
!    lp:       orbital angular momentum of left side of matrix element (bra)
!    j:        rotational quantum number of right side of matrix element (ket)
!    l:        orbital angular momentum of right side of matrix element (ket)
!    jtot:     total angular momentum
!    kp:       k quantum number of bra
!    k:        k quantum number of ket
!    lambda:   order of legendre term in expansion of potential
!    mu:       index of legendre term in expansion of potential
!    vrpm:     on return, contains primitive matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index)
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum
!  subroutines called:
!     xf3j, xf6j
!  -----------------------------------------------------------------------
implicit double precision (a-h,o-z)
logical csflag
data half, one, two, zero, four / 0.5d0, 1.d0, 2.d0, 0.d0, 4.d0/
data pi /3.1415926536d0 /
vprm = zero
xjp = jp
xj = j
xkp = kp
xk = k
xjtot = jtot
if (csflag) then
  nu = lp
  xnu = nu
end if
xlp = float(lp)
xl = float(l)
xlamda = float(lambda)
xmu = float(mu)
xnorm = (two * xjp + one) * (two * xj + one) * &
        (two * lambda + one) / (four * pi)
!  special normalization for k and/or kp = 0
!  see Eq. (46) of S. Green, j. chem. phys. 64, 3463 (1976)
if (k .eq. 0) xnorm = xnorm * half
if (kp .eq. 0) xnorm = xnorm * half
!  the desired matrix element is constructed in x
if (.not. csflag) then
!  here for cc matrix elements
  x = xf3j (xlp, xlamda, xl, zero, zero, zero)
  if (x .eq. zero) return
  x = x * xf3j (xjp, xj, xlamda, xkp, - xk, xmu)
  if (x .eq. zero) return
  x = x * xf6j (xj, xl, xjtot, xlp, xjp, xlamda)
  if (x .eq. zero) return
  iphase = jp + j + kp - jtot
  xnorm = xnorm * (two * lp + one) * (two * l + one)
else if (csflag) then
!  here for cs matrix elements
  iphase = - k - nu
  x = xf3j (xjp, xlamda, xj, -xnu, zero, xnu)
  if (x .eq. zero) return
  x = x * xf3j (xjp, xlamda, xj, -xkp, xmu, xk)
end if
vprm = ( (-1) ** iphase) * x * sqrt(xnorm)
return
end
! ----------------------------------------------------------------------
double precision function rotham1(ji, ki, jf, kf)
!
!  subroutine to compute matrix elements of the asymmmetric top hamiltionian
!  in a prolate (case Ia) basis between unsymmetrized basis functions
!  (ji,ki) and (jf,kf)
!
!  here, the body-frame quantization axis is along the C2 axis
!  (b inertial axis) of the symmetrical molecule
!
!  modification of rotham subr in hibaastp.f
!
!  author:  paul dagdigian
!  current revision date:  18-sep-2017
!  -----------------------------------------------------------------------
use mod_cosysr, only: isrcod, junkr, rspar
implicit double precision (a-h,o-z)

real(8), pointer :: arot, brot, crot, emax
arot=>rspar(1); brot=>rspar(2); crot=>rspar(3); emax=>rspar(4)

aa = brot
bb = arot
cc = crot
bpc = (bb + cc)*0.5d0
bmc = (bb - cc)*0.25d0
if (ji .ne. jf) goto 900
fjj1 = ji*(ji + 1.d0)
fk = ki
if (ki .ne. kf) goto 10
!
!     diagonal term
rotham1 = bpc*(fjj1 - fk**2) + aa*fk**2
goto 1000
!
!     off-diagonal terms
10 if (kf .ne. (ki + 2)) goto 20
rotham1 = bmc*sqrt((fjj1 - fk*(fk + 1.d0)) &
  *(fjj1 - (fk + 1.d0)*(fk + 2.d0)))
goto 1000
20 if (kf .ne. (ki - 2)) goto 900
rotham1 = bmc*sqrt((fjj1 - fk*(fk - 1.d0)) &
  *(fjj1 - (fk - 1.d0)*(fk - 2.d0)))
goto 1000
900 rotham1 = 0.d0
1000 continue
return
end
! -----------------------------------------------------------------------
subroutine syastp1 (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for C2v asymmetric top
!      + atom scattering
!
!  this basis set differs for baastp in the following way:
!    the body-frame quantization axis is along the C2 axis of a C2v symmetry
!    molecule.  This convention follows that employed in MOLSCAT
!
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
!  modified from systp subroutine
!  current revision date: 19-sep-2017 by p. dagdigian
!  -----------------------------------------------------------------------
!  variables in common bloc /cosysr/
!    isrcod:   total number of real system dependent variables
!    arot:     A rotational constant
!    brot:     B rotational constant
!    crot:     C rotational constant
!    emax:     the maximum rotational energy (in cm-1) for a channel to be
!              included in the basis
!  variables in common block /cosysi/
!    nscode:   total number of variable names which are passed to HINPUT
!              nscod must equal isrcod + isicod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    numpot:   the number of the potential used, this variable is passed
!              to the pot subroutine
!    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
!              para states will be included if iop=1 and only ortho states if
!              iop=-1
!    jmax:     the maximum rotational quantum number for the asymmetric top
!  variable in common  /cosys/
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters.  Note that the ordering
!             of the variable names in scod must correspond to the ordering
!             of the variable names in cosysi followed by the ordering of
!             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_conlam, only: nlam
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
implicit none
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: i, j, l, lc
logical existf
integer icod, ircod
character*1 dot
character*(*) fname
character*60 line, filnam, potfil, filnm1
parameter (icod=4, ircod=4)
#include "common/parbas.F90"
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

integer, pointer :: nterm, numpot, iop, jmax
real(8), pointer :: arot, brot, crot, emax
nterm=>ispar(1); numpot=>ispar(2); iop=>ispar(3); jmax=>ispar(4)
arot=>rspar(1); brot=>rspar(2); crot=>rspar(3); emax=>rspar(4)

scod(1)='NTERM'
scod(2)='NUMPOT'
scod(3)='IOP'
scod(4)='JMAX'
scod(5)='AROT'
scod(6)='BROT'
scod(7)='CROT'
scod(8)='EMAX'
scod(9)='LAMMIN'
scod(10)='LAMMAX'
scod(11)='MPROJ'
nscode = icod + ircod + 3
isicod = icod
isrcod = ircod
irpot = 1
!  set default values for asymmetric top scattering
!  (symmetric BA2 and A2BC molecules)
nterm = 4
if (iread .eq. 0) then

  mproj(1) = 0
  mproj(2) = 2
  mproj(3) = 4
  mproj(4) = 6

  lammin(1) = 1
  lammin(2) = 2
  lammin(3) = 4
  lammin(4) = 6

  lammax(1) = 8
  lammax(2) = 8
  lammax(3) = 8
  lammax(4) = 8

  iop = 1
  jmax = 3
  niout=5
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
read (8, *, err=80) iop
!  line 19
read (8, *, err=80) jmax, emax
!  line 20
read (8, *, err=80) arot, brot, crot
if(.not. readpt .or. iread .eq. 0) then
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
entry ptrastp1 (fname, readpt)
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
entry savastp1 (readpt)
!  save input parameters for asymmetric top + atom scattering
!  the order of the write statements should be identical to the read statement
!  above. for consistency with the data file written by gendat, format
!  statements should reserve the first 30 spaces for data, spaces 31-33 should
!  be left blank, and the names of the variables should be printed in spaces
!  34-80
!  line 18:
write (8, 220) iop
220 format (i4, 26x,'   iop')
!  line 20
write (8, 230) jmax, emax
230 format (i4, 3x, g12.5, 14x, 'jmax, emax')
!  line 21
write (8, 250) arot, brot, crot
250 format(3f9.4, 6x, 'arot, brot, crot')
write (8, 60) potfil
return
end
! -----------------------------------eof--------------------------------
