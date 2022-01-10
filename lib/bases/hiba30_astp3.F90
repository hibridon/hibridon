! syastp3 (savastp3/ptrastp3) defines, saves variables and reads         *
!                  potential for C2v asym top - lin mo;ecule scattering  *
! --------------------------------------------------------------------
!     This module contains (explictly) the number of terms and their
!     indices in the expansion of the PES.  Its contents should be
!     filled in the pot routine.
!
!     This module replaces lammin, lammax, mproj in hibridon.
!
module mod_asymln
implicit none
!
type lm_type
integer :: l1, m1, l2, ltot
end type lm_type
!
type(lm_type), dimension(:), allocatable :: lms
end module mod_asymln
! ----------------------------------------------------------------------
subroutine baastp3 (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, etemp, fjtemp, fktemp, fistmp, &
                  rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, &
                  jlpar, n, nmax, ntop)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of an asymmetric top molecule of C2v or Cs symmetry with a linear molecule
!
!  modification of hibaastp basis routine
!
!  here for a molecule with C2v symmetry, the body-frame z axis lies along 
!  the C2 axis of the molecule,  following green's convention.  
!  for a molecule with only Cs symmetry, the molecule-frame z axis lies 
!  along the inertial a axis of the asymmetric top 
!
!  author:  paul dagdigian
!  current revision:  11-feb-2020 by p.dagdigian
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains combined rotational quantum numbers for each
!              channel.  in terms of the rotational quantum numbers of each
!              molecule we have:  j = 10 * j1 + j2
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
!              RCUT IS NOT USED IN THIS BASIS ROUTINE
!    jtot:     total angular momentum
!              in cc calculation jtot is the total angular momentum
!              in cs calculation jtot is the l-bar quantum number
!    flaghf:   if .true., then system has half-integer spin
!              if .false., then system has integer spin
!    flagsu:   if .true., then molecule-surface collisons
!    csflag:   if .true., then coupled-states calculation
!              if .false., the close-coupled calculation
!              CSFLAG IS IGNORED IN THIS BASIS ROUTINE.  ONLY CLOSE
!              COUPLED CALCULATIONS CARRIED OUT
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true., then the molecule posseses interchange symmetry
!              (e.g. H2O, CH2, H2CO), so that only the ortho or para levels will be
!              included depending on the value of the parameter iop in common
!              /cosysi/ (see below)
!              for a moleculw with only Cs symmetry, set ihomo to .false.
!              in this case, set ihomo to .false.
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
!    b2rot:    otational constant for the linear molecule
!  variables in common block /cosysi/
!    nscode:   total number of variable names which are passed to HINPUT
!              nscod must equal isrcod + isicod + 3
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different associated legendre terms in
!              expansion of potential
!    numpot:   the number of the potential used, this variable is passed
!              to the pot subroutine
!    jmax:     the maximum rotational angular momentum for the asymmetric top
!    iop:      for a molecule of C2v symmetry, this parameter is the ortho/para 
!              label for the asymmetric top states.  para states only
!              included if iop=1, and only ortho states if iop=-1.  in this
!              case, the body-frame z axis should be along the C2 axis.
!              IMPORTANT:  for molecules with a plane of symmetry only 
!              (Cs symmetry), iop should be set to zero so that 
!              the a axis is the body-frame z axis.
!              With these considerations, the proper rotatioonal Hamiltonian
!              will be employed.
!    j2min:    minimum rotational angular momentum of molecule 2
!    j2max:    maximum rotational angular momentum of molecule 2
!    ipotsy2:  symmetry of potential.  set to 2 for homonuclear
!              molecule 2, set to 1 for heteronuclear molecule 2
!  variables in common block /coered/
!    ered:      collision energy in atomic units (hartrees)
!    rmu:       collision reduced mass in atomic units
!               (mass of electron = 1)
!  variables in common block /coconv/
!   econv:      conversion factor from cm-1 to hartrees
!   xmconv:     converson factor from amu to atomic units
!  subroutines called:
!   vastp3:     returns angular coupling coefficient for particular
!               choice of channel index
!   rotham1:    computes matrix elements of asymmetric top hamiltonian
!               here, the body-frame quantization axis is along the C2 axis
!               (usuallu the b inertial axis) of the molecule.  This subroutine
!               should be employed for molecules with C2v symmetry.  (set iop
!               equal to 1 or -1 for para and ortho nuclear spin modifications)
!   rotham:     computes matrix elements of asymmetric top hamiltonian
!               here, the body-frame quantization axis is along the a inertial 
!               axis of the molecule.  This subroutine should be used for
!               molecules with Cs symmetry.  (set iop equal to zero in this case)
! --------------------------------------------------------------------
use mod_asymln
use mod_cov2, only: nv2max, junkv => ndummy, v2
use mod_coiv2, only: iv2
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coj12, only: j12
use mod_coatpi, only: narray, isiz
use mod_coatpr, only: c
use mod_coatp1, only: ctemp
use mod_coatp2, only: chold
use mod_coatp3, only: isizh
use mod_conlam, only: nlam, nlammx, lamnum
use constants, only: econv, xmconv
implicit double precision (a-h,o-z)
logical flaghf, csflag, clist, flagsu, ihomo, bastst
#include "common/parbas.F90"
#include "common/parbasl.F90"
common /cosysi/ nscode, isicod, nterm, numpot, jmax, iop, &
  j2min, j2max, ipotsy2
common /cosysr/ isrcod, junkr, arot, brot, crot, emax, b2rot
common /coipar/ iiipar(9), iprint
common /coered/ ered, rmu
dimension j(1), l(1), is(1), jhold(1), ehold(1), &
          ishold(1), etemp(1), fjtemp(1), fktemp(1), &
          fistmp(1)
!  scratch arrays for computing asymmetric top energies and wave fns.
dimension e(narray,narray), eig(narray), vec(narray,narray), &
  sc1(narray), sc2(narray), work(1000), kp(1000), ko(1000), &
  j2rot(1000), e2rot(1000)
!
zero = 0.d0
two = 2.d0
!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (6, 5)
  write (9, 5)
5  format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
  stop
end if
if (csflag) then
  write (6, 6)
  write (9, 6)
6 format &
   ('  *** CSFLAG = .TRUE. FOR ASYMMETRIC TOP', &
    ' - LINEAR MOLECULE COLLISIONS NOT IMPLEMENTED.  ABORT ***')
  stop
end if
if (flagsu) then
  write (6, 10)
  write (9, 10)
10   format ('  *** FLAGSU = .TRUE. FOR ASYMMETRIC TOP', &
     ' + LINEAR MOLECULE COLLISIONS NOT IMPLEMENTED.  ABORT ***')
  call exit
end if
!
!  check that ipotsy2 is either 1 or 2
if (ipotsy2.lt.1 .or. ipotsy2.gt.2) then
  write(6,27) ipotsy2
  write(9,27) ipotsy2
27   format(' *** IPOTSY2 .NE. 1 .OR. 2:  ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
!
!  check that iop is either 1 or -1 if ihomo=T
if (ihomo) then
  if ( .not. (iop.eq.1 .or. iop.eq.-1) ) then
    write(6,29) 
    write(9,29)
29     format(' *** IOP .ne. 1 .or. -1 with IHOMO = true:' &
      '   ABORT ***')  
    if (bastst) then
      return
    else
      call exit
    end if
  end if
end if
!
if (bastst) then
  if (ihomo) then
    write (6,80) rmu * xmconv, arot, brot, crot, &
      b2rot, iop, ered * econv, jtot, jlpar
    write (9,80) rmu * xmconv, arot, brot, crot, &
      b2rot, iop, ered * econv, jtot, jlpar
80     format(/,' **  CC ASYMMETRIC TOP + LINEAR MOLECULE **', &
      /,'     RMU=', f9.4,'  AROT=', f7.3, '  BROT=',f7.3, &
      '  CROT=',f7.3,'  B2ROT=',f7.3/'     O/P=',i2, &
      '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
  else
    write (6,85) rmu * xmconv, arot, brot, crot, &
      b2rot, ered * econv, jtot, jlpar
    write (9,85) rmu * xmconv, arot, brot, crot, &
      b2rot, ered * econv, jtot, jlpar
85     format(/,' **  CC ASYMMETRIC TOP + LINEAR MOLECULE **', &
      /,'     RMU=', f9.4,'  AROT=', f7.3, '  BROT=',f7.3, &
      '  CROT=',f7.3,'  B2ROT=',f7.3/'     E=', f7.2, &
      '  JTOT=', i4, 2x,' JLPAR=', i2)
  end if
end if
!
!  first set up list of all j(kp,ko) states included in basis
!  for the asymmetric molecule
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
  if (iop .ne. 0) then
    e(1,1) = rotham1(ji,0,ji,0)
  else
    e(1,1) = rotham(ji,0,ji,0)
  end if
  if (isize .gt. 1) then
    do mm = 2, isize
      kk = 2*mm - 2
      if (iop .ne. 0) then
        e(mm,mm) = rotham1(ji, kk, ji, kk)
        e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
      else
        e(mm,mm) = rotham(ji, kk, ji, kk)
        e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
       end if
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
    if (iop .ne. 0) then
      e(1,1) = rotham1(ji,2,ji,2)
    else
      e(1,1) = rotham(ji,2,ji,2)
    end if
    if (isize .gt. 1) then
    do mm = 2, isize
      kk = 2*mm
      if (iop .ne. 0) then
        e(mm,mm) = rotham1(ji, kk, ji, kk)
        e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
      else
        e(mm,mm) = rotham(ji, kk, ji, kk)
        e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
      end if
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
    if (iop .ne. 0) then
      e(1,1) = rotham1(ji,1,ji,1) + rotham1(ji,1,ji,-1)
    else
      e(1,1) = rotham(ji,1,ji,1) + rotham(ji,1,ji,-1)
    end if
    if (isize .gt. 1) then
    do mm = 2, isize
      kk = 2*mm - 1
      if (iop .ne. 0) then
        e(mm,mm) = rotham1(ji, kk, ji, kk)
        e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
      else
        e(mm,mm) = rotham(ji, kk, ji, kk)
        e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
      end if
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
    if (iop .ne. 0) then
      e(1,1) = rotham1(ji,1,ji,1) - rotham1(ji,1,ji,-1)
    else
      e(1,1) = rotham(ji,1,ji,1) - rotham(ji,1,ji,-1)
    end if
    if (isize .gt. 1) then
      do mm = 2, isize
        kk = 2*mm - 1
        if (iop .ne. 0) then
          e(mm,mm) = rotham1(ji, kk, ji, kk)
          e(mm,mm-1) = rotham1(ji,kk,ji,kk-2)
        else
          e(mm,mm) = rotham(ji, kk, ji, kk)
          e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
        end if
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
!  calculate energy levels for linear molecule collision partner
!  subtract energy of j2min level
i2num = 0
do j2 = j2min, j2max, ipotsy2
  i2num = i2num + 1
  j2rot(i2num) = j2
  e2rot(i2num) = b2rot * (j2 * (j2 + 1) - j2min * (j2min + 1.d0))
end do
!
!  combine asymmetric top and linear molecule levels
nlevel = 0
do i2 = 1, i2num
  do 170 i1 = 1, nlist
    ki = fktemp(i1)
    ji = fjtemp(i1)
    isi = fistmp(i1)
!  delete state if energy > emax
    if (etemp(i1) .gt. emax) go to 170
!  check to see if (ki/isi) corresponds to an allowed ortho/para level for
!  a symmetric molecule.  the statements below correspond to BA2 and A2BC
!  type molecules like H2O, CH2, and H2CO
!
!  first compute kp and ko projection quantum numbers
    if (ji .eq. 2*(ji/2)) then
!  even j
      if (isi .eq. 1) then
!  isi = +1
        kp1 = 2 * ceiling(float(ki) * 0.5d0)
        ko1 = ji - ki
        goto 9900
      else
!  isi = -1
        kp1 = 2 * ceiling(float(ki - 2) * 0.5d0) + 1
        ko1 = ji + 1 - ki
        goto 9900
      end if
    else
!  odd j
      if (isi .eq. 1) then
!  isi = +1
        kp1 = 2 * ceiling(float(ki - 1) * 0.5d0) + 1
        ko1 = ji - ki
        goto 9900
      else
!  isi = -1
        kp1 = 2 * ceiling(float(ki - 1) * 0.5d0)
        ko1 = ji + 1 - ki
        goto 9900
      end if
    end if
9900     continue
!
!  for molecules with C2v symmetry, gather either ortho or para levels
    if (ihomo) then
      if ((kp1+ko1) .eq. 2*((kp1+ko1)/2)) then
        if (iop .eq. -1) goto 170
      else
        if (iop .eq. 1) goto 170
      end if
    end if
!
    nlevel = nlevel + 1
    ehold(nlevel) = (etemp(i1) + e2rot(i2)) / econv
    jhold(nlevel) = 10 * ji + j2rot(i2)
    ishold(nlevel) = isi * ki
    kp(nlevel) = kp1
    ko(nlevel) = ko1
!  also move e.fn coeffs
!  first, determine size of wave function expansion
    if (ki .eq. 2*(ki/2)) then
      if (isi .eq. 1) then
!  E+ block
        isize = (ji + 2)/2
      else
!  E- block
        isize = ji/2
      end if
    else
!  O+, O- blocks
      isize = (ji + 1)/2
    end if
    isizh(nlevel) = isize
    do mm = 1, isize
      isub = (nlevel - 1)*narray + mm
      isub1 = (i1 - 1)*narray + mm
      chold(isub) = ctemp(isub1)
    end do
170   continue
end do
!
!  now sort this list in terms of increasing energy
if (nlevel .gt. 1) then
  do 120 i1 = 1, nlevel - 1
    esave = ehold(i1)
    do 115 i2 = i1 + 1, nlevel
      if (ehold(i2) .lt. esave) then
!  state i2 has a lower energy than state i1, switch them
        esave = ehold(i2)
        ehold(i2) = ehold(i1)
        ehold(i1) = esave
        jsave = jhold(i2)
        jhold(i2) = jhold(i1)
        jhold(i1) = jsave
        isave = ishold(i2)
        ishold(i2) = ishold(i1)
        ishold(i1) = isave
        ksave = kp(i2)
        kp(i2) = kp(i1)
        kp(i1) = ksave
        ksave = ko(i2)
        ko(i2) = ko(i1)
        ko(i1) = ksave
        isz = isizh(i2)
        isizh(i2) = isizh(i1)
        isizh(i1) = isz
!  also move e.fn coeffs
        do mm = 1, narray
          isub1 = (i1 - 1)*narray + mm
          isub2 = (i2 - 1)*narray + mm
          scr = chold(isub2)
          chold(isub2) = chold(isub1)
          chold(isub1) = scr
        end do
      end if
115     continue
120   continue
end if
!
!  print level list if bastst = .true.
if (bastst) then
  write(6,130)
  write(9,130)
130   format(/2x,'ASYMMETRIC TOP + LINEAR MOLECULE LEVELS'/ &
    1x,65('-') / &
    '   #  J1  IS  KP  KO  J2  Eint(cm-1)  Coeffs'/1x,65('-'))
  do i = 1, nlevel
    ecm = ehold(i) * econv
    j1 = jhold(i)/10
    j2 = mod(jhold(i),10)
    write (6, 135) i, j1, ishold(i), kp(i), ko(i), j2, ecm, &
      (chold(nn),nn=(i - 1)*narray + 1, &
      (i - 1)*narray + isizh(i))
135     format (6i4, 1x, f9.3,10f8.4/34x,10f8.4/34x,10f8.4)
  end do
end if
!
!  determine number of levels which are open
nlevop = 0
do 250  i = 1, nlevel
  if (ehold(i) .le. ered) then
    nlevop = nlevop + 1
  end if
250 continue
!
!  set up channel for cc calculations
!  first calculate range of orbital angular
!  momentum quantum numbers allowed for this state
!  determine parity of molecular state [Eq. (A21) of S. Green, 
!  J. Chem. Phys. 73, 2740 (1980) and Townes and Schawlow, 
!  Microwave Spectroscopy, Eq. (3-27), p. 64.].
n = 0
do i = 1, nlevel
  j1 = jhold(i)/10
  j2 = mod(jhold(i),10)
  ki = abs(ishold(i))
  iss = sign(1, ishold(i))
  ipar = iss * (-1) ** (j1 + j2 + ki)
  j12mn = abs(j1 - j2)
  j12mx = j1 + j2
  do j12i = j12mn, j12mx
    lmin = abs(jtot - j12i)
    lmax = jtot + j12i
    do li = lmin, lmax
!  check to see if this channel has the desired parity, if so, include it
      ix = ipar * (-1) ** (li - jtot)
      if (ix .eq. jlpar) then
        n = n + 1
        if (n .gt. nmax) then
          write (6, 150) n, nmax
          write (9, 150) n, nmax
150           format(/' *** NCHANNELS=', i5, &
              ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
          stop
        end if
        eint(n) = ehold(i)
        j(n) = jhold(i)
        is(n) = ishold(i)
        l(n) = li
        j12(n) = j12i
        cent(n) = li * (li + 1)
!  move e.fn also
        isiz(n) = isizh(nlevel)
        do mm = 1, narray
          isub = (n - 1)*narray + mm
          isub1 = (i - 1)*narray + mm
          c(isub) = chold(isub1)
        end do
      end if
    end do
  end do
end do
!
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
            '; ABORT **')
    call exit
  endif
end if
!  now list channels if requested
if (bastst) then
  write (6, 305)
  write (9,305)
305   format(/2x,'CC CHANNEL BASIS'/1x,50('-')/ &
    '   #  J1  IS  J2  J12  L     EINT(CM-1)'/1x,50('-'))
  do 330  i = 1, n
    ecm = eint(i) * econv
    j1 = j(i)/10
    j2 = mod(j(i),10)
    if (bastst) then
      isize = isiz(i)
      isub = (i - 1)*narray
      if (isize .gt. 12) then
        write (6, 320) i, j1, is(i), j2, j12(i), l(i), ecm
        write (9, 320) i, j1, is(i), j2, j12(i), l(i), ecm
320         format (6i4, f12.3)
      else
        if (isize .gt. 6) then
          write (6, 3201) i, j1, is(i), j2, j12(i), l(i), ecm
          write (9, 3201) i, j1, is(i), j2, j12(i), l(i), ecm
3201            format (6i4, f12.3)
        else
          write (6, 3202) i, j1, is(i), j2, j12(i), l(i), ecm
          write (9, 3202) i, j1, is(i), j2, j12(i), l(i), ecm
3202           format (6i4, f12.3)
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
lamsum = 0
if (bastst) then
  write (6, 340)
  write (9, 340) 
340   format (/' ILAM  L1   M1   L2  LTOT    ICOL  IROW       I' &
    '       IV2       VEE')
end if
i = 0
ilam = 0
do 400 ilam = 1, nlam
  ll1 = lms(ilam)%l1
  mm1 = lms(ilam)%m1
  ll2 = lms(ilam)%l2
  lltot = lms(ilam)%ltot
  inum = 0
  do 355 icol = 1, n
    j1c = j(icol)/10
    j2c = mod(j(icol),10)
    j12c = j12(icol)
    lc = l(icol)
    do 350 irow = icol, n
      ij = ntop * (icol - 1) + irow
      j1r = j(irow)/10
      j2r = mod(j(irow),10)
      j12r = j12(irow)
      lr = l(irow)
      call vastp3(j1r,j2r,j12r,lr,j1c,j2c,j12c,lc, &
          is(irow),is(icol),irow,icol,jtot, &
          ll1,mm1,ll2,lltot,vee)
      if (abs(vee) .gt. 1.d-10) then
        i = i + 1
        if (i .le. nv2max) then
          inum = inum + 1
          v2(i) = vee
          iv2(i) = ij
          if (bastst) then
            write (6, 290) ilam, ll1, mm1, ll2, lltot, &
              icol, irow, i, iv2(i), vee
            write (9, 290) ilam, ll1, mm1, ll2, lltot, &
              icol, irow, i, iv2(i), vee
290             format (i4, 4i5, 2x, 2i6, i10, i9, e20.7)
          end if
        end if
      end if
350     continue
355   continue
  if (i .le. nv2max) lamnum(ilam) = inum
  if (bastst) then
    write (6, 370) ilam, lamnum(ilam)
    write (9, 370) ilam, lamnum(ilam)
370     format ('ILAM=',i4,' LAMNUM(ILAM) = ',i7)
  end if
  lamsum = lamsum + lamnum(ilam)
400 continue
if (i .gt. nv2max) then
  write (6, 410) i, nv2max
  write (6, 410) i, nv2max
410   format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i10, &
      ' .GT. NV2MAX=',i10,'; ABORT ***')
end if
if (bastst) then
  write (6, 420) lamsum
  write (9, 420) lamsum
420   format (' *** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', &
      i10)
end if
return
end
! --------------------------------------------------------------------
subroutine vastp3(j1r,j2r,j12r,lr,j1c,j2c,j12c,lc, &
  isr,isc,irow,icol,jtot,l1,m1,l2,ltot,vee)
!
!  subroutine to calculate v-lambda matrices for close-coupled 
!  treatment of collisions of an asymmetric top with a diatomic molecule
!  the angular dependence of the terms in the potential is given by
!  Valiron et al. [JCP 129, 134306 (2008)].  the angular expansion
!  coefficients have been muitiplied by the normalization factor in
!  the subroutine loapot
!
!  an expression for the matrix element in the space-fixed basis is given
!  by Phillips et al. [JCP 102, 6024 (1995)].
!
!  written by paul dagdigian
!  current revision date:  30-oct-2019 by p.dagdigian
!  -----------------------------------------------------------------------
!  variables in call list:
!    j1c:      rotational quantum number of #1 of left side of matrix element (bra)
!    j2c:      rotational quantum number of #2 of left side of matrix element (bra)
!    j12c:     value of j12 for right side of matrix element (bra)
!    lc:       orbital angular momentum of left side of matrix element (bra)
!    j1r:      rotational quantum number of #1 of right side of matrix element (ket)
!    j2r:      rotational quantum number of #2 of right side of matrix element (ket)
!    j12r:     value of j12 for right side of matrix element (ket)
!    lr:       orbital angular momentum of right side of matrix element (ket)
!    isr:      label for bra
!    isc:      label for ket
!    irow:     number of left-side level in channel list
!    icol:     number of right-side level in channel list
!    jtot:     total angular momentum
!    l1,m1,l2,ltot:   parameters for the given angular expansion coefficient       
!    vee:      on return, contains desired matrix element
!  subroutines called:
!    xf3j, xf6j, xf9j
!  -----------------------------------------------------------------------
use mod_coatpi, only: narray, isiz
use mod_coatpr, only: c
implicit double precision (a-h,o-z)
data one, zero, two / 1.0d0, 0.0d0, 2.0d0 /
data fourpi / 12.566370614d0 /, sqr2 / 1.414213562373095d0 /
!
vee = zero
xjtot = jtot
xj1r = j1r
xj1c = j1c
xj2r = j2r
xj2c = j2c
xj12r = j12r
xj12c = j12c
xlr = lr
xlc = lc
iepr = sign(1,isr)
if (isr .eq. 0) iepr = 1
iepc = sign(1,isc)
if (isc .eq. 0) iepc = 1
!  get quantum numbers for this angular term
xl1 = l1
xm1 = m1
xl2 = l2
xltot = ltot
!  check phase relations first
iphase = iepr * iepc * (-1) ** (j1r + j1c + m1 + l2 + ltot)
if (iphase .eq. -1) return
x1 = xf3j(xlr, xlc, xltot, zero, zero, zero)
if (x1 .eq. 0.d0) return
x2 = xf3j(xj2r, xl2, xj2c, zero, zero, zero)
if (x2 .eq. 0.d0) return
x3 = xf6j(xlr, xlc, xltot, xj12c, xj12r, xjtot)
if (x3 .eq.0.d0) return
x4 = xf9j(xj12c, xj12r, xltot, xj1c, xj1r, xl1, &
  xj2c, xj2r, xl2)
if (x4 .eq. 0.d0) return
!  first determine size and symmetry block of wave function expansions
!  ket
if (abs(isc) .eq. 2*(abs(isc)/2)) then
  if (iepc .eq. 1) then
!  E+ block
    isize = (j1c + 2)/2
    kmin = 0
  else
!  E- block
    isize = j1c/2
    kmin = 2
  end if
else
!  O+, O- blocks
    isize = (j1c + 1)/2
    kmin = 1
end if
!  bra
if (abs(isr) .eq. 2*(abs(isr)/2)) then
  if (iepr .eq. 1) then
!  E+ block
    isizp = (j1r + 2)/2
    kminp = 0
  else
!  E- block
    isizp = j1r/2
    kminp = 2
  end if
else
!  O+,O- blocks
    isizp = (j1r + 1)/2
    kminp = 1
end if
!
v = zero
do 100 mm = 1,isize
  isub = (icol - 1)*narray + mm
  kbas = kmin + (mm - 1)*2
  fkbas = kbas
do 100 nn = 1,isizp
  isubp = (irow - 1)*narray + nn
  kbasp = kminp + (nn - 1)*2
  fkbasp = kbasp
  vk1 = xf3j(xj1r,xl1,xj1c, -fkbasp, xm1, fkbas)
  vk2 = xf3j(xj1r,xl1,xj1c, fkbasp, xm1, -fkbas) 
  vk3 = xf3j(xj1r,xl1,xj1c, -fkbasp, xm1, -fkbas) 
  vk4 = xf3j(xj1r,xl1,xj1c, fkbasp, xm1, fkbas) 
!  factor of 2 at end takes care of [1 + eps*epsp*(-1)^xx] term
  vk = (vk1 + iepr * iepc * vk2 + iepc * vk3 &
    + iepr * vk4) * c(isubp) * c(isub) * two
!  special normalization for k1r and/or k1c = 0
  if (kbas .eq. 0) vk = vk / sqr2
  if (kbasp .eq. 0) vk = vk / sqr2
  v = v + vk * (-1.d0) ** kbasp / two
100 continue
!
iph = xjtot - xj1c + xj2c - xj12r - xltot
ph = (-1.d0) ** iph
sqjs = sqrt((two * xj1r + one) * (two * xj1c + one) &
  * (two * xj2r + one) * (two * xj2c + one) &
  * (two * xj12r + one) * (two * xj12c + one) &
  * (two * xlr + one) * (two * xlc + one) &
  * (two * xl2 + one) * (two * xltot + one))
!
vee = ph * sqjs * x1 * x2 * x3 * x4 * v / fourpi
if (m1 .eq. 0) vee = vee / two
!  include normalization factor depending on l1 
delm1 = 1.d0
if (m1 .eq. 0) delm1 = 2.d0
fnorm = sqrt((two * xl1 + one) * delm1 / two)
vee = vee * fnorm
!
return
end
! ----------------------------------------------------------------------
!      double precision function rotham1(ji, ki, jf, kf)
!
!  subroutine to compute matrix elements of the asymmmetric top hamiltionian
!  in a prolate (case Ia) basis between unsymmetrized basis functions
!  (ji,ki) and (jf,kf)
!
!  here, the body-frame quantization axis is along the C2 axis
!  (b inertial axis) of the symmetrical molecule
!
!  THIS ROUTINE IS IN THE FILE hibaastp1.f
!
!  author:  paul dagdigian
!  current revision date:  18-sep-2017
! ----------------------------------------------------------------------
!      double precision function rotham(ji, ki, jf, kf)
!
!  subroutine to compute matrix elements of the asymmmetric top hamiltionian
!  in a prolate (case Ia) basis between unsymmetrized basis functions
!  (ji,ki) and (jf,kf)
!
!  here, the bodyn -frame quantization axis is along the a inertila axis
!  of the molecule, as is usally employed for asymmetric tops
!
!  THIS ROUTINE IS IN THE FILE hibaastp.f
!
!  author:  paul dagdigian
!  current revision date:  16-aug-2009
!  -----------------------------------------------------------------------
!  -----------------------------------------------------------------------
subroutine syastp3 (irpot, readp, iread)
!  subroutine to read in system dependent parameters for C2v asymmetric top
!      + linear molecule scattering
!
!  current revision date: 20-jun-2019 by p. dagdigian
!  -----------------------------------------------------------------------
!  variables in common block /cosysr/
!    isrcod:   total number of real system dependent variables
!    arot:     A rotational constant
!    brot:     B rotational constant
!    crot:     C rotational constant
!    b2rot:    Rotational constant for linear molecule
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
!    jmax:     the maximum rotational quantum number for the asymmetric top
!    iop:      ortho/para label for asymmetric top states.  para states only
!              included if iop=1, and only ortho states if iop=-1
!    j2min:    minimum rotational angular momentum for molecule 2
!    j2max:    maximum rotational angular momentum for molecule 2
!    ipotsy2:  symmetry of potential.  set to 2 for homonuclear
!              molecule 2, set to 1 for heteronuclear molecule 2
!  variable in common /cosys/
!    scod:    character*8 array of dimension nscode, which contains names
!             of all system dependent parameters
!  line 16:
!    filnam:  name of file containing potential parameters
!
!  subroutines called: loapot(iunit,filnam)
!  -----------------------------------------------------------------------
#include "common/parsys_mod.F90"
use mod_conlam, only: nlam
implicit double precision (a-h,o-z)
#include "common/parsys.F90"
integer irpot
logical readp
logical airyfl, airypr, logwr, swrit, t2writ, writs, wrpart, &
        partw, xsecwr, wrxsec, noprin, chlist, ipos, flaghf, &
        csflag, flagsu, rsflag, t2test, existf, logdfl, batch, &
        readpt, ihomo, bastst, twomol, lpar
character*1 dot
character*8 scod
character*(*) fname
character*60 filnam, line, potfil, filnm1
common /cosys/ scod(lencod)
common /cosysi/ nscode, isicod, nterm, numpot, jmax, iop, &
  j2min, j2max, ipotsy2
common /cosysr/ isrcod, junkr, arot, brot, crot, emax, b2rot
#include "common/parbas.F90"
common /coskip/ nskip,iskip
common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag, &
                flaghf, flagsu, ihomo, ipos, logdfl, logwr, &
                noprin, partw, readpt, rsflag, swrit, &
                t2test, t2writ, twomol, writs, wrpart, wrxsec, &
                xsecwr,lpar(3)
#include "common/comdot.F90"
save potfil
!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
isicod = 7
scod(1) = 'NTERM'
scod(2) = 'NUMPOT'
scod(3) = 'JMAX'
scod(4) = 'IOP'
scod(5) = 'J2MIN'
scod(6) = 'J2MAX'
scod(7) = 'IPOTSY2'
!  then all the system dependent real variables
!  in the same order as in the common block /cosysr/
isrcod = 5
scod(8)='AROT'
scod(9)='BROT'
scod(10)='CROT'
scod(11)='EMAX'
scod(12)='B2ROT'
nscode = isicod + isrcod + 3
! read last few lines of the input file
if (iread .eq. 0) return
read (8, *, err=80) jmax, iop, j2min, j2max, ipotsy2
read (8, *, err=80) arot, brot, crot, emax, b2rot
read (8, *, err=80) potfil
call loapot(1,potfil)
close (8)
return
! here if read error occurs
80 write(6,90)
90 format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
return
!  -----------------------------------------------------------------------
entry ptrastp3 (fname, readp)
line = fname
readp = .true.
100 if (readp) then
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
!  -----------------------------------------------------------------------
entry savastp3 (readp)
!  save input parameters for chiral asymmetric top + atom scattering
!  the order of the write statements should be identical to the read statement
!  above. for consistency with the data file written by gendat, format
!  statements should reserve the first 30 spaces for data, spaces 31-33 should
!  be left blank, and the names of the variables should be printed in spaces
!  34-80
write (8, 230) jmax, iop, j2min, j2max, ipotsy2
230 format (5i4, 3x, 'jmax,iop,j2min,j2max,ipotsy2')
write (8, 250) arot, brot, crot, emax, b2rot
250 format(5f9.4, 6x, 'arot,brot,crot,emax,b2rot')
write (8, 60) potfil
60 format (a)
return
end
! ---------------------------------eof----------------------------------
