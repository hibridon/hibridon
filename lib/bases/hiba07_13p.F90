#include "assert.h"
module mod_hiba07_13p
real(8) :: ttrans(6,6)
contains
! sy13p (sav13p/ptr13p) defines, save variables and reads                *
!                  potential for 1S / 3P atom scattering                 *
! --------------------------------------------------
subroutine ba13p (bqs, jhold, ehold, ishold, nlevel, nlevop, &
                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential
!  for collision of a singlet and triplet atom with a structureless atom
!  authors:  brigitte pouilly and millard alexander
!  current revision date:  5-apr-2003 by mha
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
!    sc1,sc2:  scratch vectors of length at least nmax
!    sc3,sc4:  scratch vectors of length at least nmax
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
!    en(1):    asymptotic energy of j=0 fine-structure level of triplet
!    en(2):    asymptotic energy of j=1 fine-structure level of triplet
!    en(3):    asymptotic energy of j=2 fine-structure level of triplet
!    en(4):    asymptotic energy of singlet state (cm-1)
!    de(4):    morse De for 3Pi, 3Sig, 1Pi, 1Sig
!    re(4):    morse re for 3Pi, 3Sig, 1Pi, 1Sig
!    be(4):    morse beta for 3Pi, 3Sig, 1Pi, 1Sig
!    rl(4):    msv RL for 3Pi, 3Sig, 1Pi, 1Sig
!    cl(4):    msv C6 for 3Pi, 3Sig, 1Pi, 1Sig
!    cmix:     mixing coefficient of J=1 singlet and triplet levels
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different types of electronic coupling terms
!              this should be 1 here
!    nstate:   number of electronic states included
!              nstate=0:   just singlet state
!              nstate=1:   just triplet state
!              nstate=2:   both singlet and triplet state
!    ipol:     =0 if no polarization state preparation
!    ipol:     =1 if polaratization state preparation
!    npot:     potential designator (1-4 currently allowed)
!  variable in module mod_conlam
!    nlam:      the number of case(a) interaction potentials actually used
!               this is :  nlam = nlam0 + nlam1
!               if only singlet channels or only triplet channels, nlam= 2
!               if both singlet and triplet channels, nlam= 4
!  variables in common block /cosgpi/
!    nlm(1):    the total number of terms in the expansion of the v-sig
!               potential (nlamsg)
!    nlm(2):    the total number of terms in the expansion of the v-pi
!               potential (nlampi)
!    nlm(3):    the total number of terms in the expansion of the v1
!               potential (nlam1)
!    nlm(4)     the total number of terms in the expansion of the v2
!               potential (nlam2)
!  variable in common block /coconv/
!     econv:    conversion factor from cm-1 to hartrees
!     xmconv:   converson factor from amu to atomic units
!  subroutines called:
!   vlm13p:    returns angular coupling coefficient for particular
!              choice of channel index
! ------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlammx, lamnum
use mod_cosysi, only: isicod, ispar
use mod_cosysr, only: isrcod, rspar
use constants, only: econv, xmconv
use mod_parbas, only: lammin, lammax, mproj
use mod_par, only: boundc
use mod_ered, only: ered, rmu
use mod_skip, only: nskip, iskip
use mod_hitypes, only: bqs_type

implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable, target :: v2
type(bqs_type), intent(out) :: bqs
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst

dimension  jhold(1), ehold(1), sc1(1), sc2(1), sc3(1), &
          sc4(1), ishold(1)
!   econv is conversion factor from cm-1 to hartrees
!   xmconv is converson factor from amu to atomic units
integer, pointer :: nterm, nstate, ipol, npot
real(8), dimension(:), pointer :: en, de, re, be, rl, cl
real(8), pointer :: cmix
nterm=>ispar(1); nstate=>ispar(2); ipol=>ispar(3); npot=>ispar(4)
en=>rspar(1:4); de=>rspar(5:8); re=>rspar(9:12); be=>rspar(13:16)
rl=>rspar(17:20); cl=> rspar(21:24); cmix=>rspar(25)

zero = 0.d0
two = 2.d0
!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (6, 5)
  write (9, 5)
5   format (' *** FLAGHF = .TRUE. FOR 1/3 ATOM; ABORT ***' )
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (jlpar .eq. -1) then
  if (ipol .ne. 1 .and. ipol .ne. 0) then
     write (6, 6) ipol
     write (9, 6) ipol
6      format(' ** IPOL =',i3, ' SET EQUAL TO ZERO IN BA13P')
     ipol=0
  endif
else
  if (ipol .ne. 0) then
     write (6, 7) ipol
     write (9, 7) ipol
7      format(' ** IPOL =',i3, ' SET EQUAL TO ZERO IN BA13P')
     ipol=0
  endif
endif
if (csflag) then
  write (6, 8)
  write (9, 8)
8  format &
   ('  *** CSFLAG SET .FALSE. FOR 1/3 ATOM CALCULATION ***')
  csflag=.false.
end if
nsum = 0
if(nterm.ne.1) then
   write(6,9) nterm
   write(9,9) nterm
9    format(' *** NTERM = ',i3,' .NE. 1 FOR 1/3 ATOM; ABORT')
   call exit
end if
if (nstate .eq. 0) then
  nsum=2
else if (nstate .eq. 1) then
  nsum=2
else if (nstate .eq. 2) then
  nsum=4
else
  write (6, 11) nstate
  write (9, 11) nstate
11   format(' *** NSTATE = ',i3,' .NE. 0, 1, OR 2; ABORT')
  call exit
endif
if (nlam .ne. nsum) then
  if (bastst) write (6, 14) nsum
  write (9, 14) nsum
14   format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
  nlam = nsum
end if
if (clist) then
  if (flagsu) then
    write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
    write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16     format(/' **  1/3 ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4, &
      '             E=', f7.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)
  else
      write (6,20) npot, &
          rmu * xmconv, ered * econv, jtot, jlpar
      write (9,20) npot, &
          rmu * xmconv, ered * econv, jtot, jlpar
20       format(/,' **  CC 1/3 ATOM ; NPOT =',i2,' ** RMU=', f9.4, &
           '       E=',f7.2,'   JTOT=', i5, 2x,' JLPAR=',i2)
  end if
  if (.not. flagsu) write (9,30) rcut
30   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
endif
!  assign quantum numbers and energies
call bqs%init(nmax)
n=0
nlevel = 0
! here if triplet state is included
if (nstate .ge. 1) then
  jmin=0
  jmax=2
  do 120 ji=jmin, jmax
    lmin=iabs(jtot-ji)
    lmax=jtot+ji
    do 110 li=lmin, lmax
    ix = (-1) ** (li - jtot)
    if (ix .eq. jlpar) then
!  here for correct orbital angular momentum
      n = n + 1
      if (n .gt. nmax) go to 130
      bqs%lq(n) = li
      cent(n) = li * (li + 1.)
      bqs%inq(n) = 1
      bqs%jq(n) = ji
      bqs%length = n
      eint(n) = en(ji+1)/econv
    end if
110     continue
    nlevel = nlevel + 1
    ehold(nlevel) = en(ji+1)/econv
    jhold(nlevel) = ji
    ishold(nlevel) = 1
120   continue
endif
!  here if singlet state is included
if (nstate .ne. 1) then
  lmin=iabs(jtot-1)
  lmax=jtot+1
  do  125 li=lmin,lmax
    ix = (-1) ** (li - jtot)
    if (ix .eq. jlpar) then
!  here for correct orbital angular momentum
      n = n + 1
      if (n .gt. nmax) go to 130
      bqs%lq(n) = li
      cent(n) = li * (li + 1.)
      if (ipol .eq. 0) then
        bqs%inq(n) = 0
      else if (ipol.eq.1) then
        if (li .lt. jtot) bqs%inq(n)=-2
        if (li .gt. jtot) bqs%inq(n)=2
      endif
      bqs%jq(n) = 1
      bqs%length = n
      eint(n) = en(4)/econv
    end if
125   continue
  if (ipol .eq. 0) then
    nlevel=nlevel+1
    ehold(nlevel) = en(4)/econv
    jhold(nlevel) = 1
    ishold(nlevel) = 0
  else if (ipol .eq. 1) then
    do 126 ii = -2,2,4
      nlevel=nlevel+1
      ehold(nlevel) = en(4)/econv
      jhold(nlevel) = 1
      ishold(nlevel) = ii
126     continue
  endif
endif
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
if (.not.flagsu .and. rcut .gt. 0.d0 .and..not.boundc) then
  emin = 1.e+7
  do 145  i = 1, n
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
145   continue
!  now eliminate all channels with eint .ge. emin if any of the channels
!  are open asymptotically but closed at r = rcut
  if (emin.lt.ered) then
    nn = 0
    do 150 i = 1, n
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        eint(nn) = eint(i)
        bqs%inq(nn) = bqs%inq(i)
        bqs%jq(nn) = bqs%jq(i)
        cent(nn) = cent(i)
        bqs%lq(nn) = bqs%lq(i)
      end if
150     continue
!  reset number of channels
    n = nn
    bqs%length = n
  end if
end if
!  return if no channels
if (n .eq. 0) return
nlevop = nlevel
ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1

!  now list channels if requested
if (clist) then
  write (6, 255)
  write (9, 255)
255   format(/'   N   S   J    L      EINT(CM-1)')
  do 265  i = 1, n
    write (6, 260) i, bqs%inq(i), bqs%jq(i), bqs%lq(i), eint(i) * econv
    write (9, 260) i, bqs%inq(i), bqs%jq(i), bqs%lq(i), eint(i) * econv
260     format (3i4, i5, f13.3)
265   continue
end if
!  now calculate coupling matrix elements
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
do 320 il = 0, 6, 2
  lb = il
  ilam=ilam+1
  if ( ilam .le. nlam ) then  ! depending on nstate, nlam could be 2 or 4, while the il loop always iterate 4 times
    inum = 0
    ancouma => v2%get_angular_coupling_matrix(ilam)
    do icol= 1, n
      do irow = icol, n
        call vlm13p (bqs%jq(irow), bqs%lq(irow), bqs%inq(irow), bqs%jq(icol), &
                       bqs%lq(icol), bqs%inq(icol), jtot, lb, cmix, vee)
        if (vee .ne. 0) then
          i = i + 1
          inum = inum + 1
          call ancouma%set_element(irow=irow, icol=icol, vee=vee)
          if (bastst) then
            write (6, 290) ilam, lb, icol, irow, i, inum, vee
            write (9, 290) ilam, lb, icol, irow, i, inum, vee
  290             format (i4, 2i7, 2i6, i6, g17.8)
          end if
        end if
      end do
    end do
  end if
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
subroutine vlm13p (j1, l1, i1, j2, l2, i2, jtot, lb, cmix, vee)
! --------------------------------------------------------------------
!  subroutine to evaluate the angular coupling matrix element for rotationally
!  inelastic collisions of a singlet and/or triplet P atom with a
!  structureless target

!  authors:  brigitte pouilly and millard alexander
!  current revision date: 23-jul-90
! --------------------------------------------------------------------
!  variables in call list:
!  j1,l1,i1:    initial electronic orbital, orbital, and spin angular momenta
!  j2,l2,i1:    final electronic orbital, orbital, and spin angular momenta
!  jtot:     total angular momentum
!  lb:       value of expansion index:
!            lb=0:  triplet Pi potential
!            lb=2:  triplet Sigma potential
!            lb=4:  singlet Pi potential
!            lb=6:  singlet Sigma potential
!  cmix:     cosine of mixing angle of singlet and triplet states
!  v:        on return:  contains desired coupling matrix element
!  subroutines called:
!  xf3j:     evaluates 3j symbol
!  xf6j:     evaluates 6j symbol
! --------------------------------------------------------------------
use mod_hiutil, only: xf3j, xf6j
implicit double precision (a-h,o-z)
dimension vterm(4),cs(4),sn(4)
data  zero, one, two, three /0.d0, 1.d0, 2.d0, 3.d0/
if (i1 .eq. 1) then
  ind1=j1+1
else
  ind1=4
endif
if (i2 .eq. 1) then
  ind2=j2+1
else
  ind2=4
endif
do  10 i=0,6,2
  is=1
  inu=i
  if (i .gt. 2) then
    is=0
    inu=i-4
  endif
  xnu=inu
  xjtot=jtot
  xl1=l1
  xj1=j1
  xs=is
  xl2=l2
  xj2=j2

  phase=(-1)**(jtot - is)
  xnorm=three*sqrt((two*l1+one)*(two*l2+one)* &
                   (two*j1+one)*(two*j2+one))
  term1=xf3j(one,one,xnu,zero,zero,zero)
  term2=xf3j(xl1,xl2,xnu,zero,zero,zero)
  term3=xf6j(xj1,xj2,xnu,xl2,xl1,xjtot)
  term4=xf6j(xj1,xj2,xnu,one,one,xs)
  vterm(i/2+1)=phase*xnorm*term1*term2*term3*term4
10 continue
smix=sqrt(one-cmix*cmix)
do 20 i=1,4
  if (i .eq. 1 .or. i .eq. 3) then
    cs(i)=zero
    sn(i)=one
  else if (i .eq. 2) then
    cs(i)=-smix
    sn(i)=cmix
  else if (i .eq. 4) then
    cs(i)=cmix
    sn(i)=smix
  endif

20 continue
! cs is the singlet amplitude in each channel
! sn is the triplet amplitude in each channel
! note that only the j=1 channels are mixed
if (lb .le. 2) then
! here for 3Pi or 3Sigma coupling
  vee=vterm(lb/2+1)*sn(ind1)*sn(ind2)
else
  vee=vterm(lb/2+1)*cs(ind1)*cs(ind2)
endif

return
end
! -----------------------------------------
subroutine tcasea(j,jlpar)
! -----------------------------------------
!   matrix to transform from case a to case e, solely for 13P scattering
!   this matrix, tatoe, for total-j=j is returned
!   author:  thierry duhoo and millard alexander
!   latest revision date:  30-dec-1995
! -----------------------------------------
use mod_cosysr, only: isrcod, rspar
use mod_himatrix, only: mxma
use mod_hivector, only: dset
implicit double precision (a-h,o-z)
dimension tatoe(6,6), cmat(6,6)
data zero, one,two ,three/0.d0, 1.d0, 2.d0, 3.d0/
real(8), dimension(:), pointer :: en, de, re, be, rl, cl
real(8), pointer :: cmix
en=>rspar(1:4); de=>rspar(5:8); re=>rspar(9:12); be=>rspar(13:16)
rl=>rspar(17:20); cl=> rspar(21:24); cmix=>rspar(25)

if (j .lt. 2) then
! error if jtot is less than 2
  write (6, 5) jtot
5   format (' *** JTOT = ',i2,' .LT. 2 IN TCASEA; ABORT **')
  stop
endif
if (cmix .le. 0.8d0) then
  write (6, 8) cmix
8   format (' *** WARNING; CMIX =',f6.3,'; POSSIBLY TOO SMALL **')
endif
smix=sqrt(one-cmix*cmix)
! initialization of the matrix tatoe
call dset(36,zero,tatoe,1)
!
! initialize matrix which transforms from case e to spin-orbit diagonalized
! basis
call dset(36,zero,cmat,1)
xj=j
if(jlpar.gt.0) then
! here for f-labelled states
!       singlet state (assumed state 6) is not coupled
!
  tatoe(1,1)=-sqrt(1.d0/3.d0)
  tatoe(1,3)=sqrt(2.d0/3.d0)
  tatoe(2,2)=-sqrt(1.d0/2.d0)
  tatoe(2,4)=sqrt(1.d0/2.d0)
  denom=(two*xj+one)*(two*xj-one)
  tatoe(3,1)=sqrt(xj*(xj-one)/denom)
  tatoe(3,2)=sqrt((xj+one)*(xj-one)/denom)
  tatoe(3,3)=sqrt(xj*(xj-one)/(two*denom))
  tatoe(3,4)=sqrt((xj+one)*(xj-one)/denom)
  tatoe(3,5)=sqrt((xj+one)*(xj+two)/(two*denom))
  denom=(two*xj+three)*(two*xj-one)
  tatoe(4,1)=-sqrt((two*xj*(xj+one))/(three*denom))
  tatoe(4,2)=-sqrt(three/(two*denom))
  tatoe(4,3)=-sqrt((xj*(xj+one))/(three*denom))
  tatoe(4,4)=-sqrt(three/(two*denom))
  tatoe(4,5)=sqrt((three*(xj+two)*(xj-one))/(denom))
  denom=(two*xj+three)*(two*xj+one)
  tatoe(5,1)=sqrt(((xj+two)*(xj+one))/(denom))
  tatoe(5,2)=-sqrt((xj*(xj+two))/(denom))
  tatoe(5,3)=sqrt(((xj+two)*(xj+one))/(two*denom))
  tatoe(5,4)=-sqrt((xj*(xj+two))/(denom))
  tatoe(5,5)=sqrt((xj*(xj-one))/(two*denom))
  tatoe(6,6)=one
! now matrix which transforms from case e to spin-orbit diagonalized
! basis
  cmat(1,1)=one
  cmat(2,2)=cmix
  cmat(2,6)=-smix
  cmat(6,2)=smix
  cmat(3,3)=one
  cmat(4,4)=one
  cmat(5,5)=one
  cmat(6,6)=cmix
! here for e-labelled states
else
  denom=two*(two*xj+one)
  tatoe(1,1)=-sqrt((xj+one)/denom)
  tatoe(1,2)=sqrt((two*xj)/denom)
  tatoe(1,3)=sqrt((xj+one)/denom)
  tatoe(2,1)=-sqrt((xj)/denom)
  tatoe(2,2)=-sqrt((two*(xj+one))/denom)
  tatoe(2,3)=sqrt((xj)/denom)
  tatoe(3,1)=sqrt((xj-one)/denom)
  tatoe(3,3)=sqrt((xj-one)/denom)
  tatoe(3,4)=sqrt(two*(xj+two)/denom)
  tatoe(4,1)=-sqrt((xj+two)/denom)
  tatoe(4,3)=-sqrt((xj+two)/denom)
  tatoe(4,4)=sqrt((two*(xj-one))/denom)
  denom=two*xj+one
  tatoe(5,5)=sqrt(xj/denom)
  tatoe(5,6)=sqrt((xj+one)/denom)
  tatoe(6,5)=-sqrt((xj+one)/denom)
  tatoe(6,6)=sqrt(xj/denom)
! now matrix which transforms from case e to spin-orbit diagonalized
! basis
  cmat(1,1)=cmix
  cmat(1,5)=-smix
  cmat(5,1)=smix
  cmat(2,2)=cmix
  cmat(2,6)=-smix
  cmat(6,2)=smix
  cmat(3,3)=one
  cmat(4,4)=one
  cmat(5,5)=cmix
  cmat(6,6)=cmix
endif
call mxma(cmat,1,6,tatoe,1,6,ttrans,1,6,6,6,6)
!      call mxoutd(94, cmat, 6,6,0,.false.)
!      call mxoutd(94, tatoe, 6,6,0,.false.)
!      call mxoutd(94, ttrans, 6,6,0,.false.)
!      call exit
return
end
! -----------------------------------------------------------------------
subroutine sy13p (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for collisions of
! atom in singlet and/or triplet P electronic state with closed shell atom
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  current revision date: 17-jun-1992 by mha
!  -----------------------------------------------------------------------
!  variable in common cosysi/
!    nscode:  total number of system dependent parameters
!             nscode = isicod + isrcod +3
!    isicod:  number of integer system dependent parameters
!    nterm:    number of different electronic coupling terms
!              this should be 1 here
!    nstate:   number of electronic states included
!              nstate=0:   just singlet state
!              nstate=1:   just triplet state
!              nstate=2:   both singlet and triplet state
!  variable in common cosys/
!    scod:    character*8 array of dimension lcode, which contains names
!             of all system dependent parameters
!  line 16:
!    filnam:  name of file containing potential parameters
!
!  subroutines called: loapot(iunit,filnam)
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, rspar
use funit, only: FUNIT_INP
use mod_parbas, only: lammin, lammax, mproj
use mod_skip, only: nskip, iskip
use mod_hiutil, only: gennam, get_token
implicit none
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
integer :: i, j, l, lc
real(8) :: rms
logical existf

character*1 dot
character*(*) fname
character*60 filnam, line, potfil, filnm1
#include "common/comdot.F90"
save potfil
integer, pointer, save :: nterm, nstate, ipol, npot
real(8), dimension(:), pointer, save :: en, de, re, be, rl, cl
real(8), pointer, save :: cmix, alphg, rgaus, agaus, demor, remor, bemor, dissmor
nterm=>ispar(1); nstate=>ispar(2); ipol=>ispar(3); npot=>ispar(4)
en=>rspar(1:4); de=>rspar(5:8); re=>rspar(9:12); be=>rspar(13:16)
rl=>rspar(17:20); cl=>rspar(21:24); cmix=>rspar(25); alphg=>rspar(26)
rgaus=>rspar(27); agaus=>rspar(28); demor=>rspar(29); remor=>rspar(30)
bemor=>rspar(31); dissmor=>rspar(32) 

!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
isicod = 4
scod(1) = 'NTERM'
scod(2) = 'NSTATE'
scod(3) = 'IPOL'
scod(4) = 'NPOT'
!  then all the system dependent real variables
!  in the same order as in the common block /cosysr/
isrcod = 32
scod(5) = 'E-3P0'
scod(6) = 'E-3P1'
scod(7) = 'E-3P2'
scod(8) = 'E-1P1'
scod(9) = 'DE(3PI)'
scod(10) = 'DE(3SG)'
scod(11) = 'DE(1PI)'
scod(12) = 'DE(1SG)'
scod(13) = 'RE(3PI)'
scod(14) = 'RE(3SG)'
scod(15) = 'RE(1PI)'
scod(16) = 'RE(1SG)'
scod(17) = 'BE(3PI)'
scod(18) = 'BE(3SG)'
scod(19) = 'BE(1PI)'
scod(20) = 'BE(1SG)'
scod(21) = 'RL(3PI)'
scod(22) = 'RL(3SG)'
scod(23) = 'RL(1PI)'
scod(24) = 'RL(1SG)'
scod(25) = 'C(3PI)'
scod(26) = 'C(3SG)'
scod(27) = 'C(1PI)'
scod(28) = 'C(1SG)'
scod(29) = 'CMIX'
scod(30) = 'ALPHG'
scod(31) = 'RGAUS'
scod(32) = 'AGAUS'
scod(33) = 'DEMOR'
scod(34) = 'REMOR'
scod(35) = 'BEMOR'
scod(36) = 'DISSMOR'
nscode = isicod + isrcod
!  set default values for 1/3 P atom scattering
nterm = 1
nstate = 2
ipol=1
lammin(1)= 1
lammax(1) = 4
mproj(1)=0
niout=2
indout(1)=0
indout(2)=1
potfil = ' '
if(iread.eq.0) return
read (8, *, err=888) nstate, ipol, npot
read (8, *, err=888) (en(i), i=1,3)
read (8, *, err=888) en(4), cmix
read (8, *, err=888)  de(1), re(1), be(1), rl(1), cl(1)
read (8, *, err=888)  de(2), re(2), be(2), rl(2), cl(2)
read (8, *, err=888)  de(3), re(3), be(3), rl(3), cl(3)
read (8, *, err=888)  de(4), re(4), be(4), rl(4), cl(4)
read (8, *, err=888)  demor, remor, bemor, dissmor
read (8, *, err=888)  rgaus, agaus, alphg
rms=0.d0
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
entry ptr13p (fname,readpt)
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
entry sav13p (readpt)
!  save input parameters for singlet-sigma + atom scattering
write (FUNIT_INP, 300) nstate, ipol, npot
300 format(3i4, 18x,'   nstate, ipol, npot')
write (FUNIT_INP, 310) (en(i), i=1,3)
310 format(3(1pg12.4),4x,'   E-3P0, 3P1, 3P2')
write (FUNIT_INP, 315) en(4), cmix
315 format (2(1pg12.4),16x,'   E-1P1, CMIX')
write (FUNIT_INP, 320) de(1), re(1), be(1), rl(1), cl(1)
320 format (5(1pg12.4),' 3 Pi Parameters')
write (FUNIT_INP, 325) de(2), re(2), be(2), rl(2), cl(2)
325 format (5(1pg12.4),' 3 Sig Parameters')
write (FUNIT_INP, 330) de(3), re(3), be(3), rl(3), cl(3)
330 format (5(1pg12.4),' 1 Pi Parameters')
write (FUNIT_INP, 335) de(4), re(4), be(4), rl(4), cl(4)
335 format (5(1pg12.4),' 1 Sig Parameters')
write (FUNIT_INP, 336) demor, remor, bemor, dissmor
336 format(4(1pg12.4),12x,' 1 Sigma Morse')
write (FUNIT_INP, 337) rgaus, agaus, alphg
337 format (3(1pg12.4),24x,' Gaussian Coupling')
write (FUNIT_INP, 285) potfil
return
end
end module mod_hiba07_13p
