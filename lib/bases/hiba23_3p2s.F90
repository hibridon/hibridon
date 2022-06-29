#include "assert.h"
module mod_hiba23_3p2s
contains
! sy3p2s (sav3p2s/ptr3p2s) defines, saves variables and reads            *
!                  potentials for 3P atom + 2S atom                      *
! --------------------------------------------------------------------
subroutine ba3p2s (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, &
     bastst, ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine coupling potential (electrostatic + spin-orbit)
!  for collision involving the 3P state of atom in a p^4 electronic
!  configuration with an atom in a 2S state
!  (based on the ba1d3p basis routine)
!
!  THIS IS A REVISION OF DAGDIGIAN'S ORIGINAL BASIS ROUTINE, USING THE
!  TRANSFORMAITON WORKED OUT BY SINGER, FREED, AND BAND, FREED, JCP 79,
!  6060 (1983).  UNLIKE THE ROUTINE hiba3p2s.f~origPJD, THE SPIN-ORBIT
!  COUPLING IS PARAMETERIZED WITH JUST TWO PARAMETERS and APERP, WORKED
!  OUT IN MHA'S NOTES.  APAR AND APERP ARE ASSUMED TO BE CONSTANTS.
!
!  j12 array is in module mod_coj12 to pass to other subrs
!
!  author:  paul dagdigian
!  current revision date:  14-jan-2019 by p.dagdigian
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains electronic angular momentum quantum number
!              for each channel (j = 0, 1, 2 for the 3P spin-orbit levels
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
!    isc1:     scratch vector of length at least nmax
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
!    isrcod:   number of real system dependent variables (none here)
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   total number of integer system dependent variables
!    nterm:    number of different types of electronic coupling terms
!              this should be 1 here
!  variables in common block /coered/
!    ered:      collision energy in atomic units (hartrees)
!    rmu:       collision reduced mass in atomic units
!               (mass of electron = 1)
!  variable in module mod_conlam
!    nlam:      the number of case(a) interaction potentials actually used
!               this is :  nlam = nlam0 + nlam1
!               if only 1D channels, nlam = 3 (1Sigma+, 1Pi, 1Delta PE curves)
!               if only 3P channels, nlam = 6 (3Sigma- and 3Pi PE curves
!               and two spin-orbit matrix elements Axy and Axz, and the
!               R=inf values of Axy and Azy)
!               if both 1D and 3P channels, nlam= 18 (4 PE curves mentioned above,
!               2 spin-orbit 3P matrix elements, and 5 spin-orbit matrix
!               elements describing 1D-3P coupling, and the R-inf values of
!               7 spin-orbit matrix elements)

!  variable in common block /coconv/
!     econv:    conversion factor from cm-1 to hartrees
!     xmconv:   converson factor from amu to atomic units
!  subroutines called:
!   vlm13p:    returns angular coupling coefficient for particular
!              choice of channel index
!
!  subroutines called:
!  xf3j:     evaluates 3j symbol
!  xf9j:     evaluates 9j symbol
! ------------------------------------------------------------
use mod_coeig2, only: t12, t32
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coj12, only: j12
use mod_covvl, only: vvl
use mod_conlam, only: nlam
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use constants, only: econv, xmconv
implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst
#include "common/parbas.F90"
#include "common/parbasl.F90"
!  common blocks cojtot, coja, and coel used to transmit to ground subroutine
common /cojtot/ jjtot, jjlpar
common /coja/  jja(9)
common /coel/  ll(9)
common /coered/ ered, rmu
common /coskip/ nskip, iskip
!  arrays in argument list
dimension j(1), l(1), is(1),jhold(1),ehold(1),ishold(1)
!  matrices for transformation between atomic and molecular BF functions
dimension c12(5,5), c32(3,3)
!  quantum numbers for BF atomic and basis functions
dimension &
  fjo32(3), fj32(3), flam32(3), sigm32(3), stot32(3), &
  fjo12(5), fj12(5), flam12(5), sigm12(5), stot12(5)
data flo, so, flh, sh, fjh /1.d0, 1.d0, 0.d0, 0.5d0, 0.5d0/
!  atomic BF basis functions (omega=3/2)
data fjo32 /2.d0, 2.d0, 1.d0/, fj32 /2.5d0, 1.5d0, 1.5d0/
!  molecular BF basis functions (omega=3/2)
data flam32 /1.d0, 0.d0, 1.d0/, sigm32 /0.5d0, 1.5d0, 0.5d0/, &
  stot32 /1.5d0, 1.5d0, 0.5d0/
!  atomic BF basis functions (omega=1/2)
data fjo12 /2.d0, 2.d0, 1.d0, 1.d0, 0.d0/, &
  fj12 /2.5d0, 1.5d0, 1.5d0, 0.5d0, 0.5d0/
!  molecular BF basis functions (omega=1/2)
data flam12 /-1.d0, 1.d0, 0.d0, 1.d0, 0.d0/, &
  sigm12 /1.5d0, -0.5d0, 0.5d0, -0.5d0, 0.5d0/, &
  stot12 /1.5d0, 1.5d0, 1.5d0, 0.5d0, 0.5d0/
!  scratch arrays for computing asymmetric top energies and wave fns.
dimension en0(4), en12(5), en32(3), vec(5,5), work(288)
integer, pointer :: nterm, nstate
real(8), pointer :: en1d
nterm=>ispar(1); nstate=>ispar(2)
en1d=>rspar(1)

zero = 0.d0
two = 2.d0
!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (6, 5)
  write (9, 5)
5   format (' *** FLAGHF = .TRUE. FOR 3P/2S ATOM; ABORT ***' )
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (csflag) then
  write (6, 8)
  write (9, 8)
8  format &
   ('  *** CSFLAG SET .FALSE. FOR 3P/2S ATOM CALCULATION ***')
  csflag=.false.
end if
if (bastst) then
  write (6, 14) nlam
  write (9, 14) nlam
14   format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2)
end if
if (bastst) then
  if (flagsu) then
    write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
    write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16     format(/' **  3P/2S ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4, &
      ' E=', f7.2, /, '         JTOT=', i5, 2x,' JLPAR=',i2)
  else
    write (6,20) npot, &
        rmu * xmconv, ered * econv, jtot, jlpar
    write (9,20) npot, &
        rmu * xmconv, ered * econv, jtot, jlpar
20     format(/,' **  CC 3P/2S ATOM ; NPOT =',i2,' ** RMU=', f9.4, &
        ' E=', f10.2, /, '         JTOT=', i5, 2x,' JLPAR=',i2)
  end if
  if (.not. flagsu) write (9,30) rcut
30   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
endif
!
!  form list of spin-orbit levels and determine the energetically open levels
!  order of the j levels for p^4 orbital occupancy
!
!  get atomic spin-orbit parameter (apar at large R)
r = 100.d0
call pot(vv0 ,r)
a  = vvl(5) * econv
nlevel = 0
do jj = 2, 0, -1
  nlevel = nlevel + 1
  jhold(nlevel) = jj
  ishold(nlevel) = 3
  eso = - 0.5d0 * a * (jj * (jj + 1.d0) - 4.d0)
!  set zero of energy as E(j=2)
  ehold(nlevel) = (eso + 0.5d0 * a * 2.d0) / econv
end do
!
!  form list of energetically open levels and sort list
!  to put closed channels at end
nlevop = 0
if (nlevel .gt. 1) then
  do 180 i = 1, nlevel - 1
    if (ehold(i) .le. ered) then
      nlevop = nlevop + 1
    else
      do 175 ii = i + 1, nlevel
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
          goto 180
        endif
175       continue
    endif
180   continue
  if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
else
!  here for only one level
  if (ehold(1) .le. ered) then
    nlevop = 1
  else
    nlevop = 0
  endif
endif
if (nlevop .le. 0) then
  write (9,185)
  write (6,185)
185   format('*** NO OPEN LEVELS IN BA3P2S; ABORT')	
  if (bastst) return
  call exit
endif
!
!  set up CC channel basis
n = 0
do 120 jj = 1, nlevel
  ji = jhold(jj)
!   determine j12 values - vector coupling of j(O) + j(H) [=1/2]
  fj12mn = abs(ji - fjh)
  fj12mx = ji + fjh
  j12min = fj12mn
  j12max = fj12mx
  do 119 j12i = j12min, j12max
    lmin = abs(jtot-j12i)
    lmax = jtot + j12i + 1
    do 110 li=lmin,lmax
!  parity in symmetry-adapted basis
      ix = (-1) ** li
      if (ix .eq. jlpar) then
!  here for correct orbital angular momentum
        n = n + 1
        if (n .gt. nmax) go to 130
        l(n) = li
        cent(n) = li * (li + 1.)
        is(n) = 3
        j(n) = ji
        j12(n) = j12i
!  next two variables for transmission to ground subroutine
        jja(n) = j(n)
        ll(n) = l(n)
        eint(n) = ehold(jj)
      end if
110     continue
119   continue
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
!
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
ntop = max(n, nlevop)
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
!        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
!     :       ntop = ntop + 1
!++++  INCREASE SIZE OF NTOP
!      ntop = ntop + 4
!
!  now list channels if requested
if (bastst) then
  write (6,1255)
  write (9,1255)
1255   format(/' ASYMPTOTIC ATOMIC ENERGY LEVELS' &
    /'   N   S   J      EINT(CM-1)')
  do 1265 i = 1, nlevel
    write (6,2260) i, ishold(i), jhold(i), ehold(i) * econv
    write (9,2260) i, ishold(i), jhold(i), ehold(i) * econv
2260     format (3i4, f13.3)
1265   continue
  write (6, 255)
  write (9, 255)
255   format(/' CHANNEL BASIS FUNCTIONS' &
      /'   N   S   J   J12   L      EINT(CM-1)')
  do 265  i = 1, n
    write (6, 260) i, is(i), j(i), j12(i), l(i), eint(i) * econv
    write (9, 260) i, is(i), j(i), j12(i), l(i), eint(i) * econv
260     format (3i4, 2i5, f13.3)
265   continue
end if
!
!  now calculate coupling matrix elements
if (bastst) then
  write (6, 280)
  write (9, 280)
280   format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
end if
!
!  compute transformation between jo, jh, j and Lambda, S, Sigma
!  BF basis sets
!
!  omega = 1/2
omega = 0.5d0
do ia = 1, 5
  do im = 1, 5
    iph = flo - stot12(im) + 0.5d0
    t12(im,ia) = sqrt((2.d0*stot12(im)+1.d0)*(2.d0*fjh+1.d0) &
      *(2.d0*fjo12(ia)+1.d0)*(2.d0*flo+1.d0)) &
      *(-1.d0)**iph * sqrt(2.d0*fj12(ia)+1.d0) &
      *xf3j(flo,stot12(im),fj12(ia),flam12(im), &
          sigm12(im),-omega) &
      *xf9j(flo,so,fjo12(ia), flh,sh,fjh, &
          flo,stot12(im),fj12(ia))
  enddo
enddo
!
!  omega = 3/2
omega = 1.5d0
do ia = 1, 3
  do im = 1, 3
    iph = flo - stot32(im) + 1.5d0
    t32(im,ia) = sqrt((2.d0*stot32(im)+1.d0)*(2.d0*fjh+1.d0) &
      *(2.d0*fjo32(ia)+1.d0)*(2.d0*flo+1.d0)) &
      *(-1.d0)**iph * sqrt(2.d0*fj32(ia)+1.d0) &
      *xf3j(flo,stot32(im),fj32(ia),flam32(im), &
          sigm32(im),-omega) &
      *xf9j(flo,so,fjo32(ia), flh,sh,fjh, &
          flo,stot32(im),fj32(ia))
  enddo
enddo
!
!  i counts v2 elements
!  inum counts v2 elements for given lambda
!  ilam counts number of v2 matrices
i = 0
ilam=0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 320 lb = 1, nlam
  ilam = ilam + 1
  ancouma => v2%get_angular_coupling_matrix(ilam)
  inum = 0
  do icol= 1, n
    do irow = icol, n
      call vlm3p2s (j(irow), j12(irow), l(irow), &
        j(icol), j12(icol), l(icol), jtot, lb, vee)
      if (vee .ne. 0) then
        i = i + 1
        inum = inum + 1
        call ancouma%set_element(irow=irow, icol=icol, vee=vee)
        if (bastst) then
          write (6, 290) ilam, lb, icol, irow, i, vee
          write (6, 290) ilam, lb, icol, irow, i, vee
290             format (i4, 2i7, 2i6, g17.8)
        endif
      end if
    end do
  end do
if (bastst) then
  write (6, 315) ilam, ancouma%get_num_nonzero_elements()
  write (9, 315) ilam, ancouma%get_num_nonzero_elements()
315   format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
end if
320 continue
if (bastst) then
  write (6, 360) i
  write (9, 360) i
360   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS', &
           i6)
end if
return
end
! --------------------------------------------------------------------
subroutine vlm3p2s (j1, j12_1, l1, j2, j12_2, l2, jtot, lb, vee)
! --------------------------------------------------------------------
!  subroutine to evaluate the angular coupling matrix element for
!  inelastic collisions of a 3P atom with a 2S atom

!  author:  paul dagdigian
!  current revision date: 18-dec-2014
! --------------------------------------------------------------------
!  variables in call list:
!  j1,j12_1,l1:    initial electronic and orbital angular momenta
!  j2,j12_2,l2:    final electronic and orbital angular momenta
!  jtot:     total angular momentum
!  lb:       value of expansion index:
!            lb=1:  X2Pi potential
!            lb=2:  4Pi potential
!            lb=3:  2Sigma- potential
!            lb=4:  4Sigma- potential
!            lb=5:  apar spin-orbit matrix element (assumed independent of R)
!            lb=6:  aperp spin-orbit matrix element (assumed independent of R)
!
!  vee:      on return:  contains desired coupling matrix element
!  subroutines called:
!  xf3j:     evaluates 3j symbol
! --------------------------------------------------------------------
use mod_coeig2, only:  t12, t32
implicit double precision (a-h,o-z)
dimension v12(5,5),v32(3,3),a12(5,5),a32(3,3)
vee = 0.d0
xl1 = l1
xl2 = l2
xj1 = j1
xj2 = j2
xj12_1 = j12_1 + 0.5d0
xj12_2 = j12_2 + 0.5d0
xjtot = jtot + 0.5d0
xl12 = sqrt((2.d0*xl1 + 1.d0) * (2.d0*xl2 + 1.d0))
xfj12 = xf3j(xl1, xj12_1, xjtot, 0.d0, 0.5d0, -0.5d0) &
  * xf3j(xl2, xj12_2, xjtot, 0.d0, 0.5d0, -0.5d0)
xfj32 = xf3j(xl1, xj12_1, xjtot, 0.d0, 1.5d0, -1.5d0) &
  * xf3j(xl2, xj12_2, xjtot, 0.d0, 1.5d0, -1.5d0)
xfj52 = xf3j(xl1, xj12_1, xjtot, 0.d0, 2.5d0, -2.5d0) &
  * xf3j(xl2, xj12_2, xjtot, 0.d0, 2.5d0, -2.5d0)
!  clear body-frame Lambda,S,Sigma pot/Hso matrices for Omega = 1/2, 3/2
do i = 1,5
  do j = 1,5
    v12(i,j) = 0.d0
  enddo
enddo
do i = 1,3
  do j = 1,3
    v32(i,j) = 0.d0
  enddo
enddo
! get basis #'s for omega = 1/2 and 3/2 BF components of atomic basis functions
i11 = 1
i13 = 1
if (j1.eq.2 .and. j12_1.eq.1) then
  i11 = 2
  i13 = 2
endif
if (j1.eq.1 .and. j12_1.eq.1) then
  i11 = 3
  i13 = 3
endif
if (j1.eq.1 .and. j12_1.eq.0) then
  i11 = 4
  i13 = 0
endif
if (j1.eq.0) then
  i11 = 5
  i13 = 0
endif
!
i21 = 1
i23 = 1
if (j2.eq.2 .and. j12_2.eq.1) then
  i21 = 2
  i23 = 2
endif
if (j2.eq.1 .and. j12_2.eq.1) then
  i21 = 3
  i23 = 3
endif
if (j2.eq.1 .and. j12_2.eq.0) then
  i21 = 4
  i23 = 0
endif
if (j2.eq.0) then
  i21 = 5
  i23 = 0
endif
!
sign = 1.d0
goto (10,20,30,40,50,60),lb
goto 1000
!
!  X2Pi PE curve - mat elements below in Lambda,S,Sigma,Omega basis
10 v12(4,4) = 1.d0
v32(3,3) = 1.d0
!  compute matrix element in j1,j2,j12,Omega basis
call dgemm('t','n',5,5,5,1.d0,t12,5,v12,5,0.d0,a12,5)
call dgemm('n','n',5,5,5,1.d0,a12,5,t12,5,0.d0,v12,5)
vee12 = v12(i11,i21)
vee32 = 0.d0
if (i13.ne.0 .and. i23.ne.0) then
  call dgemm('t','n',3,3,3,1.d0,t32,3,v32,3,0.d0,a32,3)
  call dgemm('n','n',3,3,3,1.d0,a32,3,t32,3,0.d0,v32,3)
  vee32 = v32(i13,i23)
endif
vee = 2.d0 * xl12 * (vee12 * xfj12 + vee32 * xfj32)
goto 1000
!
!  4Pi PE curve
20 v12(1,1) = 1.d0
v12(2,2) = 1.d0
v32(1,1) = 1.d0
vee52 = 1.d0
!  compute matrix element in j1,j2,j12,Omega basis
call dgemm('t','n',5,5,5,1.d0,t12,5,v12,5,0.d0,a12,5)
call dgemm('n','n',5,5,5,1.d0,a12,5,t12,5,0.d0,v12,5)
vee12 = v12(i11,i21)
vee32 = 0.d0
if (i13.ne.0 .and. i23.ne.0) then
  call dgemm('t','n',3,3,3,1.d0,t32,3,v32,3,0.d0,a32,3)
  call dgemm('n','n',3,3,3,1.d0,a32,3,t32,3,0.d0,v32,3)
  vee32 = v32(i13,i23)
endif
vee = 2.d0 * xl12 * (vee12 * xfj12 + vee32 * xfj32 &
  + vee52 * xfj52)
goto 1000
!
!  2Sigma- curve
30 v12(5,5) = 1.d0
!  compute matrix element in j1,j2,j12,Omega basis
call dgemm('t','n',5,5,5,1.d0,t12,5,v12,5,0.d0,a12,5)
call dgemm('n','n',5,5,5,1.d0,a12,5,t12,5,0.d0,v12,5)
vee12 = v12(i11,i21)
vee = 2.d0 * xl12 * vee12 * xfj12
goto 1000
!
!  4Sigma- curve
40 v12(3,3) = 1.d0
v32(2,2) = 1.d0
!  compute matrix element in j1,j2,j12,Omega basis
call dgemm('t','n',5,5,5,1.d0,t12,5,v12,5,0.d0,a12,5)
call dgemm('n','n',5,5,5,1.d0,a12,5,t12,5,0.d0,v12,5)
vee12 = v12(i11,i21)
vee32 = 0.d0
if (i13.ne.0 .and. i23.ne.0) then
  call dgemm('t','n',3,3,3,1.d0,t32,3,v32,3,0.d0,a32,3)
  call dgemm('n','n',3,3,3,1.d0,a32,3,t32,3,0.d0,v32,3)
  vee32 = v32(i13,i23)
endif
vee = 2.d0 * xl12 * (vee12 * xfj12 + vee32 * xfj32)
goto 1000
!
!  apar spin-orbit matrix element
50 vee = 0.d0
goto 1000
!
!  aperp spin-orbit matrix element
60 vee = 0.d0
1000 return
end
! ------------------------------eof-----------------------------------
! -----------------------------------------------------------------------
subroutine sy3p2s (irpot, readp, iread)
!  subroutine to read in system dependent parameters for collisions of
!  atom in 3P state with an atom in 2S state
!  current revision date: 1-oct-2018 by p.dagdigian
!  -----------------------------------------------------------------------
!  NOTE:  no real variables for this basis type
!  variable in common cosysi/
!    nscode:   total number of system dependent parameters
!              nscode = isicod + isrcod + 3
!    isicod:   number of integer system dependent parameters
!    nterm:    number of different types of electronic coupling terms
!              this should be 1 here
!    nvib:     initial vibrational state of A state, in photodissociation calculation
!  variable in common  /cosys/
!    scod:    character*8 array of dimension lcode, which contains names
!             of all system dependent parameters
!  line 16:
!    filnam:  name of file containing potential parameters
!
!  subroutines called: loapot(iunit,filnam)
!  -----------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_conlam, only: nlam
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod
implicit double precision (a-h,o-z)
integer irpot
logical readp
logical airyfl, airypr, logwr, swrit, t2writ, writs, wrpart, &
        partw, xsecwr, wrxsec, noprin, chlist, ipos, flaghf, &
        csflag, flagsu, rsflag, t2test, existf, logdfl, batch, &
        readpt, ihomo, bastst, twomol, lpar

character*1 dot
character*(*) fname
character*60 filnam, line, potfil, filnm1
#include "common/parbas.F90"
common /coskip/ nskip,iskip
common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag, &
                flaghf, flagsu, ihomo, ipos, logdfl, logwr, &
                noprin, partw, readpt, rsflag, swrit, &
                t2test, t2writ, twomol, writs, wrpart, wrxsec, &
                xsecwr,lpar(3)
#include "common/comdot.F90"
save potfil
integer, pointer :: nterm, nvib
nterm=>ispar(1); nvib=>ispar(2)
!  number and names of system dependent parameters
!  first all the system dependent integer variables
!  in the same order as in the common block /cosysi/
!  variable names must be in uppercase of maximum length 6 characters
isicod = 2
scod(1) = 'NTERM'
scod(2) = 'NVIB'
!  then all the system dependent real variables
!  in the same order as in the common block /cosysr/
isrcod = 0
nscode = isicod + isrcod
!  set default values for 3P + 2S atom scattering
nterm = 1
lammin(1)= 1
lammax(1) = 6
mproj(1)=0
nlam = 6
niout=1
indout(1)=3
potfil = ' '
if(iread.eq.0) return
read (8, *, err=888) nvib
line=' '
if(.not.readp.or.iread.eq.0) then
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
entry ptr3p2s (fname,readp)
line = fname
readp = .true.
286 if (readp) then
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
! --------------------------------------------------------------
entry sav3p2s (readp)
!  save input parameters for 3P atom + 2S atom scattering
write (8, 1285) nvib
1285 format(i4, 29x,'nvib')
write (8, 285) potfil
return
end
end module mod_hiba23_3p2s
