#include "assert.h"
module mod_hiba08_2mol
contains
! sy2mol (sav2mol/ptr2mol) defines, save variables and reads             *
!                  potential for 2 singlet sigma molecule scattering     *
subroutine ba2mol (j, l, is, jhold, ehold, ishold, nlevel, nlevop, &
                  j1, j2, j12, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  nch, nmax, ntop, v2)
! --------------------------------------------------------------------
!  subroutine to determine angular coupling potential for collision
!  of two rigid heternuclear diatomic molecules
!  authors:  millard alexander and peter vohralik
!  current revision date:  7-apr-03 by mha
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains combined rotational quantum numbers for each
!              channel.  in terms of the rotational quantum numbers of each
!              molecule we have:  j = 10 * j1 + j2
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       additional quantum index of each channel:  on return this is
!              set equal to j12 (the vector resultant of j1 and j2) times nsym
!    jhold:    on return contains rotational quantum numbers for each
!              rotational level, irrespective of projection degeneracy and the
!              degeneracy associated with different values of j12
!              note that jhold = 10 * j1hold + j2hold
!    ehold:    on return contains energy in hartrees of each rotational
!              level
!    ishold:   on return this is set equal to nsym
!    nlevel:   on return contains number of energetically distinct
!              rotational levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              rotational levels used in channel basis which are open
!              asymptotically
!    j1:       rotational quantum number of molecule 1 in each channel
!    j2:       rotational quantum number of molecule 2 in each channel
!    j12:      the vector resultant of the two molecular rotational angular
!              momenta:  j12 = j1 + j2
!              note:  j1, j2, and j12 are not preserved in calling program!!!
!    sc4:      scratch vector (not used here)
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
!              the homonuclear option is not specifically implemented here
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              parity=(-1)**jlpar (by definition parity=(-1)**(j1+j1+l)
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to subroutines
!              propag and soutpt.  ntop is set in basis only if nu = numin
!              otherwise it is unchanged from the value supplied by the
!              calling program
!  variables in common block /cosysr/
!    isrcod:   number of real system dependent variables
!    brot:     rotational constant in cm-1
!    drot:     centrifugal distortion constant in cm-1
!    hrot:     next centrifugal distortion constant in cm-1
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!    isicod:   numbr of integer system dependent variables
!    nterm:    number of different lambda terms included in potential
!     nterm = 1, include lam1=1  lam2=1  lam=2   dipole-dipole
!     nterm = 2, include dipole-dipole and
!                dipole-quadrupole (lam1=1  lam2=2  lam=3)
!     nterm = 3, include dipole-dipole, dipole-quadrupole, and
!                short-range term (lam1=0  lam2=1  lam=1)
!    nsym:     interchange symmetry of included channels
!  variables in common block /cotwo/
!    numj:     number of j1-j2 values
!    nj1j2:    specific j1-j2 values (up to a maximum of 50)
!              N.B. this dimension is set here
!  variables in common block /coered/
!    ered:      collision energy in atomic units (hartrees)
!    rmu:       collision reduced mass in atomic units
!               (mass of electron = 1)
!  variable in module mod_conlam
!               nlam is set equal to nterm; see above
!  subroutines called:
!   vlmlml:    returns molecule-molecule angular coupling coefficient for
!              particular choice of channel index
! --------------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_conlam, only: nlam
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use constants, only: econv, xmconv
use mod_par, only: iprint
#include "common/parbasl.F90"

implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
logical ihomo, flaghf, csflag, clist, flagsu, bastst
#include "common/parbas.F90"
common /cotwo/ numj,nj1j2(50)
common /coselb/ ibasty
common /coered/ ered, rmu
dimension j(1), l(1), is(1), jhold(1), ehold(1), j12(1), j1(1), &
          j2(1), sc4(1), ishold(1)
data ione, itwo, ithree / 1, 2, 3 /
integer, pointer :: nterm, nsym
real(8), pointer :: brot, drot, hrot
nterm=>ispar(1); nsym=>ispar(2)
brot=>rspar(1); drot=>rspar(2); hrot=>rspar(3)

!  check for consistency in the values of flaghf and csflag
if (flaghf) then
  write (9, 5)
5   format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***' )
  if (bastst) then
    return
  else
    call exit
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
    call exit
  end if
end if
!  check for consistency in values of numj, nterm
if (numj .gt. 50) then
  write (6, 15) numj
  write (9, 15) numj
15   format (' *** NUMBER OF J1-J2 PAIRS=',i3, &
          ' .GT. 50; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (nterm .gt. 3) then
  write (6, 20) nterm
  write (9, 20) nterm
20   format (' *** NUMBER OF ANISOTROPIC TERMS IN POTENTIAL=', i3, &
          ' .GT. 3; ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
if (csflag) then
  write (6, 25)
  write (9, 25)
25   format (' *** CSFLAG=.TRUE.; NOT ALLOWED IN HF-HF BASIS;', &
          ' ABORT ***')
  if (bastst) then
    return
  else
    call exit
  end if
end if
write (9,  90) ered*econv, jtot, nsym, jlpar, nlam
90 format(/2x, ' **  CC  HF-HF  **  E=', f9. 3, 2x, 'JTOT=', i4, &
   '  SYM=', i2, '  PAR=', i2, '  N-LAM=', i2)
write (9, 95) rcut
95   format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=', &
          f8.2)
nlam = nterm
!  assign quantum numbers and energies for rotational levels
nlevop = 0
nch = 0
do 200  i = 1, numj
  jj1 = nj1j2(i) / 10
  jj2  = mod(nj1j2(i), 10)
  jmax = jj1 + jj2
  jmin = iabs(jj1 - jj2)
  ijpar=(-1)**(jj1+jj2)
  do 180  jj12 = jmin, jmax
    lmin = iabs(jtot - jj12)
    lmax = jtot + jj12
    if (lmin .gt. lmax) go to 180
    do 170   ll = lmin, lmax
      ix = ( - 1) ** (jj1 + jj2 + ll)
      if (ix .eq. jlpar) then
        if (jj1 .eq. jj2) then
          ix = ( - 1) ** (jj1 + jj2 + jj12 + ll)
          if (ix .ne. nsym) go to 170
        end if
        nch = nch + 1
        if (nch .gt. nmax) go to 220
        j1(nch) = jj1
        j2(nch) = jj2
        j(nch) = 10 * jj1 + jj2
        j12(nch) = jj12
        l(nch) = ll
        is(nch)=ijpar
!  now calculate the diagonal matrix elements of the hamiltonian
!  first the internal rotational energy
        fj1 = jj1 * (jj1 + ione)
        fj2 = jj2 * (jj2 + ione)
        eint(nch) = (brot * (fj1 + fj2) - &
                     drot * (fj1 * fj1 + fj2 * fj2) + &
                     hrot * (fj1 ** 3 + fj2 ** 3) ) / econv
!  now the centrifugal potential
        cent(nch) =  ll * (ll + ione)
      end if
170     continue
180   continue
!  form list of all energetically open rotational levels included in the
!  calculations and their energies
  jhold(i)  = 10 * jj1 + jj2
  fj1       = jj1 * ( jj1 + ione )
  fj2       = jj2 * ( jj2 + ione )
  ehold(i)  =   (brot/econv) * ( fj1 + fj2 ) &
               - (drot/econv) * ( fj1**2 + fj2**2 ) &
               + (hrot/econv) * ( fj1**3 + fj2**3 )
  ishold(i) = ijpar
!-------------------------------------------------------------------
!***  PFV'S addition
  if(ered-ehold(i).gt.0.) then
    nlevop = nlevop + 1
    if(i.ne.nlevop) then
      write(9,271)
271       format( &
    ' ROTATIONAL LEVELS ARE NOT ORDERED WITH ALL THE '/ &
         ' OPEN CHANNELS BEFORE THE CLOSED ONES . PFV. ')
      call exit
     end if
   end if
!***
200 continue
220 if (nch .gt. nmax) then
  write (9, 230) nch, nmax
  write (6, 230) nch, nmax
230   format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF', &
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
  do 260  i = 1, nch
    if (eint(i) .le. ered) then
!  here if channel is open asymptotically
      if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut) &
          .gt. (ered - eint(i)) ) then
!  here if channel is open asymptotically but closed at r = rcut
        if (eint(i) .lt. emin) emin = eint(i)
!  emin now contains the lowest channel energy for which this
!  condition is met
      end if
    end if
260   continue
!  now eliminate all channels with eint .ge. emin if any of the channels
!  are open asymptotically but closed at r = rcut
  if (emin .lt. ered) then
    nn = 0
    do 265 i = 1, nch
      if (eint(i) .lt. emin) then
!  here if this channel is to be included
        nn = nn + 1
        eint(nn) = eint(i)
        cent(nn) = cent(i)
        j(nn) = j(i)
        j1(nn) = j1(i)
        j2(nn) = j2(i)
        j12(nn) = j12(i)
        l(nn) = l(i)
      end if
265     continue
!  reset number of channels
    nch = nn
  end if
end if
!  return if no channels
if (nch .eq. 0) return
if (nu .eq. numin) then
  ntop = max(nn, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
  if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) &
       ntop = ntop + 1
else
  if (nn.gt.ntop) then
    write (6, 270) nn, nu, ntop
    write (9, 270) nn, nu, ntop
270     format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3, &
            '; ABORT **',/, &
    '     CHECK RCUT')
    call exit
  endif
end if
!  now list channels if requested

if (clist) then
  if (bastst) write (6,280)
  write(9,280)
280   format(' **  LIST OF ROTATIONAL LEVELS :' &
   /,'   N   J1  J2   E-(CM-1)',/)
  do 300 i = 1 , numj
    jj1 = nj1j2(i) / 10
    jj2  = mod(nj1j2(i), 10)
    if (bastst) write(6,290) i , jj1, jj2, ehold(i)*econv
    write(9,290) i , jj1, jj2, ehold(i)*econv
290     format(i4,i5,i4,f9.3)
300   continue
end if
nlevel = numj
! i counts v2 elements
! inum counts v2 elements for given lambda
! ilam counts number of v2 matrices
! ij is address of given v2 element in present v2 matrix
i  =  0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do 360  ilam = 1, nterm
  ancouma => v2%get_angular_coupling_matrix(ilam)
  if (bastst) write (6, 315) ilam
315   format (' *** ILAM =', i2)
  if (bastst .and. iprint.gt.1) then
     write (9, 310) jlpar, nsym
     write (6, 310) jlpar, nsym
310     format (' CHECK OF VLMTWO:  JLPAR = ',i2,', NSYM = ',i2,/, &
      ' IROW  ICOL  I    IJ  J1RW J2RW J12R LROW', &
      ' J1CL J2CL J12C LCOL JTOT    VEE')
  endif
  inum = 0
  do icol  =  1, nch
    do irow  =  icol,  nch
      !  determine and store angular coupling matrix elements
      call vlmtwo (ilam, jlpar, nsym, j1(irow), j2(irow), &
                   j12(irow), l(irow), j1(icol), j2(icol), &
                   j12(icol), l(icol), jtot, vee)
      if (vee .ne. 0) then
        i = i + 1
        inum = inum + 1
        if (bastst .and. iprint.gt.1) then
          write (6, 320) irow, icol, &
                   i, j1(irow), &
                   j2(irow), j12(irow), l(irow), j1(icol), &
                   j2(icol), j12(icol), l(icol), jtot, vee
          write (9, 320) irow, icol, &
                   i, j1(irow), &
                   j2(irow), j12(irow), l(irow), j1(icol), &
                   j2(icol), j12(icol), l(icol), jtot, vee
320           format (i4, 11i5,1pe16.7)
        endif
        call ancouma%set_element(irow=irow, icol=icol, vee=vee)
      end if
    end do
  end do
  if (bastst .and. iprint.ge.1) then
    write (6, 355) ilam,inum
    write (9, 355) ilam,inum
355     format(' LAMBDA=',i2, &
           ', NUMBER OF NONZERO MATRIX ELEMENTS = ',i5)
  end if
360 continue
if (clist) then
  if (bastst) write(6, 380)
  write(9, 380)
380   format (/' **  LIST OF CHANNELS:', &
   /,'   N   J1  J2  J12  L   E-(CM-1)   E-(HARTREE)',/)
endif
do 400 i = 1 , nch
  if (clist) then
    if (bastst) write (6,390) i, j1(i), j2(i), j12(i), l(i), &
                  eint(i) * econv, eint(i)
    write (9,390) i, j1(i), j2(i), j12(i), l(i), &
                  eint(i) * econv, eint(i)
390     format(i4, i5, i4, i4, i5, f9.3, 1pe15.5)
  endif
  j(i) = 10 * j1(i) + j2(i)
  is(i) = nsym * j12(i)
400 continue
return
end
! --------------------------------------------------------------------
subroutine vlmtwo(ilam,npar,nsym,j1,j2,j12,l,j1p,j2p,j12p,lp,j,vs)
implicit double precision (a-h,o-z)
! ilam = 1: dipole-dipole
! ilam = 2: dipole-quadrupole
! ilam = 3: short-range
if (ilam .eq. 1) then
  lb1=1
  lb2=1
  lb=2
else if (ilam .eq. 2) then
  lb1=1
  lb2=2
  lb=3
else if (ilam .eq. 3) then
  lb1=0
  lb2=1
  lb=1
endif
veelam=f2mol(lb1,lb2,lb,j1,j2,j12,l,j1p,j2p,j12p,lp,j)
if(lb1.eq.lb2) go to 21
veelam=veelam+((-1)**(lb1+lb2))* &
  f2mol(lb2,lb1,lb,j1,j2,j12,l,j1p,j2p,j12p,lp,j)
21 if(j1p.ne.j2p) go to 22
vbar=veelam
go to 23
22 vbar=f2mol(lb1,lb2,lb,j1,j2,j12,l,j2p,j1p,j12p,lp,j)
if(lb1.eq.lb2) go to 23
vbar=vbar+((-1)**(lb1+lb2))* &
  f2mol(lb2,lb1,lb,j1,j2,j12,l,j2p,j1p,j12p,lp,j)
23 d1=2.
if(j1.eq.j2) d1=d1+2
d1=1./sqrt(d1)
d2=2.
if(j1p.eq.j2p) d2=d2+2
d2=1./sqrt(d2)
vs=2.*d1*d2*(veelam+nsym*npar*((-1)**j12p)*vbar)
return
end
! --------------------------------------------------------------------
function f2mol(lb1,lb2,lb,j1,j2,j12,l,j1p,j2p,j12p,lp,j)
! this routine evaluates the non-symmetrized potential matrix
! elements for a given lamda1(lb1),lamda2(lb2), and lamda(lb)
! coefficient index in the expansion of the intermolecular
! potential.  see eq.(16) of alexander and depristo, 'symmetry
! considerations in the quantum treatment of collisions
! between two diatomic molecules'.
! note that (4*pi)**3=1984.40171
implicit double precision (a-h,o-z)
data zero /0.d0/
f2mol=zero
if(iabs(j1-j1p).gt.lb1) return
if((j1+j1p).lt.lb1) return
if(iabs(j2-j2p).gt.lb2) return
if((j2+j2p).lt.lb2) return
if((-1)**(lb1+lb2+lb).lt.0.or.(-1)**(l+lp+lb).lt.0) return
if((-1)**(j1+j1p+lb1).lt.0.or.(-1)**(j2+j2p+lb2).lt.0) return
b1=f6j(l,lp,lb,j12p,j12,j)
if(b1.eq.zero) return
a=(2*j1+1)*(2*j2+1)*(2*j12+1)*(2*l+1)*(2*lb1+1)
ap=(2*j1p+1)*(2*j2p+1)*(2*j12p+1)*(2*lp+1)*(2*lb2+1)
factor=(-1)**(j+j1+j2+j12p)*sqrt((a*ap))*(2*lb+1)
b1=b1*f3j0(j1,j1p,lb1)*f3j0(j2,j2p,lb2)*f3j0(l,lp,lb)
if(b1.eq.zero)  return
! now we evalute a sum of coefficients to get the rest of the
! expansion of the 9j symbol.
kmin=max0(iabs(j12-lb2),iabs(j2p-j1))
kmin=max0(kmin,iabs(j12p-lb1))
kmax=min0((j12+lb2),(j2p+j1))
kmax=min0(kmax,(j12p+lb1))
sum=zero
do 10 k=kmin,kmax
b2=f6j(lb1,lb2,lb,j12,j12p,k)
if(b2.eq.zero) go to 10
b3=f6j(j1,j1p,lb1,j12p,k,j2p)
if(b3.eq.zero) go to 10
b4=f6j(j2,j2p,lb2,k,j12,j1)
if(b4.eq.zero) go  to 10
sum=sum+(2*k+1)*b2*b3*b4
10 continue
f2mol=sum*b1*factor
return
end
! --------------------------------------------------------------
subroutine sy2mol (irpot, readpt, iread)
!  subroutine to read in system dependent parameters for hf-hf collisions
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  authors:  millard alexander and peter vohralik
!  current revision date: 19-aug-1991
!  -----------------------------------------------------------------------
!  variables in common block /cotwo/
!    numj:     number of j1-j2 values
!    nj1j2:    specific j1-j2 values (up to a maximum of 50)
!              N.B. this dimension is set in hiba2mol
!  variables in common block /cosysr/
!    isrcod:   number of real system dependent variables
!    brot:     rotational constant in cm-1
!    drot:     centrifugal distortion constant in cm-1
!    hrot:     next centrifugal distortion constant in cm-1
!  variables in common block /cosysi/
!    nscode:   total number of system dependent variables
!              nscode = isicod + isrcod
!    isicod:   numbr of integer system dependent variables
!    nterm:    number of different lambda terms included in potential
!              see below
!    nsym:     interchange symmetry of included channels
!  variable in module mod_conlam
!    nlam:     the total number of angular coupling terms
!              nlam is set equal to nterm; see above
!     nterm = 1, include lam1=1  lam2=1  lam=2   dipole-dipole
!     nterm = 2, include dipole-dipole and
!                dipole-quadrupole (lam1=1  lam2=2  lam=3)
!     nterm = 3, include dipole-dipole, dipole-quadrupole, and
!                short-range term (lam1=0  lam2=1  lam=1)
!  line 13:
!    nterm:    total number of anistropic terms in potential
!  line 14:
!    brot, drot, hrot:    rotational constants of the molecule in cm-1
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
integer :: i, ihigh, ij, iline, ilow, itop, j, l, lc, nj1j2, numj
logical existf
character*1 dot
character*(*) fname
character*60 filnam, line, potfil, filnm1
common /cotwo/ numj,nj1j2(3)
save potfil
#include "common/comdot.F90"

integer, pointer :: nterm, nsym
real(8), pointer :: brot, drot, hrot
nterm=>ispar(1); nsym=>ispar(2)
brot=>rspar(1); drot=>rspar(2); hrot=>rspar(3)

!     number and names of system dependent parameters
isicod = 2
isrcod = 3
nscode=isicod+isrcod
scod(1) = 'NTERM'
scod(2) = 'NSYM'
scod(3) = 'BROT'
scod(4) = 'DROT'
scod(5) = 'HROT'
!  set default values for hf-hf collisions
nterm = 1
nsym = 1
numj = 3
brot =  20.559741
drot = 2.1204498e-3
hrot = 0.168e-6
nj1j2(1) = 00
nj1j2(2) = 11
nj1j2(3) = 02
if(iread.eq.0) return
!  line 1
iline=1
read (8, *, err=888) nterm, nsym, numj
nlam = nterm
if (numj .le. 0) then
  write (6, 80) numj
80   format(' ** NUMJ = ',i3,'; .LE. O IN SY2MOL; ABORT')
  call exit
endif
!  line 2
iline=2
read (8, *, err=888) brot, drot, hrot
if (numj .le. 20) then
  iline=3
  read (8, *, err=888) (nj1j2(i), i=1,numj)
else
  ilow=1
  ihigh=20
  do  100 ij = 1, numj, 20
    iline=iline+1
    itop=min(ihigh,numj)
    read (8, *, err=888) (nj1j2(i), i=ilow,itop)
    ilow=itop+1
    ihigh=ihigh+20
100   continue
endif
!  line 15 name of file containing potential parameters
iline=iline+1
line=' '
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  return
endif
read (8, 285, end=286) line
285 format (a)
goto 286
! here if read error occurs
888 write(6,1000) iline
1000 format(/'   *** ERROR DURING READ OF SY2MOL INPUT; LINE: ',i2)
return
!
entry ptr2mol (fname,readpt)
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
entry sav2mol (readpt)
!  save input parameters for hf-hf scattering
!  line 13:
write (8, 300) nterm, nsym,numj
300 format (3i4,18x,'   nterm, nsym, numj')
!  line 14
write (8, 320) brot, drot, hrot
320 format(f8.4, 2g11.4, '   brot, drot, hrot')
!  line 15
if (numj .le. 20) then
  iline=3
  write (8, 330) (nj1j2(i)/10,mod(nj1j2(i),10), i=1,numj)
330 format (1x,20(2i1,'  '),t65,'nj1j2')
else
  ilow=1
  ihigh=20
  do  350 ij = 1, numj, 20
    iline=iline+1
    itop=min(ihigh,numj)
    write (8, 330) (nj1j2(i)/10,mod(nj1j2(i),10), &
                    i=ilow,itop)
    ilow=itop+1
    ihigh=ihigh+20
350   continue
endif
return
end
end module mod_hiba08_2mol
