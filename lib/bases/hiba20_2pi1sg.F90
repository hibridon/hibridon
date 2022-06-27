#include "assert.h"
! sy2pi1sg (sav2pi1sg/ptr2pi1sg) defines, saves variables and            *
!                  reads potential for doublet pi-singlet sigma* system  *
!     Basis subroutine for the collision of a doublet-Pi molecule with
!     a singlet-Sigma molecule.
!
!     Author: Qianli Ma
!
!     ------------------------------------------------------------------
!     This module contains (explictly) the number of terms and their
!     indices in the expansion of the PES.  Its contents should be
!     filled in the pot routine.
!     This module replaces lammin, lammax, mproj in hibridon.
module mod_2pi1sg
implicit none
!
type lm_type
integer :: l1, l2, l
logical :: is_diag
end type lm_type
!
type lev_type
integer :: tj1, fi1, eps1, tj2
real(8) :: c32, c12, e
end type lev_type
!
type chn_type
integer :: ilev, tj12, tl
end type chn_type
!
integer :: nv
type(lm_type), dimension(:), allocatable :: lms
type(lev_type), dimension(:), allocatable :: levs
type(chn_type), dimension(:), allocatable :: chns
end module mod_2pi1sg

module mod_hiba20_2pi1sg
contains

!     ------------------------------------------------------------------
double precision function v2pi1sg(tjtot, tj1p, eps1p, c32p, c12p, &
     tj2p, tj12p, tlp, tj1, eps1, c32, c12, tj2, tj12, tl, &
     tlam1, tlam2, tlam, isdiag)
implicit none
!
!     The subroutine calculate the coupling matrix elements, shown in
!     Eq. (27) in the notes of Q. Ma.
!
!     If (omeg1p .eq. omeg1), the coefficient before B is calculated;
!     otherwise the coefficient before F is calculated.
!
integer :: tjtot, tj1p, eps1p, tj2p, tj12p, tlp, tj1, eps1, &
     tj2, tj12, tl, tlam1, tlam2, tlam
real(8) :: c12p, c32p, c12, c32
logical :: isdiag
integer :: iphase
real(8) :: phase, pref, threej, sixj, ninej
real(8) :: tf3j, tf6j, tf9j, tf3jm0
real(8), parameter :: machep=epsilon(0d0)
!
v2pi1sg = 0d0
!
iphase = eps1p * eps1 * (-1) ** ((tj1p + tj1 + tlam1) / 2)
if (iphase .eq. 1) return
!
threej = tf3jm0(tj2p, tlam2, tj2) * tf3jm0(tlp, tlam, tl)
if (dabs(threej) .lt. machep) return
!     omega-dependent part
if (isdiag) then
   threej = threej * &
        (c12p * c12 * tf3j(tj1p, tlam1, tj1, -1, 0, 1) &
        - c32p * c32 * tf3j(tj1p, tlam1, tj1, -3, 0, 3))
else
   threej = threej * dble(eps1) * &
        (c12p * c32 * tf3j(tj1p, tlam1, tj1, -1, 4, -3) &
        - c32p * c12 * tf3j(tj1p, tlam1, tj1, -3, 4, -1))
end if
if (dabs(threej) .lt. machep) return
!
sixj = tf6j(tj12, tl, tjtot, tlp, tj12p, tlam)
if (dabs(sixj) .lt. machep) return
ninej = tf9j(tj1, tj2, tj12, tj1p, tj2p, tj12p, &
     tlam1, tlam2, tlam)
if (dabs(ninej) .lt. machep) return
!
!     Again 1 is added to compensate the dropped half-integer part
iphase = (tjtot + tlam1 - tlam2 + tj1 - tj2 + tj12p - tl - tlp &
     - 1) / 2
if (mod(iphase, 2) .eq. 0) then
   phase = 1d0
else
   phase = -1d0
end if
!
pref = (tj1p + 1d0) * (tj2p + 1d0) * (tj12p + 1d0) * (tlp + 1d0) &
     * (tj1 + 1d0) * (tj2 + 1d0) * (tj12 + 1d0) * (tl + 1d0) &
     * (tlam + 1d0)
pref = dsqrt(pref)
!
v2pi1sg = phase * pref * threej * sixj * ninej
return
end

!     ------------------------------------------------------------------
subroutine ba2pi1sg(jchn, lchn, ischn, jlev, elev, islev, nlev, &
     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, bastst, &
     ihomo, nu, numin, jlpar, twomol, nchn, nmax, ntop, v2)
use mod_2pi1sg
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cchn => cent
use mod_coeint, only: echn => eint
use mod_coj12, only: j12chn => j12
use mod_conlam, only: nlam
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar 
use mod_hibasutil, only: raise
use mod_par, only: iprint
implicit none
!
!     The following arrays store the parameters of channels and levels.
integer :: jchn(*), lchn(*), ischn(*), jlev(*), islev(*)
real(8) :: elev(*)
!     The following parameters store the number of channels, opened
!     levels and levels
integer :: nchn, nlevop, nlev
!     Various parameters used through out hibridon
integer :: jtot, jlpar, nu, numin
logical :: flaghf, flagsu, csflag, clist, bastst, ihomo, twomol
!     Limits of channels
integer :: nmax, ntop
!
real(8) :: rcut
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
!
common /coered/ ered, rmu
real(8) :: ered, rmu
!
integer :: i, ilev, iv, irow, icol, inum, i1, i2
real(8) :: vee
real(8), parameter :: machep=epsilon(0d0)
!
integer, pointer :: j1max, npar, j2min, j2max, iptsy2
real(8), pointer :: brot, aso, p, q, drot
j1max=>ispar(1); npar=>ispar(2); j2min=>ispar(3); j2max=>ispar(4); iptsy2=>ispar(5)
brot=>rspar(1); aso=>rspar(2); p=>rspar(3); q=>rspar(4); drot=>rspar(5)

if (.not. flaghf) &
     call raise('FLAGHF = .FALSE. FOR DOUBLET SYSTEM')
if (ihomo) call raise('HOMONUCLEAR 2PI NOT IMPLEMENTED')
if (flagsu) &
     call raise('FLAGSU = .TRUE. FOR MOL-MOL COLLISION')
if (csflag) &
     call raise('CS CALCULATION NOT IMPLEMENTED')
if (rcut .ge. 0d0) &
     call raise('RCUT NOT IMPLEMENTED, PLEASE SET IT ' &
     // 'NEGATIVE.')
if (.not. twomol) &
     call raise('TWOMOL IS FALSE FOR MOL-MOL COLLISION.')
!
call genlev_2pi1sg(2 * j1max + 1, 2 * j2min, 2 * j2max, &
     2 * iptsy2, npar, brot, aso, p, q, drot)
call sortlev_2pi1sg()
if (bastst) call prtlev_2pi1sg()
call genchn_2pi1sg(2 * jtot + 1, jlpar)
if (bastst .and. iprint .ge. 1) call prtchn_2pi1sg()
!
!     Copy level and channel info to hibridon arrays
nlev = size(levs)
nchn = size(chns)
ntop = max(nchn, nlev)
if (nchn .gt. nmax) &
     call raise('TOO MANY CHANNELS.')
if (nchn .eq. 0) return
do ilev = 1, nlev
   jlev(ilev) = levs(ilev)%tj1 / 2 * 10 + levs(ilev)%tj2 / 2
   islev(ilev) = levs(ilev)%eps1 * levs(ilev)%fi1
   elev(ilev) = levs(ilev)%e
end do
do i = 1, nchn
   jchn(i) = jlev(chns(i)%ilev)
   ischn(i) = islev(chns(i)%ilev)
   j12chn(i) = chns(i)%tj12 / 2
   lchn(i) = chns(i)%tl / 2
   cchn(i) = dble((lchn(i) * (lchn(i) + 1)))
   echn(i) = elev(chns(i)%ilev)
end do
!
!     Count number of asymtotically open levels
nlevop = 0
do ilev = 1, nlev
   if (elev(ilev) .le. ered) then
      nlevop = nlevop + 1
   else
      exit
   end if
end do
!
!     Calculate coupling matrix elements
nlam = nv
i = 0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do iv = 1, nv
   ancouma => v2%get_angular_coupling_matrix(iv)
   inum = 0
   do icol = 1, nchn
      i1 = chns(icol)%ilev
      do irow = icol, nchn
         i2 = chns(irow)%ilev
         vee = v2pi1sg(2*jtot+1, levs(i1)%tj1, levs(i1)%eps1, &
              levs(i1)%c32, levs(i1)%c12, levs(i1)%tj2, &
              chns(icol)%tj12, chns(icol)%tl, levs(i2)%tj1, &
              levs(i2)%eps1, levs(i2)%c32, levs(i2)%c12, &
              levs(i2)%tj2, chns(irow)%tj12, chns(irow)%tl, &
              2*lms(iv)%l1, 2*lms(iv)%l2, 2*lms(iv)%l, &
              lms(iv)%is_diag)
         if (dabs(vee) .lt. machep) cycle
         i = i + 1
         inum = inum + 1
         call ancouma%set_element(irow=irow, icol=icol, vee=vee)
         if (bastst .and. iprint .ge. 2) write (6, 346) iv, &
              lms(iv)%l1, lms(iv)%l2, lms(iv)%l, icol, irow, &
              i, vee
346          format(i4, 3i3, 2i4, i5, f17.10)
      end do
   end do
   if (bastst) write (6, 347) iv, lms(iv)%l1, lms(iv)%l2, &
        lms(iv)%l, ancouma%get_num_nonzero_elements()
347    format ('ILAM=', i3, ' LAM1=', i3, ' LAM2=', i3, &
        ' LAM=', i3, ' LAMNUM(ILAM) = ', i6)
end do
if (bastst) write (6, 350) i
350 format ('TOTAL NUMBER OF NON-ZERO V2 TERMS: ', i6)
return
end subroutine ba2pi1sg
!     ------------------------------------------------------------------
subroutine genchn_2pi1sg(tjtot, jlpar)
use mod_2pi1sg
implicit none
integer, intent(in) :: tjtot, jlpar
integer :: nchn, ilev, tj12, tl
!     Count number of channels
nchn = 0
do ilev = 1, size(levs)
   do tj12 = iabs(levs(ilev)%tj1 - levs(ilev)%tj2), &
        levs(ilev)%tj1 + levs(ilev)%tj2, 2
      do tl = iabs(tjtot - tj12), tjtot + tj12, 2
         nchn = nchn + 1
      end do
   end do
end do
nchn = nchn / 2
!
if (allocated(chns)) deallocate(chns)
allocate(chns(nchn))
nchn = 0
do ilev = 1, size(levs)
   do tj12 = iabs(levs(ilev)%tj1 - levs(ilev)%tj2), &
        levs(ilev)%tj1 + levs(ilev)%tj2, 2
      do tl = iabs(tjtot - tj12), tjtot + tj12, 2
         if (jlpar .ne. levs(ilev)%eps1 * &
              (-1) ** ((levs(ilev)%tj1 + levs(ilev)%tj2 &
              - tjtot + tl) / 2)) cycle
         nchn = nchn + 1
         chns(nchn)%ilev = ilev
         chns(nchn)%tj12 = tj12
         chns(nchn)%tl = tl
      end do
   end do
end do
return
end subroutine genchn_2pi1sg
!     ------------------------------------------------------------------
subroutine sortlev_2pi1sg()
use mod_2pi1sg
implicit none
type(lev_type) :: temp
integer :: i, j
real(8) :: esave
do i = 1, size(levs) - 1
   esave = levs(i)%e
   do j = i + 1, size(levs)
      if (levs(j)%e .lt. esave) then
         esave = levs(j)%e
         temp = levs(j)
         levs(j) = levs(i)
         levs(i) = temp
      end if
   end do
end do
esave = levs(1)%e
do i = 1, size(levs)
   levs(i)%e = levs(i)%e - esave
end do
end subroutine sortlev_2pi1sg
!     ------------------------------------------------------------------
subroutine genlev_2pi1sg(tj1max, tj2min, tj2max, tipotsy2, npar, &
     brot, aso, p, q, drot)
use mod_2pi1sg
use mod_hibasutil, only: raise
use constants, only: econv
implicit none
integer, intent(in) :: tj1max, tj2min, tj2max, tipotsy2, npar
real(8), intent(in) :: brot, aso, p, q, drot
integer :: tj1, tj2, eps1, fi1, ilev, nlev
real(8) :: x, o11, o22, o12, tho, roteng
!
nlev = 2 * tj1max * ((tj2max - tj2min) / tipotsy2 + 1)
if (iabs(npar) .eq. 1) nlev = nlev / 2
if (allocated(levs)) deallocate(levs)
allocate(levs(nlev))
ilev = 0
do tj1 = 1, tj1max, 2
   x = dble((tj1 + 1) / 2)
   do eps1 = -1, 1, 2
      if (npar .eq. 1 .and. eps1 .eq. -1) cycle
      if (npar .eq. -1 .and. eps1 .eq. 1) cycle
!     o11 is the matrix element for omega=3/2
      o11 = 0.5d0 * aso + (x ** 2 - 2d0) * brot &
           + 0.5d0 * (x ** 2 - 1d0) * q
!     o22 is the matrix element for omega=1/2
      o22 = -0.5d0 * aso + x ** 2 * brot &
           + 0.5d0 * (1d0 - eps1 * x) * p &
           + 0.5d0 * (1d0 - eps1 * x) ** 2 * q
!     o12 is the off-diagonal matrix element
      o12 = -dsqrt(x ** 2 - 1d0) * brot &
           - 0.25d0 * dsqrt(x ** 2 - 1d0) * p &
           - 0.5d0 * (1d0 - eps1 * x) * dsqrt(x ** 2 - 1d0) * q
!     tho is the mixing angle between omega=3/2 and omega=1/2 states
      tho = 0.5d0 * datan(2 * o12 / (o11 - o22))
      do fi1 = 1, 2
         if (tj1 .eq. 1 .and. fi1 .eq. 1) cycle
         if (fi1 .eq. 1) then
            roteng = o11 * dcos(tho) ** 2 + o22 * dsin(tho) ** 2 &
                 + o12 * dsin(2 * tho)
         else
            roteng = o11 * dsin(tho) ** 2 + o22 * dcos(tho) ** 2 &
                 - o12 * dsin(2 * tho)
         end if
         do tj2 = tj2min, tj2max, tipotsy2
            ilev = ilev + 1
            levs(ilev)%tj1 = tj1
            levs(ilev)%fi1 = fi1
            levs(ilev)%eps1 = eps1
            levs(ilev)%tj2 = tj2
            levs(ilev)%e = (roteng + drot * (tj2 * (tj2 + 2) / 4)) &
                 / econv
            if (fi1 .eq. 1) then
               levs(ilev)%c32 = dcos(tho)
               levs(ilev)%c12 = dsin(tho)
            else
               levs(ilev)%c32 = dsin(tho)
               levs(ilev)%c12 = -dcos(tho)
            end if
         end do
      end do
   end do
end do
if (ilev .ne. nlev) call raise('NUMBER OF CHANNELS ' &
     // 'INCONSISTENT WITH PREDICTION, CHECK PARAMETERS')
return
end subroutine genlev_2pi1sg
!     ------------------------------------------------------------------
subroutine prtlev_2pi1sg()
use mod_2pi1sg
use constants, only: econv
implicit none
integer :: i
write (6, 10)
10 format (/, ' ** LEVEL LIST')
do i = 1, size(levs)
   write (6, 20) i, levs(i)%tj1, levs(i)%fi1, levs(i)%eps1, &
        levs(i)%tj2 / 2, levs(i)%e * econv, levs(i)%c32, &
        levs(i)%c12
end do
20 format (i3, ' J1=', i2, '/2 F', i1, ' EPS1=', i2, '  J2=', i2, &
     '  E=', f9.4, '  C3/2=', f7.4, ' C1/2=', f7.4)
return
end subroutine prtlev_2pi1sg
!     ------------------------------------------------------------------
subroutine prtchn_2pi1sg()
use mod_2pi1sg
use constants, only: econv
implicit none
integer :: i, ilev
character :: ifchar
write (6, 10)
10 format (/, ' ** CHANNEL LIST')
do i = 1, size(chns)
   ilev = chns(i)%ilev
   if (levs(ilev)%eps1 .eq. 1) then
      ifchar = 'e'
   else
      ifchar = 'f'
   end if
   write (6, 280) i, levs(ilev)%tj1, levs(ilev)%fi1, ifchar, &
        levs(ilev)%tj2 / 2, chns(i)%tj12, chns(i)%tl / 2, &
        levs(ilev)%c32, levs(ilev)%c12, levs(ilev)%e * econv
280    format (i4, 1x, 'J1=', i2, '/2 F', i1, a, 2x, &
        'J2=', i2, 2x, 'J12=', i2, '/2  L=', i3, 2x, &
        'C3/2=', f7.4, 1x, 'C1/2=', f7.4, 2x, &
        'E=', f9.4)
end do
write (6, *)
return
end subroutine prtchn_2pi1sg
!     ------------------------------------------------------------------
!  ---------------------------------------------------------------------
subroutine sy2pi1sg(irpot, readpt, iread)
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use mod_hibasutil, only: raise
use funit, only: FUNIT_INP
implicit none
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
character*(*) fname
!     NUMBER OF BASIS-SPECIFIC VARIABLES, MODIFY ACCORDINGLY.
integer icod, ircod
parameter (icod=5, ircod=5)

character*40 potfil
save potfil

integer, pointer :: j1max, npar, j2min, j2max, iptsy2
real(8), pointer :: brot, aso, p, q, drot
j1max=>ispar(1); npar=>ispar(2); j2min=>ispar(3); j2max=>ispar(4); iptsy2=>ispar(5)
brot=>rspar(1); aso=>rspar(2); p=>rspar(3); q=>rspar(4); drot=>rspar(5)

!     DEFINE THE NAMES HERE
scod(1)='J1MAX'
scod(2)='NPAR'
scod(3)='J2MIN'
scod(4)='J2MAX'
scod(5)='IPTSY2'
scod(6)='BROT'
scod(7)='ASO'
scod(8)='P'
scod(9)='Q'
scod(10)='DROT'
nscode = icod+ircod
isicod = icod
isrcod = ircod
!     KEEP THE FOLLOWING LINE
if (iread .eq. 0) return
!     READ THE LAST FEW LINES OF THE INPUT FILE
read (8, *, err=80) j1max, npar
read (8, *, err=80) j2min, j2max, iptsy2
read (8, *, err=80) brot, aso, p, q
read (8, *, err=80) drot
read (8, *, err=80) potfil
call loapot(10, potfil)
close (8)
return
80 call raise('ERROR READING FROM INPUT FILE.')
return
entry ptr2pi1sg(fname, readpt)
return
entry sav2pi1sg(readpt)
!     WRITE THE LAST FEW LINES OF THE INPUT FILE.
write (FUNIT_INP, 230) j1max, npar
230 format (2i4, 22x, '   j1max, npar')
write (FUNIT_INP,231) j2min, j2max, iptsy2
231 format (3i4, 18x,'   j2min, j2max, iptsy2')
write (FUNIT_INP, 250) brot, aso, p, q
250 format (4(f10.4, 1x), 'brot, aso, p, q' )
write (FUNIT_INP, 251) drot
251 format (f12.6, 18x,'   drot')
write (FUNIT_INP, *) potfil
return
end
end module mod_hiba20_2pi1sg
