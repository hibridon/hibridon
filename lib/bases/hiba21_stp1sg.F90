#include "assert.h"
! systp1sg (savstp1sg/ptrstp1sg) defines, saves variables and reads      *
!                  potential for symmetric top - singlet sigma system    *
!     Basis subroutine for the collision of a symmetric top with
!     a singlet-Sigma molecule
!
!     system dependent parameters:
!       ipotsy - cylindrical symmetry of the potential.  presently
!          only ipotsy=3 supported (e.g. NH3, CD3)
!       iop  - bitwise flag to set rotational basis (see the .html help
!       file for details)
!          CH3:  ortho levels, iop=2; para levels, iop=4 or 8
!          CD3:  A1 levels, iop=1; A1 levels, iop=2; E levels, iop=4 or 8
!          NH3:  ortho levels, iop=3; para levels, iop=12
!          ND3:  A1 levels, 17; A2 levels, 34; E levels, iop=12 (added ny p.j.d.)
!       j1max - maximum rotational momentum in channel expansion of sym. top
!       e1max - maximum energy of state to be included in sym. top expansion
!       j2min, j2max - minimum and maximum rotational angular momentum of
!          the linear molecule
!       brot, crot - rotational constants of the symmetric top
!       drot - rotational constant of the linear molecule
!       delta - inversion splitting of the symmetric top (assumed to
!          be the same for all levels
!       see the .html help file for the angular expansion of the potential used
!
!     Author: Qianli Ma
!     written:  Jul-2013
!
!     extended to include the 3 nuclear symmetry permutations of ND3
!     (20-nov-2014, p. dagdigian)
!
!     current revision:  24-nov-2014 (p. dagdigian)
!     ------------------------------------------------------------------
!     This module contains (explictly) the number of terms and their
!     indices in the expansion of the PES.  Its contents should be
!     filled in the pot routine.
!     This module replaces lammin, lammax, mproj in hibridon.
module mod_stp1sg
implicit none
!
type lm_type
integer :: l1, l2, l, mu1
end type lm_type
!
type lev_type
integer :: j1, k1, eps1, j2
real(8) :: e
end type lev_type
!
type chn_type
integer :: ilev, j12, l
end type chn_type
!
integer :: nv
type(lm_type), dimension(:), allocatable :: lms
type(lev_type), dimension(:), allocatable :: levs
type(chn_type), dimension(:), allocatable :: chns
end module mod_stp1sg
!     ------------------------------------------------------------------
module mod_hiba21_stp1sg
   contains
!     ------------------------------------------------------------------
   real(8) function vstp1sg(j1p, lp, j1, l, j2p, j2, j12p, j12, &
   jtot, kp, k, lam1, mu1, lam2, lam, epsp, eps)
implicit none
integer :: j1p, kp, epsp, j2p, j12p, lp, j1, k, eps, j2, j12, l, &
   jtot, lam1, mu1, lam2, lam
integer :: iphase, tj1p, tj1, tlam1, tkp, tmu1, tk
real(8) :: pref, threej, sixj, ninej, tf3jm0, tf3j, f6j, f9j
real(8), parameter :: machep = epsilon(0d0)
!
vstp1sg = 0d0
iphase = epsp * eps * (-1) ** (j1p + j1 + lam2 + lam + mu1)
if (iphase .eq. -1) return
threej = tf3jm0(2 * j2p, 2 * lam2, 2 * j2)
if (dabs(threej) .lt. machep) return
threej = threej * tf3jm0(2 * lp, 2 * lam, 2 * l)
if (dabs(threej) .lt. machep) return
!
tj1p = 2 * j1p
tlam1 = 2 * lam1
tj1 = 2 * j1
tkp = 2 * kp
tmu1 = 2 * mu1
tk = 2 * k
threej = threej * (tf3j(tj1p, tlam1, tj1, -tkp, tmu1, tk) &
   + epsp * eps * tf3j(tj1p, tlam1, tj1, tkp, tmu1, -tk) &
   + epsp * tf3j(tj1p, tlam1, tj1, tkp, tmu1, tk) &
   + eps * tf3j(tj1p, tlam1, tj1, -tkp, tmu1, -tk))
if (dabs(threej) .lt. machep) return
!
sixj = f6j(j12, l, jtot, lp, j12p, lam)
if (dabs(sixj) .lt. machep) return
ninej = f9j(j1, j2, j12, j1p, j2p, j12p, lam1, lam2, lam)
if (dabs(ninej) .lt. machep) return
!
iphase = (-1) ** (jtot + lam1 - lam2 + j1 - j2 + j12p - l - lp + k &
   + mu1)
if (mu1 .eq. 0) then
 pref = 0.5d0
else
 pref = 1d0
end if
if (kp .eq. 0) pref = pref * dsqrt(0.5d0)
if (k .eq. 0) pref = pref * dsqrt(0.5d0)
pref = pref * dsqrt( &
     dble((tj1p + 1)) &
   * dble((2 * j2p + 1)) &
   * dble((2 * j12p + 1) &
   * dble((2 * lp + 1)) &
   * dble((tj1 + 1)) &
   * dble((2 * j2 + 1)) &
   * dble((2 * j12 + 1)) &
   * (dble(2 * l + 1)) &
   * dble((2 * lam + 1))) &
)
!
vstp1sg = pref * dble(iphase) * threej * sixj * ninej
return
end function vstp1sg
 
subroutine bastp1sg(jchn, lchn, ischn, jlev, elev, islev, nlev, &
     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, bastst, &
     ihomo, nu, numin, jlpar, twomol, nchn, nmax, ntop, v2)
use mod_stp1sg
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cchn => cent
use mod_coeint, only: echn => eint
use mod_coj12, only: j12chn => j12
use mod_conlam, only: nlam
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use mod_hibasutil, only: raise
use mod_par, only: iprint
use mod_ered, only: ered, rmu
implicit none
!
!     The following arrays store the parameters of channels and levels.
!     Note that j12 is stored in a common block
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
integer :: i, ilev, iv, irow, icol, inum, i1, i2
real(8) :: vee
real(8), parameter :: machep=epsilon(0d0)
!
integer, pointer :: ipotsy, iop, ipotsy2, j1max, j2min, j2max
real(8), pointer :: brot, crot, delta, e1max, drot
ipotsy=>ispar(1); iop=>ispar(2); ipotsy2=>ispar(3); j1max=>ispar(4); j2min=>ispar(5); j2max=>ispar(6)
brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); e1max=>rspar(4); drot=>rspar(5); 

if (flaghf) &
     call raise('FLAGHF = .TRUE. FOR SINGLET SYSTEM')
if (ihomo) call raise('IHOMO NOT USED IN THIS BASIS,' &
     //' PLEASE SET IT FALSE')
if (flagsu) &
     call raise('FLAGSU = .TRUE. FOR MOL-MOL COLLISION')
if (csflag) &
     call raise('CS CALCULATION NOT IMPLEMENTED')
!$$$      if (rcut .ge. 0d0 .and. .not. bastst)
!$$$     $     call raise('RCUT NOT IMPLEMENTED, PLEASE SET IT '
!$$$     $     // 'NEGATIVE.')
if (.not. twomol) &
     call raise('TWOMOL IS FALSE FOR MOL-MOL COLLISION.')
!
call genlev_stp1sg(ipotsy, iop, ipotsy2, j1max, j2min, j2max, &
     brot, crot, delta, drot, e1max)
call sortlev_stp1sg()
if (bastst) call prtlev_stp1sg()
call genchn_stp1sg(jtot, jlpar)
if (bastst .and. iprint .ge. 1) call prtchn_stp1sg()
!
!     Copy level and channel info to hibridon arrays
nlev = size(levs)
nchn = size(chns)
ntop = max(nchn, nlev)
if (nchn .gt. nmax) &
     call raise('TOO MANY CHANNELS.')
if (nchn .eq. 0) return
do ilev = 1, nlev
   jlev(ilev) = levs(ilev)%j1 * 10 + levs(ilev)%j2
   islev(ilev) = -levs(ilev)%eps1 * (-1) ** levs(ilev)%j1 &
        * levs(ilev)%k1
   elev(ilev) = levs(ilev)%e
end do
do i = 1, nchn
   jchn(i) = jlev(chns(i)%ilev)
   ischn(i) = islev(chns(i)%ilev)
   j12chn(i) = chns(i)%j12
   lchn(i) = chns(i)%l
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
         vee = vstp1sg(levs(i1)%j1, chns(icol)%l, levs(i2)%j1, &
              chns(irow)%l, levs(i1)%j2, levs(i2)%j2, &
              chns(icol)%j12, chns(irow)%j12, jtot, &
              levs(i1)%k1, levs(i2)%k1, lms(iv)%l1, lms(iv)%mu1, &
              lms(iv)%l2, lms(iv)%l, levs(i1)%eps1, levs(i2)%eps1)
         if (dabs(vee) .lt. machep) cycle
         i = i + 1
         inum = inum + 1
         call ancouma%set_element(irow=irow, icol=icol, vee=vee)
         if (bastst .and. iprint .ge. 2) write (6, 346) iv, &
              lms(iv)%l1, lms(iv)%l2, lms(iv)%l, lms(iv)%mu1, &
              icol, irow, i, vee
346          format(i4, 4i3, 2i4, i5, f17.10)
      end do
   end do
   if (bastst) write (6, 347) iv, lms(iv)%l1, lms(iv)%l2, &
        lms(iv)%l, lms(iv)%mu1, ancouma%get_num_nonzero_elements()
347    format (' ILAM=', i3, '  L1=', i3, '  L2=', i3, &
        '  L=', i3, '  MU=', i3, '  LAMNUM(ILAM) = ', i6)
end do
if (bastst) write (6, 350) i
350 format ('TOTAL NUMBER OF NON-ZERO V2 TERMS: ', i6)
return
end subroutine bastp1sg
!     ------------------------------------------------------------------
subroutine genchn_stp1sg(jtot, jlpar)
use mod_stp1sg
implicit none
integer, intent(in) :: jtot, jlpar
integer :: nchn, ichn, ilev, j12, l
!     First pass (nchn = -1): Count the number of channels
!     Second pass (nchn is the number of channels): Save channel info
nchn = -1
do while (.true.)
   ichn = 0
   if (nchn .ne. -1) then
      if (allocated(chns)) deallocate(chns)
      allocate(chns(nchn))
   end if
   do ilev = 1, size(levs)
      do j12 = iabs(levs(ilev)%j1 - levs(ilev)%j2), &
           levs(ilev)%j1 + levs(ilev)%j2
         do l = iabs(jtot - j12), jtot + j12
            if (jlpar .ne. -levs(ilev)%eps1 * &
                 (-1) ** (levs(ilev)%k1 + levs(ilev)%j1 &
                 + levs(ilev)%j2 + l - jtot)) cycle
            ichn = ichn + 1
            if (nchn .eq. -1) cycle
            chns(ichn)%ilev = ilev
            chns(ichn)%j12 = j12
            chns(ichn)%l = l
         end do
      end do
   end do
   if (nchn .eq. -1) then
      nchn = ichn
   else
      exit
   end if
end do
return
end subroutine genchn_stp1sg
!     ------------------------------------------------------------------
subroutine sortlev_stp1sg()
use mod_stp1sg
use constants, only: econv
implicit none
type(lev_type) :: temp
type(lev_type), dimension(:), allocatable :: levs1
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
return
end subroutine sortlev_stp1sg
!     ------------------------------------------------------------------
subroutine genlev_stp1sg(ipotsy, iop, ipotsy2, j1max, j2min, &
     j2max, brot, crot, delta, drot, e1max)
use mod_stp1sg
use mod_hibasutil, only: raise
use constants, only: econv
implicit none
integer, intent(in) :: ipotsy, iop, ipotsy2, j1max, j2min, j2max
real(8), intent(in) :: brot, crot, delta, drot, e1max
integer :: j1, k1, eps1, j2, ilev, nlev, igroup
real(8) :: roteng
!
!     Only C3v/D3h symmetric top implemented
if (ipotsy .ne. 3) call raise('ONLY C3v/D3h SYMMETRIC ' &
     //'TOP IMPLEMENTED')
!     First pass (nlev .eq. -1): Count the number of levels
!     Second pass (nlev is the number of channels): Save level info
nlev = -1
do while (.true.)
   if (nlev .ne. -1) then
      if (allocated(levs)) deallocate(levs)
      allocate(levs(nlev))
   end if
   ilev = 0
   do j1 = 0, j1max
      do k1 = 0, j1
         do eps1 = -1, 1, 2
            if (k1 .eq. 0 .and. eps1 .eq. -1) cycle
            if (mod(k1, 3) .eq. 0) then
               igroup = 0
            else
               igroup = 2
            end if
            if ((k1.ne.0) .and. (mod(k1,3).eq.0) .and. &
               ((iop.eq.17) .or. (iop.eq.34))) goto 25
            if (eps1 * (-1) ** j1 .eq. 1) then
               igroup = igroup + 1
               if (k1.eq.0 .and. iop.eq.17) goto 25
            else
               if (k1.eq.0 .and. iop.eq.34) goto 25
            end if
            if (.not. btest(iop, igroup)) cycle
25             roteng = brot * j1 * (j1 + 1) &
                 + (crot - brot) * k1 ** 2
            if (eps1 * (-1) ** j1 .eq. 1) then
               if (iop.eq.17 .and. k1.eq.0) then
                  roteng = roteng
               else
                  roteng = roteng + delta
               end if
            else
               if (iop.eq.17 .and. k1.eq.0) then
                  roteng = roteng + delta
               end if
            end if
            if (roteng .gt. e1max) cycle
            if (nlev .eq. -1) then
               ilev = ilev + (j2max - j2min) / ipotsy2 + 1
               cycle
            end if
            do j2 = j2min, j2max, ipotsy2
               ilev = ilev + 1
               levs(ilev)%j1 = j1
               levs(ilev)%k1 = k1
               levs(ilev)%eps1 = eps1
               levs(ilev)%j2 = j2
               levs(ilev)%e = (roteng + drot * j2 * (j2 + 1)) &
                    / econv
            end do
         end do
      end do
   end do
   if (nlev .eq. -1) then
      nlev = ilev
   else
      exit
   end if
end do
return
end subroutine genlev_stp1sg
!     ------------------------------------------------------------------
subroutine prtlev_stp1sg()
use mod_stp1sg
use constants, only: econv
implicit none
integer :: i
write (6, 10)
10 format (/' ** CC SYMMETRIC TOP-LINEAR MOLECULE' &
   //10x, 'SORTED LEVEL LIST'/ &
   '    N   J1   K1  EPS1  J2   ENT(CM-1)')
do i = 1, size(levs)
   write (6, 20) i, levs(i)%j1, levs(i)%k1, levs(i)%eps1, &
        levs(i)%j2, levs(i)%e * econv
end do
20 format(5i5,f11.3)
return
end subroutine prtlev_stp1sg
!     ------------------------------------------------------------------
subroutine prtchn_stp1sg()
use mod_stp1sg
use constants, only: econv, xmconv
implicit none
integer :: i, ilev
write (6, 10)
10 format (/, ' ** CHANNEL LIST')
do i = 1, size(chns)
   ilev = chns(i)%ilev
   write (6, 280) i, levs(ilev)%j1, levs(ilev)%k1, &
        levs(ilev)%eps1, levs(ilev)%j2, chns(i)%j12, &
        chns(i)%l, levs(ilev)%e * econv
280    format (i4, 1x, 'J1=', i2, ' K1=', i2, ' EPS1=', i2, &
        2x, 'J2=', i2, 2x, 'J12=', i2, 1x, 'L=', i3, 2x, &
        'E=', f9.4)
end do
write (6, *)
return
end subroutine prtchn_stp1sg
!     ------------------------------------------------------------------
!     THIS SUBROUTINE GOVERNS THE INPUT/OUTPUT OF THE BASIS ROUTINE.
!     ONLY IREAD IS USED: RETURN DIRECTLY IF ZERO.
subroutine systp1sg(irpot, readpt, iread)
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use mod_hibasutil, only: raise
use funit, only: FUNIT_INP
implicit none
!
integer, intent(out) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
character*(*) fname
!     NUMBER OF BASIS-SPECIFIC VARIABLES, MODIFY ACCORDINGLY.
integer icod, ircod
parameter (icod=6, ircod=5)
character*40 potfil
save potfil

integer, pointer :: ipotsy, iop, ipotsy2, j1max, j2min, j2max
real(8), pointer :: brot, crot, delta, e1max, drot
ipotsy=>ispar(1); iop=>ispar(2); ipotsy2=>ispar(3); j1max=>ispar(4); j2min=>ispar(5); j2max=>ispar(6)
brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); e1max=>rspar(4); drot=>rspar(5); 

!     DEFINE THE NAMES HERE
scod(1)='IPOTSY'
scod(2)='IOP'
scod(3)='IPOTSY2'
scod(4)='J1MAX'
scod(5)='J2MIN'
scod(6)='J2MAX'
scod(7)='BROT'
scod(8)='CROT'
scod(9)='DELTA'
scod(10)='E1MAX'
scod(11)='DROT'
nscode = icod+ircod
isicod = icod
isrcod = ircod
!     KEEP THE FOLLOWING LINE
if (iread .eq. 0) return
!     READ THE LAST FEW LINES OF THE INPUT FILE
read (8, *, err=80) ipotsy, iop
read (8, *, err=80) j1max, e1max
read (8, *, err=80) j2min, j2max, ipotsy2
read (8, *, err=80) brot, crot, delta
read (8, *, err=80) drot
read (8, *, err=80) potfil
call loapot(10, potfil)
close (8)
return
80 call raise('error read from input file.')
return
!     ------------------------------------------------------------------
entry ptrstp1sg(fname, readpt)
return
!     ------------------------------------------------------------------
entry savstp1sg(readpt)
!     WRITE THE LAST FEW LINES OF THE INPUT FILE.
write (FUNIT_INP, 220) ipotsy, iop
220 format (2i4, 25x,'ipotsy, iop')
write (FUNIT_INP, 230) j1max, e1max
230 format (i4, f8.2, 21x, 'j1max, e1max')
write (FUNIT_INP, 231) j2min, j2max, ipotsy2
231 format (3i4, 21x,'j2min, j2max, ipotsy2')
write (FUNIT_INP, 250) brot, crot, delta
250 format (3f10.4, 3x, 'brot, crot, delta, emax' )
write (FUNIT_INP, 251) drot
251 format (f10.4, 23x,'drot')
write (FUNIT_INP, *) potfil
return
end
!     ------------------------------------------------------------------
end module mod_hiba21_stp1sg
