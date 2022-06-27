!     pot_ohh2_2013.f
!     authors: Qianli Ma
!
!     Pot routine for a collision between OH and H2.
!
!     PES to be read from data file specified in the input file.
!
!     This basis routine should have ibasty = 100, so hibridon will
!     treat it as a mol-mol problem with j12.
!
!     DUMMY SUBROUTINES FOR ALL USER-DEFINED BASIS/POT.
#include "assert.h"
#include "common/ground.F90"
!
!     ------------------------------------------------------------------
!     THE FOLLOWING SOUBROUTINE WILL BE THE MAIN FUNCTION FOR MAKEPOT.
!     NOTE THAT VVL IS IN HARTREES.
subroutine driver
!
#include "pot_ohh2_bausr_common.F90"
character*40 filenm
double precision r, vv0
integer i
!
print *, 'Input the name of potdata file, without _b/_f:'
read (5, *, end=99) filenm
call loapot(10, filenm)
10 print *, 'R (bohr), Ctrl+D to exit:'
read (5, *, end=99) r
call pot(vv0, r)
write (6, *) 'B-coefficients'
write (6, 20) (lam1b(i), lam2b(i), lamb(i), &
     vvl(i) * econv, i=1, nvb)
write (6, *) 'F-coefficients'
write (6, 20) (lam1f(i), lam2f(i), lamf(i), &
     vvl(nvb + i) * econv, i=1, nvf)
20 format (3(3(i2, 1x), 1pe16.8, 2x))
goto 10
99 end
!     ------------------------------------------------------------------
!     DATA FILE, IF REQUIRED, CAN BE LOADED WITH THIS SOUBROUTINE.
subroutine loapot(iunit, filnam)
use mod_hibasutil, only: raise
!
#include "pot_ohh2_bausr_common.F90"
#include "common/parpot.F90"
!
character*(*) filnam
integer iunit, ir, iv
character*255 datfl
!     
!     WHEN HIBRIDON LOADS, A STRING CONTAINING ONLY ONE SPACE WILL BE
!     PASSED TO THIS SUBROUTINE.  MAKE SURE NOTHING WILL BE DONE.
if (filnam .eq. ' ') return
!
potnam = 'MA-DAGDIGIAN-KLOS-ALEXANDER OH--H2 MRCI PES'
!
!     LOAD DATA TO A COMMON BLOCK ONLY USED IN THIS FILE, AND DO PRE-
!     PROCESSING IF NECESSARY.
!     
!     Read B-Coefficients
call datfln(trim(filnam) // '_b', datfl)
open (unit=iunit, file=datfl, status='old')
read (iunit, *) nr
if (nr .gt. MAX_NR) call raise('too many R in data file.')
read (iunit, *) (rr(ir), ir=1, nr)
read (iunit, *) nvb
if (nvb .gt. MAX_NVB) &
     call raise('too many vlm terms in data file.')
do iv = 1, nvb
   read (iunit, *) lam1b(iv), lam2b(iv), lamb(iv)
   read (iunit, *) (bcoef(ir, iv), ir=1, nr)
end do
close (unit=iunit)
do iv = 1, nvb
   call dscal(nr, 1d0 / econv, bcoef(1, iv), 1)
end do
!     
!     Read F-Coefficients
call datfln(trim(filnam) // '_f', datfl)
open (unit=iunit, file=datfl, status='old')
read (iunit, *) nr
if (nr .gt. MAX_NR) call raise('too many R in data file.')
read (iunit, *) (rr(ir), ir=1, nr)
read (iunit, *) nvf
if (nvf .gt. MAX_NVF) &
     call raise('too many vlm terms in data file.')
do iv = 1, nvf
   read (iunit, *) lam1f(iv), lam2f(iv), lamf(iv)
   read (iunit, *) (fcoef(ir, iv), ir=1, nr)
end do
close (unit=iunit)
do iv = 1, nvf
   call dscal(nr, 1d0 / econv, fcoef(1, iv), 1)
end do
!
!     spline parameters are prepared here; pot routine only call the
!     accompanying function seval.
do iv = 1, nvb
   call spline(nr, rr, bcoef(1, iv), splb_b(1, iv), splb_c(1, iv), &
        splb_d(1, iv))
end do
do iv = 1, nvf
   call spline(nr, rr, fcoef(1, iv), splf_b(1, iv), splf_c(1, iv), &
        splf_d(1, iv))
end do
return
end
!     ------------------------------------------------------------------
!     A POT ROUTINE RETURNS VVL ARRAY (COVVL BLOCK) FOR A GIVEN R.
!     ALWAYS SET VV0 = 0 TO AVOID CONFUSION.  THE ISOTROPIC TERM CAN BE
!     EASILY TREAT AS ONE TERM IN THE POTENTIAL EXPANSION.  NOTE THAT
!     VVL SHOULD BE IN HARTREES.
subroutine pot(vv0, r_raw)
!
#include "pot_ohh2_bausr_common.F90"
real(8), intent(in) :: r_raw
real(8), intent(out) :: vv0
real(8) seval, r
integer :: iv
!
if (r_raw .lt. 3.5d0) then
   r = 3.5d0
else
   r = r_raw
end if
vv0 = 0d0
do iv = 1, nvb
   vvl(iv) = seval(nr, r, rr, bcoef(1, iv), &
        splb_b(1, iv), splb_c(1, iv), splb_d(1, iv))
end do
do iv = 1, nvf
   vvl(iv + nvb) = seval(nr, r, rr, fcoef(1, iv), &
        splf_b(1, iv), splf_c(1, iv), splf_d(1, iv))
end do
return
end
!     ------------------------------------------------------------------
!     THIS SUBROUTINE GOVERNS THE INPUT/OUTPUT OF THE BASIS ROUTINE.
!     ONLY IREAD IS USED: RETURN DIRECTLY IF ZERO.
subroutine syusr(irpot, readpt, iread)
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only : isrcod, junkr, rspar
use mod_hibasutil, only: raise
use funit, only: FUNIT_INP
!
#include "pot_ohh2_bausr_common.F90"
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
80 call raise('error read from input file.')
return
!     ------------------------------------------------------------------
entry ptrusr(fname, readpt)
return
!     ------------------------------------------------------------------
entry savusr(readpt)
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
!     ------------------------------------------------------------------
!     THE BASIS ROUTINE.  JHOLD, EHOLD, ISHOLD, NLEVEL, NLEVOP ARE TO
!     BE RETURNED FOR LEVEL INFORMATION; J, L, IS, EINT, CENT, N,
!     NTOP, AND POSSIBLY J12 ARE TO BE RETURNED FOR CHANNEL INFORMATION;
!     ALSO CALCULATE NON-ZERO COUPLING MATRIX ELEMENTS AND RETURN IN V2,
!     IV2, LAMNUM.
subroutine bausr(j, l, is, jhold, ehold, ishold, nlevel, nlevop, &
     sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, csflag, &
     clist, bastst, ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coj12, only: j12
use mod_conlam, only: nlam
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, junkr, rspar
use mod_hibasutil, only: raise
use mod_par, only: iprint
!
#include "pot_ohh2_bausr_common.F90"
integer, intent(out) :: j(:)
integer, intent(out) :: l(:)
integer, intent(out) :: is(:)
integer, intent(out), dimension(:) :: jhold
real(8), intent(out), dimension(:) :: ehold
integer, intent(out), dimension(:) :: ishold
integer, intent(out) :: nlevel
integer, intent(out) :: nlevop
real(8), intent(out), dimension(:), target :: sc1  ! k
real(8), intent(out), dimension(:), target :: sc2  ! ieps
real(8), intent(out), dimension(:), target :: sc3  ! jtemp
real(8), intent(out), dimension(:), target :: sc4  ! ktemp
real(8), intent(in) :: rcut
integer, intent(in) :: jtot
logical, intent(in) :: flaghf
logical, intent(in) :: flagsu
logical, intent(in) :: csflag
logical, intent(in) :: clist
logical, intent(in) :: bastst
logical, intent(in) :: ihomo
integer, intent(in) :: nu
integer, intent(in) :: numin
integer, intent(in) :: jlpar
integer, intent(out) :: n
integer, intent(in) :: nmax
integer, intent(out) :: ntop
type(ancou_type), intent(out), allocatable, target :: v2
type(ancouma_type), pointer :: ancouma
!     sc1--sc4 are scratches, whose type can be arbitary
!     In this pot routine, sc1 is the mixing angle between omega=3/2 and
!     omega=1/2 states for each level; sc2 is that mixing angle for each
!     channel; sc3 and sc4 are coefficients for omega=3/2 and omega=1/2
!     states, respectively, for each channel.
common /coered/ ered, rmu
double precision ered, rmu
!
integer nlist, ji1, eps1, fi1, ji2, li, ji1p, eps1p, fi1p, &
     ji2p, lip, jsave, isave, ipar, j12min, ji12, ji12p, &
     lpar, i, ilam, ivx, iv, icol, irow, inum, &
     i1, i2
double precision roteng, esave, vee, s1save, c12p, c32p, c12, c32
double precision v2pisg
double precision x, o11, o12, o22, tho
character(3) :: strfi
integer, pointer :: j1max, npar, j2min, j2max, iptsy2
real(8), pointer :: brot, aso, p, q, drot
j1max=>ispar(1); npar=>ispar(2); j2min=>ispar(3); j2max=>ispar(4); iptsy2=>ispar(5)
brot=>rspar(1); aso=>rspar(2); p=>rspar(3); q=>rspar(4); drot=>rspar(5)
!
if (.not. flaghf) call &
     raise('FLAGHF = .FALSE. FOR DOUBLET SYSTEM')
if (ihomo) call raise ('HOMONUCLEAR 2PI NOT IMPLEMENTED')
!     TWOMOL should be true, but this parameter is not passed here
!$$$      if (.not. twomol) call
!$$$     $     raise('TWOMOL = .FALSE. FOR MOL-MOL COLLISION')
if (flagsu) call raise('FLAGSU = .TRUE. FOR MOL-MOL COLLISION')
if (csflag) call raise('CS CALCULATION NOT IMPLEMENTED')
!     
!     Generate level list
nlist = 0
do ji1 = 0, j1max
   x = dble(ji1 + 1)
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
291       format (2i4, 4f10.3)
      do fi1 = 1, 2
         if (ji1 .eq. 0 .and. fi1 .eq. 1) cycle
         if (fi1 .eq. 1) then
            roteng = o11 * dcos(tho) ** 2 + o22 * dsin(tho) ** 2 &
                 + o12 * dsin(2 * tho)
         else
            roteng = o11 * dsin(tho) ** 2 + o22 * dcos(tho) ** 2 &
                 - o12 * dsin(2 * tho)
         end if
         do ji2 = j2min, j2max, iptsy2
            nlist = nlist + 1
            jhold(nlist) = 10 * ji1 + ji2
            ishold(nlist) = eps1 * fi1
            ehold(nlist) = (roteng + drot * ji2 * (ji2 + 1)) &
                 / econv
!     sc1 stores, for this basis, the mixing angle
            sc1(nlist) = tho
         end do
      end do
   end do
end do
!     
!     Sort the level list
do i1 = 1, nlist - 1
   esave = ehold(i1)
   do i2 = i1 + 1, nlist
      if (ehold(i2) .lt. esave) then
         esave = ehold(i2)
         ehold(i2) = ehold(i1)
         ehold(i1) = esave
         jsave = jhold(i2)
         jhold(i2) = jhold(i1)
         jhold(i1) = jsave
         isave = ishold(i2)
         ishold(i2) = ishold(i1)
         ishold(i1) = isave
         s1save = sc1(i2)
         sc1(i2) = sc1(i1)
         sc1(i1) = s1save
      end if
   end do
end do
!     
!     Set the lowest level as reference
esave = ehold(1)
do i1 = 1, nlist
   ehold(i1) = ehold(i1) - esave
end do
!     
if (bastst) call pr_lev_ohh2(nlist, jhold, ishold, ehold, sc1)
!     
!     Create a channel list from the level list
n = 0
do nlevel = 1, nlist
   ji1 = jhold(nlevel) / 10
   ji2 = mod(jhold(nlevel), 10)
   eps1 = isign(1, ishold(nlevel))
   ipar = eps1 * (-1) ** ji1
!     Since j1 is actually a half-integer, care must be
!     taken when calculate the bounds of j12 and l
   if (ji1 .ge. ji2) then
      j12min = ji1 - ji2
   else
      j12min = ji2 - ji1 - 1
   end if
   do ji12 = j12min, ji1 + ji2
      do li = iabs(jtot - ji12), jtot + ji12 + 1
         lpar = (-1) ** (li - jtot + ji2)
         if (ipar * lpar .ne. jlpar) cycle
         n = n + 1
         if (n .gt. nmax) call raise('too many channels.')
         j(n) = jhold(nlevel)
         is(n) = ishold(nlevel)
         j12(n) = ji12
         eint(n) = ehold(nlevel)
         l(n) = li
         cent(n) = li * (li + 1)
         sc2(n) = sc1(nlevel)
         if (iabs(ishold(nlevel)) .eq. 1) then
            sc3(n) = dcos(sc2(n))
            sc4(n) = dsin(sc2(n))
         else
            sc3(n) = dsin(sc2(n))
            sc4(n) = -dcos(sc2(n))
         end if
      end do
   end do
end do
nlevel = nlist
!     
if (bastst .and. iprint .ge. 1) then
   write (6, *)
   write (6, 275) jtot * 2 + 1, jlpar
275    format (' ** CHANNEL LIST FOR JTOT=', i3, '/2, JLPAR=', i2)
   do i = 1, n
      ji1 = j(i) / 10
      ji2 = mod(j(i), 10)
      select case (is(i))
         case (-2)
            strfi = 'F2f'
         case (-1)
            strfi = 'F1f'
         case (1)
            strfi = 'F1e'
         case (2)
            strfi = 'F2e'
      end select
      write (6, 280) i, 2 * ji1 + 1, strfi, ji2, 2 * j12(i) + 1, &
           l(i), sc3(i), sc4(i), eint(i) * econv
280       format (i4, 1x, 'J1=', i2, '/2 ', a, 2x, &
           'J2=', i2, 2x, 'J12=', i2, '/2  L=', i3, 2x, &
           'C3/2=', f7.4, 1x, 'C1/2=', f7.4, 2x, &
           'E=', f9.4)
   end do
   write (6, *)
end if
!     
!     Determine the number of asymtotically open levels
nlevop = 0
do i = 1, nlist
   if (ehold(i) .le. ered) nlevop = nlevop + 1
end do
!     
if (n .eq. 0) return
if (nu .eq. numin) then
   ntop = max(n, nlevop)
!  ntop is the maximum row dimension of all matrices passed in the
!  call list of subroutines propag and soutpt.
!  for fps make sure this is an odd number, for faster bank access.
!  this has no effect on vax or cray
   if (mod(ntop, 2) .eq. 0 .and. ntop .lt. nmax) &
        ntop = ntop + 1
else
   if (n .gt. ntop) call raise('nch greater than ntop.')
end if
!
!     Calculate coupling matrix elements
nlam = nvb + nvf
i = 0
ilam = 0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do ivx = 1, nvb + nvf
   ilam = ilam + 1
   ancouma => v2%get_angular_coupling_matrix(ilam)
   inum = 0
   do icol = 1, n
      do irow = icol, n
         ji1p = j(icol) / 10
         fi1p = iabs(is(icol))
         eps1p = isign(1, is(icol))
         ji2p = mod(j(icol), 10)
         ji12p = j12(icol)
         lip = l(icol)
         c32p = sc3(icol)
         c12p = sc4(icol)
!     
         ji1 = j(irow) / 10
         fi1 = iabs(is(irow))
         eps1 = isign(1, is(irow))
         ji2 = mod(j(irow), 10)
         ji12 = j12(irow)
         li = l(irow)
         c32 = sc3(irow)
         c12 = sc4(irow)
!     
         if (ivx .le. nvb) then
            iv = ivx
            vee = v2pisg(jtot, ji1p, eps1p, c12p, c32p, ji2p, &
                 ji12p, lip, ji1, eps1, c12, c32, ji2, ji12, &
                 li, lam1b(iv), lam2b(iv), lamb(iv), .true.)
         else
            iv = ivx - nvb
            vee = v2pisg(jtot, ji1p, eps1p, c12p, c32p, ji2p, &
                 ji12p, lip, ji1, eps1, c12, c32, ji2, ji12, &
                 li, lam1f(iv), lam2f(iv), lamf(iv), .false.)
         end if
         if (dabs(vee) .lt. machep) cycle
         i = i + 1
         inum = inum + 1
         call ancouma%set_element(irow=irow, icol=icol, vee=vee)
         if (bastst .and. iprint .ge. 2) then
            if (ivx .le. nvb) then
               write (6, 345) ilam, lam1b(iv), lam2b(iv), &
                    lamb(iv), icol, irow, i, vee
            else
               write (6, 345) ilam, lam1f(iv), lam2f(iv), &
                    lamf(iv), icol, irow, i, vee
            end if
         end if
345          format (i4, 3i3, 2i4, i5, f17.10)
      end do
   end do
   if (bastst .and. ivx .le. nvb) write (6, 347) &
        ilam, lam1b(iv), lam2b(iv), lamb(iv), ancouma%get_num_nonzero_elements()
   if (bastst .and. ivx .gt. nvb) write (6, 347) &
        ilam, lam1f(iv), lam2f(iv), lamf(iv), ancouma%get_num_nonzero_elements()
347    format ('ILAM=', i3, ' LAM1=', i3, ' LAM2=', i3, &
        ' LAM=', i3, ' LAMNUM(ILAM) = ', i6)
end do
if (bastst .and. iprint .ge. 2) then
   call v2%print(unit=6)
end if
return
end
!     ------------------------------------------------------------------
subroutine pr_lev_ohh2(n, js, iss, es, thos)
implicit none
!
integer n, js(*), iss(*)
double precision es(*), thos(*)
integer i, j1, j2, fi1, eps1
double precision ecm, econv, tho, c12, c32
parameter (econv=219474.6315343234)
write (6, 125)
125 format (/, 10x, &
     'SORTED LEVEL LIST', /, '     N    J1  F1/2  EPS1    J2', &
     '  EINT(CM-1)       C-1/2       C-3/2')
do i = 1, n
   j2 = mod(js(i), 10)
   j1 = js(i) / 10
   fi1 = iabs(iss(i))
   eps1 = isign(1, iss(i))
   ecm = es(i) * econv
   tho = thos(i)
   if (fi1 .eq. 1) then
      c32 = dcos(tho)
      c12 = dsin(tho)
   else
      c32 = dsin(tho)
      c12 = -dcos(tho)
   end if
   write (6, 126) i, dble(j1) + 0.5d0, fi1, eps1, j2, &
        ecm, c12, c32
126    format (i6, f6.1, 3i6, 3f12.3)
end do
return
end

!     ------------------------------------------------------------------
double precision function v2pisg(jtot, j1p, eps1p, c12p, c32p, &
     j2p, j12p, lp, j1, eps1, c12, c32, j2, j12, l, &
     lam1, lam2, lam, isdiag)
implicit none
!     
!     The subroutine calculate the coupling matrix elements, shown in
!     Eq. (27) in the notes of Q. Ma.
!     
!     If (omeg1p .eq. omeg1), the coefficient before B is calculated;
!     otherwise the coefficient before F is calculated.
!     
integer :: jtot, j1p, eps1p, j2p, j12p, lp, j1, eps1, &
     j2, j12, l, lam1, lam2, lam
real(8) :: c12p, c32p, c12, c32
logical :: isdiag
integer :: iphase
real(8) :: phase, pref, threej, sixj, ninej
real(8) :: xj1p, xj2p, xj12p, xlp, &
     xj1, xj2, xj12, xl, xjtot, xlam1, xlam2, xlam
real(8) :: xf3jm0, xf3j, xf6j, xf9j
real(8), parameter :: machep=epsilon(0d0)
!     
!     One is added since j1p and j1 should be half integers
iphase = eps1p * eps1 * (-1) ** (j1p + j1 + lam1 + 1)
if (iphase .eq. 1) then
   v2pisg = 0d0
   return
end if
!     
xj1p = dble(j1p) + 0.5d0
xj2p = dble(j2p)
xj12p = dble(j12p) + 0.5d0
xlp = dble(lp)
xj1 = dble(j1) + 0.5d0
xj2 = dble(j2)
xj12 = dble(j12) + 0.5d0
xl = dble(l)
xjtot = dble(jtot) + 0.5d0
xlam1 = dble(lam1)
xlam2 = dble(lam2)
xlam = dble(lam)
!     
!     
v2pisg = 0d0
!
threej = xf3jm0(xj2p, xlam2, xj2) * xf3jm0(xlp, xlam, xl)
if (dabs(threej) .lt. machep) return
!     omega-dependent part
if (isdiag) then
   threej = threej * &
        (c12p * c12 * xf3j(xj1p, xlam1, xj1, -0.5d0, 0d0, 0.5d0) &
        - c32p * c32 * xf3j(xj1p, xlam1, xj1, -1.5d0, 0d0, 1.5d0))
else
   threej = threej * dble(eps1) * &
        (c12p * c32 * xf3j(xj1p, xlam1, xj1, -0.5d0, 2d0, -1.5d0) &
        - c32p * c12 * xf3j(xj1p, xlam1, xj1, &
        -1.5d0, 2d0, -0.5d0))
end if
if (dabs(threej) .lt. machep) return
!     
sixj = xf6j(xj12, xl, xjtot, xlp, xj12p, xlam)
if (dabs(sixj) .lt. machep) return
ninej = xf9j(xj1, xj2, xj12, xj1p, xj2p, xj12p, &
     xlam1, xlam2, xlam)
if (dabs(ninej) .lt. machep) return
!     
!     Again 1 is added to compensate the dropped half-integer part
iphase = jtot + lam1 - lam2 + j1 - j2 + j12p - l - lp  + 1
if (mod(iphase, 2) .eq. 0) then
   phase = 1d0
else
   phase = -1d0
end if
!     
pref = (2d0 * xj1p + 1d0) * (2d0 * xj2p + 1d0) &
     * (2d0 * xj12p + 1d0) * (2d0 * xlp + 1d0) &
     * (2d0 * xj1 + 1d0) * (2d0 * xj2 + 1d0) * (2d0 * xj12 + 1d0) &
     * (2d0 * xl + 1d0) * (2d0 * xlam + 1d0)
pref = dsqrt(pref)
!     
v2pisg = phase * pref * threej * sixj * ninej
return
end

subroutine datfln(filenm, fullnm)
character (len=*) :: filenm, fullnm
fullnm = 'potdata/' // trim(filenm)
return
end subroutine datfln
