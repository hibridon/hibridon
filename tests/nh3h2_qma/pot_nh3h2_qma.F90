!     pot_nh3h2_qma.f
!     authors: Qianli Ma
!
!     refit of Grenoble NH3-H2 ab initio points
!     (see Maret at el., MNRAS 399, 425 (2009)
!     and data provided by Claire Rist - Nov-2011)
!
!     Pot routine for a collision between a symmetric top and a linear
!     molecule.
!
!     PES to be read from data file specified in the input file.
!
!     This basis routine should have ibasty = 100, so hibridon will
!     treat it as a mol-mol problem with j12.
!
!     DUMMY SUBROUTINES FOR ALL USER-DEFINED BASIS/POT.
#include "assert.h"
#include "unused.h"
#include "common/ground.F90"
!
!     ------------------------------------------------------------------
!     THE FOLLOWING SOUBROUTINE WILL BE THE MAIN FUNCTION FOR MAKEPOT.
!     NOTE THAT VVL IS IN HARTREES.

module mod_nh3h2_qma
   character*40 :: potfil

contains

!     ------------------------------------------------------------------
subroutine pr_lev_nh3h2(n, js, ks, iepss, es)
implicit none
!
integer n, js(*), ks(*), iepss(*)
double precision es(*)
integer i, j1, j2, isym
double precision ecm, econv
parameter (econv=219474.6315343234)
write (6, 125)
125 format (/, 10x, &
     'SORTED LEVEL LIST', /, '   N   J   K  EPS INV J2   ', &
     'EINT(CM-1)')
do i = 1, n
   j2 = mod(js(i), 10)
   j1 = js(i) / 10
   isym = -iepss(i) * (-1) ** j1
   ecm = es(i) * econv
   write (6, 126) i, j1, ks(i), iepss(i), isym, j2, ecm
126    format (6i4, f10.3)
end do
return
end

end module mod_nh3h2_qma

subroutine driver
use mod_covvl, only: vvl
use mod_hibasutil, only: raise
use constants, only: econv
use mod_hipot, only: loapot, pot
implicit none
!
integer MAX_NR, MAX_NV
parameter (MAX_NR=300, MAX_NV=300)

common /stpln1/ nr, nv, rr, v_pot
common /stpln2/ lb1, mu1, lb2, mu2
common /stplnb/ spl_b
common /stplnc/ spl_c
common /stplnd/ spl_d
integer nr, nv
integer lb1(MAX_NV), mu1(MAX_NV), lb2(MAX_NV), mu2(MAX_NV)
double precision rr(MAX_NR), v_pot(MAX_NR, MAX_NV), &
     spl_b(MAX_NR, MAX_NV), spl_c(MAX_NR, MAX_NV), &
     spl_d(MAX_NR, MAX_NV)

character*40 filenm
double precision r, vv0
integer i
!
print *, 'Please input the name of potdata file:'
read (5, *, end=99) filenm
call loapot(10, filenm)
10 print *, 'R (bohr), Ctrl+D to exit:'
read (5, *, end=99) r
call pot(vv0, r)
write (6, 20) (vvl(i) * econv, i=1, nv)
20 format (10(1pe16.8))
goto 10
99 end
!     ------------------------------------------------------------------
!     DATA FILE, IF REQUIRED, CAN BE LOADED WITH THIS SOUBROUTINE.
subroutine loapot(iunit, filnam)
use mod_hibasutil, only: raise
use constants, only: econv
use mod_hiblas, only: dscal
use mod_hipotutil, only: spline, datfln
implicit none
integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)    
!
integer MAX_NR, MAX_NV
parameter (MAX_NR=300, MAX_NV=300)

common /stpln1/ nr, nv, rr, v_pot
common /stpln2/ lb1, mu1, lb2, mu2
common /stplnb/ spl_b
common /stplnc/ spl_c
common /stplnd/ spl_d
integer nr, nv
integer lb1(MAX_NV), mu1(MAX_NV), lb2(MAX_NV), mu2(MAX_NV)
double precision rr(MAX_NR), v_pot(MAX_NR, MAX_NV), &
     spl_b(MAX_NR, MAX_NV), spl_c(MAX_NR, MAX_NV), &
     spl_d(MAX_NR, MAX_NV)
!
integer ir, iv
character*255 datfl
!
!     WHEN HIBRIDON LOADS, A STRING CONTAINING ONLY ONE SPACE WILL BE
!     PASSED TO THIS SUBROUTINE.  MAKE SURE NOTHING WILL BE DONE.
if (filnam .eq. ' ') return
!
!     LOAD DATA TO A COMMON BLOCK ONLY USED IN THIS FILE, AND DO PRE-
!     PROCESSING IF NECESSARY.
call datfln(filnam, datfl)
open (unit=iunit, file=datfl, status='old')
read (iunit, *) nr
if (nr .gt. MAX_NR) call raise('too many R in data file.')
read (iunit, *) (rr(ir), ir=1, nr)
read (iunit, *) nv
if (nv .gt. MAX_NV) call raise('too many vlm terms in data file.')
do iv = 1, nv
   read (iunit, *) lb1(iv), mu1(iv), lb2(iv), mu2(iv)
   read (iunit, *) (v_pot(ir, iv), ir=1, nr)
end do
close (unit=iunit)
do iv = 1, nv
   call dscal(nr, 1d0 / econv, v_pot(1, iv), 1)
end do
!
!     spline parameters are prepared here; pot routine only call the
!     accompanying function seval.
do iv = 1, nv
   call spline(nr, rr, v_pot(1, iv), spl_b(1, iv), spl_c(1, iv), &
        spl_d(1, iv))
end do
return
end
!     ------------------------------------------------------------------
!     A POT ROUTINE RETURNS VVL ARRAY (COVVL BLOCK) FOR A GIVEN R.
!     ALWAYS SET VV0 = 0 TO AVOID CONFUSION.  THE ISOTROPIC TERM CAN BE
!     EASILY TREAT AS ONE TERM IN THE POTENTIAL EXPANSION.  NOTE THAT
!     VVL SHOULD BE IN HARTREES.
subroutine pot(vv0, r)
use mod_covvl, only: vvl
use mod_hipotutil, only: seval
implicit none
real(8), intent(out) :: vv0
real(8), intent(in) :: r  ! intermolecular distance
!
integer MAX_NR, MAX_NV
parameter (MAX_NR=300, MAX_NV=300)

common /stpln1/ nr, nv, rr, v_pot
common /stpln2/ lb1, mu1, lb2, mu2
common /stplnb/ spl_b
common /stplnc/ spl_c
common /stplnd/ spl_d
integer nr, nv
integer lb1(MAX_NV), mu1(MAX_NV), lb2(MAX_NV), mu2(MAX_NV)
double precision rr(MAX_NR), v_pot(MAX_NR, MAX_NV), &
     spl_b(MAX_NR, MAX_NV), spl_c(MAX_NR, MAX_NV), &
     spl_d(MAX_NR, MAX_NV)

integer iv
!
vv0 = 0d0
do iv = 1, nv
   vvl(iv) = seval(nr, r, rr, v_pot(1, iv), &
        spl_b(1, iv), spl_c(1, iv), spl_d(1, iv))
end do
return
end
!     ------------------------------------------------------------------
!     THIS SUBROUTINE GOVERNS THE INPUT/OUTPUT OF THE BASIS ROUTINE.
!     ONLY IREAD IS USED: RETURN DIRECTLY IF ZERO.
subroutine syusr(irpot, readpt, iread)
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysr, only: isrcod, rspar
use mod_hibasutil, only: raise
use mod_nh3h2_qma, only: potfil
use mod_hipot, only: loapot
implicit none
!
integer MAX_NR, MAX_NV
parameter (MAX_NR=300, MAX_NV=300)

common /stpln1/ nr, nv, rr, v_pot
common /stpln2/ lb1, mu1, lb2, mu2
common /stplnb/ spl_b
common /stplnc/ spl_c
common /stplnd/ spl_d
integer nr, nv
integer lb1(MAX_NV), mu1(MAX_NV), lb2(MAX_NV), mu2(MAX_NV)
double precision rr(MAX_NR), v_pot(MAX_NR, MAX_NV), &
     spl_b(MAX_NR, MAX_NV), spl_c(MAX_NR, MAX_NV), &
     spl_d(MAX_NR, MAX_NV)

integer, intent(inout) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
!     NUMBER OF BASIS-SPECIFIC VARIABLES, MODIFY ACCORDINGLY.
integer icod, ircod
parameter (icod=5, ircod=5)


integer, pointer :: nterm, ipotsy, iop, ninv, jmax, ipotsy2, j2max, j2min
real(8), pointer :: brot, crot, delta, emax, drot
nterm=>ispar(1); ipotsy=>ispar(2); iop=>ispar(3); ninv=>ispar(4); jmax=>ispar(5); ipotsy2=>ispar(6)
j2max=>ispar(7); j2min=>ispar(8)

brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); emax=>rspar(4); drot=>rspar(5)
UNUSED_DUMMY(irpot)
UNUSED_DUMMY(readpt)
!     DEFINE THE NAMES HERE
scod(1)='NTERM'
scod(2)='IPOTSY'
scod(3)='IOP'
scod(4)='NINV'
scod(5)='JMAX'
scod(6)='IPOTSY2'
scod(7)='J2MAX'
scod(8)='J2MIN'
scod(9)='BROT'
scod(10)='CROT'
scod(11)='DELTA'
scod(12)='EMAX'
scod(13)='DROT'
nscode = icod+ircod+3
isicod = icod
isrcod = ircod
!     KEEP THE FOLLOWING LINE
if (iread .eq. 0) return
!     READ THE LAST FEW LINES OF THE INPUT FILE
read (8, *, err=80) ipotsy, iop, ninv, ipotsy2
read (8, *, err=80) jmax
read (8, *, err=80) j2min, j2max
read (8, *, err=80) brot, crot, delta, emax
read (8, *, err=80) drot
read (8, *, err=80) potfil
call loapot(10, potfil)
close (8)
return
80 call raise('error read from input file.')
return
end subroutine
!     ------------------------------------------------------------------
subroutine ptrusr(fname, readpt)
implicit none
character*(*), intent(inout) :: fname
logical, intent(inout) :: readpt
UNUSED_DUMMY(fname)
UNUSED_DUMMY(readpt)
return
end subroutine
!     ------------------------------------------------------------------
subroutine savusr()
use mod_cosysi, only: ispar
use mod_cosysr, only: rspar
use funit, only: FUNIT_INP
use mod_nh3h2_qma, only: potfil
implicit none
integer, pointer :: nterm, ipotsy, iop, ninv, jmax, ipotsy2, j2max, j2min
real(8), pointer :: brot, crot, delta, emax, drot
nterm=>ispar(1); ipotsy=>ispar(2); iop=>ispar(3); ninv=>ispar(4); jmax=>ispar(5); ipotsy2=>ispar(6)
j2max=>ispar(7); j2min=>ispar(8)

brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); emax=>rspar(4); drot=>rspar(5)
!     WRITE THE LAST FEW LINES OF THE INPUT FILE.
write (FUNIT_INP, 220) ipotsy, iop, ninv, ipotsy2
220 format (4i4, 14x,'   ipotsy, iop, ninv, ipotsy2')
write (FUNIT_INP, 230) jmax
230 format (i4, 26x, '   jmax')
write (FUNIT_INP,231) j2min, j2max
231 format (2i4, 22x,'   j2min, j2max')
write (FUNIT_INP, 250) brot, crot, delta, emax
250 format (3f8.4, f8.2, ' brot, crot, delta, emax' )
write (FUNIT_INP, 251) drot
251 format (f12.6, 18x,'   drot')
write (FUNIT_INP, *) potfil
return
end
!     ------------------------------------------------------------------
!     THE BASIS ROUTINE.  JHOLD, EHOLD, ISHOLD, NLEVEL, NLEVOP ARE TO
!     BE RETURNED FOR LEVEL INFORMATION; J, L, IS, IEPS, EINT, CENT, N,
!     NTOP, AND POSSIBLY J12 ARE TO BE RETURNED FOR CHANNEL INFORMATION;
!     ALSO CALCULATE NON-ZERO COUPLING MATRIX ELEMENTS AND RETURN IN V2,
!     IV2, LAMNUM.
subroutine bausr(bqs, jhold, ehold, ishold, nlevel, nlevop, &
     sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, csflag, &
     clist, bastst, ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
use mod_assert, only: fassert
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coamat, only: ietmp ! ietmp(1)
use mod_conlam, only: nlam
use mod_cosysi, only: ispar
use mod_cosysr, only: rspar
use mod_hibasutil, only: vlmstpln, raise
use constants, only: econv
use, intrinsic :: ISO_C_BINDING   ! for C_LOC and C_F_POINTER
use mod_par, only: iprint
use mod_ered, only: ered
use mod_hitypes, only: bqs_type
use mod_nh3h2_qma, only: pr_lev_nh3h2
implicit none
type(bqs_type), intent(out) :: bqs
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
integer, dimension(:), pointer :: k
integer, dimension(:), pointer :: ieps
integer, dimension(:), pointer :: jtemp
integer, dimension(:), pointer :: ktemp
type(ancouma_type), pointer :: ancouma
integer MAX_NR, MAX_NV
parameter (MAX_NR=300, MAX_NV=300)
common /stpln1/ nr, nv, rr, v_pot
common /stpln2/ lb1, mu1, lb2, mu2
common /stplnb/ spl_b
common /stplnc/ spl_c
common /stplnd/ spl_d
integer nr, nv
integer lb1(MAX_NV), mu1(MAX_NV), lb2(MAX_NV), mu2(MAX_NV)
double precision rr(MAX_NR), v_pot(MAX_NR, MAX_NV), &
     spl_b(MAX_NR, MAX_NV), spl_c(MAX_NR, MAX_NV), &
     spl_d(MAX_NR, MAX_NV)
!
integer nlist, ki, ji, numeps, iep, iepsil, isym, ji2, i1, i2, &
     jsave, ksave, isave, ipar, ji12, &
     li, lpar, i, ilam, iv, icol, &
     irow, jc, jr, j2c, j2r, kc, kr, j12c, j12r, lc, lr, inum
double precision roteng, esave, vee
!
integer, pointer :: nterm, ipotsy, iop, ninv, jmax, ipotsy2, j2max, j2min
real(8), pointer :: brot, crot, delta, emax, drot
UNUSED_DUMMY(rcut)
UNUSED_DUMMY(clist)
UNUSED_DUMMY(ihomo)
nterm=>ispar(1); ipotsy=>ispar(2); iop=>ispar(3); ninv=>ispar(4); jmax=>ispar(5); ipotsy2=>ispar(6)
j2max=>ispar(7); j2min=>ispar(8)
brot=>rspar(1); crot=>rspar(2); delta=>rspar(3); emax=>rspar(4); drot=>rspar(5)

call C_F_POINTER (C_LOC(sc1), k, [nmax])
call C_F_POINTER (C_LOC(sc2), ieps, [nmax])
call C_F_POINTER (C_LOC(sc3), jtemp, [nmax])
call C_F_POINTER (C_LOC(sc4), ktemp, [nmax])

if (flaghf) call raise('FLAGHF = .TRUE. FOR SINGLET SYSTEM')
if (flagsu) call raise('FLAGSU = .TRUE. FOR MOL-MOL COLLISION')
if (csflag) call raise('CS calculation not implemented')
nlist = 0
do ki = 0, jmax
   do ji = ki, jmax
      numeps = 2
      if (ki .eq. 0) numeps = 1
      do iep = 1, numeps
         if (iop .eq. -1 .and. mod(ki, ipotsy) .ne. 0) cycle
         if (iop .eq. 1 .and. mod(ki, ipotsy) .eq. 0) cycle
         roteng = brot * ji * (ji + 1) + (crot - brot) * ki ** 2
         if (roteng .gt. emax) cycle
         iepsil = 1 - 2 * (iep - 1)
         isym = -iepsil * (-1) ** ji
         if (isym .eq. -1) roteng = roteng + delta
         if ((ninv .ne. 2) .and. (ninv .ne. isym)) cycle
         do ji2 = j2min, j2max, ipotsy2
            nlist = nlist + 1
            jtemp(nlist) = 10 * ji + ji2
            ktemp(nlist) = ki
            ietmp(nlist) = iepsil
            ehold(nlist) = (roteng + drot * ji2 * (ji2 + 1)) &
                 / econv
         end do
      end do
   end do
end do
do i1 = 1, nlist - 1
   esave = ehold(i1)
   do i2 = i1 + 1, nlist
      if (ehold(i2) .lt. esave) then
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
   end do
end do
do i = 1, nlist
   jhold(i) = jtemp(i)
   ishold(i) = -ietmp(i) * (-1) ** (jhold(i) / 10) * ktemp(i)
end do
if (bastst) call pr_lev_nh3h2(nlist, jtemp, ktemp, ietmp, ehold)
!
!     now set up channel and level list for scattering calculation
call bqs%init(nmax)
n = 0
do nlevel = 1, nlist
   ki = ktemp(nlevel)
   ji = jtemp(nlevel) / 10
   ji2 = mod(jtemp(nlevel), 10)
   iepsil = ietmp(nlevel)
   ipar = -iepsil * (-1) ** (ji + ki)
   do ji12 = iabs(ji - ji2), ji + ji2
      do li = iabs(jtot - ji12), jtot + ji12
         lpar = (-1) ** (li - jtot + ji2)
         if (ipar * lpar .ne. jlpar) cycle
         n = n + 1
         if (n .gt. nmax) call raise('too many channels.')
         bqs%jq(n) = jtemp(nlevel)
         bqs%inq(n) = -iepsil * (-1) ** ji * ki
         ieps(n) = iepsil
         bqs%j12(n) = ji12
         eint(n) = ehold(nlevel)
         bqs%lq(n) = li
         bqs%length = n
         cent(n) = li * (li + 1)
         if (bastst .and. iprint .ge. 2) &
              write (6, 280) n, bqs%jq(n), bqs%inq(n), ieps(n), bqs%j12(n), &
              bqs%lq(n), eint(n) * econv
280          format (6i6, f12.3)
      end do
   end do
end do
nlevel = nlist
!  also determine number of levels which are open
nlevop = 0
do i = 1, nlist
   if (ehold(i) .le. ered) nlevop = nlevop + 1
end do
!  return if no channels
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
nlam = nv
i = 0
ilam = 0
ASSERT(.not. allocated(v2))
v2 = ancou_type(nlam=nlam, num_channels=ntop)
do iv = 1, nv
   ilam = ilam + 1
   ancouma => v2%get_angular_coupling_matrix(ilam)
   inum = 0
   do icol = 1, n
      do irow = icol, n
         jc = bqs%jq(icol) / 10
         jr = bqs%jq(irow) / 10
         j2c = mod(bqs%jq(icol), 10)
         j2r = mod(bqs%jq(irow), 10)
         kc = iabs(bqs%inq(icol))
         kr = iabs(bqs%inq(irow))
         j12c = bqs%j12(icol)
         j12r = bqs%j12(irow)
         lc = bqs%lq(icol)
         lr = bqs%lq(irow)
         call vlmstpln(jr, lr, jc, lc, j2r, j2c, j12r, j12c, &
              jtot, kr, kc, lb1(iv), mu1(iv), lb2(iv), mu2(iv), &
              ieps(irow), ieps(icol), vee, csflag)
         if (vee .eq. 0d0) cycle
         i = i + 1
         inum = inum + 1
         call ancouma%set_element(irow=irow, icol=icol, vee=vee)
         if (bastst .and. iprint .ge. 2) &
              write (6, 345) ilam, lb1(iv), icol, irow, i, vee
345          format (i4, 2i7, 2i6, f17.10)
      end do
   end do
   if (bastst) write (6, 347) ilam, lb1(iv), mu1(iv), lb2(iv), &
         mu2(iv), ancouma%get_num_nonzero_elements()
347    format ('ILAM=', i3, ' LAM=', i3, ' MU=', i3, &
        ' LAM2=', i3, ' MU2=', i3, ' LAMNUM(ILAM) = ', i6)
end do
if (bastst .and. iprint .ge. 2) then
   call v2%print(unit=6)
end if
return
end
