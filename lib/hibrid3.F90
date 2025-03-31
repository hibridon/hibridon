#include "assert.h"
#include "unused.h"
module mod_hibrid3
  use mod_assert, only: fassert
contains
!************************************************************************
!                         hibridon 3  library                           *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!   4. potmin     function, determines minimum of potential             *
!   4a. testpt     to print out potential as a function of r or theta   *
!   4b. testptn   testpot for basistyp 1 and 4                          *
!   4c. testpt20  testpot for basistyp 20, added by L. Ma oct. 2012     *
!   5. propag     controls propagation                                  *
!   6. rles       solves a * x = c                                      *
!      logdb      log derivativ propagator                             *
!   7. runlog     log-derivative propagator from r to r = rend          *
!   8. scairy     scaled airy functions and derivatives                 *
!   9. smatop     subroutine  to compute s-matrix                       *
!  10. smatrx           see 9.                                          *
!  11. expand     expands smatrix from open to full basis               *
!     the subroutines below have been moved to hibound.f                *
!  (12. bound      susan gregurick's bound state program)               *
!  (13. bound_wavfn)                                                    *
!  (14. h_basis)                                                        *
!                                                                       *
!************************************************************************
! -----------------------------------------------------------------------

! -------------------------------------------------------------------------
subroutine druckq(w,n,m,string)
implicit double precision (a-h,o-z)
character*(*) string
dimension w(n,n)
print*,' '
print*,string
print*,' '
do 10 i=1,m
10 write(6,20) (w(i,j),j=1,m)
20 format(1x,10f12.6)
return
end
! -------------------------------------------------------------------------
function potmin()
!  current revision date: 25-sept-87
implicit double precision (a-h,o-z)
r = 4.0d0
dr = 0.5d0
call pot(vv0,r)
10 r = r+dr
elast = vv0
call pot(vv0,r)
if(elast-vv0) 40,50,10
40 dr = -dr*0.5d0
if(abs(dr).gt.0.01d0) goto 10
50 potmin = r
return
end
! -------------------------------------------------------------------------
subroutine propag (z, w, zmat, &
                   bqs, &
                   ien, nerg, en, eshift, rstart, rendld, spac, &
                   tolhi, rendai, rincr, fstfac, tb, tbm, &
                   ipos, prlogd, noprin, airyfl, prairy, &
                   nch, nopen, nmax, v2, isteps, nsteps)
! ------------------------------------------------------------------------
!  subroutine to:
!    1.  propagate the log-derivative matrix from rstart to rendld
!        using the log-derivative propagator of manolopoulos
!    2.  propagate the log-derivative matrix from rendld to rendai
!        using the airy propagator of alexander and manolopoulos
!    3.  obtain the s-matrix and t-matrix squared at r=rendai
!  author:  millard alexander
!  revision date: 28-may-1993
!  revision date: 13-nov-1996 by mby adding bound calculation
!  latest revision (algorithm): 22-apr-1997 by mha
!  modified call to bound to use q.ma's revised version of this
!    subroutine:  28-jun-2013, p.dagdigian
! ----------------------------------------------------------------------------
!  definition of variables in call list:
!    z:           matrix of maximum dimension nmax x nmax
!                 on return:  the upper left nopen x nopen block of z
!                          contains the modulus squared of the t-matrix
!    w:           on return:  the upper-left nopen x nopen block of sr
!                          contains the real part of the s-matrix
!    zmat:        on return: the upper-left nopen x nopen block of si
!                         contains the imaginary part of the s-matrix
!    bqs:         on entry: contain the rotational angular momenta, orbital
!                 angular momenta, and other quantum index for each channel
!                 on return: the first nopen elements contain the rotational
!                 angular momenta, and other quantum index for each open chann
!    nerg:        the number of total energies at which calculation is to be
!                 performed
!    energ:       the array of total energies at which calculation is to be
!                 performed
!    ien:         the current ordinal total energy, i.e. if ien = 2,
!                 current calculation corresponds to the 2nd total energy
!    en:          current total energy in hartree
!    eshift:      2 m [en - energ(1)] / h-bar**2 in atomic units
!    rstart:      starting point (bohr) for logd integration
!    rendld:      ending point for logd integration
!    spac:        step size for logd integration
!    tolhi:       error parameter to determine step sizes in airy integration
!    rendai:      outer limit for airy integration
!    rincr:       point at which increase in airy steps can occur
!    fstfac:       factor by which logd step size is multiplied to give
!                 initial airy step size
!                 e.g.  drnow-first = fstfac * spac
!    tb:          cpu time (sec) used in determination of channel basis
!    tbm:         wall clock time (sec) used in determination of channel basis
!  logical variables:
!    ipos         if .true., then 132 column printer
!                 if .false., then 80 column printer
!    prlogd        if .true., then lower triangle of the log-derivative matrix
!                 is printed out at end of logd and airy integration
!    noprin       if .true., then all printing is suppressed
!    iprint:      if .true., then print out of step-by-step information
!    airyfl:      if .true., then airy propagation will occur
!                 if .false., then no airy propagation will occur, the
!                 integration will stop at r=renld
!    prairy:      if .true., then interval size, interval midpoint, and maximu
!                 estimated diagonal and off-diagonal corrections are printed
!                 out in airy propagation
!    nopen        on return:  number of energetically open channels
!    nch          number of coupled equations
!    nmax         leading dimension of matrices z, w, zmat, as well a
!                 all vectors

!  variables in common block /copmat/
!    rtmn,rtmx: minimum and maximum turning points
!    iflag:     variable used in determination of turning points (not used her
!           iflag = 0 if all channels are in classically forbidden region
!           iflag = 1 if some channels are open
!           iflag = 2 if all asymptotically open channels are open at r
!  ------------------------------------------------------------------
use mod_ancou, only: ancou_type
use mod_hibrid2, only: mxoutd
use funit, only: FUNIT_TRANS_MAT, FUNIT_QUAD_MAT
use mod_phot, only: boundf
use mod_pmat, only: rtmn, rtmx, iflag
use mod_cputim, only: cpuld, cpuai, cpupot, cpusmt, cpupht
use mod_hiutil, only: mtime, gettim
use mod_hibrid1, only: airprp
use mod_hitypes, only: bqs_type
use mod_hiblas, only: dcopy
implicit none
!   square matrices
real(8), intent(out) :: z(nmax, nch)
real(8), intent(out) :: w(nmax, nch)
real(8), intent(out) :: zmat(nmax, nch)
type(bqs_type), intent(inout) :: bqs
integer, intent(in) :: ien
integer, intent(in) :: nerg
real(8), intent(in) :: en
real(8), intent(in) :: eshift
real(8), intent(in) :: rstart
real(8), intent(in) :: rendld  ! end radius for log derivative method
real(8), intent(in) :: spac
real(8), intent(in) :: tolhi
real(8), intent(in) :: rendai
real(8), intent(in) :: rincr
real(8), intent(in) :: fstfac
real(8), intent(in) :: tb
real(8), intent(out) :: tbm
logical, intent(in) :: ipos
logical, intent(in) :: prlogd
logical, intent(in) :: noprin
logical, intent(in) :: airyfl
logical, intent(in) :: prairy
integer, intent(in) :: nch
integer, intent(out) :: nopen
integer, intent(in) :: nmax
type(ancou_type), intent(in) :: v2
integer, intent(inout) :: isteps
integer, intent(inout) :: nsteps

logical :: twoen
logical ::  first
character*10 :: tbs,tps,tds,tws,tairys, twfs
character*10 :: time

real(8) :: r
real(8) :: t1, t2
real(8) :: t11, t22
real(8) :: tb1, tb2, tbx
real(8) :: ttx, tty
real(8) :: drnow
real(8) :: cpup
integer :: icol
integer :: itwo
integer :: nrow
real(8) :: tairy
real(8) :: td, tdm, tp, tpm, tw, twm, twf, twfm

!  vectors
!  prec is precision of single precision floating pt number
real(8), parameter :: prec = 1.d+11
UNUSED_DUMMY(tbm)
tbm = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.

call mtime(t11,t22)
first = .true.
r = rstart
iflag = 0
if (nerg .gt. 1) then
  twoen = .true.
  itwo = ien - 1
  rewind FUNIT_TRANS_MAT
  rewind FUNIT_QUAD_MAT
else
  twoen = .false.
  itwo = -1
end if
!  return with message if the number of channels equals zero
if (nch .le. 0) then
  nopen = 0
  write (6, 80)
  write (9, 80)
80   format (' *** NCH.LE.0 IN PROPAG, RETURN WITHOUT PROPAGATION')
  return
end if

if (boundf) then
  call mtime(tb1,tbx)
! revise call to bound (28-jun-2013, p.dagdigian)
!        call bound(z,w,zmat,amat,bmat,nch,nmax)
  call bound(nch, nmax, v2)
! return to flow after bound calculation
  call mtime(tb2,tbx)
  call gettim(tb2-tb1,time)
  write (6, 105) time
105   format (' *** CPU TIME FOR BOUND = ',a)
  return
end if

!****************************************************************

!  now integrate coupled equations from r=rstart to r=rendld
!  using manolopoulos log-derivative integrator
cpup=cpupot

call mtime(ttx,tty)
call runlog (z, &
             r, rendld, spac, eshift, itwo, twoen, &
             td, tdm, tp, tpm, twf, twfm, prlogd, noprin, &
             ipos, nch, nmax, v2, isteps, nsteps)

!  on return from runlog, z contains the log-derivative matrix at r = rendld
!  branch to airy integration if desired. integrate coupled equations
!  from r=rendld to r=rendai using alexander-manolopoulos linear reference
!  potential integrator
call mtime(t1,t2)
cpuld=cpuld+t1-t11-cpupot+cpup-twf
cpup=cpupot
cpupht=cpupht+twf
if (airyfl) then
  drnow = spac * fstfac
!  symmetrize (fill in upper triangle) of logd matrix
  if (nch .gt. 1) then
    do 110 icol = 1, nch -1
!  nrow is the number of subdiagonal elements in column icol
    nrow = nch - icol
          call dcopy (nrow,z(icol+1,icol),1,z(icol,icol+1),nmax)
110     continue
  endif

  call airprp (z, &
               r, rendai, drnow, en, &
               tolhi, rincr, eshift, nch, nmax, itwo, prairy, &
               twoen,noprin, v2)
!  on return from airprp, z contains the log-derivative matrix at r = rendai
end if
call mtime(t11,t22)
tairy = t11 - t1
! note that we don't need to subtract delta-time for calculation
! of ground state wavefunction here, since it is done only in logd ste
tp=tp+cpupot-cpup
tairy=tairy-cpupot+cpup
cpuai=cpuai+tairy
if (prlogd .and. airyfl) then
  write (9, 260) r
260   format(/' ** LOG-DERIVATIVE MATRIX AFTER AIRPRP; R = ', 1pe15.8)
  call mxoutd (9, z, nch, nmax, 0, ipos)
end if
!  now calculate s-matrix and t-matrix squared
call smatrx (z, w, zmat, &
             bqs, r, prec, tw, twm, nopen, nch, nmax, &
             prlogd,ipos)

! convert to time string
  call mtime(t1,t2)
  cpusmt=cpusmt+t1-t11
  if (.not.noprin) then
  call gettim(td,tds)
  call gettim(tairy,tairys)
  call gettim(tb,tbs)
  call gettim(tp,tps)
  call gettim(tw,tws)
  call gettim(twf,twfs)
  write (6, 290) tbs, tps, tds, tairys, twfs, tws
290   format (' ** CPU TIMES:  BASIS=',a,'  POT=',a, &
          '  LOGD=',(a),'  AIRY= ',a,/14x, &
          '  PSI0=',a, '  SMAT=',a)
  write (9, 300)
300   format(/' ** TIMING (CPU/ELAPSED)')
  write (9, 310) tbs,tps,tds,tws,tairys,twfs
310   format('   BAS=',a,' POT=',a,' LOGD=',a, &
         ' SMAT=',a,/ &
         '  AIRY=',a,' PSI0=',a)
  write (9, 320) rtmn,rtmx
320   format(' ** TURNING POINTS:  MIN=', f7.3,'  MAX=',f7.3,' BOHR')
end if

return
end
! ------------------------------------------------------------------------
subroutine rles (a, c, n, m, nmax)
!  subroutine to solve linear equations a * x = c
!  written by:  millard alexander
!  current revision date: 6-apr-2004
!  --------------------------------------------------------------------
!  variables in call list:
!    a:       on input: contains n x n coefficient matrix
!    c:       on input: contains n x m right hand side
!             on return: contains n x m solution matrix
!                    (original c matrix is destroyed)
!    kpvt:    scratch vector of length n
!    n:       size of matrix a, number of rows in matrix c
!    m:       number of columns in matrix c
!    nmax:    maximum row dimension of matrices
!  subroutines called:
!    sgefa, sgesl:  standard linpack routines
!                   note that in some versions izero should be eliminated
!                   from call to sgesl
!  --------------------------------------------------------------------
implicit double precision (a-h,o-z)
integer ierr, izero, m, n, nmax
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
integer :: icol, iptc
#endif
integer :: kpvt(n)
dimension a(1), c(1)
data izero /0/
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
!  factor the matrix
!           return ierr=2 if a is singular
!           return ierr=0 otherwise
#endif
#if defined(HIB_UNIX)  && !defined(HIB_UNIX_IBM) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
call sgefa (a, nmax, n, kpvt, ierr)
#endif
#if defined(HIB_UNIX_IBM)
call dgef(a,nmax,n,kpvt)
ierr=0
#endif
#if defined(HIB_UNIX_CONVEX)
call dgefa (a, nmax, n, kpvt, ierr)
#endif
#if defined(HIB_UNIX)  && !defined(HIB_UNIX_IBM) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
if (ierr .eq. 2) then
!  here if singular matrix
  write (9, 150)
  write (6, 150)
150   format (' *** SINGULAR MATRIX IN SGEFA; ABORT')
  stop
end if
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
!  now solve the linear equations (one-by-one)
iptc = 1
do 200  icol = 1, m
!  iptc points to top of column icol for matrix c stored in packed column fo
#endif
#if defined(HIB_UNIX)  && !defined(HIB_UNIX_IBM) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
  call sgesl (a, nmax, n, kpvt, c(iptc), izero)
#endif
#if defined(HIB_UNIX_IBM)
  call dges (a, nmax, n, kpvt, c(iptc), izero)
#endif
#if defined(HIB_UNIX_CONVEX)
  call dgesl (a, nmax, n, kpvt, c(iptc), izero)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
  iptc = iptc + nmax
200 continue
!  solution matrix now in c
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
call dgetrf(n,n,a,nmax,kpvt,ierr)
if (ierr.ne.0) then
    write (6,210) ierr
    write (9,210) ierr
210     format ('*** IERR = ',i2,' IN DGETRF; ABORT')
    stop
endif
call dgetrs('N',n,m,a,nmax,kpvt,c,nmax,ierr)
if (ierr.ne.0) then
    write (6,220) ierr
    write (9,220) ierr
220     format ('*** IERR = ',i2,' IN DGETRS; ABORT')
    stop
endif
#endif
return
end
! -------------------------------------------------------------------
! -------------------------------------------------------------------
subroutine logdb (z, nmax, nch, rmin, rmax, nsteps, &
                  eshift, iread, iwrite, tl, tp, twf, v2)
!     routine to initialise the log derivative matrix, y, at r = rmin,
!     and propagate it from rmin to rmax using the method described in
!     d.e.manolopoulos, j.chem.phys., 85, 6425 (1986)
!     references to this paper appear below in comments
!     author:  david manolopoulos and millard alexander
!     current revision date (the propagator): 28-nov-2007
!     revised on 30-mar-2012 by q. ma for stream I/O of wfu files
!     current revision: 20-apr-2012 by q. ma
!**********************************************************************
!***   this integrator is not released for general public use      ****
!***   all use must be by specific prior arrangement with:         ****
!***     millard alexander, department of chemistry,               ****
!***     university of maryland, college park, md, 20742           ****
!***     tel: 1.301.405.1823; email: mha@umd.edu                   ****
!***   no part of this program may be copied or used for other     ****
!***   purposes without the author's permission.                   ****
!**********************************************************************
!  ------------------------------------------------------------------
!     variables in call list:
!     z             array of dimension (nmax,nch)
!                   on return z contains the log derivative matrix
!                   at r = rmax
!     nmax          leading dimension of arrays z and w
!     nch           number of coupled equations
!                   = actual order of matrices z and w
!     rmin,rmax     integration range limits
!                   rmin is assumed to lie inside the classically
!                   forbidden region for all channels
!     nsteps        number of sectors partitioning integration range
!                   this version uses a constant step size throughout
!                   the integration range, h = (rmax-rmin)/(2*nsteps)
!           note:  if nsteps=0, then this routine only initializes the logd
!                  matrix; no propagation is done
!     eshift        energy shift for diagonal elements of the coupling
!                   matrix at subsequent energies:
!                   new energy = first energy + eshift
!     iread,iwrite  logical variables for i/o of energy-independent
!                   information from/to unit 11
!     tl            elapsed time used in logd integration
!                   exclusive of calls to potential and ground state
!                   wavefunction (in photodissociation calculations)
!     tp            elapsed time used in calls to potential
!     twf           elapsed time used in determination of ground
!                   state wavefunction (this should be zero for scattering)
!                   exclusive of calls to potential
!                   this timing information comes from repeated calls
!                   of the form call mtime(elapsed) where "elapsed"
!                   is the current elapsed time in seconds
!     blas routines daxpy and dscal are used in o(n*n) loops
!     symmetry of the matrices z and w is not exploited in these loops,
!     and blas routines are not used for o(n) loops
!  ------------------------------------------------------------------
use mod_coqvec, only: mxphot, nphoto, q ! q is an output of this subroutine
use mod_ancou, only: ancou_type
use mod_wave, only: irec, ifil, nchwfu, nrlogd, iendwv, get_wfu_logd_rec_length

use funit
use mod_phot, only: photof, wavefn, writs
use mod_hiutil, only: daxpy_wrapper
use mod_hiutil, only: mtime
use mod_himatrix, only: mxma
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
use mod_himatrix, only: syminv
#endif
use mod_hivector, only: dset, matmov
use mod_hibrid1, only: potmat
use mod_hiblas, only: dscal

implicit double precision (a-h,o-z)
real(8), intent(out) :: z(nmax*nch)
integer, intent(in) :: nmax
integer, intent(in) :: nch
real(8), intent(inout) :: rmin
real(8), intent(inout) :: rmax
integer, intent(inout) :: nsteps
real(8), intent(in) :: eshift
logical, intent(in) :: iread
logical, intent(in) :: iwrite
real(8), intent(out) ::  tl
real(8), intent(out) ::  tp
real(8), intent(out) ::  twf
type(ancou_type), intent(in) :: v2

!     wref          scratch array of dimension nch
!                   used as workspace for the reference potential
!     z1,z2         scratch arrays of dimension nch
!                   used as workspace for the homogeneous
!                   half-sector propagators
!     scr1,scr2     scratch arrays of dimension nch
!                   used as workspace in calls to smxinv
real(8) :: wref(nch)
real(8) :: z1(nch)
real(8) :: z2(nch)
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
real(8) :: scr1(nch)
real(8) :: scr2(nch)
#endif


integer ich, icode, icol, idiag, ierr, ij, irow, istep, kstep, &
        ncol, ndiag, nrow
!      real arg, d3, d6, eight, eshift, fac, h, half, hi, one, r, rmax,
!     :     rmin, sqrt, tan, tanh, t0, t0w, t1, t1w, tf, tfw, th, tl,
!     :     tn, tp,  zero, wdiag
!      real w, z
!      real scr1, scr2, wref, z1, z2
!      external mtime, potmat, daxpy, smxinv, dscal
!     matrices z and w are stored column by column as one-dimensi
data zero,  one,  two, three,six,  eight &
    / 0.d0, 1.d0, 2.d0, 3.d0, 6.d0, 8.d0 /
data izero, ione /0, 1/
!     amat,bmat     scratch matrices of dimension (nmax,nch)
!                   used as workspace for the dissociation overlap
!     amat and bmat are stored column by column as one-dimensi
real(8),allocatable :: amat(:)
real(8),allocatable :: bmat(:)

!     w             scratch array of dimension (nmax,nch)
!                   used as workspace for the potential matrix
real(8), allocatable :: w(:)

integer :: nqmax
real(8) :: simpwt  ! simpson quadrature weight
nqmax = 0
!

!     make sure that rmin, rmax and nsteps are the same if second
!     energy calculation
if (iwrite) write (FUNIT_QUAD_MAT) rmin, rmax, nsteps
if (iread ) read  (FUNIT_QUAD_MAT) rmin, rmax, nsteps
simpwt = 0.d0
if (nsteps .ne. 0) then
  h = (rmax - rmin) / (2 * nsteps)
  hi = one / h
  simpwt= h/three
  d3 = h * h / three
  d6 = - h * h / six
  half = h / two
  if (photof) then
     nqmax=nch*nphoto
!          initialize q vector
#if defined(HIB_UNIX) || defined(HIB_CRAY)
     call dset(nqmax, zero, q, 1)
#endif
   end if 
end if

allocate(w(nmax*nch))
!     row, column and diagonal increments for matrices z and w
nrow = 1
ncol = nmax
ndiag = nmax + 1  ! stride between 2 consecutive diagonal elements (nmax to skip a row, then add 1 to reach the next column)
!     obtain coupling matrix, w, at r = rmin
!     diagonal elements must be shifted at subsequent energies
tp = zero
tpw = zero
twf = zero
twfw = zero
call mtime(tf,tfw)
if (iread) then
   icol = 1
   do 5  ich = 1, nch
      read (FUNIT_QUAD_MAT) (w(ij), ij = icol, icol + nch - 1)
      icol = icol + ncol
5    continue
   idiag = 1
   do 10  ich = 1, nch
      w(idiag) = w(idiag) - eshift
      idiag = idiag + ndiag
10    continue
if (photof) read (FUNIT_QUAD_MAT) (q(i), i=1, nqmax)
else
   istep = 0
   r = rmin
   call mtime(t0,t0w)
   call potmat(w, r, nch, nmax, v2)
   call mtime(t1,t1w)
   tp = tp + t1 - t0
   tpw = tpw + t1w - t0w
!     determine ground state wavefunction times dipole moment at beginning of
!     first sector
   if (photof) then
     call mtime(t0,t0w)

     call ground(q, r, nch, nphoto, mxphot)
!     wt is weight for simpson's rule quadrature, isimp is alternator for
!     simpson's rule quadrature
     wt = simpwt
     isimp=1
!     multiply psi(0) mu(0) by simpson's rule wt at first point
#if defined(HIB_UNIX) || defined(HIB_CRAY)
     call dscal (nqmax, wt, q, 1)
#endif
     call mtime(t1,t1w)
     twf=twf+t1-t0
     twfw=twfw+t1w-t0w
   endif
!
   if (iwrite) then
      icol = 1
      do 15  ich = 1, nch
         write (FUNIT_QUAD_MAT) (w(ij), ij = icol,icol+nch-1)
         icol = icol + ncol
15       continue
      if (photof) write (FUNIT_QUAD_MAT) (q(i), i=1, nqmax)
   endif
endif
!     use diagonal approximation to wkb initial value for log
!     derivative matrix  (eqn 16)
!     rmin is assumed to lie inside the classically forbidden
!     region in all channels  (w(ii) > 0)
!     first zero out z matrix
icol = 1
do 20  ich = 1, nch
#if defined(HIB_UNIX) || defined(HIB_CRAY) || defined(HIB_MAC)
   call dset (nch, zero, z(icol), nrow)
#endif
   icol = icol + ncol
20 continue

idiag = 1
do  30 ich = 1, nch
   wdiag = w(idiag)
   if (wdiag .le. 0) then
     write (9, 25) ich, wdiag
     write (6, 25) ich, wdiag
25      format (' *** WAVEVECTOR > 0 IN CHANNEL', i5, &
        '     IN LOGD INITIALIZATION;  ABORT ***', &
       /'     V - E =',1pe15.7)
     call exit
   end if
   z(idiag) = sqrt(wdiag)
   idiag = idiag + ndiag
30 continue

if (nsteps .le. 0) then
  if (photof) then
    write (9, 35) nsteps
    write (6, 35) nsteps
35     format (' *** WARNING: IN LODG PHOTOF=TRUE AND NSTEPS=',i3, &
            ' SHOULD BE .GT. 0 ***')
    goto 261
  endif
endif

if (nsteps /= 0) then
  !     propagate z matrix from rmin to rmax

  !     with a constant step size it is convenient to propagate
  !     the matrix z = h*y rather than the log derivative matrix, y
  !     (eqns 10, 12, 13 and 14 are simply multiplied through by h)
  fac = h
  icol = 1
  do ich = 1, nch

#if defined(HIB_UNIX) || defined(HIB_CRAY)
     call dscal (nch, fac, z(icol), nrow)
#endif
     icol = icol + ncol
  end do

  allocate(amat(nmax*nch))
  allocate(bmat(nmax*nch))

  do 250  kstep = 1, nsteps
  !     apply quadrature contribution at beginning of sector,
  !     r = a
  !     after this loop z contains z(a)+(h^2/3)w(a)
  !     todo : handle the case where w is lower triangular
  !     (uninitialized values in upper triangle), see issue 49
     fac = d3
     icol = 1
     do  50 ich = 1, nch
        call daxpy_wrapper (nch, fac, w(icol), nrow, z(icol), nrow)
        icol = icol + ncol
  50    continue
  !     the reference potential for the sector is the diagonal
  !     of the coupling matrix evaluated at the centre of the
  !     sector, r = c  (eqn 15)
     if (iread) then
        read (FUNIT_QUAD_MAT) (wref(ich), ich = 1, nch)
        do 60  ich = 1, nch
           wref(ich) = wref(ich) - eshift
  60       continue
     else
        istep = istep + 1
        r = rmin + istep * h
        call mtime(t0,t0w)

        call potmat(w, r, nch, nmax, v2)
        call mtime(t1,t1w)
        tp = tp + t1 - t0
        tpw = tpw + t1w - t0w
        idiag = 1
        do 70  ich = 1, nch
           wref(ich) = w(idiag)
           idiag = idiag + ndiag
  70       continue
        if (iwrite) then
           write (FUNIT_QUAD_MAT) (wref(ich), ich = 1, nch)
        endif
     endif
  !     adjust quadrature contribution at r = a to account for
  !     sector reference potential  (eqn 11)
  !     after this loop z contains z(a)+(h^2/3)w1(a)
     fac = - d3
     idiag = 1
     do  80 ich = 1, nch
        z(idiag) = z(idiag) + fac * wref(ich)
        idiag = idiag + ndiag
  80    continue
  !     evaluate homogeneous half sector propagators  (eqn 10)
  !     z(i) = h * y(i), i=1,2,3,4
     do 85  ich = 1, nch
        arg = half * sqrt(abs(wref(ich)))
           if (wref(ich) .lt. zero) then
              tn = tan(arg)
              z1(ich) = arg / tn - arg * tn
              z2(ich) = arg / tn + arg * tn
           else
              th = tanh(arg)
              z1(ich) = arg / th + arg * th
              z2(ich) = arg / th - arg * th
           endif
  !           z3(ich) = z2(ich)
  !           z4(ich) = z1(ich)
  85    continue
  !     propagate z matrix across the first half sector
  !     from r = a to r = c  (eqn 14)
     idiag = 1
     do  90 ich = 1, nch
        z(idiag) = z(idiag) + z1(ich)
        idiag = idiag + ndiag
  90    continue

  !     z now contains z(a)+z1(a,c)
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
     call smxinv(z, nmax, nch, scr1, scr2, ierr)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
    call syminv(z,nmax,nch,ierr)
#endif
     if (ierr .ne. 0) then
        write (6, 9000) kstep, ierr
        write (9, 9000) kstep, ierr
        call exit
     endif
  !     z now contains [z(a)+z1(a,c)]^-1
     icol = 1
     do  95 ich = 1, nch
        fac = z2(ich)
        call dscal (nch, fac, z(icol), nrow)
        icol = icol + ncol
  95    continue
  !     z now contains [z(a)+z1(a,c)]^-1 z2(a,c)
  !     if photodissociation calculation or wavefunction desired:
  !     save this matrix, which is hg(a,m), in amat
     if (photof .or. wavefn) &
       call matmov(z, amat, nch, nch, nmax, nmax)
     irow = 1
     do  110 ich = 1, nch
        fac = - z2(ich)
        call dscal (nch, fac, z(irow), ncol)
        irow = irow + nrow
  110    continue
  !     z now contains - z3(a,c) [z(a)+z1(a,c)]^-1 z2(a,c)
     idiag = 1
     do  120 ich = 1, nch
        z(idiag) = z(idiag) + z1(ich)
        idiag = idiag + ndiag
  120    continue

  !     evaluate quadrature contribution at sector mid-point,
  !     r = c  (eqn 12)     (first energy calculation only)
     if (iread) then
        icol = 1
        do  125 ich = 1, nch
           read (FUNIT_QUAD_MAT) (w(ij), ij = icol,icol+nch-1)
           icol = icol + ncol
  125       continue
     else
        fac = d6
        icol = 1
        do  130 ich = 1, nch
           call dscal (nch, fac, w(icol), nrow)
           icol = icol + ncol
  130       continue
        idiag = 1
        do  140 ich = 1, nch
           w(idiag) = one
           idiag = idiag + ndiag
  140       continue

#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
        call smxinv(w, nmax, nch, scr1, scr2, ierr)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
        call syminv(w,nmax,nch,ierr)
#endif
        if (ierr .eq. 2) then
           icode = 2
           write (9,9000) kstep, icode
           call exit
        endif
        idiag = 1

        do  150 ich = 1, nch
           w(idiag) = w(idiag) - one
           idiag = idiag + ndiag
  150       continue
        if (iwrite) then
           icol = 1
           do  155 ich = 1, nch
              write (FUNIT_QUAD_MAT) (w(ij), ij = icol,icol+nch-1)
              icol = icol + ncol
  155          continue
        endif
     endif
  !     apply quadrature contribution at sector mid-point
  !     corrections to z4(a,c) and z1(c,b) are applied
  !     simultaneously  (eqn 13)
     fac = eight
     icol = 1
     do  160 ich = 1, nch
        call daxpy_wrapper (nch, fac, w(icol), nrow, z(icol), nrow)
        icol = icol + ncol
  160    continue

  !     propagate z matrix across the second half sector
  !     from r = c to r = b   (eqn 14)
     idiag = 1
     do  170 ich = 1, nch
        z(idiag) = z(idiag) + z1(ich)
        idiag = idiag + ndiag
  170    continue
  !     at this point z contains z(c) + z1(c,b)
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
     call smxinv(z, nmax, nch, scr1, scr2, ierr)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
     call syminv(z,nmax,nch,ierr)
#endif
     if (ierr .eq. 2) then
        icode = 3
        write (9,9000) kstep, icode
        call exit
     endif
  !     at this point z contains [z(c) + z1(c,b)]^-1
     icol = 1
     do 175 ich = 1, nch
        fac = z2(ich)
        call dscal (nch, fac, z(icol), nrow)
        icol = icol + ncol
  175    continue
  !     z now contains [z(c)+z1(c,b)]^-1 z2(c,b)
  !     if photodissociation calculation or wavefunction desired:
     if (photof .or. wavefn) then
  !     first save this matrix, which is g(c,b), in bmat
       call matmov(z, bmat, nch, nch, nmax, nmax)
  !     use bmat and w as temporary storage here
       call mxma (amat, 1, nmax, bmat, 1, nmax, w, 1, nmax, &
                  nch, nch, nch)
  !     w now contains the matrix g(a,m)g(m,b)=g(a,b)
  !     if wavefunction desired, save this matrix
      if (wavefn .and. writs) then
         irec = irec + 1
  !     nrlogd is the number of LOGD records - used to seek the wfu file
         nrlogd = nrlogd + 1
         write (ifil, err=950) r - h, r, (w(i), i=1, nch)
         icol = 1
         do ich = 1, nch
            write (ifil, err=950) (w(icol - 1 + i), i=1, nch)
            icol = icol + nmax
         end do
         write (ifil, err=950) 'ENDWFUR', char(mod(irec, 256))
         iendwv = iendwv + get_wfu_logd_rec_length(nchwfu, 0)
      endif
    endif

  !     accumulate overlap matrix with ground state
  !     if photodissociation calculation
    if (photof) then
  !     premultiply g(a,b) by wt psi(a) mu(a)
  !     use bmat as temporary storage here
       call mxma(q,nch,1,w,1,nmax,bmat,nch,1,nphoto,nch,nch)
  !     bmat now contains [...+wt*psi(a)mu(a)] g(a,b)
  !       stored as successive columns with each initial
  !       state corresponding to a column
  !     now determine psi(b)mu(b), save this in q
       rnew = rmin + (istep+1) * h
       if (iread) then
         read (FUNIT_QUAD_MAT) (q(i), i=1, nqmax)
       else
         call mtime(t0,t0w)
         call ground(q, rnew, nch, nphoto, mxphot)
  ! recalculate simpson's rule wt for this point
         wt=(3.d0+isimp)*simpwt
         isimp=-isimp
  !     multiply psi(b) mu(b) by simpson's rule wt at r=rnew
         call dscal (nqmax, wt, q, 1)
         call mtime(t1,t1w)
         twf=twf+t1-t0
         twfw=twfw+t1w-t0w
         if (iwrite) write (FUNIT_QUAD_MAT) (q(i), i=1, nqmax)
       endif
  !     add wt*psi(b)mu(b) to bmat and resave
       fac=one
       call daxpy_wrapper(nqmax, fac, bmat, 1, q, 1)
     endif
  !     now premultiply [z(c)+z1(c,b)]^-1 z2(c,b) by z2(c,b)
     irow=1
     do  190 ich = 1, nch
       fac = - z2(ich)
       call dscal (nch, fac, z(irow), ncol)
       irow = irow + nrow
  190    continue
  !      z now contains - z2(c,b) [z(c)+z1(c,b)]^-1 z2(c,b)
     idiag = 1
     do  195 ich = 1, nch
        z(idiag) = z(idiag) + z1(ich)
        idiag = idiag + ndiag
  195    continue
  !     apply reference potential adjustment to quadrature
  !     contribution at r = b  (eqn 11)
     fac = - d3
     idiag = 1

     do  200 ich = 1, nch
        z(idiag) = z(idiag) + fac * wref(ich)
        idiag = idiag + ndiag
  200    continue
  !     obtain coupling matrix, w, at end of sector
     if (iread) then
        icol = 1
        do  205 ich = 1, nch
           read (FUNIT_QUAD_MAT) (w(ij), ij = icol, icol + nch - 1)
           icol = icol + ncol
  205       continue
        idiag = 1
        do  210 ich = 1, nch
           w(idiag) = w(idiag) - eshift
           idiag = idiag + ndiag
  210       continue
     else
        istep = istep + 1
        r = rmin + istep * h
        call mtime(t0,t0w)
        call potmat(w, r, nch, nmax, v2)
        call mtime(t1,t1w)
        tp = tp + t1 - t0
        tpw = tpw + t1w - t0w
        if (iwrite) then
           icol = 1
           do  215 ich = 1, nch
              write (FUNIT_QUAD_MAT) (w(ij), ij = icol,icol+nch-1)
              icol = icol + ncol
  215          continue
        endif
     endif
  !     apply quadrature contribution at r = b  (eqn 12)
     fac = d3
     icol = 1
     do 220  ich = 1, nch
        call daxpy_wrapper (nch, fac, w(icol), nrow, z(icol), nrow)
        icol = icol + ncol
  220    continue
  !     propagation loop ends here
  250 continue

  deallocate(w)
  deallocate(bmat)
  deallocate(amat)

  !     recover the log derivative matrix, y, at r = rmax
  fac = hi
  icol = 1
  do ich = 1, nch
     call dscal (nch, fac, z(icol), nrow)
     icol = icol + ncol
  end do

end if

261 call mtime(tl,tlw)
tl = tl - tf - tp - twf
tlw = tlw - tfw - tpw - twfw

return
!
950 write (0, *) ' *** ERROR WRITING WFU FILE (LOGD). ABORT'
call exit()
return
!
9000 format(' *** MATRIX INVERSION ERROR IN LOGDB AT KSTEP =', &
 i4,  /'            ERROR OCCURRED AT ROW =', i2, &
 '; ABORT ***')
end
! ----------------------------------------------------------------------
subroutine runlog (z, &
                   r, rend, &
                   spac, eshift, itwo, twoen, tl, tlw, tp, tpw, &
                   twf, twfw, prlogd, noprin, ipos, nch, nmax, v2, isteps, nsteps)
!     log-derivative propagator from r to r = rend
!     the logd code is based on the improved log-derivative method
!     for reference see  d.e.manolopoulos, j.chem.phys., 85, 6425 (1986)
!     author: david manolopoulos
!     current revision (algorithm): 5-may-1997 by mha
!  ------------------------------------------------------------------
!     variables in call list:
!     z             on return: contains the log-derivative matrix
!                   at r = rend
!     wref,z1,z2    scratch vectors of length nch
!     scr1, scr2,
!     r             on entry: r = initial interparticle separation
!                   on return: r = final interparticle separation
!     rend          ending point for log-derivative integration
!            note:  if rend = r, then the log-derivative propagator is not
!                   used, but the z matrix is still initialized
!
!     spac          step size for integration
!     eshift        energy shift for diagonal elements of the coupling
!                   matrix at subsequent energies -
!     twoen         if .true., if calculation is to be performed at more
!                   than one energy
!     itwo          if twoen = .true. and itwo = 0,  then subroutine is being
!                   called at first energy of multiple energy calculation, so
!                   energy independent information will be written to unit 11
!                   if twoen = .true. and itwo > 0,  then subroutine is being
!                   called at subsequent energy of multiple energy calculation
!                   so energy independent information will be read from unit 1
!     tl,tlw        cpu and wall clock times for log-derivative
!                   integration exclusive of potential evaluation
!     tp,tpw        cpu and wall clock times for potential evaluation
!     twf,twfw      cpu and wall clock times for evaluation of ground
!                   state wf (for photodissociation calculation)
!                   this timing information is provided by repeated
!                   calls of the form call mtime(cpu,wall), where cpu
!                   and wall refer to the current cpu and wall clock
!                   times in seconds
!     prlogd         if .true., then lower triangle of the z
!                   matrix is printed out at end of log-derivative
!                   integration
!     noprin        if .true., then all printing is suppressed
!     ipos          if .true., then 132 column printer
!                   if .false., then 80 column printer
!     nch           number of coupled equations
!     nmax          leading dimension of arrays z and w
!  ------------------------------------------------------------------
use mod_coqvec, only: nphoto, q
use mod_ancou, only: ancou_type
use mod_hibrid2, only: mxoutd, mxoutr
use mod_phot, only: photof
implicit double precision (a-h, o-z)
type(ancou_type), intent(in) :: v2
real(8), intent(out) :: z(nmax*nch)
real(8), intent(inout) :: r
real(8), intent(in) :: rend
real(8), intent(in) :: spac
real(8), intent(in) :: eshift
integer, intent(in) :: itwo
logical, intent(in) :: twoen
real(8), intent(out) :: tl
real(8), intent(out) :: tlw
real(8), intent(out) :: tp
real(8), intent(out) :: tpw
real(8), intent(out) :: twf
real(8), intent(out) :: twfw
logical, intent(in) :: prlogd
logical, intent(in) :: noprin
logical, intent(in) :: ipos
integer, intent(in) :: nch
integer, intent(in) :: nmax
integer, intent(inout) :: isteps
integer, intent(inout) :: nsteps

!      real eshift, r, rend, rmax, rmin, spac, tl, tlw, tp, tpw
!      real z
!  internal logical variables
logical iread, iwrite, print

!  z, w, amat, and bmat are stored column by column in one dimensional arrays

UNUSED_DUMMY(tlw)
tlw = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(tpw)
tpw = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED_DUMMY(twfw)
twfw = 0.0  ! silences warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.

iwrite = twoen .and. (itwo.eq.0)
iread =  twoen .and. (itwo.gt.0)
print = .not. noprin
!  set up for log-derivative integration
rmin = r
rmax = rend
if (rmax .eq. rmin) then
!  here if no logd propagation, just initialization of matrix
  nsteps = 0
else
  if (isteps .eq. 0) then
    nsteps = int ((rmax - rmin) / spac) + 1
    isteps=1
  endif
  if (print) then
     write (9, 10) rmin, rmax, spac, nsteps, nch
10      format(/' ** LOGD PROPAGATION:', /, &
           '         RSTART =', f8.3,'  REND =', f8.3, &
       '  SPAC =', f8.4,'  NSTEP =', i5, '  NCH =', i5)
     write (6, 20) rmin, rmax, spac, nsteps
20      format (' ** LOGD:  RSTART =',f7.3, &
             '   REND =', f7.3,'   SPAC = ', f7.3, &
             '   NSTEP =', i5)
  endif
endif
!  logdb initialises the log-derivative matrix at r = rmin
!  and propagates it from rmin to rmax
call logdb (z, &
            nmax, nch, &
            rmin, rmax, nsteps, eshift, iread, iwrite, &
            tl, tp, twf, v2)

!  on return: z(i,j) i=1,nch j=1,nch  now contains the log-derivative matrix
!  at r = rmax (the final interparticle separation)
r = rmax
!  print out log-derivative matrix at end of logd integration
!  ( if desired )
if (prlogd .and. print) then

  write (9, 40) r
40   format(/' ** LOG-DERIVATIVE MATRIX AFTER LOGDB; R = ', 1pe15.8)
  call mxoutd (9, z, nch, nmax, 0, ipos)
  if (photof) then
    write (9, 45) r
!  print out <psi0 | mu matrix at the end of logd-integration
!  if desired
45     format (/' ** <PSI0|MU| VECTOR AFTER LOGDB; R =',1pe15.8)
    isym=0
    call mxoutr(9, q, nch, nphoto, nch, isym, ipos)
  endif
end if

return
end
! -------------------------------------------------------------------------
subroutine scairy (z, scai, scbi, scaip, scbip, zeta)
!     scaled airy functions and derivatives
!     this program written by d.e. manolopoulos (sept. 1986)
!     current revision date: sept-1986
!  ----------------------------------------------------------------
!     for  z .lt. (-5.0d0)
!     ai(z) = scai*cos(zeta) + scbi*sin(zeta)
!     bi(z) = scbi*cos(zeta) - scai*sin(zeta)
!     ai'(z) = scaip*cos(zeta) + scbip*sin(zeta)
!     bi'(z) = scbip*cos(zeta) - scaip*sin(zeta)
!     where zeta = (2/3)*(-z)**(3/2) + pi/4
!     for  (-5.0d0) .le. z .le. (+0.0d0)
!     ai(z) = scai
!     bi(z) = scbi
!     ai'(z) = scaip
!     bi'(z) = scbip
!     and zeta = 0
!     for  z > 0.d0
!     ai(z) = scai*exp(-zeta)
!     bi(z) = scbi*exp(+zeta)
!     ai'(z) = scaip*exp(-zeta)
!     bi'(z) = scbip*exp(+zeta)
!     where zeta = (2/3)*(+z)**(3/2)
!  ----------------------------------------------------------------
!     evaluation of the functions is based on a number of
!     chebyshev expansions
!
!     this version is suitable for machines with full word precision
!  ----------------------------------------------------------------
implicit double precision (a-h,o-z)
double precision a, b, c, c1, c2, df, dg, ex, exp1z, exp2z, f, &
                 g, pib4, root4z, rootz, rt3, scai, scaip, scbi, &
                 scbip, t, t2, xeps, y, z, zcube, zeta, zsq
data c1   / 3.55028053887817239d-01 /
data c2   / 2.58819403792806798d-01 /
data rt3  / 1.73205080756887729d+00 /
data pib4 / 7.85398163397448310d-01 /
zeta = 0.0d0
xeps = 0.0d0
! ----------------------------------------------------------------------------
!  here if near enough origin to use 3 term power series
if ( dabs(z) .le. 0.025d0) then
  zsq = z * z
  zcube = zsq * z
!  evaluate power series ( three terms is sufficient for abs(x) < 0.025)
  df = 1.d0 + zcube / 6.d0 + zcube * zcube / 180.d0
  dg = z * (1.d0 + zcube / 12.d0 + zcube * zcube / 504.d0)
  scai = c1 * df - c2 * dg
  scbi = rt3 * (c1 * df + c2 * dg)
!  now for derivatives
  df = zsq / 2.d0 + zsq * zcube / 30.d0
  dg = 1.d0 + zcube / 3.d0 + zcube * zcube / 72.d0
  scaip = c1 * df - c2 * dg
  scbip = rt3 * (c1 * df + c2 * dg)
!  scale the functions by exp(zeta) if z .gt. 0
  if (z .gt. 0.d0) then
    rootz = dsqrt(z)
    zeta = 2.0d0 * z * rootz / 3.0d0
    ex = exp(zeta)
    scai = scai * ex
    scaip = scaip * ex
    scbi = scbi / ex
    scbip = scbip / ex
  end if
  return
end if
if (z.lt.(+9.0d0)) go to 10
rootz = dsqrt(z)
root4z = dsqrt(rootz)
zeta = 2.0d0*z*rootz/3.0d0
t = 36.0d0/zeta - 1.0d0
y = ((((((((( +1.16537795324979200d-15*t &
  -1.16414171455572480d-14)*t +1.25420655508401920d-13)*t &
  -1.55860414100340659d-12)*t +2.21045776110011276d-11)*t &
  -3.67472827517194031d-10)*t +7.44830865396606612d-09)*t &
  -1.95743559326380581d-07)*t +7.44672431969805149d-06)*t &
  -5.28651881409929932d-04)*t +2.81558489585006298d-01
scai = y/root4z
y = ((((((((((( +4.50165999254528000d-15*t &
  +1.56232018374502400d-14)*t +5.26240712559918080d-14)*t &
  +2.97814898856618752d-13)*t +1.97577620975625677d-12)*t &
  +1.53678944110742706d-11)*t +1.45409933537455235d-10)*t &
  +1.71547326972380087d-09)*t +2.61898617129147064d-08)*t &
  +5.49497993491833009d-07)*t +1.76719804365109334d-05)*t &
  +1.12212109935874117d-03)*t +5.65294557558522063d-01
scbi = y/root4z
y = ((((((((( +1.20954638924697600d-15*t &
  -1.21281218539020800d-14)*t +1.31303723724964224d-13)*t &
  -1.64152781754533677d-12)*t +2.34672185025709461d-11)*t &
  -3.94507329122119338d-10)*t +8.13125005420910243d-09)*t &
  -2.19736365932356533d-07)*t +8.83993515227257822d-06)*t &
  -7.43456339972080231d-04)*t -2.82847316336379200d-01
scaip = y*root4z
y = ((((((((((( -4.59170437029478400d-15*t &
  -1.59840960512122880d-14)*t -5.41258863340784640d-14)*t &
  -3.07414589507261184d-13)*t -2.04866616770522650d-12)*t &
  -1.60321415915690897d-11)*t -1.52922073861488292d-10)*t &
  -1.82445639488695332d-09)*t -2.83250890588806503d-08)*t &
  -6.11130377639012647d-07)*t -2.07842147963678572d-05)*t &
  -1.56350017663858255d-03)*t +5.62646283094843014d-01
scbip = y*root4z
return
10 if (z.lt.(+4.5d0)) go to 20
rootz = dsqrt(z)
zeta = 2.0d0*z*rootz/3.0d0
exp1z = dexp(zeta-2.5d0*z)
exp2z = dexp(zeta-2.625d0*z)
t = 4.0d0*z/9.0d0 - 3.0d0
y = ((((((((((((((((((((( +9.69081960415394529d-11*t &
  +3.24436136050920784d-10)*t -3.57419513430644674d-09)*t &
  -3.84461320827974687d-09)*t +8.88116699085949212d-08)*t &
  -6.26105174374717557d-08)*t -1.69051051004298110d-06)*t &
  +3.80731416363041759d-06)*t +2.43840529113057777d-05)*t &
  -9.74379632673654766d-05)*t -2.45324254437931970d-04)*t &
  +1.69517926953312785d-03)*t +1.19638433540225211d-03)*t &
  -2.15255594590357451d-02)*t +9.33777073522844198d-03)*t &
  +1.98716159257796883d-01)*t -2.54001858882057718d-01)*t &
  -1.27148775197878180d+00)*t +2.52046376168394778d+00)*t &
  +5.04987271423387057d+00)*t -1.33120978544419281d+01)*t &
  -9.34903846550381088d+00)*t +3.10330812950257837d+01
scai = y*exp1z
y = (((((((((((((((((((((((( +3.79210935744593920d-14*t &
  -4.16346635040194560d-14)*t -3.63110681886588928d-13)*t &
  +1.38932592029414195d-12)*t -4.00489068810888806d-12)*t &
  +1.39019501834951721d-11)*t -4.50877182237241508d-11)*t &
  +1.38942309844733264d-10)*t -3.92503498108710093d-10)*t &
  +1.20125005161756928d-09)*t -3.14234550677825531d-09)*t &
  +1.03100587323694771d-08)*t -2.35240060783126760d-08)*t &
  +8.98525670958611253d-08)*t -1.57273011181242048d-07)*t &
  +7.77696763289738864d-07)*t -8.40211181188135235d-07)*t &
  +6.34887361301864569d-06)*t -2.73464023289055762d-06)*t &
  +4.54606729925166230d-05)*t +2.20459155042947089d-06)*t &
  +2.58823388957588056d-04)*t +7.31023768389466446d-05)*t &
  +1.01013806904596356d-03)*t +2.64794416332118755d-04)*t &
  +1.97499785553709145d-03
scbi = y/exp1z
y = ((((((((((((((((((((( -4.40679918437492851d-10*t &
  +1.30954945449348301d-10)*t +1.30052079376596751d-08)*t &
  -2.21315827945437064d-08)*t -2.56850909380644963d-07)*t &
  +8.66960855365698346d-07)*t +3.75622307499741911d-06)*t &
  -2.15396233361107222d-05)*t -3.55804094667597110d-05)*t &
  +3.95317852914037711d-04)*t +5.03369361986934094d-05)*t &
  -5.54634417403436820d-03)*t +5.29658186908372832d-03)*t &
  +5.91311623537658225d-02)*t -1.09446664596286554d-01)*t &
  -4.63589435529194219d-01)*t +1.25323269822030972d+00)*t &
  +2.50138108959469254d+00)*t -9.12668774193995449d+00)*t &
  -8.14385732036876466d+00)*t +4.00134082550833019d+01)*t &
  +1.15396202931444799d+01)*t -8.17378314444550419d+01
scaip = y*exp1z
y = (((((((((((((((((((((((( -1.12976379481423872d-13*t &
  +2.84163275199873024d-13)*t +9.21367859618119680d-14)*t &
  -6.47465116933029888d-13)*t +5.66210442158931968d-13)*t &
  -3.03158042458901709d-12)*t +1.32640217809876419d-11)*t &
  -3.03558223041639219d-11)*t +5.32290407073565901d-11)*t &
  +1.67561690905544950d-11)*t -3.35234276365918044d-10)*t &
  +2.92807773020050397d-09)*t -8.76900994127464369d-09)*t &
  +4.69138029321003869d-08)*t -1.00929917942876779d-07)*t &
  +5.40401934648687824d-07)*t -8.19977129258456927d-07)*t &
  +5.13367651438974580d-06)*t -4.77800617725922708d-06)*t &
  +4.02415391117897098d-05)*t -1.74571192912274417d-05)*t &
  +2.45332091645215217d-04)*t -2.22916383050374016d-05)*t &
  +1.02535993549737948d-03)*t +5.94033287658300975d-05)*t &
  +2.17420627539345627d-03
scbip = y/exp2z
return
20 if (z.le.(+0.0d0)) go to 40
rootz = dsqrt(z)
zeta = 2.0d0*z*rootz/3.0d0
exp1z = dexp(zeta-1.5d0*z)
exp2z = dexp(zeta-1.375d0*z)
t = 4.0d0*z/9.0d0 - 1.0d0
if (z.lt.(+xeps)) go to 30
y = ((((((((((((((((((((((( +4.97635854909020570d-12*t &
  -3.25024150273916928d-11)*t -5.15773946723072737d-11)*t &
  +8.66802872160017711d-10)*t -9.51292671519803048d-10)*t &
  -1.33268133924677102d-08)*t +4.37061406144179625d-08)*t &
  +1.18943714086308365d-07)*t -8.66980482244589319d-07)*t &
  -2.46768077494905499d-08)*t +1.10610939830483627d-05)*t &
  -1.80475663535516462d-05)*t -9.22213518989192294d-05)*t &
  +3.15767712665407001d-04)*t +4.08626419412850994d-04)*t &
  -3.12704269924340764d-03)*t +6.27899244118607949d-04)*t &
  +1.99062142478229001d-02)*t -2.27427058211322122d-02)*t &
  -7.94869698136278246d-02)*t +1.54261999158247445d-01)*t &
  +1.75618463128730757d-01)*t -5.05223670654169859d-01)*t &
  -1.49695902416050331d-01)*t +6.91290454439828966d-01
scai = y*exp1z
y = (((((((((((((((((((((((((((-8.01144609907912212d-11*t &
  +2.67566208080291037d-10)*t +1.74416971406971503d-10)*t &
  -3.12642164666800066d-09)*t +1.22114569059570056d-08)*t &
  -2.93647730218878800d-08)*t +1.76951994785830839d-08)*t &
  +2.13143266932123830d-07)*t -1.15569603602267288d-06)*t &
  +3.34394065752949896d-06)*t -5.20143492253259528d-06)*t &
  -3.21937890029830155d-06)*t +5.00360593064643409d-05)*t &
  -1.77449408434194908d-04)*t +3.86357389967150628d-04)*t &
  -4.53337922165622921d-04)*t -2.60866378774883161d-04)*t &
  +3.01355585350049504d-03)*t -8.39047077309199055d-03)*t &
  +1.63240267627966090d-02)*t -1.90830727084112485d-02)*t &
  +1.65592661387959142d-02)*t +1.76101803014184860d-02)*t &
  -3.36652019472526494d-02)*t +1.23831258886916327d-01)*t &
  -6.48342330363017516d-02)*t +2.20310550882807725d-01)*t &
  -1.03883014957365224d-02)*t +2.06857611342460346d-01
scbi = y/exp2z
30 y = ((((((((((((((((((((((( -2.31635825886515692d-11*t &
  +8.43840142802870600d-11)*t +3.68028065271203758d-10)*t &
  -2.61043232825754937d-09)*t -4.65110871930215858d-10)*t &
  +4.46164842334855713d-08)*t -9.24599436690579710d-08)*t &
  -4.55809882095931368d-07)*t +2.21024501804834447d-06)*t &
  +1.50251398952558802d-06)*t -2.91830008657289876d-05)*t &
  +3.51391100964982453d-05)*t +2.37966767002002741d-04)*t &
  -7.00969870295148024d-04)*t -9.84923358717942729d-04)*t &
  +6.68935321740601810d-03)*t -1.66398286740112083d-03)*t &
  -3.83618654865390504d-02)*t +4.80463615092658847d-02)*t &
  +1.28359791076466449d-01)*t -2.80267155846714091d-01)*t &
  -2.06049815358004057d-01)*t +7.63522843530878467d-01)*t &
  +6.47699892977822355d-02)*t -8.32940737409625965d-01
scaip = y*exp2z
y = (((((((((((((((((((((((((((+2.69330665471830131d-10*t &
  -1.25313111217921013d-09)*t +1.45057587508619405d-09)*t &
  +5.82827351134571594d-09)*t -3.96093412314305685d-08)*t &
  +1.37346521367521144d-07)*t -2.78927594518121271d-07)*t &
  +2.96531845420687661d-08)*t +2.27734981888044076d-06)*t &
  -1.02295902888535994d-05)*t +2.65515218319523965d-05)*t &
  -3.86457370206378782d-05)*t -1.52212232476268640d-05)*t &
  +2.84765225803690646d-04)*t -9.65798046252914453d-04)*t &
  +2.04618065580453522d-03)*t -2.68702422147972510d-03)*t &
  +8.36839039610090712d-04)*t +6.87131161447866570d-03)*t &
  -2.10563741100004648d-02)*t +4.13290131622517073d-02)*t &
  -5.03310394511775398d-02)*t +5.95467795825179773d-02)*t &
  -1.64213101223235839d-02)*t +5.02536006477020710d-02)*t &
  +5.75601787687195966d-02)*t +1.33220031651076020d-01)*t &
  +7.76356357899154668d-02)*t +2.11213324176049168d-01
scbip = y/exp1z
return
40 if (z.lt.(-5.0d0)) go to 60
t = z/5.0d0
t = -t*t*t
t = 2.0d0*t - 1.0d0
t2 = 2.0d0*t
if (z.gt.(-xeps)) go to 50
a = +1.63586492025000000d-18
b = t2*a -1.14937368283025000d-16
c = t2*b-a +7.06090635856696000d-15
a = t2*c-b -3.75504581033290114d-13
b = t2*a-c +1.70874975807662448d-11
c = t2*b-a -6.56273599013291800d-10
a = t2*c-b +2.09250023300659871d-08
b = t2*a-c -5.42780372893997236d-07
c = t2*b-a +1.11655763472468469d-05
a = t2*c-b -1.76193215080912647d-04
b = t2*a-c +2.03792657403144947d-03
c = t2*b-a -1.61616260941907957d-02
a = t2*c-b +7.87369695059018748d-02
b = t2*a-c -1.88090320218915726d-01
c = t2*b-a +8.83593328666433903d-02
a = t2*c-b +9.46330439565858235d-02
f = t*a-c +7.60869994141726643d-02
a = +1.23340698467000000d-19
b = t2*a -9.05440546731800000d-18
c = t2*b-a +5.83052348377146000d-16
a = t2*c-b -3.26253073273305810d-14
b = t2*a-c +1.56911825099665634d-12
c = t2*b-a -6.40386375393414830d-11
a = t2*c-b +2.18414557202733054d-09
b = t2*a-c -6.11127835033401880d-08
c = t2*b-a +1.37095478225289560d-06
a = t2*c-b -2.39464595313812449d-05
b = t2*a-c +3.13306256975299299d-04
c = t2*b-a -2.90953380590207648d-03
a = t2*c-b +1.76972907074092250d-02
b = t2*a-c -6.17055677164122241d-02
c = t2*b-a +9.52472833367213949d-02
a = t2*c-b -4.32381694223484894d-02
g = t*a-c +3.76828717701544063d-02
scai = f - g*z
scbi = rt3*(f + g*z)
50 a = -2.51308436743000000d-18
b = t2*a +1.65543326242034000d-16
c = t2*b-a -9.49237123028142500d-15
a = t2*c-b +4.68795260455788096d-13
b = t2*a-c -1.96942895842729954d-11
c = t2*b-a +6.93493715818491929d-10
a = t2*c-b -2.01076965264476206d-08
b = t2*a-c +4.69655735896232104d-07
c = t2*b-a -8.59527033121202608d-06
a = t2*c-b +1.18871496270269531d-04
b = t2*a-c -1.18244097697332692d-03
c = t2*b-a +7.87645202148185146d-03
a = t2*c-b -3.14174372672396468d-02
b = t2*a-c +6.20464642445295805d-02
c = t2*b-a -4.83824291776351778d-02
f = t*c-b +2.64808460123486707d-02
a = +5.89382778069400000d-18
b = t2*a -4.04811810887971000d-16
c = t2*b-a +2.42680453287673090d-14
a = t2*c-b -1.25683910148099294d-12
b = t2*a-c +5.55607745069567295d-11
c = t2*b-a -2.06683376304577072d-09
a = t2*c-b +6.35924425685425485d-08
b = t2*a-c -1.58422527393619013d-06
c = t2*b-a +3.11007119112993551d-05
a = t2*c-b -4.64189437787271433d-04
b = t2*a-c +5.00970025411579034d-03
c = t2*b-a -3.62166342717373453d-02
a = t2*c-b +1.53114671641953510d-01
b = t2*a-c -2.69270807740667256d-01
c = t2*b-a -9.61843661149152853d-02
a = t2*c-b +2.07099372879297732d-01
g = t*a-c +9.79943887874547828d-02
scaip = z*z*f - g
scbip = rt3*(z*z*f + g)
return
60 rootz = dsqrt(-z)
root4z = -dsqrt(rootz)
zeta = 2.0d0*(-z)*rootz/3.0d0
t = -250.0d0/(z*z*z) - 1.0d0
a = ((((((((((((( -4.50071772808806400d-15*t &
  +1.11777933477806080d-14)*t -1.39959545848483840d-14)*t &
  +4.93110187870320640d-14)*t -2.02193307034590720d-13)*t &
  +7.53585452663569920d-13)*t -3.14632365928501299d-12)*t &
  +1.52351450024952975d-11)*t -8.75801572233507014d-11)*t &
  +6.27349413509555121d-10)*t -6.02183526555303242d-09)*t &
  +8.70043536788235270d-08)*t -2.32935044050984079d-06)*t &
  +1.83605337367638430d-04)*t -5.64003555099413391d-01
scbi = a/root4z
b = (((((((((((((((( -4.12972759036723200d-15*t &
  +8.36512465551360000d-15)*t -2.05945081774080000d-16)*t &
  +6.23733840790323200d-15)*t -5.81333983959859200d-14)*t &
  +1.52893566095288320d-13)*t -4.11064788026333184d-13)*t &
  +1.33820884559538637d-12)*t -4.74293914921785574d-12)*t &
  +1.84868021228605050d-11)*t -8.15686769476673166d-11)*t &
  +4.19373390376196942d-10)*t -2.61584084406303574d-09)*t &
  +2.10021454539364698d-08)*t -2.37847770210509358d-07)*t &
  +4.43114636962516363d-06)*t -1.83241371436579068d-04)*t &
  +3.89918976811026487d-02
scai = (b/zeta)/root4z
a = ((((((((((((( -4.58484390222233600d-15*t &
  +1.13969221615738880d-14)*t -1.43160328250060800d-14)*t &
  +5.04734978526300160d-14)*t -2.07055957015081472d-13)*t &
  +7.73043520694004480d-13)*t -3.23454581960357018d-12)*t &
  +1.57043540332660220d-11)*t -9.06023827679991573d-11)*t &
  +6.52303613917050367d-10)*t -6.30993998756281944d-09)*t &
  +9.23711460831703303d-08)*t -2.54030284953639173d-06)*t &
  +2.17448385781448409d-04)*t +5.64409671680379110d-01
scaip = a*root4z
b = (((((((((((((((( +4.19612197958451200d-15*t &
  -8.50454708509081600d-15)*t +2.31421341122560000d-16)*t &
  -6.39683104557465600d-15)*t +5.92509321833062400d-14)*t &
  -1.56008660983891968d-13)*t +4.20106807813331968d-13)*t &
  -1.36926896339755520d-12)*t +4.86000800286762854d-12)*t &
  -1.89780061819570625d-11)*t +8.39314701970122041d-11)*t &
  -4.32843814802265754d-10)*t +2.71124934991469715d-09)*t &
  -2.19026888712002973d-08)*t +2.50504395196083566d-07)*t &
  -4.75245434337472120d-06)*t +2.05252791097940732d-04)*t &
  -5.46414841607309762d-02
scbip = (b/zeta)*root4z
zeta = zeta + pib4
return
end
! ----------------------------------------------------------------------
subroutine psiasy(fj,fn,sr,si,psir,psii,nopen,nmax)
! subroutine to determine real and imaginary part of asymptotic wavefunction o
! derivative of these
!  asmptotically, in the case of inelastic scattering, the wavefunction is
!  exp[-i(kr-l pi/2)] - S exp[i(kr-l pi/2)]
!  whereas in the case of photodissociation,
!  exp[-i(kr-l pi/2)] S - exp[i(kr-l pi/2)]
!  this is equivalent to, in the case of inelastic scattering
!  - yl (1-Sr) + jl Si + i [-jl(1+Sr)+yl Si]
!  and, for photodissociation,
!   yl (1-Sr) + jl Si + i [-jl(1+Sr)-yl Si]
!  written by:  millard alexander
!  current revision date:  16-jun-1990
! ---------------------------------------------------------------------
!  variables in call list:
!    fj             contains (for wavefunction calculation) the normalized
!                   ricatti bessel function jl
!                   contains (for derivative calculation) the derivative with
!                   respect to r of the normalized ricatti bessel function jl
!    fn             contains (for wavefunction calculation) the normalized
!                   ricatti bessel function yl
!                   contains (for derivative calculation) the derivative with
!                   respect to r of the normalized ricatti bessel function yl
!    sr, si         matrices of order nmax x nmax which contain
!                   on input: real and imaginary parts of s-matrix
!                   on return: real and imaginary parts of asymptotic
!                   wavefunction
!    psir           on return contains nopen x nopen real part of
!                   asymptotic wavefunction (or derivative)
!    psii           on return contains nopen x nopen imag part of asymptotic
!                   wavefunction (or derivative)
!
!    nopen          number of open channels
!    nmax           row dimension of matrices
! ----------------------------------------------------------------------------
use mod_phot, only: photof
use mod_hiutil, only: daxpy_wrapper
use mod_hivector, only: matmov
use mod_hiblas, only: dscal
implicit none
real(8), intent(in) :: fj(nopen)
real(8), intent(in) :: fn(nopen)
real(8), intent(in) :: sr(nmax,nmax)
real(8), intent(inout) :: si(nmax,nmax)
real(8), intent(out) :: psir(nmax,nmax)
real(8), intent(out) :: psii(nmax,nmax)
integer, intent(in) :: nopen
integer, intent(in) :: nmax

real(8), parameter :: one=1.d0
real(8), parameter :: onemin=-1.d0
integer :: icol, irow
real(8) :: fac
!    unit           scratch vector
real(8) :: unit(nopen)
  !   put unit vector into array unit
  do icol = 1, nopen
    unit(icol) = one
  end do
  ! first we want to calculate real part of wavefunction at infinity
  ! that is   yl(kr) (Sr-1) + jl(kr) Si for scattering or
  !         - yl(kr) (Sr-1) + jl(kr) Si for photodissociation
  ! first move Sreal into psii
  call matmov (sr, psii, nopen, nopen, nmax, nmax)
  ! now subtract unit matrix
  call daxpy_wrapper (nopen, onemin, unit, 1, psii(1, 1), nmax + 1)
  ! now premultiply by diagonal matrix -yl(kr) for photodissociation or
  ! +yl(kr) for scattering
  do irow = 1, nopen
    fac=one*fn(irow)
    if (photof) fac=-fac
    call dscal(nopen, fac, psii(irow,1), nmax)
  end do
  ! now store simag in psir
  call matmov(si, psir, nopen, nopen, nmax, nmax)
  ! premultiply by diagonal matrix jl(kr)
  do irow = 1, nopen
    call dscal(nopen, fj(irow), psir(irow,1), nmax)
  end do
  ! now evaluate +/- yl(kr) (Sr-1) + jl(kR) Si, save in psir
  do icol = 1, nopen
    call daxpy_wrapper(nopen, one, psii(1, icol), 1, psir(1,icol), 1)
  end do
  ! psir now contains real part of asymptotic scattering wavefunction
  ! now compute imaginary part of asymptotic wavefunction
  ! that is - jl(kr) (1+Sr) + yl(kr) Si for scattering or
  !         - jl(kr) (1+Sr) - yl(kr) Si for photodissociation
  ! now move Sreal into psii
  call matmov (sr, psii, nopen, nopen, nmax, nmax)
  ! now add unit matrix
  call daxpy_wrapper (nopen, one, unit, 1, psii(1, 1), nmax + 1)
  ! now premultiply by diagonal matrix -jl(kr)
  do irow = 1, nopen
    fac=-fj(irow)
    call dscal(nopen, fac, psii(irow,1), nmax)
  end do
  ! replace real part of s matrix by real part of asymptotic wavefunction
  call matmov(psir,sr,nopen, nopen, nmax, nmax)
  ! premultiply Simag by diagonal matrix yl(kr) for scattering or by
  ! -yl(kr) for photodissociation
  do irow = 1, nopen
    fac=fn(irow)
    if (photof) fac=-fac
    call dscal(nopen, fac, si(irow,1), nmax)
  end do
  ! now evaluate - jl(kr) (1+Sr) +/- yl(kR) Si, save in psii
  do icol = 1, nopen
    call daxpy_wrapper(nopen, one, si(1, icol), 1, psii(1,icol), 1)
  end do
  ! replace imaginary part of s matrix by imaginary part of
  ! asymptotic wavefunction
  call matmov(psii,si,nopen, nopen, nmax, nmax)
return
end
! ----------------------------------------------------------------------
subroutine smatop (tmod, sr, si, scmat, lq, r, prec, &
                   nopen, nmax,kwrit,ipos)
!  subroutine to compute s-matrix and, alternatively, the asymptotic
!  scattering wavefunction
!  asmptotically, in the case of inelastic scattering, the wavefunction is
!  exp[-i(kr-l pi/2)] - s exp[i(kr-l pi/2)]
!  this is equivalent to, in the case of inelastic scattering
!  - yl (1-sr) + jl si + i [-jl(1+sr)+yl si]
!  written by:  millard alexander
!  current revision date (algorithm):  5-mar-1997
!  revised on 30-mar-2012 by q. ma for stream I/O of wfu files
!  revised on 21-dec-2020 by p. dagdigian to check to overflow/underflow
!    in ricatti-bessel funcstions and abort run
!  current revision: 21-dec-2020 by p. dagdigian
! ---------------------------------------------------------------------
!  variables in call list:
!    tmod:    on input:  tmod contains the log-derivative matrix at r
!             if closed channels are present, it is assumed that prior
!             to entry into smat, all the closed channel components in the
!             log-derivative matrix have been eliminated, so that tmod
!             contains the open-open block of the log-derivative matrix,
!             packed into the lower nopen x nopen submatrix
!             on return:  tmod contains the modulus squared of the t-matrix
!               if the logical variable flagsu, contained in modurle mod_surf is
!               .true., then the calculation is assumed to be that of a molecule
!               colliding with a surface, in which case tmod contains the
!               modulus squared of the s-matrix
!             if the flag photof=.true., then tmod contains the photodissociat
!             transition probabilities into each final state, with each row
!             corresponding to a different initial state
!             if wavefn .eq. true, then tmod contains real part of the
!               derivative of the asymptotic wavefunction
!    sr:      on return contains the real part of the s-matrix
!             if the flag wavefn=.true., then sr contains the real part
!             of the asymptotic wavefunction, defined above
!    si:      on return contains the imaginary part of the s-matrix
!             if the flag wavefn=.true., then si contains the imaginary part
!             of the asymptotic wavefunction, defined above
!    scmat:   scratch matrix
!             if wavefn .eq. true, then scmat contains imaginary part of the
!               derivative of the asymptotic wavefunction
!    sr:      on return contains the real part of the s-matrix
!    lq:      array of orbital angular momenta for the open channels
!    fj,fn,   scratch vectors dimensioned at least nopen in calling program
!    fpj,fpn,
!    pk
!    derj, dern
!    r:       interparticle separation in bohr
!    prec:    precision of single precision floating point number
!    nopen:   number of open channels, this must have been determined before
!             call to smatopen
!    nmax:    maximum row dimensions of matrices
!    kwrit:   if true, k matrix is printed out
!    ipos:    if true, 132 line printer

!  subroutines called:
!    vsmul:     scalar times a vector
!    cbesn,cbesj  ricatti-bessel functions (from b.r. johnson)
!    rles:      linear equation solver
!    rgmmul:    sub-matrix multiply
!    daxpy:     linpack blas routine
!    vmul:      vector times a vector
! --------------------------------------------------------------------
use mod_coqvec, only: nphoto, q
use mod_coeint, only: eint
use mod_cosysi, only: ispar
! temporary storage for smatrices
use mod_cotq1, only: srsave => dpsir ! srsave(100)
use mod_cotq2, only: sisave => tq2 ! sisave(100)
use mod_hibrid2, only: mxoutd, mxoutr
use mod_par, only: prsmat, jlpar! spac=>scat_spac
use mod_wave, only: irec, ifil, ipos2, ipos3, nrlogd, iendwv, ipos2_location
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_phot, only: photof, wavefn
use mod_surf, only: flagsu
use mod_hiutil, only: gennam
use mod_hiutil, only: daxpy_wrapper
use mod_hibrid1, only: cbesj, cbesn
use mod_himatrix, only: transp
use mod_hivector, only: matmov, vsmul, vmul
use mod_hiiolib1, only: openf
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
use mod_himatrix, only: mxma
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
use mod_himatrix, only: syminv
#endif
use mod_hiblas, only: dscal, dcopy, dgemm

implicit double precision (a-h,o-z)
real(8), dimension(nmax, nmax), intent(inout) :: tmod
real(8), dimension(nmax, nmax), intent(out) :: sr
real(8), dimension(nmax, nmax), intent(out) :: si
real(8), dimension(nmax, nmax), intent(out) :: scmat
integer, dimension(nopen), intent(in) :: lq
real(8), intent(in) :: r
real(8), intent(in) :: prec
integer, intent(in) :: nopen
integer, intent(in) :: nmax
logical, intent(in) :: kwrit
logical, intent(in) :: ipos

! ricatti-bessel functions
real(8), dimension(nopen) :: fj  
real(8), dimension(nopen) :: fn
real(8), dimension(nopen) :: fpj
real(8), dimension(nopen) :: fpn

real(8), dimension(nopen) :: pk  ! wavevectors for open channels
real(8), dimension(nopen) :: derj
real(8), dimension(nopen) :: dern

integer isw, i, icol, l
#if defined(HIB_UNIX_IBM)
character*1 forma
#endif
character*40 flxfil
!     The following three variables are used to determine the (machine
!     dependent) size of built-in types
integer int_t
double precision dble_t
character char_t
data isw / 0 /
integer, pointer :: ipol
ipol=>ispar(3)
one = 1.0d0
izero=0
onemin = -1.d0
twomin = -2.d0

! init arrays
 sr = 0d0 ; si = 0d0 ; scmat = 0d0
 
!  calculate asymptotic wavevectors of each channel
do  30   i = 1, nopen
  p2 = 2 * rmu * (ered - eint(i))
  pk(i) = sqrt(p2)
30 continue
!  calculate K-matrix
! if photodissocation calculation, save log-derivative matrix in si
if (photof) call matmov(tmod,si,nopen,nopen,nmax,nmax)
do  60   i = 1, nopen
  l = lq(i)
  p = pk(i)
  call cbesj (l, p, r, cj, cpj)
  call cbesn (l, p, r, cn, cpn)
!
!  check for overflow/underflow in ricatti-bessel functions
  if (isnan(cj)) then
    write(6,445)
    write(9,445)
445     format(/' *** OVERFLOW/UNDERFLOW IN RICATTI-BESSEL FUNCTIONS.', &
      'ABORT ***'/)
    stop
  end if
!
  fj(i) = cj
  fn(i) = cn
  fpj(i) = cpj
  fpn(i) = cpn
  if (cj .eq. 0. .or. cn .eq. 0) write (9 ,50) i
  if (abs(cj) .lt. (100./prec) .or. abs(cn) .gt. (prec/100.)) &
      then
    write (6, 50) i
    write (9, 50) i
  end if
50   format &
      (' *** R-END TOO SMALL FOR BESSEL FUNCTIONS IN CHANNEL', &
       i3, '; ABORT ***')
60 continue
! if wavefunction desired, then save in record 2 of direct access file
! 1. number of open channels,
! 2. total number of records on wavefunction
! 3. asymptotic interparticle separation
! 4. wavevectors of open channels
! 5. bessel functions and derivatives for open channels
if (wavefn) then
!
!     As of the 12.1.3 version of ifort, a big hole is likely to be
!     created here in the wfu file.  Most likely to be a bug of ifort.
!
!     If the bug is fixed in the future, use the following line and
!     remove ALL references to iendwv
!$$$         inquire (ifil, pos=ipos2)
!
   ipos2 = iendwv
   write (ifil, err=950, pos=ipos2_location) ipos2, ipos3, nrlogd
   write (ifil, err=950, pos=ipos2) irec, nopen, nphoto, &
        r, (pk(i), i=1, nopen), (fj(i), i=1, nopen), &
        (fpj(i), i=1, nopen), (fn(i), i=1, nopen), &
        (fpn(i), i=1, nopen)
   iendwv = iendwv + 3 * int(sizeof(int_t), kind(int_t)) + &
        (5 * nopen + 1) * int(sizeof(dble_t), kind(int_t))
endif
! save derivatives
call dcopy(nopen,fpj,1,derj,1)
call dcopy(nopen,fpn,1,dern,1)
do  70   icol = 1, nopen
  cj = fj(icol)
  cn = - fn(icol)
  call vsmul (tmod(1,icol), 1, cn, sr(1,icol), 1, nopen)
  sr(icol,icol) = sr(icol,icol) + fpn(icol)
  call vsmul (tmod(1,icol), 1, cj, tmod(1,icol), 1, nopen)
  tmod(icol,icol) = tmod(icol,icol) - fpj(icol)
70 continue
!  the above loops refer to eq.(18) of bob johnson's contribution to
!  nrcc workshop "algorithms and computer codes ..." vol. 1
!  sr is equal to - ( (y(xn)n(xn) - n'(xn) )
!  tmod is equal to y(xn) j(xn) - j'(xn)
!  solve linear equations for k-matrix (negative of r-matrix)
!  fpj is used as scratch vector here
call rles (sr, tmod, nopen, nopen, nmax)
!  since r-matrix is negative of k-matrix, r = -tmod at this point
isym=1
if (kwrit) then
  write (9, 76)
76   format(/,' ** K MATRIX')
  call mxoutd (9, tmod, nopen, nmax, isym, ipos)
endif
if (.not. photof) then
! here for scattering
!  calculate r**2, store in sr
#if defined(HIB_NONE)
call rgmmul (isw, nopen, nopen, nopen, tmod, 1, nmax, &
             tmod, 1, nmax, sr, 1, nmax)
#endif
#if defined(HIB_CRAY)
call mxma (tmod, 1, nmax, tmod, 1, nmax, sr, 1, nmax, &
            nopen, nopen, nopen)
#endif
#if defined(HIB_UNIX)
call dgemm('n','n',nopen,nopen,nopen,1.d0,tmod,nmax, &
            tmod,nmax,0.d0,sr,nmax)
#endif
!  determine imaginary part of s-matrix
!  also put unit vector into array fpn
do 80  icol = 1, nopen
  fpn(icol) = one
  call vsmul (tmod(1,icol), 1, twomin, si(1,icol), 1, &
              nopen)
80 continue
call daxpy_wrapper (nopen, one, fpn, 1, sr(1,1), nmax + 1)
!  solve linear equations for s-imaginary
!  see eq.(67) of r.g. gordon, meth. comp. phys. 10 (1971) 81
!  fpj is used as scratch vector here
call rles (sr, si, nopen, nopen, nmax)
!  determine real part of s-matrix
!  see eq.(68) of r.g. gordon, meth. comp. phys. 10 (1971) 81
#if defined(HIB_NONE)
call rgmmul (isw, nopen, nopen, nopen, tmod, 1, nmax, &
             si, 1, nmax, sr, 1, nmax)
#endif
#if defined(HIB_CRAY) || defined(HIB_UNIX) || defined(HIB_MAC)
call mxma (tmod, 1, nmax, si, 1, nmax, sr, 1, nmax, &
            nopen, nopen, nopen)
#endif
!  the matrix sr now contains s-real - del, i.e. the negative of the
!  t-matrix
!  now calculate t-squared
do 90  icol = 1, nopen
  call vmul (sr(:,icol), 1, sr(:,icol), 1,tmod(:,icol), 1, &
             nopen)
!  tmod now contains (s-real - del)**2
  call vmul (si(:,icol), 1, si(:,icol), 1, scmat(:,icol), 1, &
             nopen)
!  scmat now contains s-imag **2
!  now add s-imag**2 to (s-real - del)**2
  call daxpy_wrapper (nopen, one, scmat(1,icol), 1, tmod(1,icol), 1)
90 continue
!  now add 1 back to the diagonal elements of sr to get true s-real
call daxpy_wrapper (nopen, one, fpn, 1, sr(1,1), nmax + 1)
! TO RENDER NEXT LOOP OPERATIVE, CHANGE .EQ.10 TO .EQ.1
if (ibasty.eq.7 .and. jlpar.eq.-1 .and. ipol.eq.10) then
! special for singlet-triplet mixing:
! rotate transition probabilities to correspond to state assignment:
! psi-parallel = (Jtot+1)^1/2 |el=J+1> - Jtot^1/2 |el=J-1>
! psi-perpendicular = Jtot^1/2 |el=J+1> + (Jtot+1)^1/2 |el=J-1>
! with assignment that original column five is singlet with el=J-1 and
! column 6 is singlet with el=J+1
! use orthogonal plane rotation
    xjtot=lq(5)+1
    cs=sqrt(xjtot+1.d0)
    sn=sqrt(xjtot)
    xnorm=sqrt(2.d0*xjtot+1.d0)
    sn=sn/xnorm
    cs=cs/xnorm
    cs2=cs*cs
    sn2=sn*sn
    cssn=cs*sn
    call drot(4,sr(1,5),1,sr(1,6),1,cs,sn)
    call drot(4,si(1,5),1,si(1,6),1,cs,sn)
    call dcopy(4,sr(1,5),1,sr(5,1),nmax)
    call dcopy(4,sr(1,6),1,sr(6,1),nmax)
    call dcopy(4,si(1,5),1,si(5,1),nmax)
    call dcopy(4,si(1,6),1,si(6,1),nmax)
    do 95 i=1, 4
      tmod(5,i)=sr(5,i)*sr(5,i)+si(5,i)*si(5,i)
      tmod(6,i)=sr(6,i)*sr(6,i)+si(6,i)*si(6,i)
95     continue
    call dcopy(4,tmod(5,1),nmax,tmod(1,5),1)
    call dcopy(4,tmod(6,1),nmax,tmod(1,6),1)

endif

!  if molecule-surface collisions, then the diagonal elements of
!  tmod should be s-real**2 + s-imag**2
if (flagsu) then
  do 100  icol = 1, nopen
    tmod(icol,icol) = scmat(icol,icol) + sr(icol,icol) ** 2
100   continue
end if
if (.not. wavefn) return
! if wavefunction desired, then save
! real and imaginary part of s-matrix in record 2 of direct access file

do icol=1, nopen
   write (ifil, err=950) (sr(i, icol), i=1, nopen)
end do
do icol=1, nopen
   write (ifil, err=950) (si(i, icol), i=1, nopen)
end do
    if (kwrit .and. photof) then
      write (9,121)
121         format(/,'** REAL PART OF S MATRIX')
        call mxoutd (9, sr, nopen, nmax, isym, ipos)
      write (9,122)
122         format(/,'** IMAGINARY PART OF S MATRIX')
        call mxoutd (9, si, nopen, nmax, isym, ipos)
    endif
write (ifil, err=950) 'ENDWFUR', char(2)
iendwv = iendwv + 8 * sizeof(char_t) &
     + (2 * nopen ** 2) * sizeof(dble_t)
! save smatrix temporarily
    call matmov(sr,srsave,nopen,nopen,nmax,nmax)
    call matmov(si,sisave,nopen,nopen,nmax,nmax)
! here if wavefunction wanted
  call psiasy(fj,fn,sr,si,tmod,scmat,nopen,nmax)
! on return, sr contains real part of asymptotic wavefunction, si contains
! imaginary part of asymptotic wavefunction
! also determine real and imaginary part of derivative of
! asymptotic wavefunction
 call psiasy(derj,dern,srsave,sisave,tmod, &
             scmat,nopen,nmax)
endif
! here for photodissociation, at this point:
!              tmod contains K matrix
!              si contains log-derivative matrix
if (photof) then
! copy K matrix into sr
  call matmov(tmod,sr,nopen,nopen,nmax,nmax)
! invert K matrix
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
  call smxinv(sr,nmax,nopen,srsave,sisave,ierr)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
  call syminv(sr,nmax,nopen,ierr)
#endif
! save lhs for determination of real part of transition amplitudes
  do  300   icol = 1, nopen
    cj=fj(icol)
    do 290 irow=1, nopen
      sr(irow,icol)=sr(irow,icol)+tmod(irow,icol)
! sr now contains K + K^-1
      scmat(irow,icol)=si(irow,icol)*cj
290     continue
    scmat(icol,icol)=scmat(icol,icol)-derj(icol)
300   continue
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
    call mxma(scmat,1,nmax,sr,1,nmax,sisave,1,nmax, &
                nopen,nopen,nopen)
#endif
#if defined(HIB_UNIX_IBM)
  forma='N'
  call dgemul(scmat,nmax,forma,sr,nmax,forma, &
                sisave,nmax,nopen,nopen,nopen)
#endif
! solve linear equations for real part of transition amplitudes
! fpj is used as scratch here
  call rles(sisave,q,nopen,nphoto,nmax)
  call dscal(nopen*nphoto,onemin,q,1)
! q now contains real part of transition amplitudes
! store real part of transition amplitudes in sr
  call matmov(q,sr,nopen,nopen,nopen,nmax)
! now calculate imaginary part of transition amplitudes
  call dscal(nopen*nphoto,onemin,q,1)
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
    call mxma(tmod,1,nmax,q,1,nopen,srsave,1,nmax, &
                nopen,nopen,nphoto)
#endif
#if defined(HIB_UNIX_IBM)
  forma='N'
  call dgemul(tmod,nmax,forma,q,nopen,forma, &
                srsave,nmax,nopen,nopen,nphoto)
#endif
! store imaginary part of transition amplitudes in si
  call matmov(srsave,si,nopen,nopen,nmax,nmax)
! transpose transition amplitudes for output
  call transp(sr, nopen, nmax)
  call transp(si,nopen,nmax)
  if (kwrit) then
    write (9, 360)
360     format (/,'** REAL PART OF TRANSITION AMPLITUDES')
    isym=0
    call mxoutr(9, sr, nphoto, nopen, nmax, isym, ipos)
    write (9, 365)
365     format (/,'** IMAGINARY PART OF TRANSITION AMPLITUDES')
    call mxoutr(9, si, nphoto, nopen, nmax, isym, ipos)
  endif
! save transition amplitudes
  if (wavefn) then

!          call dbwi(nphoto,1,ifil,REC_LAST_USED)
     do jrow = 1, nphoto
        write (ifil, err=950) (sr(jrow, jcol), jcol=1, nopen)
     end do
! NB there were 2*nphoto!
     do jrow = 1, nphoto
        write (ifil, err=950) (si(jrow, jcol), jcol=1, nopen)
     end do

     write (ifil, err=950) 'ENDWFUR', char(2)
     iendwv = iendwv + 8 * sizeof(char_t) &
          + (2 * nopen ** 2) * sizeof(dble_t)
  endif
!  now determine transition probabilities by squaring
!  and normalize to unit total probability
  do 370 nst=1, nphoto
! use dern as scratch here
    dern(nst) = &
      ddot(nopen,sr(nst,1),nmax,sr(nst,1),nmax) &
      + ddot(nopen,si(nst,1),nmax, si(nst,1),nmax)
    dern(nst) = dern(nst)/rmu
370   continue
  write (9, 371)
371   format (/,'** SUM OF PHOTOFRAGMENT FLUXES (AU)')
  isym=0
  call mxoutr(9, dern, 1, nphoto, 1, isym, ipos)
  do 380 nst = 1, nphoto
    fac=one/rmu
    do 375 ncol = 1, nopen
      scmat(nst,ncol) = (sr(nst,ncol)*sr(nst,ncol)+ &
            si(nst,ncol)*si(nst,ncol))*fac
375       continue
380   continue
  t1=scmat(1,1)+scmat(1,2)
  t2=scmat(1,3)+scmat(1,4) &
           +scmat(1,5)+scmat(1,6)

! just for HCL problem
  if (ibasty .eq. 10) then
    call gennam(flxfil,'hclflux',1,'flx',lenft)
    call openf(3,flxfil,'sf',0)
    write (3, 381) ered*219474.6,dern(1),t1,t2,t2/t1,t1/dern(1)
    write (6, 381) ered*219474.6,dern(1),t1,t2,t2/t1,t1/dern(1)
    close (3)
381     format(f10.2,3g14.5,2f8.4)
  endif
!  scmat now contains transition probabilities into all final
!  states (columns) for each initial state (rows)
  write (9,395)
395   format (/,'** PHOTOFRAGMENT FLUXES (AU)', &
          ' (COLUMNS ARE FINAL STATES')
!        write (6, 396) ered*219474.6,dern(1),(scmat(1,i), i=1,nopen)
!** for cnne flux calc *** mby
!         call openf(4,'resenergy.dat','sf',0)
!         write (4,396) ered*219474.6d0,dern(1)
!         close(4)
  call mxoutr(9, scmat, nphoto, nopen, nmax, isym, ipos)
! diagnostic for two state problem; leave in for now 1/21/92
!        write (23,400) spac, scmat(1,1), scmat(1,2)
! 400   format (f19.11,2(f25.12))
  do 410 nst = 1, nphoto
    fac=one/(rmu*dern(nst))
    do 405 ncol = 1, nopen
      scmat(nst,ncol) = (sr(nst,ncol)*sr(nst,ncol)+ &
            si(nst,ncol)*si(nst,ncol))*fac
405       continue
410   continue
!  scmat now contains normalized photofragment fluxes into all final
!  states (columns) for each initial state (rows)
  write (9,420)
420   format (/,'** NORMALIZED PHOTOFRAGMENT FLUXES', &
          ' (COLUMNS ARE FINAL STATES')
  call mxoutr(9, scmat, nphoto, nopen, nmax, isym, ipos)
! determine real and imaginary parts of chi (save these in sr and si)
! determine real and imaginary parts of derivatives (save these in tmod
! and scmat
  if (wavefn.or.prsmat) then
! retranspose transition amplitudes
  call transp(sr, nopen, nmax)
  call transp(si,nopen,nmax)
    do  450 icol=1,nphoto
    do  450 irow=1,nopen
      srr=sr(irow,icol)
      sii=si(irow,icol)
      si(irow,icol)=-fn(irow)*sii+fj(irow)*srr
      sr(irow,icol)=-fn(irow)*srr-fj(irow)*sii
      scmat(irow,icol)=-fpn(irow)*sii+derj(irow)*srr
      tmod(irow,icol)=-fpn(irow)*srr-derj(irow)*sii
450     continue
  endif
endif
return
!
950 write (0, *) '*** ERROR WRITING WFU FILE (SMATOP). ABORT'
call exit()
end
! -----------------------------------------------------------------------
subroutine smatrx (z, sr, si, &
                   bqs, r, prec, ts, tsw, nopen, nch, nmax, &
                   kwrit,ipos)
! -----------------------------------------------------------------------
!  subroutine to:
!                1. eliminate closed-channel components in the log-
!                   derivative matrix
!                2. obtain the s-matrix
!                3. if wavefunction desired, store asymptotic wavefunction
!  current revision date (algorithm): 9-feb-1992
!   --------------------------------------------------------------------------
!  variables in call list
!    z:       on input:  z contains the log-derivative matrix at r
!             on return:  the upper left nopen x nopen block of z
!                         contains the modulus squared of the t-matrix
!    sr:      on return:  the upper-left nopen x nopen block of sr
!                         contains the real part of the s-matrix
!    si:      on return: the upper-left nopen x nopen block of si
!                        contains the imaginary part of the s-matrix
!    amat,    scratch matrices
!     bmat
!    bqs:     rotational angular momenta, orbital angular momenta, and
!             additional quantum index for each channel
!    isc1,sc1,
!    sc2,sc3,   scratch vectors of dimension at least equal to the number of
!    sc4,sc5:   channels
!    sc6, sc7
!
!    r:       interparticle separation
!    prec:    precision of single precision floating-point number
!    ts,tsw:  on return: contain cpu and wall clock time in seconds spent in
!             determination of s-matrix
!    nopen    on return:  number of energetically open channels
!    nch      on entry:  number of channels
!    nmax     on entry:  maximum row dimension of matrices
!    kwrit    if true, k matrix is printed out


!  ---------------------------------------------------------------------------
use mod_coqvec, only: nphoto, q
use mod_coeint, only: eint
use mod_hibrid2, only: mxoutd, mxoutr
use mod_wave, only: ifil, ipos2, ipos3, nrlogd, iendwv, ipos2_location
use mod_ered, only: ered
use mod_phot, only: photof, wavefn
use mod_hiutil, only: mtime
use mod_hitypes, only: bqs_type
use mod_hiblas, only: dcopy
implicit double precision (a-h,o-z)
real(8), intent(inout) :: z(nmax,nmax)
real(8), intent(out) :: sr(nmax,nmax)
real(8), intent(out) :: si(nmax,nmax)
type(bqs_type), intent(inout) :: bqs
real(8), intent(in) :: r
real(8), intent(in) :: prec
real(8), intent(out) :: ts
real(8), intent(out) :: tsw
integer, intent(out) :: nopen
integer, intent(in) :: nch
integer, intent(in) :: nmax
logical, intent(in) :: kwrit
logical, intent(in) :: ipos


real(8), allocatable :: amat(:,:)
integer :: isc1(nch)

data izero /0/
!     The following variables are used to determine the (machine
!     dependent) size of built-in types
double precision dble_t
character char_t
!  if kwrit (prlogd) = .true. and photodissociation calculation, print out
!  <psi|mu matrix at end of airprp
if (kwrit .and. photof) then
    write (9, 20)
20     format (/,' ** GAMMA2 VECTOR AFTER AIRPRP')
    isym=0
    call mxoutr(9, q, nch, nphoto, nch, isym, ipos)
  endif
!  first eliminate all closed channel components in log-derivative matrix
call mtime(t1,t2)
nopen = 0
do   50  i = 1, nch
  if (eint(i) .le. ered) then
!  here if this channel is open
    nopen = nopen + 1
    isc1(nopen) = i
    eint(nopen) = eint(i)
    bqs%jq(nopen) = bqs%jq(i)
    bqs%inq(nopen) = bqs%inq(i)
    bqs%lq(nopen) = bqs%lq(i)
  end if
50 continue
bqs%length = nopen
allocate(amat(nmax, nmax))
if (nopen .lt. nch) then
!  now pack the log-derivative matrix into a matrix of size nopen x nopen
!  keeping only the open-channel components
!  if photodissociation calculation, pack gamma2 also
  do  120  icol = 1, nopen
  ic = isc1(icol)
    do  100  irow = 1, nopen
      ir = isc1(irow)
      amat(irow,icol) = z(ir,ic)
100     continue
    call dcopy (nopen, amat(1, icol), 1, z(1, icol), 1)
120   continue
  do 140 icol =1, nphoto
    npoint=0
    do 130 irow = 1, nopen
      ir = isc1(irow)
      amat(irow,icol)=q(ir+npoint)
130     continue
    npoint=npoint+nch
140   continue
  npoint=1
  do  150 icol=1, nphoto
    call dcopy(nch, amat(1,icol), 1, q(npoint), 1)
    npoint=npoint+nopen
150   continue
  if (kwrit .and. photof) then
    write (9, 165)
165     format (/,' ** PACKED GAMMA2 VECTOR BEFORE SMATOP')
    isym=0
    call mxoutr(9, q, nopen, nphoto, nopen, isym, ipos)
  endif
endif
!  now determine s-matrix and modulus squared t-matrix
!  isc1, sc1, sc2, sc3, and sc4 are all used as scratch arrays
!  scmat is used as scratch matrix here
!  this uses new smat routine involving just open channels
call smatop (z, sr, si, amat, bqs%lq, r, prec, nopen, nmax, kwrit,ipos)
if (wavefn) then
! if wavefunction desired, then
! sr and si contain open channel portion of asymptotic wavefunction
! and z and amat contain derivative (real and imag) of asymptotic wfn
! now save channel packing list and
! real and imaginary part of wavefunction
! in record 3 of direct access file
!
!     Please refer to subroutine smatop for the info regarding the
!     following commented statement
!$$$         inquire (ifil, pos=ipos3)
   ipos3 = iendwv
   write (ifil, pos=ipos2_location) ipos2, ipos3, nrlogd
   write (ifil, err=950, pos=ipos3) (isc1(i), i=1, nopen)
   do icol = 1, nopen
      write (ifil, err=950) (sr(i, icol), i=1, nopen)
   end do
   do icol = 1, nopen
      write (ifil, err=950) (si(i, icol), i=1, nopen)
   end do
   do icol = 1, nopen
      write (ifil, err=950) (z(i, icol), i=1, nopen)
   end do
   do icol = 1, nopen
      write (ifil, err=950) (amat(i, icol), i=1, nopen)
   end do
   write (ifil, err=950) 'ENDWFUR', char(3)
   iendwv = iendwv + 8 * sizeof(char_t) &
        + (4 * nopen ** 2 + nopen) * sizeof(dble_t)
endif
deallocate(amat)
call mtime(t11,t22)
ts=t11 - t1
tsw=t22 -t2
return
!
950 write (0, *) '*** ERROR WRITING WFU FILE (SMATRX). ABORT.'
call exit()
end
!  ---------------------------------------------------------------------------
subroutine expand(ncol,nopen,nch,nmax,ipack,sr,si,bmat)
! expands the first ncol columns of sr and si from nopen*nopen
! to nch*nch inserting zeros for closed channel components
! NB if photodissocation calculation, only 1st column of sr and si
! are expanded to include zeros, but 1st column remains 1st column
! author:  millard alexander
! revision date:  30-dec-1995
! variables in call list
!  ---------------------------------------------------------------------------
! nopen:    number of open channels
! nch:      full number of channels
! ipack:    channel packing list; ipack(n) is the location of the
!           nth open channel in the full channel list
! nmax:     maximum row dimension of matrices
! sr:       on input:  real part of nopen x nopen s matrix
!           on return: real part of nch x nch s matrix
! si:       on input:  imaginary part of nopen x nopen s matrix
!           on return: imaginary part of nch x nch s matrix
! bmat      scratch matrix
!  ---------------------------------------------------------------------------
use mod_par, only: photof
use mod_hivector, only: dset, matmov
implicit double precision (a-h,o-z)
dimension sr(nopen,nopen),si(nopen,nopen), &
          bmat(nmax,nmax),ipack(15)
zero=0.d0
if (nch .eq. nopen) return
if (nch .lt. nopen) then
  write (6, 100) nch, nopen
100   format(/' **** NCH = ',i3,' .LT. NOPEN =',i3, &
         ' IN EXPAND; ABORT')
  call exit
endif
if (.not. photof) then
  do 200 icol=1, nch
    call dset(nch, zero, bmat(1,icol), 1)
200   continue
  do 230 icol=1, ncol
! ic is index of icolth open channel in the full channel list
    ic=ipack(icol)
    do 225 irow=1, nopen
      ir=ipack(irow)
      bmat(ir,ic)=sr(irow,icol)
225     continue
230   continue
  call matmov(bmat,sr,nch,nch,nmax,nmax)
  do 250 icol=1, ncol
! ic is index of icolth open channel in the full channel list
    ic=ipack(icol)
    do 245 irow=1, nopen
      ir=ipack(irow)
      bmat(ir,ic)=si(irow,icol)
245     continue
250   continue
  call matmov(bmat,si,nch,nch,nmax,nmax)
else
! here for photodissociation, expansion just applies to 1st column
  call dset(nch, zero, bmat, 1)
  do 325 irow=1, nopen
! ir is index of irowth open channel in the full channel list
    ir=ipack(irow)
    bmat(ir,1)=sr(irow,1)
325   continue
  call matmov(bmat,sr,nch,nch,nmax,nmax)
  do 345 irow=1, nopen
    ir=ipack(irow)
    bmat(ir,1)=si(irow,1)
345   continue
  call matmov(bmat,si,nch,nch,nmax,nmax)
endif
return
end
end module mod_hibrid3
