!************************************************************************
!                                                                       *
!                         hibridon 1  library                           *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!   2. airprp        airy zeroth-order propagator                       *
!   2a. gndloc       computes ground state wf and transform to local    *
!                    basis                                              *
!   3. airymp        returns the moduli and phases of the airy          *
!                    functions and their derivatives                    *
!   4. cbesj         centrifugal scattering function and its            *
!                    derivative.(first kind)                            *
!   5. cbesn         centrifugal scattering function and its            *
!                    derivative.(second kind)                           *
!   6. corr          determines approximate values for diagonal and     *
!                    off-diagonal correction terms in airy propagator   *
!   7. dtrans        computes  b * a * b-transpose                      *
!   8. difs (compar,compt2)         compares two s-matrices             *
!   NOTES:                                                              *
!      hypxsc moved to a separate file (hihypxsc.f)                     *
!      difcrs(ampli/sphn/plm) moved to a separate file (hidifcrs.f)     *
!                                                                       *
!************************************************************************
! NB cstart aix-unix uses fortran rather than essl routines
!************************************************************************
subroutine airprp (z, &
   xf, rend, drnow, en, &
   tolai, rincr, eshift, nch, nmax, itwo, iprint, twoen, noprin)
!  airy zeroth-order propagator from r=xf to r=rend
!  for reference see m. alexander, "hybrid quantum scattering algorithms ...",
!                    j. chem. phys. 81, 4510 (1984)
!                and m. alexander and d. manolopoulos, "a stable linear
!                    reference potential algorithm for solution ..."
!                    j. chem. phys. 86, 2044 (1987)
!**********************************************************************
!***   this integrator is not released for general public use      ****
!***   all use must be by specific prior arrangement with:         ****
!***     millard alexander, department of chemistry,               ****
!***     university of maryland, college park, md, 20742-2021      ****
!***     tel: 1.301.405.1823; email: mha@umd.edu                   ****
!***   no part of this program may be copied or used for other     ****
!***   purposes without the author's permission.                   ****
!**********************************************************************
!  author:  millard alexander
!  current revision date (the propagator): 30-dec-1995
!  revised on 30-mar-2012 by q. ma for stream I/O of wfu files
!  current revision: 8-oct-2012 by q. ma
! ----------------------------------------------------------------------------
!  definition of variables in call list:
!   z:               matrix of maximum dimension nmax*nmax
!                    on entry z contains the initial z-matrix at r=xf
!                    on return z contains the z-matrix at r=rend
!   w, tmat, vecnow
!    , vecnew:       scratch matrices of dimension at least nch*nch
!  eigold, eignew
!    , hp, y1, y2
!    , cc, y4, gam1,
!      gam2  :       scratch vectors of dimension at least nch
!  xf:               on entry: contains initial value of interparticle distanc
!                    on exit:  contains final value of interparticle distance
!                              this is equal to rend if normal termination
!                              otherwise an error message is printed
!  drnow:            on entry:  contains initial interval size
!                    on exit:  contains final interval size
!  en:               collision energy in atomic units
!  tolai:            parameter to determine step sizes
!                    if tolai .lt. 1, then estimated errors are used to
!                    determine next step sizes following the procedure outline
!                    in m.h. alexander, "hybrid quantum scattering algorithms
!                    if tolai .ge. 1, then step sizes are controlled by the
!                    algorithm:  drnext = tolai * drnow
!  rincr:            step size increase will occur only if rnext > rincr
!  variables in common block /cophot/
!     photof        true if photodissociation calculation
!                   false if scattering calculation
!     wavefn        true if G(a,b) transformation matrices are saved
!                   to be used later in computing the wavefunction
!  variable in common block /cowave/
!     irec          record number of last written G(a,b) matrix
!     ifil          local unit number for G(a,b) file
!  variables in common block /coered/
!    ered:      collision energy in atomic units (hartrees)
!    rmu:       collision reduced mass in atomic units (mass of electron = 1)
!  logical variables:
!     iprint:       if .true., then print out of step-by-step information
!     twoen:        if .true., then
!           itwo.eq.0, transformation matrices and relevant information
!             for propagation is written into file 10 for second energy
!           itwo.gt.0, transformation matrices and relevant information
!             for propagation is read from file 10 for second energy
!             calculation at energy=en+eshift
!     noprin:       if .true., then most printing is suppressed
! ----------------------------------------------------------------------------
use mod_coqvec, only: nphoto, q
use mod_cosc10, only: sc10
use mod_hibrid3, only: outmat, potent
implicit double precision (a-h, o-z)
!  matrix dimensions (row dimension = nmax, matrices stored column by column)
real(8), dimension(nmax*nmax), intent(inout) :: z

real(8), intent(inout) :: xf
real(8), intent(in) :: rend
real(8), intent(inout) :: drnow
real(8), intent(in) :: en
real(8), intent(in) :: tolai
real(8), intent(in) :: rincr
real(8), intent(in) :: eshift
integer, intent(in) :: nch
integer, intent(in) :: nmax
integer, intent(inout) :: itwo
logical, intent(in) :: iprint
logical, intent(in) :: twoen
logical, intent(in) :: noprin
integer i, icol, iend, ierr, ipt, izero, kstep, maxstp, &
        ncol, npt, nskip
logical photof, wavefn, boundf, wrsmat

common /cophot/ photof, wavefn, boundf, wrsmat
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
common /coered/ ered, rmu
common /coipar/ ipar(9),jprint
common /coselb/ ibasty
#if defined(HIB_UNIX_IBM)
character*1 forma, formb
#endif
data izero, ione, zero, one /0, 1, 0.d0, 1.d0/
!  powr is the power at which step sizes increase
!  only for tolai < 1
!  this used to be an input; now it is not
data powr /3.d0/
!     The following variables are for size-determination of (machine
!     dependent) built-in types
integer int_t
double precision dble_t
character char_t
real(8) :: w(nmax*nmax)
real(8), dimension(nmax*nmax) :: tmat
real(8), dimension(nmax*nmax) :: vecnow
real(8), dimension(nmax*nmax) :: vecnew
!  vectors dimensioned nch
real(8), dimension(nch) :: eigold
real(8), dimension(nch) :: eignow
real(8), dimension(nch) :: hp
real(8), dimension(nch) :: y1
real(8), dimension(nch) :: y2
real(8), dimension(nch) :: cc
real(8), dimension(nch) :: y4
real(8), dimension(nch) :: gam1
real(8), dimension(nch) :: gam2
! ----------------------------------------------------------------------------
if (.not.twoen) itwo = -1
if (itwo .gt. 0) go to 60
if (photof) then
  nqmat=nphoto*nch

  write (6, 20)
  write (9, 20)
20   format (' GAMMA2 INITIALIZED AT BEGINNING OF AIRPRP')
  call dset(nqmat,zero,q,1)
endif
spcmx = zero
spcmn = rend - xf
rmin = xf
!  determine local wavevectors at rmin to use in estimating second derivatives
!  hp and y1 are used as scratch vectors here
call wavevc (w, eigold, hp, y1, rmin, nch, nmax)
!  local wavevectors at rmin are returned in eigold
drfir = drnow
drmid = drnow * 0.5
rlast = xf
rold = xf
rnow = rlast + drmid
rnext = rlast + drnow
!  define local basis at rnow and carry out transformations
!  vecnew is used as scratch matrix and y1 is used as scratch vector here
call potent (w, vecnow, vecnew, eignow, hp, y1, &
             rnow, drnow, en, xlarge, nch, nmax)
!  vecnow is transformation from free basis into local basis
!  in first interval
!  e.g. p1=vecnow  ; see eq.(23) of
!  m.h. alexander, "hybrid quantum scattering algorithms ..."
!  store vecnow in tmat
call matmov (vecnow, tmat, nch, nch, nmax, nmax)
!  determine approximate values for diagonal and off-diagonal
!  correction terms
call corr (eignow, eigold, hp, drnow, drmid, xlarge, cdiag, &
           coff, nch)
maxstp = ( (rend-xf) / drnow ) * 5
xf = rend
if (iprint) then
  write (9, 40)
40   format(/' ** AIRY PROPAGATION (NO DERIVATIVES):')
  write (9, 50)
50   format('   STEP   RNOW', 5x, 5hdrnow, 5x, 5hcdiag, 6x, 4hcoff)
  if (jprint .ge. 2) write (9, 55)
55   format ('   ALSO ADIABATIC ENERGIES (HARTREE)')
end if
60 iend = 0
if (itwo .lt. 0) go to 70
!  write or read relevant information
call outmat (tmat, eigold, hp, eshift, drnow, rnow, &
             nch, nmax, itwo)
!  start airy propagation
! ----------------------------------------------------------------------------
70 do 200  kstep = 1, maxstp
!  transform log-deriv matrix from local basis in last interval to
!  local basis in present interval.  see eq.(23) of
!  m.h. alexander, "hybrid quantum scattering algorithms ..."
!  w is used as scratch matrix here, and y1 is scratch array
!  if photodissociation calculation, transform gamma2 into local
!  interval
  if (photof) then
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
    call mxma(tmat,1,nmax,q,1,nphoto,w,1,nmax,nch,nch,nphoto)
#endif
#if defined(HIB_UNIX_IBM)
    forma='N'
    call dgemul(tmat,nmax,forma,q,nmax,forma, &
                w,nmax,nch,nch,nphoto)
#endif
  call dcopy(nqmat,w,1,q,1)
  endif
!  determine ground state wavefunction and derivative and then
!  transform these into local basis (photodissociation)
!  y1 and y4 are used as scratch vectors here
if (photof) call gndloc(vecnow,w,rnow,drnow,nch,nmax)
!  transform logderivative matrix into current interval
call dtrans ( z, tmat, w, y1, xlarge, nch, nmax, izero)
!  tmat is no longer needed
!  solve for log-derivative matrix at right-hand side of
!  present interval.  this uses new algorithm of manalopoulos and alexander
!  namely
!               (n)    (n)      -1   (n)      (n)
!     z    = - y    [ y    + z ]    y     +  y
!      n+1      2      1      n      2        4
!  where y  , y  , and y   are the (diagonal) elements of log-derivative
!         1    2        4
!  propagator defined in alexander and manolopoulos
!  determine these diagonal matrices
!  eqs. (38)-(44) of m. alexander and d. manolopoulos, "a stable linear
!                    reference potential algorithm for solution ..."
  call spropn ( rnow, drnow, eigold, hp, y1, y4, y2, &
                   gam1, gam2, nch)
!  set up matrix to be inverted
!  nskip is spacing between diagonal elements of matrix stored column by colum
  nskip = nmax + 1
  call daxpy_wrapper (nch, one, y1, 1, z, nskip)
!  invert (y  +  z )
!           1     n
!  hp and cc are used as scratch arrays here
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
  call smxinv (z, nmax, nch, hp, cc, ierr)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
  call syminv(z,nmax,nch,ierr)
#endif
  if (ierr .ne. 0) then
    write (9, 80) kstep
    write (6, 80) kstep
80     format ('  *** INSTABILITY IN SMXINV IN AIRPRP, KSTEP=', i3, &
          ' ABORT ***')
    call exit
  end if
!  z now contains Z(0,rlast,rnext) in DEM's notation
!  if photodissociation calculation,
  if (photof) then
!  add gamma1(rlast,rnext) to gamma2(0,rlast) to form zeta(0,rlast,nrext)
    call vadd(1,q,1,gam1,1,nqmat)
!  then premultiply by Z(0,rlast,rnext)
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
    call mxma(z,1,nmax,q,1,nphoto,tmat,1,nmax,nch,nch,nphoto)
#endif
#if defined(HIB_UNIX_IBM)
    forma='N'
    call dgemul(z,nmax,forma,q,nmax,forma, &
                tmat,nmax,nch,nch,nphoto)
#endif
! if wavefunction desired, temporarily save local mu(a,b) propagator
! this is now in the first column of tmat
  if (wavefn) call dcopy(nch,tmat,1,sc10,1)
!  premultiply by y3 and add to existing q
!  cc is used as scratch array here
    ind=1
    jnd=1
    do 85  i = 1, nphoto
      call vmul(y2,1,tmat(jnd),1,cc,1,nch)
      call vadd(1,gam2(ind),1,cc,1,nch)
      call dcopy(nch,gam2(ind),1,q(ind),1)
      ind=ind+nch
      jnd=jnd+nmax
85     continue
! q now contains gamma2(0,rnext) in local basis
  endif

!                            -1
!  evaluate  - y  ( y  + z )    y
!               2    1    n      2
!  in the next loop evaluate the full, rather than lower triangle
!  changed for photodissociation calculation, old commands kept with
!  c in first column
npt = 1
do 90  i = 1, nch
  ncol = nch
!        ncol = nch - i + 1
  fact = y2(i)
  call dscal (ncol, fact, z(npt), 1)
  npt = npt + nmax
!        npt = npt + nskip
90 continue
!                            -1
!  z now contains  ( y  + z )    y  , this is G(n-1,n) in the local basis
!                     1    n      2
if (wavefn) then
! save this matrix as well as transformation matrix
! into local interval and local propagators
!
   irec = irec + 1
   write (ifil, err=950) -rlast, drnow
!     Adiabatic energies
   write (ifil, err=950) (eigold(i), i=1, nch)
!     The following information will not be written if wrsmat set to F
   if (wrsmat) then
      icol = 1
      do ich = 1, nch
         write (ifil, err=950) (z(icol - 1 + i), i=1, nch)
         icol = icol + nmax
      end do
      icol = 1
      do ich = 1, nch
         write (ifil, err=950) (vecnow(icol - 1 + i), i=1, nch)
         icol = icol + nmax
      end do
!
      write (ifil, err=950) (y1(i), i=1, nch), (y2(i), i=1, nch), &
           (y4(i), i=1, nch), (gam1(i), i=1, nch), &
           (sc10(i), i=1, nch)
      lrairy = (2 * nchwfu ** 2 + 6 * nchwfu + 2) &
           * sizeof(dble_t) + 8 * sizeof(char_t)
   else
      lrairy = (nchwfu + 2) * sizeof(dble_t) + 8 * sizeof(char_t)
   end if
!
   write (ifil, err=950) 'ENDWFUR', char(mod(irec, 256))
   iendwv = iendwv + lrairy
end if
!
do 110  i = 1, nch
  fact = - y2(i)
  call dscal (i, fact, z(i), nmax)
110 continue
!  add on  y
!           4
npt = 1
call daxpy_wrapper (nch, one, y4, 1, z, nskip)
!  now fill in upper half of z matrix
ipt = 2
if (nch .gt. 1) then
  do  120  icol = 1, nch - 1
!       ncol is number of subdiagonal elements in column icol
!       ipt points to the first subdiagonal element in column icol
!       (ipt + nmax - 1) points to the first superdiagonal element in row icol
    ncol = nch - icol
    call dcopy (ncol, z(ipt), 1, z(ipt + nmax -1), nmax)
    ipt = ipt + nskip
120   continue
end if
if (itwo .gt. 0) go to 160
!  obligatory write of step information if deviations from linear
!  potential are unusually large
!  this is only done if tolai .lt. 1, in which case the largest correction
!  is used to estimate the next step
if (tolai .lt. 1.d0) then
  cmax = max (cdiag, coff)
  if (cmax .gt. (5.d0* tolai)) then
    write (9,125)
    write (6,125)
125     format &
    (' ** ESTIMATED CORRECTIONS LARGER THAN 5*TOLAI IN AIRPRP')
    if (kstep .eq. 1) then
      write (9, 130)
      write (6, 130)
130       format ('    THE INITIAL VALUE OF DRNOW (SPAC*FSTFAC) IS', &
              ' PROBABLY TOO LARGE')
    else
      write (9, 140)
      write (6, 140)
140       format &
      ('   CHECK FOR DISCONTINUITIES OR UNPHYSICAL OSCILLATIONS', &
     /,'   IN YOUR POTENTIAL')
    end if
    if (.not. iprint) then
      write (9, 50)
      write (9,150) kstep, rnow, drnow, cdiag, coff
    end if
  end if
end if
!     write out information about step just completed
if (iprint) then
  write (9,150) kstep, rnow, drnow, cdiag, coff
150   format (i6, 4e10.3)
end if
!     get set for next step
160 if (iend .eq. 1) go to 250
if (itwo .gt. 0) go to 180
!  if tolai .lt. 1, predict next step size from largest correction
if (tolai .lt. 1.) then
!  note that the following statement is slightly different from  eq. (30)
!  of m.h. alexander, "hybrid quantum scattering algorithms ...  and that
!  the step-size algorithm is only approximately  related to any real
!  estimate of the error coff and cdiag should be approximately tolai, so
!  from eq. (27):
!  drnow(at n+1) = (12 tolai/kbar(n+1)w(n+1)-tilda')**(1/3)
!  which is approximately = (12 tolai/kbar(n)w(n)-tilda')**(1/3)
!                         = ((12 coff/kbar w-tilda') (tolai/coff))**(1/3)
!                         = drnow(at n) (tolai/coff)**(1/3)
!  or from eq. (29):
!                   drnow = drnow (tolai/cdiag)**(1/3)
!  then, using the larger error and allowing pow to vary:
! (note 1/23/92 powr fixed at 3 in data statement)
!  step size increase occurs only if rnext > rincr
    if (rnext .gt. rincr) &
            drnow = drnow * (tolai/cmax) ** (1.d0 / powr)
  else
!  if tolai .ge. 1, then
!  minimum step size is first interval width
    if (kstep .eq. 1) spcmn = drnow
!  and next step size is tolai * present step size
! only if rnext > powr
    if (rnext .gt. rincr) drnow = tolai * drnow
  end if
!  drnow is step size in next interval
rlast = rnext
rnext = rnext + drnow
if (rnext .lt. rend) go to 170
iend = 1
rnext = rend
drnow = rnext - rlast
170 rnew = rlast + 0.5d0 * drnow
if (kstep .gt. 1 .and. iend .ne. 1) then
  if (tolai .lt. 1) then
    if (drnow .lt. spcmn) spcmn = drnow
  end if
  if (drnow .gt. spcmx) spcmx = drnow
end if
drmid = rnew - rnow
!  restore eigenvalues
call dcopy (nch, eignow, 1, eigold, 1)
!  define local basis at rnew and carry out transformations
!  tmat is used as scratch matrix and y1 is used as scratch vector here
call potent (w, vecnew, tmat, eignow, hp, y1, &
             rnew, drnow, en, xlarge, nch, nmax)
!  determine matrix to transform log-deriv matrix into new interval
!  see eq. (22) of m.h. alexander, "hybrid quantum scattering algorithms ..."
!  on return from subroutine 'steppr':
!    the matrix pn [eq.(22)] is stored in
!    and vecnew has been transfered to vecnow; i.e. vecnow contains
!    the matrix tn [eq.(22)]
call steppr (vecnow, vecnew, tmat, nmax, nch)
!  restore radius values
rnow = rnew
!  determine approximate values for diagonal and off-diagonal
!  correction terms
call corr (eignow, eigold, hp, drnow, drmid, xlarge, cdiag, &
           coff, nch)
if (itwo .lt. 0) go to 200
if (iend .eq. 1) rnow = - rnow
!  write or read relevant information
180 call outmat (tmat, eigold, hp, eshift, drnow, rnow, &
             nch, nmax, itwo)
if (itwo .eq. 0) go to 200
!  negative rnow is cue for last step in second energy calculation
if (rnow .gt. 0.d0) go to 200
rnow = - rnow
iend = 1
!     go back to start new step
200 continue
!  the following statement is reached only if the integration has
!  not reached the asymptotic region in maxstp steps
write (9,210) maxstp, rnext
210 format (' *** AIRY PROPAGATION NOT FINISHED IN', i4, &
        ' STEPS:  R-FIN SET TO', f8.4,' ***',/)
xf = rnext
250 continue
if (itwo .lt. 0) go to 260
call outmat (vecnow, eigold, hp, eshift, drnow, xf, nch, &
             nmax, itwo)
!  transform log-deriv matrix into free basis.  transformation matrix is
!  just vecnow-transpose; see eq.(24) of m.h. alexander, "hybrid quantum
!  scattering algorithms ..."
260 call transp (vecnow, nch, nmax)
! if photodissociation calculation, also transform gamma2 to free basis
if (photof) then
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
    call mxma(vecnow,1,nmax,q,1,nphoto,w,1,nmax,nch,nch,nphoto)
#endif
#if defined(HIB_UNIX_IBM)
    forma='N'
    call dgemul(vecnow,nmax,forma,q,nmax,forma, &
                w,nmax,nch,nch,nphoto)
#endif
  ind=1
  jnd=1
  do 265 i=1,nphoto
    call dcopy(nch,w(jnd),1,q(ind),1)
    ind=ind+nch
    jnd=jnd+nmax
265   continue
endif
call dtrans (z, vecnow, w, hp, xlarge, nch, nmax, izero)
! if 2s-2p collisions, restore asymptotic case (e) energies
if (ibasty .eq. 10) call energ22
if (noprin) go to 320
if (itwo .lt. 0) write (9,280)
if (itwo .eq. 0) write (9,290)
if (itwo .gt. 0) write (9,300)
280 format (' ** AIRY PROPAGATION - FIRST ENERGY;', &
        ' TRANSFORMATION MATRICES NOT WRITTEN')
290 format (' ** AIRY PROPAGATION - FIRST ENERGY;', &
        ' TRANSFORMATION MATRICES WRITTEN')
300 format (' ** AIRY PROPAGATION - SECOND ENERGY;', &
        ' TRANSFORMATION MATRICES READ')
write (9,305) rmin, rend, tolai, kstep
write (9,310) spcmn, spcmx, rincr
305 format ('         RBEGIN =', f7.3, '  REND =', f7.3, &
        '     TOLAI =', 1pe8.1, '  NINTERVAL =', i3)
310 format ('         DR-MIN =', f7.3, '  DR-MAX =', f8.3, &
        '  R-INCR =', f7.3)
write (6, 315) rmin, rend, rincr, spcmn, spcmx, kstep
315 format (' ** AIRY:  RSTART =' ,f7.3,'  REND =',f7.3, &
        '   RINCR =',f7.3, &
        '   DRMIN =',f7.3, '   DRMAX =',f7.3,'   NSTEP =', i4)
320 continue
return
!
950 write (0, *) ' *** ERROR WRITING WFU FILE (AIRY). ABORT.'
call exit()
end
! ----------------------------------------------------------------------
subroutine gndloc (vecnow, scr, rnow, drnow, nch, nmax)
! ----------------------------------------------------------------------
!  this subroutine first determines the ground state wavefunction at
!    rnow +/- 0.5 drnow/sqrt(3), to estimate the function and its
!    derivative to use in airy propagation.  the chosen points are the
!    nodes of a two-point gauss-legendre quadrature
!  the function and its derivative are then transformed into the local basis
!  author:  millard alexander
!  current revision date: 20-sep-1994
! ---------------------------------------------------------------------
!  variables in call list:
!    vecnow:   contains Tn matrix (transformation into local basis)
!    scr:      scratch matrix
!    rnow:     midpoint of the current interval
!    drnow:    width of the current interval
!    nch:      number of channels
!    nmax:     maximum row dimension of matrices and maximum dimension of
!              vectors
! ----------------------------------------------------------------------
use mod_coqvec2, only: q => q2
use mod_cotq1, only: tmat => dpsir ! tmat(80)
implicit double precision (a-h,o-z)
!  square matrices (of row dimension nmax)
dimension vecnow(80), scr(80)
data one, onemn, half, sq3 /1.d0, -1.d0, 0.5d0, 1.732050807d0/
ra = rnow - half * drnow / sq3
rb = rnow + half * drnow / sq3
fact = sq3 / drnow
nphoto=1
mxphot=nch*nphoto
!  save vecnow in matrix tmat
call matmov(vecnow,tmat,nch,nch,nmax,nch)
call ground(scr, ra, nch, nphoto, mxphot)
call ground(q, rb, nch, nphoto, mxphot)
call dcopy(nch,q,1,q(nch+1),1)
call daxpy_wrapper(nch,one, scr,1,q,1)
call daxpy_wrapper(nch,onemn, scr,1,q(nch+1),1)
call dscal(nch,half,q,1)
call dscal(nch,fact,q(nch+1),1)
!  activate the next statement to kill derivatives of ground state
!  function
!      call dset(nch,0.d0, q(nch+1),1)
!  average ground state wavefunction now stored in first column of q
!  derivative stored in second column of q
!  now transform these into local basis
call mxma(vecnow,1,nmax,q,1,nch,scr,1,nch,nch,nch,2)
call dcopy(2*nch,scr,1,q,1)
return
end
! -----------------------------------------------------------------------
subroutine cbesn(l,p,x,cn,cnp)
!
!     centrifugal scattering function and its derivative.(second kind)
!     (defined j. comp. phys. vol 13,page 447)
!     this routine returns cn = ricatti-bessel function yl(p x) / sqrt(p)
!                          chp = d cn / dx = sqrt(p) d yl(p,x) / d(px)
!     or, in smp notation, if
!     u = p x
!     cn : yl_hat[u] / sqrt(p)
!     cnp : d[cn,x] : sqrt(p) d[yl_hat[u],u]
!     here we have used the definitions in abramowitz and stegun
!     this subroutine is unchanged from the original version written
!     by b.r. johnson
!     obtained from nrcc 1979
!
!     current revision date: 24-sept-87
!
implicit double precision (a-h,o-z)
ps = sqrt(p)
r = p*x
if(l-1) 1,2,3
1 cn = -cos(r)/ps
cnp = ps*sin(r)
return
2 cn = -(cos(r)/r+sin(r))/ps
cnp = -ps*cos(r)-cn/x
return
3 fnm1 = -cos(r)/r
fn = (fnm1 - sin(r))/r
tnp1dr = 1.0d0/r
tr = 2.0d0/r
do  100   i = 2,l
tnp1dr = tnp1dr + tr
fnp1 = tnp1dr*fn-fnm1
fnm1 = fn
100 fn = fnp1
cn = fn*ps*x
cnp = ps*(r*fnm1-dble(l)*fn)
return
end
! -----------------------------------------------------------------------
subroutine cbesj(l,p,x,cj,cjp)
!
!     centrifugal scattering function and its derivative.(first kind)
!     (defined j. comp. phys. vol 13,page 447)
!     this routine returns cj = ricatti-bessel function jl(p x) / sqrt(p)
!                          cjp = d cj / dx = sqrt(p) d jl(p,x) / d(px)
!     or, in smp notation, if
!     u = p x
!     cj : jl_hat[u] / sqrt(p)
!     cjp : d[cj,x] : sqrt(p) d[jl_hat[u],u]
!     here we have used the definitions in abramowitz and stegun
!     this subroutine is unchanged from the original version written
!     by b.r. johnson
!     obtained from nrcc 1979
!
!     current revision date: 24-sept-87
!
implicit double precision (a-h,o-z)
data zero,one,two /0.d0,1.d0,2.d0/
ps = sqrt(p)
r = p*x
if(l.ne.0)  go to 100
cj = sin(r)/ps
cjp = ps*cos(r)
return
100 if(r.ge.dble(l))  go to 500
!
!     millers method (recur down)
!
lim = l
test = 1.d8
tl1 = 2*lim + 1
tl1dr = tl1/r
tdr = two/r
fn = one
fnm1 = zero
do  401   i = 1,1000
fnp1 = tl1dr*fn-fnm1
tl1dr = tl1dr + tdr
fnm1 = fn
fn = fnp1
if(abs(fn).lt.test)  go to 401
max = i
go to 403
401 continue
write(9,404)
404 format(' *** no convergence in downward recurrence for cbesj;', &
       '  abort ***')
stop
403 n = lim + max
maxx = n - l
tnp = 2*n + 3
coef = tnp/r
fnp1 = zero
fn = one
do  405   i = 1,maxx
coef = coef - tdr
fnm1 = coef*fn - fnp1
fnp1 = fn
405 fn = fnm1
fnp1 = fnp1/fn
fn = one
ratio = fnp1
do  410   i = 1,l
coef = coef - tdr
fnm1 = coef*fn - fnp1
fnp1 = fn
410 fn = fnm1
cj = ps*x/((fn-r*fnp1)*cos(r) + r*fn*sin(r))
cjp = cj*(dble(l+1)/x - p*ratio)
return
!
!     recur up
!
500 if(l.gt.1)  go to 3
cj = (sin(r)/r-cos(r))/ps
cjp = ps*sin(r)-cj/x
return
3 fnm1 = sin(r)/r
fn = (fnm1-cos(r))/r
tnp1dr = one/r
tr = two/r
do  600   i = 2,l
tnp1dr = tnp1dr + tr
fnp1 = tnp1dr*fn-fnm1
fnm1 = fn
600 fn = fnp1
cj = fn*ps*x
cjp = ps*(r*fnm1-dble(l)*fn)
return
end
! -----------------------------------------------------------------------
subroutine corr (eignow, eigold, hp, drnow, drmid, xlarge, &
                 cdiag, coff, nch)
!  subroutine to determine approximate values for diagonal and off-diagonal
!  correction terms in airy propagator
!  also copies new eigenvalues from array eignow into array eigold
!  author:  millard alexander
!  current revision date: 27-sept-87
!
!  -------------------------------------------------------------------------
!  variables in call list:
!    eignow:    on entry: vector containing eigenvalues of wavevector matrix
!               in current interval
!    eigold:    on entry: vector containing eigenvalues of wavevector matrix
!               in previous interval
!               on return: vector containing eigenvalues of wavevector matrix
!               in current interval
!    hp:        vector containing diagonal elements of derivative of
!               transformed hamiltonian matrix in current interval
!               this is the same as the negative of the diagonal elements of
!               the wn-tilde-prime matrix
!    drnow:     width of current interval
!    drmid:     distance between mid-point of current interval and mid_point o
!               previous interval
!    xlarge:    largest off-diagonal elemtent in transformed wavevector matrix
!               in current interval
!    cdiag:     on return:  contains estimate of error due to neglected
!               diagonal elements of wn-tilde-double prime matrix
!               see eq.(29) of m.h. alexander, "hybrid quantum scattering
!               algorithms"
!    coff:     on return:  contains estimate of error due to neglected
!               off-diagonal elements of wn-tilde-prime matrix
!               see eq.(26) of m.h. alexander, "hybrid quantum scattering
!               algorithms"
!    nch:       number of channels
!  ----------------------------------------------------------------------
implicit double precision (a-h,o-z)
data zero,one,two /0.d0, 1.d0, 2.d0/
!      real  cay, cdiag, coff, drmid, drnow, factor, w2p, xlarge
!      real eignow, eigold, hp
!      real sqrt
integer i, nch
!  arrays, must be dimensioned at least nch
dimension eignow(1), eigold(1), hp(1)
factor = two / (drmid**2)
cay = zero
cdiag =zero
do 30  i = 1 , nch
!  ----------------------------------------------------------------------
!  estimate second derivative of wavevector by power series expansion
!                                                        2   2    2
!       w(r ) = w(r ) + (r  - r ) (dw/dr) + 0.5 (r  - r )  (d w/dr )
!          2       1      2    1         r        2    1            r
!                                         1                          1
!  which can be rearranged to give [since drmid = r - r  and hp = - (dw/dr) ]
!                                                  1   2
!         2    2                                                  2
!       (d w/dr ) = - 2 [ w(r ) - w(r ) + drmid * hp(r ) ] / drmid
!                r           1       2
!                 1
!  ----------------------------------------------------------------------
  w2p =  - factor * (eignow(i) - eigold(i) + drmid * hp(i))
  cdiag = cdiag + abs(w2p)
  cay = cay + sqrt (abs(eignow(i)))
30 continue
cay = cay / dble(nch)
cdiag = cdiag / dble(nch)
!  cay now contains average wavevector magnitude
!  cdiag now contains average magnitude of the second derivative of the
!    wavevector array
!  now calculate estimate of error
cdiag = (drnow**3) * cdiag / 12.d0
coff = cay * xlarge * (drnow**3) / 12.d0
!  now copy new eigenvalue array into eigold
call dcopy (nch, eignow, 1, eigold, 1)
return
end
! -----------------------------------------------------------------------
subroutine dtrans (a, b, c, diag, xlarge, n, nmax, ifind)
! -----------------------------------------------------------------------
!  to compute b * a * b-transpose
!  author:  millard alexander
!  current revision date:  28-dec-2003
! -----------------------------------------------------------------------
!  variables in call list:
!    a:       on return contains b * a * b-transpose
!    b:       on entry contains multiplicand matrix, this matrix is destroyed
!    c:       scratch matrix
!    diag:    if ifind .gt. 0, then on return the vector diag contains the
!             diagonal elements of b * a * b-transpose
!    xlarge:  if ifind .gt. 0, then on return the maximum (in magnitude)
!             off-diagonal element of b * a * b-transpose is returned as xlarg
!    n:       order of matrices and length of vector diag
!    nmax:    maximum row dimension of a, b, c in calling program
!    ifind:   integer variable, if ifind .gt. 0, then diag and xlarge are
!             computed as described above
!  all matrices are stored in packed column form
!  subroutines called:
!  rgmmul:      generalized matrix multiply ( a * b = c or a * b-transpose = c
!  dcopy:       linpack blas
!  maxmgv:      find maximum (absolute value) element in a vector
! -----------------------------------------------------------------------
implicit double precision (a-h,o-z)
#if defined(HIB_UNIX_IBM)
character*1 forma, formb
#endif
integer ic, icol, ifind, ipt, n, ncol, nmax, isw
dimension a(1), b(1), c(1), diag(1)
isw = 0
#if defined(HIB_NONE)
call rgmmul (isw, n, n, n, b, 1, nmax, a, 1, nmax, c, 1, nmax)
call rgmmul (isw, n, n, n, c, 1, nmax, b, nmax, 1, a, 1, nmax)
#endif
!#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
#if defined(HIB_NONE)
 call mxma (b,1,nmax,a,1,nmax,c,1,nmax,n,n,n)
 call mxma (c,1,nmax,b,nmax,1,a,1,nmax,n,n,n)
#endif
#if defined(HIB_UNIX_IBM)
 forma='N'
 formb='T'
 call dgemul (b,nmax,forma,a,nmax,forma,c,nmax,n,n,n)
 call dgemul (c,nmax,forma,b,nmax,formb,a,nmax,n,n,n)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_IBM)
 call dgemm('n','n',n,n,n,1.d0,b,nmax,a,nmax,0d0,c,nmax)
 call dgemm('n','t',n,n,n,1.d0,c,nmax,b,nmax,0d0,a,nmax)
#endif

!  a now contains desired product
if (ifind .le. 0) return
!  store diagonal elements in vector diag
!  nmax + 1 is the skip distrance for diagonal elements of the a-matrix
call dcopy (n, a(1), nmax + 1, diag(1), 1)
!  now look for maximum (in absolute value) off-diagonal element
!  since matrix is symmetrical, only the lower triangle is searched
xlarge = 0.d0
ipt = 2
do 15  icol = 1, n - 1
!  ipt points to first sub-diagonal element in column icol
!  ncol is the number of sub-diagonal elements in column icol
  ncol = n - icol
  call maxmgv (a(ipt), 1, zabs, ic, ncol)
  if (zabs .gt. xlarge) xlarge = zabs
  ipt = ipt + nmax + 1
15 continue
return
end
!  -------------------------------------------------------------
subroutine difs(fname1,fname2,ienerg,iprint,acc,accmx,thrs, &
  imx,jmx,ityp)
!  -------------------------------------------------------------
!   reads and compares two s-matrix files
!   fname1,fname: job names
!   ienerg: energy number
!   all s-matrices in first file are read
!   second file is searched for corresponding jtot
!   authors:  joachim werner and millard alexander
!   current revision date 19-mar-1996
!  ------------------------------------------------------------
use mod_codim, only: mmax
use mod_coamat, only: simag2 ! simag2(1)
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_cojhld, only: jout1 => jhold ! jout1(1)
use mod_coehld, only: jout2 => eholdint ! jout2(1)
use mod_coinhl, only: jlev => inhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
use mod_coisc2, only: ipack1 => isc2 ! ipack1(1)
use mod_coisc3, only: ipack2 => isc3 ! ipack2(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_cosc2, only: jpack1 => sc2int ! jpack1(1)
use mod_cosc3, only: jpack2 => sc3int ! jpack2(1)
use mod_cosc4, only: lpack1 => sc4int ! lpack1(1)
use mod_cosc5, only: lpack2 => sc5int ! lpack2(1)
use mod_coz, only: sreal1 => z_as_vec ! sreal1(1)
use mod_cow, only: sreal2 => w_as_vec ! sreal2(1)
use mod_cozmat, only: simag1 => zmat_as_vec ! simag1(1)
use mod_hibrid5, only: sread

implicit double precision (a-h,o-z)
character*(*) fname1,fname2
character*20 cdate1,cdate2
character*40 xnam1,xnam2
character*48 potnam1,potnam2,label1,label2

logical existf,csflg1,csflg2,flghf1,flghf2,flgsu1,flgsu2
logical twoml1,twoml2,nucr1,nucr2
#include "common/parpot.F90"
!
zero=0
acc=0
accmx=0
imx=0
jmx=0
ityp=0
call gennam(xnam1,fname1,ienerg,'smt',lenx1)
inquire (file=xnam1(1:lenx1), exist=existf)
if (.not. existf) then
  write (6, 20) xnam1(1:lenx1)
20   format(/' FILE ',(a),' NOT FOUND')
  return
end if
call gennam(xnam2,fname2,ienerg,'smt',lenx2)
inquire (file=xnam2(1:lenx2), exist=existf)
if (.not. existf) then
  write (6, 20) xnam2(1:lenx2)
  return
end if
call openf(1,xnam1(1:lenx1),'tu',0)
call openf(2,xnam2(1:lenx2),'tu',0)
call rdhead(1,cdate1, &
     ered1,rmu1,csflg1,flghf1,flgsu1, &
   twoml1,nucr1,jfirst,jfinal,jtotd,numin,numax,nud, &
   nlevel,nlevop,nnout1,jlev,inlev,elev,jout1)
label1=label
potnam1=potnam
call rdhead(2,cdate2, &
   ered2,rmu2,csflg2,flghf2,flgsu2, &
   twoml2,nucr2,jfirst,jfinal,jtotd,numin,numax,nud, &
   nlevel,nlevop,nnout2,jlev,inlev,elev,jout2)
label2=label
potnam2=potnam
if (thrs .gt. zero) then
  label = 'modulus of S matrix'
  length = 19
else
  label = 'S matrices'
  length = 10
end if
if(iprint.ge.1) then
  write(6,25) label(1:length), 'FIRST' ,xnam1(1:lenx1), &
              cdate1,label1,potnam1,'SECOND',xnam2(1:lenx2), &
              cdate2,label2,potnam2
25   format(/' COMPARE ',(a),/ &
    1x,(a),' FILE: ',(a)/ &
   ' WRITTEN:   ',(a)/ &
   ' LABEL:     ',(a)/ &
   ' POT NAME:  ',(a)// &
    1x,(a),' FILE: ',(a)/ &
   ' WRITTEN:   ',(a)/ &
   ' LABEL:     ',(a)/ &
   ' POT NAME:  ',(a))
end if
idif=0
if (rmu1.ne.rmu2) idif=idif+1
if (ered1.ne.ered2) idif=idif+1
if (nnout1.ne.nnout2) idif=idif+1
if (csflg1.neqv.csflg2) idif=idif+1
if (flghf1.neqv.flghf2) idif=idif+1
if (flgsu1.neqv.flgsu2) idif=idif+1
if (twoml1.neqv.twoml2) idif=idif+1
do 26 i=1,iabs(nnout1)
26 if (jout1(i).ne.jout2(i)) idif=idif+1
!
30 nopen1 = 0
call sread (0,sreal1, simag1, jtot1, jlpar1, nu1, &
                  jq, lq, inq, ipack1, jpack1, lpack1, &
                  1, mmax, nopen1, lengt1, ierr)
if(ierr.eq.-1) goto 200
if(ierr.lt.-1) then
  write(6,35) xnam1
35   format(' ERROR READING FILE ',(a))
  goto 200
end if
nopen2 = 0
call sread (0,sreal2, simag2, jtot2, jlpar2, nu2, &
                  jq, lq, inq, ipack2, jpack2, lpack2, &
                  2, mmax, nopen2, lengt2, ierr)
if(ierr.eq.-1) goto 200
if(ierr.lt.-1) then
  write(6,35) xnam2
  goto 200
end if
if(nu1.ne.nu2) idif=idif+1
if(jtot1.ne.jtot2) idif=idif+1
if(jlpar1.ne.jlpar2) idif=idif+1
if(lengt1.ne.lengt2) idif=idif+1
do 60 i=1,lengt1
if(jpack1(i).ne.jpack2(i)) idif=idif+1
if(lpack1(i).ne.lpack2(i)) idif=idif+1
60 if(ipack1(i).ne.ipack2(i)) idif=idif+1
if(idif.ne.0) then
  write(6,70) jtot1,jtot2
70   format(/' PARAMETERS NOT EQUAL FOR JTOT1=',i3,'  JTOT2=',i3)
  goto 200
end if
ncol=lengt1
if(nnout1.lt.0) then
  if(nopen1.ne.nopen2) then
    write(6,80) nopen1,nopen2
80     format(/' NOPEN1.NE.NOPEN2 FOR NNOUT.LT.0:', &
     '  NOPEN1=',i3,'  NOPEN2=',i3,'  NNOUT=',i2)
    call closf(1)
    call closf(2)
    return
  end if
  ncol=nopen1
end if
if(nnout1.gt.0) nopen1=0
if (thrs .lt. zero) then
!  here for comparison of s matrices
  if(iprint.ge.1) write(6,90) 'real',jtot1,jlpar1,nu1
90   format(/' comparing ',a,' part of S matrices for jtot=',i3, &
       '  jlpar=',i2,'  NU=',i3)
  call compar(sreal1,sreal2,mmax,lengt1,nopen1, &
    erabs,errel,ermabs,ermrel,thrs,n,iprint,im,jm)
  acc=max(acc,abs(errel))
  if(abs(ermrel).gt.accmx) then
     accmx=abs(ermrel)
     imx=im
     jmx=jm
     ityp=1
  end if
  ij=(jm-1)*mmax+im
  if(iprint.ge.1) then
    write(6,100) erabs,ermabs,errel,ermrel,im,jm, &
                 sreal1(ij),im,jm,sreal2(ij)
    write(6,90) 'IMAGINARY',jtot1,jlpar1,nu1
  end if
  call compar(simag1,simag2,mmax,lengt1,nopen1, &
    erabs,errel,ermabs,ermrel,thrs,n,iprint,im,jm)
  acc=max(acc,abs(errel))
  if(abs(ermrel).gt.accmx) then
     accmx=abs(ermrel)
     imx=im
     jmx=jm
     ityp=2
  end if
  ij=(jm-1)*mmax+im
  if(iprint.gt.1) write(6,95)
95   format()
  if(iprint.ge.1) &
    write(6,100) erabs,ermabs,errel,ermrel,im,jm, &
                 simag1(ij),im,jm,simag2(ij), thrs
100   format(' average absolute difference: ',f11.8/ &
         ' largest absolute difference: ',f11.8/ &
         ' average relative difference: ',f10.2,'%'/ &
         ' largest relative difference: ',f10.2,'%'/ &
  ' S1(',i2,',',i2,') =',g12.4,'  S2(',i2,',',i2,') =',g12.4,/ &
  ' inspection threshold is ',1pg8.1)
else
!  here for comparison of modulus of s-matrices
  if(iprint.ge.1) write(6,105) jtot1,jlpar1,nu1
105   format(/' moduli of S matrices for jtot=',i3, &
       '  jlpar=',i2,'  NU=',i3)
  call compt2(sreal1,simag1,sreal2,simag2,mmax,lengt1,nopen1, &
       erabs,errel,ermabs,ermrel,thrs,n,iprint,im,jm)
  acc=max(acc,abs(errel))
  if(abs(ermrel).gt.accmx) then
     accmx=abs(ermrel)
     imx=im
     jmx=jm
  end if
  ij=(jm-1)*mmax+im
  if (iprint .ge. 1) &
    write(6,190) erabs,ermabs,errel, &
    ermrel,im,jm,sqrt(simag1(ij)**2 + sreal1(ij)**2),im,jm, &
    sqrt(simag2(ij)**2 + sreal2(ij)**2),thrs
190  format(/ &
  ' average absolute difference: ',f11.8/ &
  ' largest absolute difference: ',f11.8/ &
  ' average relative difference: ',f10.2,'%'/ &
  ' largest relative differencE: ',f10.2,'%'/ &
  ' mod(S1)(',i2,',',i2,') =',g12.4, &
  ' mod(S2)(',i2,',',i2,') =',g12.4,/ &
  ' inspection threshold is ',1pg8.1)

end if
goto 30
200 call closf(1)
call closf(2)
return
end
subroutine compar(s1,s2,ndim,l1,n1,erabs,errel,ermabs,ermrel, &
  thrs,n,iprint,im,jm)
implicit double precision (a-h,o-z)
!   author: h.-j. werner
!   current revision date 23-sept-87
dimension s1(ndim,1),s2(ndim,1)
erabs=0
errel=0
ermabs=0
ermrel=0
ncol=l1
if(n1.ne.0) ncol=n1
n=0
m=0
do 100 i=1,ncol
je=i
if(n1.ne.0) je=l1
jee=min0(6,je)
if(iprint.gt.1) then
  write(6,10)
  write(6,10) (s1(i,j),j=1,jee)
  write(6,10) (s2(i,j),j=1,jee)
10   format(1x,6g12.4)
end if
do 100 j=1,je
n=n+1
er=abs(s1(i,j)-s2(i,j))
erabs=erabs+er
if(er.gt.ermabs) ermabs=er
if(abs(s1(i,j)).lt.thrs) goto 100
m=m+1
er=2.0d0*er/(abs(s1(i,j))+abs(s2(i,j)))
if(er.gt.ermrel) then
  ermrel=er
  im=i
  jm=j
end if
errel=errel+er
100 continue
n=max0(n,1)
m=max0(m,1)
erabs=erabs/n
errel=100.d0*errel/m
ermrel=100.d0*ermrel
return
end
! -----------------------------------------------------------------------
subroutine compt2(s1r,s1i,s2r,s2i,ndim,l1,n1,erabs,errel, &
  ermabs,ermrel,thrs,n,iprint,im,jm)
implicit double precision (a-h,o-z)
dimension s1r(ndim,1),s2r(ndim,1)
dimension s1i(ndim,1),s2i(ndim,1)
erabs=0
errel=0
ermabs=0
ermrel=0
ncol=l1
if(n1.ne.0) ncol=n1
n=0
m=0
do 100 i=1,ncol
je=i
if(n1.ne.0) je=l1
jee=min0(6,je)
if(iprint.gt.1) then
  write(6,10)
  write(6,10) (sqrt(s1r(i,j)**2+ s1i(i,j)**2),j=1,jee)
  write(6,10) (sqrt(s2r(i,j)**2+ s2i(i,j)**2),j=1,jee)
10   format(1x,6g12.4)
end if
do 100 j=1,je
n=n+1
s1mod = sqrt((s1r(i,j))**2 + s1i(i,j)**2)
s2mod = sqrt((s2r(i,j))**2 + s2i(i,j)**2)
if(s1mod.gt.thrs .and. s2mod.gt.thrs) then
  er=abs(s1mod - s2mod)
  erabs=erabs+er
  if(er.gt.ermabs) ermabs=er
  m=m+1
  er=2.0d0*er/(s1mod + s2mod)
  if(er.gt.ermrel) then
    ermrel=er
    im=i
    jm=j
  end if
  errel=errel+er
end if
100 continue
n=max0(n,1)
m=max0(m,1)
erabs=erabs/n
errel=100.d0*errel/m
ermrel=100.d0*ermrel
return
end
!  -------------------------------------------------------------
