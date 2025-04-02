#include "assert.h"
#include "unused.h"
module mod_hibrid1
use mod_assert, only: fassert
contains
!************************************************************************
!                                                                       *
!                         hibridon 1  library                           *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!   1. outmat        writes or reads transformation matrix              *
!   3. potmat        determines w(r) matrix                             *
!   2. potent        sets up wavevector matrix, derivativ etc.          *
!   3. spropn        this subroutine calculates the diagonal matrices   *
!                    to propagate the log-derivative matrix through     *
!                    the current interval                               *
!   4. steppr        determines matrix to transform log-deriv matrix    *
!                    into new interval                                  *
!   7. wavevc        sets up wavevector matrix and diagonalizes it      *
!   2. airprp        airy zeroth-order propagator                       *
!   2a. gndloc       computes ground state wf and transform to local    *
!                    basis                                              *
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

subroutine outmat (tmat, eigold, hp, eshift, drnow, rnow, &
                   n, nmax, itwo)
!  subroutine to either write or read transformation matrix and
!  relevant information from file FUNIT_TRANS_MAT
!  called from spropn
!  author:  millard alexander
!  current revision date: 14-feb-91
!  -------------------------------------------------------------------------
!  variables in call list:
!    tmat:     n x n matrix to contain transformation matrix
!    eigold:   array of dimension n which contains local wavevectors
!    hp:       array of dimension n which contains derivatives of hamiltonian
!              matrix.  this is just the negative of the derivatives of the
!              wavevector matrix
!    eshift:   amount local wavevectors will be shifted in second energy
!              calculation:  2         2
!                           k (new) = k (old) + eshift
!    drnow:    width of current interval
!    rnow:     midpoint of current interval
!    n:        number of channels
!    nmax:     maximum row dimension of matrix tmat
!    itwo:     if = 0, then subroutine called at first energy of multiple
!              energy calculation, so transformation matrix and relevant
!              information will be written
!              if > 0, then subroutine called at subsequent energy of multiple
!              energy calculation, so transformation matrix and relevant
!              information will be read
!  ------------------------------------------------------------------------
use funit, only: FUNIT_TRANS_MAT
implicit none
real(8), intent(inout) :: tmat(nmax*n)
real(8), intent(inout) :: eigold(n)
real(8), intent(inout) :: hp(:)
real(8), intent(in) :: eshift
real(8), intent(inout) :: drnow
real(8), intent(inout) :: rnow
integer, intent(in) :: n
integer, intent(in) :: nmax
integer, intent(in) :: itwo

integer :: i
logical :: isecnd
integer :: nsq 

isecnd = .false.
if (itwo .gt. 0) isecnd = .true.
!  if first energy calculation, isecnd = .false.
!    in which case logical unit FUNIT_TRANS_MAT will be written
!  if subsequent energy calculation, isecnd = .true.
!    in which case logical unit FUNIT_TRANS_MAT will be written
!  read/write rnow, drnow, diagonal elements of transformed dw/dr matrix,
!  and diagonal elements of transformed w matrix
nsq = n * nmax
if (isecnd) then
  read (FUNIT_TRANS_MAT) rnow, drnow, (hp(i) , i = 1, n), &
        (eigold(i) , i = 1, n), (tmat(i), i=1, nsq)
else
  write (FUNIT_TRANS_MAT) rnow, drnow, (hp(i) , i = 1, n), &
        (eigold(i) , i = 1, n), (tmat(i), i=1, nsq)
endif
!  now shift energies (if subsequent energy)
if (isecnd) then
  do  i = 1, n
    eigold(i) = eigold(i) + eshift
  end do
end if
return
end

!#define ANCOUMA_READ_USE_ASSOCIATE

#define ANCOUMA_READ_METHOD_NORMAL 0
#define ANCOUMA_READ_METHOD_INLINE_LEVEL1 1
#define ANCOUMA_READ_METHOD_INLINE_LEVEL2 2
#define ANCOUMA_READ_METHOD_INLINE_LEVEL3 3

#define ANCOUMA_READ_METHOD ANCOUMA_READ_METHOD_INLINE_LEVEL3


#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL1)

#define ANCOUMA_GET_ELEMENT(ancouma, element_index, ij, vee) \
ij = ancouma%v2i%get_element(iv2_element) ; \
vee = ancouma%v2d%get_element(iv2_element)

#endif


! ------------------------------------------------------------------------
subroutine potmat (w, r, nch, nmax, v2)
!  determine negative of lower triangle of w(r) matrix.  see eq. (3) of
!  m.h. alexander, "hybrid quantum scattering algorithms ..."
!  author:  millard alexander
!  latest revision date:  24-apr-1997 by mha
!  -------------------------------------------------------------------
!  variables in call list:
!    w:        matrix of maximum row dimension nmax
!              on output contains negative of lower triangle of w(r)
!    r:        value of interparticle distance at which -w(r) is to
!              be evaluated
!    nch:      actual size of matrix w
!    nmax:     maximum row and column dimension of matrix w
!   lamnum:     number of non-zero v2 matrix elements for each lambda

!  variables in common block /copmat/
!    rtmn,rtmx: minimum and maximum turning points (not used here)
!    iflag:     variable used in determination of turning points (not used her
!           iflag = 0 if all channels are in classically forbidden region
!           iflag = 1 if some channels are open
!           iflag = 2 if all asymptotically open channels are open at r
!  subroutines called:
!    pot:      returns r-dependence of each angular term in the potential
!    vsmul:    multiplies vector by scalar and stores result in another
!              vector
!  -------------------------------------------------------------------
use mod_ancou, only: ancou_type, ancouma_type
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_covvl, only: vvl
use mod_grovec, only: igrovec_type_block, dgrovec_type_block
use mod_hiba10_22p, only: trans22
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_pmat, only: rtmn, rtmx, iflag
use mod_cputim, only: cpupot
use mod_hivector, only: dset
use constants, only: zero, two
use mod_hiblas, only: dscal, dcopy
use mod_hipot, only: pot
use mod_hiutil, only: get_cpu_time
implicit none
real(8), dimension(*), intent(out) :: w
real(8), intent(in) :: r
integer, intent(in) :: nch
integer, intent(in) :: nmax
type(ancou_type), intent(in) :: v2

integer i, ilam
integer :: iv2_element
integer :: ij
integer :: num_nz_elements
real(8) :: vee

type(ancouma_type), pointer :: ancouma

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL2) || (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL3)
integer :: block_size
integer :: block_index
integer :: el_index_in_block
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL3)
integer :: num_full_blocks
integer :: num_remaining_nz_elements
type(igrovec_type_block), pointer :: blocki
type(dgrovec_type_block), pointer :: blockd
#endif


integer :: icol, icolpt, ioff, ipt, irowpt, iwpt, ncol, nmaxp1
real(8) :: r2, twormu, vv0, wmax, wmin


!ABER only for testing potential-matrix
!     if (r.le.3.41) icount=0
!ABER
!  calculate coefficients of each angular term
cpupot = cpupot - get_cpu_time()
call pot( vv0, r)

!  vv0 is the isotropic term in the potential
!  the coefficients for each angular term in the coupling potential
!  [ vvl(i) for i = 1, nlam ] are returned in common block vvl contains
!  multiply all vvl terms by twice the reduced mass
twormu = two * rmu
call dscal(v2%nlam, twormu, vvl, 1)
!  now loop over angular coupling matrix to calculate the potential matrix
!    w(ij) = 2 * rmu * vv0 + sum [ 2 * rmu * vvl(ilam) * v2(ij,ilam) ]
!  first zero out lower triangle of potential matrix
#if defined(ISSUE49_IS_FIXED)
iwpt = 1
do 20 icol = 1, nch
  ncol = nch - icol + 1
  call dset(ncol, zero, w(iwpt), 1)
  iwpt = iwpt + nmax + 1
20    continue
#else
! although w is a lower triangular matrix, initialize both upper and
! lower part of it as some code is still treating it as full and therefore
! reading uninitialized values
iwpt = 1
ncol = nch
do icol = 1, nch
  call dset(ncol, zero, w(iwpt), 1)
  iwpt = iwpt + nmax
end do
#endif
ioff = 0

!ABER
!       call druckq(w,nmax,nch,'potential matrix',icount)
!ABER
if (v2%nlam .gt. 0) then
  do ilam = 1, v2%nlam
      ! write(6,*) 'ilam=', ilam, 'v2%get_angular_coupling_matrix(ilam)%get_num_nonzero_elements()=', v2%get_angular_coupling_matrix(ilam)%get_num_nonzero_elements()
      !ancouma => v2%ancouma(ilam)
#if defined(ANCOUMA_READ_USE_ASSOCIATE)
      associate( ancouma => v2%get_angular_coupling_matrix(ilam) )
#else
      ancouma => v2%get_angular_coupling_matrix(ilam)
#endif
      num_nz_elements = ancouma%get_num_nonzero_elements()
      ASSERT(num_nz_elements >= 0)

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_NORMAL )
      do iv2_element = 1, num_nz_elements
        call ancouma%get_element(iv2_element, ij, vee)  
        w(ij) = w(ij) + vee * vvl(ilam)
      end do
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL1)
      do iv2_element = 1, num_nz_elements
        ANCOUMA_GET_ELEMENT(ancouma, iv2_element, ij, vee)
        w(ij) = w(ij) + vee * vvl(ilam)
      end do
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL2)
      block_size = ancouma%v2i%block_size
      do iv2_element = 1, num_nz_elements
        block_index =  (iv2_element-1) / block_size
        el_index_in_block = iv2_element - 1 - (block_index * block_size)
        ij = ancouma%v2i%blocks(block_index)%p(el_index_in_block)
        vee = ancouma%v2d%blocks(block_index)%p(el_index_in_block)

        w(ij) = w(ij) + vee * vvl(ilam)
      end do
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL3)
      !1023 -> 0
      !1024 -> 1
      !1025 -> 1 
      block_size = ancouma%v2i%block_size
      num_full_blocks = ((num_nz_elements - 1) / block_size)
      iv2_element = 0
      do block_index = 0, num_full_blocks-1
        blocki => ancouma%v2i%blocks(block_index)
        blockd => ancouma%v2d%blocks(block_index)
        do el_index_in_block = 0, block_size-1
          ij = blocki%p(el_index_in_block)
          vee = blockd%p(el_index_in_block)
          iv2_element = iv2_element + 1
          w(ij) = w(ij) + vee * vvl(ilam)
        end do
      end do
      num_remaining_nz_elements = num_nz_elements - iv2_element
      if (num_remaining_nz_elements > 0) then
        blocki => ancouma%v2i%blocks(num_full_blocks)
        blockd => ancouma%v2d%blocks(num_full_blocks)
        do el_index_in_block = 0, num_remaining_nz_elements-1
          ij = blocki%p(el_index_in_block)
          vee = blockd%p(el_index_in_block)
          iv2_element = iv2_element + 1
          w(ij) = w(ij) + vee * vvl(ilam)
        end do
      end if
#endif

#if defined(ANCOUMA_READ_USE_ASSOCIATE)
      end associate
#endif

  end do
endif
!ABER
!       call druckq(w,nmax,nch,'potential matrix',icount)
!ABER
!  now add on isotropic term plus centrifugal barrier and subtract wavevectors
!  from diagonal terms
r2 = 1.d0 / ( r * r )
ipt = 1
do i = 1, nch
  w(ipt) = w(ipt) + twormu * ( vv0  - (ered - eint(i)) ) &
           + r2 * cent(i)
  ipt = ipt + nmax + 1
end do

!ABER
!       call druckq(w,nmax,nch,'potential matrix',icount)
!ABER
!  look for innermost and outermost turning points
wmax = - 1.e+20
wmin = - wmax
ipt = 1
do i = 1, nch
  if (w(ipt) .lt. wmin) wmin = w(ipt)
  !  ignore closed channels
  if (eint(i) .lt. ered) then
    if (w(ipt) .gt. wmax) wmax = w(ipt)
  end if
  ipt = ipt+ nmax + 1
end do
if (iflag .ge. 2) go to 70
if (wmin .gt. zero) go to 90
if (iflag - 1) 60, 65, 90
60 iflag = 1
rtmn = r
go to 90
65 if (wmax .gt. zero) go to 90
iflag = 2
rtmx = r
!  check that no centrifugal barrier is present
70 if (wmin .lt. zero) go to 75
rtmn = r
75 if (wmax .lt. zero) go to 90
rtmx = r
90 continue

cpupot = cpupot + get_cpu_time()
! here for 2s-2p scattering
! fill in upper triangle of w matrix
!  first fill in upper half of original matrix
if (ibasty .eq. 10) then

  nmaxp1 = nmax + 1
  icolpt = 2
  irowpt = nmaxp1
  do icol = 1, nch - 1
!  icolpt points to first sub-diagonal element in column icol
!  irowpt points to first super-diagonal element in row icol
!  ncol is number of subdiagonal elements in column icol
    ncol = nch - icol
    call dcopy (ncol, w(icolpt), 1, w(irowpt), nmax)
    icolpt = icolpt + nmaxp1
    irowpt = irowpt + nmaxp1
  end do
!  transform w matrix into case e basis
  call trans22(w,nch, nmax)
endif

return
end
! ----------------------------------------------------------------------
subroutine potent (w, vecnow, scmat, eignow, hp, scr, &
   rnow, drnow, xlarge, nch, nmax, v2)
! ----------------------------------------------------------------------
!  this subroutine first sets up the wave-vector matrices:
!    w = w[rnow + 0.5 drnow/sqrt(3)] and w = w[rnow - 0.5 drnow/sqrt(3)]
!     b                                   a
!  then diagonalizes the average; i.e. 0.5 (w  + w )
!                                            b    a
!  the radial derivative of the wavevector matrix is calculated by finite
!  difference, using the nodes of a two-point gauss-legendre quadrature
!              1/2
!   d(w)/dr = 3    (w  - w ) / drnow
!                    b    a
!  this is then transformed into the local basis
!  author:  millard alexander
!  current revision date: 23-feb-2004
! ---------------------------------------------------------------------
!  variables in call list:
!    w:        on return:  contains transform of dh/dr
!                          this is the same as the negative of the
!                          wn-tilde-prime matrix
!    vecnow:   on return:  contains matrix of eigenvectors
!    scmat:    scratch matrix
!    eignow:   on return:  contains eigenvalues of wavevector matrix
!    hp:       on return: contains diagonal elements of transformed dh/dh
!                         this is the same as the negative of the diagonal
!                         elements of the wn-tilde-prime matrix
!    scr:      scratch vector
!    rnow:     midpoint of the current interval
!    drnow:    width of the current interval
!    xlarge:   on return contains largest off-diagonal element in
!              wn-tilde-prime matrix
!    nch:      number of channels
!    nmax:     maximum row dimension of matrices and maximum dimension of
!              vectors
! ----------------------------------------------------------------------
   use mod_ancou, only: ancou_type
   use mod_himatrix, only: transp
   use mod_hiblas, only: dscal, dcopy, daxpy_wrapper, dsyevr_wrapper
   implicit none
!  square matrices (of row dimension nmax)
real(8), dimension(nmax*nmax), intent(out) :: w
real(8), dimension(nmax*nmax), intent(out) :: vecnow
real(8), dimension(nmax*nmax), intent(out) :: scmat
!  vectors dimensioned at least nch
real(8), dimension(nmax), intent(out) :: eignow
real(8), dimension(nmax), intent(out) :: hp
real(8), dimension(nmax), intent(out) :: scr
real(8), intent(in) :: rnow
real(8), intent(in) :: drnow
real(8), intent(out) :: xlarge
integer, intent(in) :: nch
integer, intent(in) :: nmax
type(ancou_type), intent(in) :: v2

!      real eignow, hp, scmat, scr, vecnow, w
!      real drnow, en, fact, half, one, ra, rb, rnow, sq3, xlarge, xmin1
integer icol, ierr, ipt, nrow
!  local arrays (for lapack dsyevr)
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
integer, dimension(2*nch) :: isuppz
integer, dimension(10*nch) :: iwork
real(8), dimension(57*nch) :: work
#endif
data ione / 1 /
integer :: ione
data one,xmin1,half,sq3 /1.d0,-1.d0,0.5d0,1.732050807568877d0/
real(8) :: one, xmin1, half, sq3
integer :: nmaxp1, nmaxm1
real(8) :: ra, rb
real(8) :: fact, abstol, vl, vu
integer :: lwork, liwork, lsup, il, iu, m, I
nmaxp1 = nmax + 1
nmaxm1 = nmax - 1
ra = rnow - half * drnow / sq3
rb = rnow + half * drnow / sq3
!  scmat is used to store the wavevector matrix at rb
call potmat (w, ra, nch, nmax, v2)
call potmat (scmat, rb, nch, nmax, v2)
!  since potmat returns negative of lower triangle of w(r) matrix (eq.(3) of
!  m.h. alexander, "hybrid quantum scattering algorithms ..."),
!  next loop changes its sign
ipt = 1
do 100 icol = 1, nch
!  nrow is the number of (diagonal plus subdiagonal) elements in column icol
!  ipt points to the diagonal element in column icol for a matrix stored in
!  packed column form
  nrow = nch - icol + 1
  call dscal (nrow, xmin1, w(ipt), 1)
  call dscal (nrow, xmin1, scmat(ipt), 1)
  ipt = ipt + nmaxp1
100 continue
!  next loop stores average wavevector matrix in scmat and derivative of
!  hamiltonian matrix, in free basis, in w
fact =  - sq3 / drnow
!  the additional minus sign in the preceding expression is introduced because
!  dh/dr =-dw/dr;  see eq.(9) of
!  m.h. alexander, "hybrid quantum scattering algorithms ..."
ipt = 1
do 105 icol = 1, nch
!  nrow is the number of (diagonal plus subdiagonal) elements in column icol
!  ipt points to the diagonal element in column icol for a matrix stored in
!  packed column form
!  hp and scr are used as scratch vectors here
  nrow = nch - icol + 1
  call dcopy (nrow, scmat(ipt), 1, scr, 1)
  call daxpy_wrapper (nrow, one, w(ipt), 1, scmat(ipt), 1)
  call daxpy_wrapper (nrow, xmin1, w(ipt), 1, scr, 1)
  call dscal (nrow, half, scmat(ipt), 1)
  call dscal (nrow, fact, scr, 1)
  call dcopy (nrow, scr, 1, w(ipt), 1)
  ipt = ipt + nmaxp1
105 continue
!  next loop fills in upper triangles of w and scmat
if (nch .gt. 1) then
  ipt = 2
  do 110 icol = 1, nch -1
!  ipt points to the first subdiagonal element in column icol
!  nrow is the number of subdiagonal elements in column icol
    nrow = nch - icol
    call dcopy (nrow, w(ipt), 1, w(ipt + nmaxm1), nmax)
    call dcopy (nrow, scmat(ipt), 1, scmat(ipt + nmaxm1), nmax)
    ipt = ipt + nmaxp1
110   continue
end if
! ----------------------------------------------------------------------
!  diagonalize scmat at rnow and transpose matrix of eigenvectors
!  after transposition, the vecnow matrix is identical to the tn matrix
!  of eq.(6) of m.h. alexander, "hybrid quantum scattering algorithms ..."
!  now call eispack eigenvalue and eigenvector routine (hp is used as
!  a scratch vector here)
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
call rs (nmax, nch, scmat, eignow, ione, vecnow, scr, hp, ierr)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
lwork=57*nch
liwork=10*nch
abstol=1.d-16
lsup=2*nch
vl = 0.0
vu = 0.0

call dsyevr_wrapper('V','A','L',nch,scmat,nmax,vl,vu,il,iu,abstol,m, &
   eignow,vecnow,nmax,isuppz,work,lwork,iwork,liwork,ierr)
#endif
if (ierr .ne. 0) then
  write (6, 115) ierr
  write (9, 115) ierr
115   format (' *** IERR =',i3,' IN AIRPRP/POTENT/RS;  ABORT ***')
  write (9, 120) (eignow (i), i=1, nch)
120   format (' EIGENVALUES ARE:',/,8(1pe16.8) )
  call exit
end if
call transp (vecnow, nch, nmax)
!  transform the derivative into the local basis
!  subroutine dtrans returns the negative of the wn-tilde-prime matrix;
!  eq.(9) of m.h. alexander, "hybrid quantum scattering algorithms ..."
call dtrans (w, vecnow, scmat, hp, xlarge, nch, nmax, ione)
return
end
! -----------------------------------------------------------------------
subroutine spropn (rnow, width, eignow, hp, y1, y4, y2, &
                   gam1, gam2, nch)
!-----------------------------------------------------------------------------
!  this subroutine calculates the diagonal matrices to propagate the
!  log-derivative matrix through the current interval
!  also calculated are the ihomogeneous propagators (explained below)
!-----------------------------------------------------------------------------
!  variables in call list:
!    rnow:       midpoint of the current interval
!    width:      width of the current interval
!    eignow:     array containing the wavevectors
!                these are defined by eq. (6) of m.alexander,
!                j. chem. phys. 81,4510 (1984)
!    hp:         array containing the negative of diagonal elements of the
!                derivative of the wavevector matrix at the center of the
!                current interval [see eq. (9) of m.alexander,
!                j. chem. phys. 81,4510 (1984)
!                this array thus contains the derivative of the diagonal
!                elements of the transformed hamiltonian matrix
!    y1, y2, y4: on entry, contain the desired diagonal elements of the
!                homogeneous propagator
!    gam1, gam2: on return, if photof .true, contain the desired diagonal
!                elements of the ihomogeneous propagators
!                otherwise gam1 and gam2 are returned as zero
!    nch:        the number of channels, this equals the dimensions of the
!                eignow, hp, y1, y2, y2, gam1, and gam2 arrays
!-----------------------------------------------------------------------------
!  the key equations, reproduced below, are taken from
!  m. alexander and d. manolopoulos, "a stable linear reference potential
!  algorithm for solution ..."
!  each uncoupled equation can be written as:
!         2    2
!     [ d / dr + eignow - hp * r ] f(r) = 0
!     where r is the distance from the midpoint of the current interval
!  the linearly indepedent solutions are the airy functions ai(x) and bi(x)
!  where  x = alpha (r + beta)
!                   1/3
!  with   alpha = hp   , and beta = (-eignow) / hp
!  the three diagonal elements of the cauchy propagator necessary to propagate
!  the log-derivative matrix are:
!    b = pi [ ai(x ) bi(x ) - ai(x )bi(x ) ] / alpha
!                 1      2        2     1
!    a = pi [ - ai'(x ) bi(x ) + ai(x ) bi'(x ) ]
!                    1      2        2       1
!    d = pi [ ai(x ) bi'(x ) - ai'(x ) bi(x ) ]
!                 1       2         2      1
!    where x  = alpha ( beta + width / 2) and
!           2
!          x  = alpha ( beta - width / 2)
!           1
!  here "width" denotes the width of the interval
!  the diagonal elements of the "imbedding type" propagator are given in terms
!  of the diagonal elements of the cauchy propagator by:
!     y = a/b     y = y = 1/b    and   y = d/b
!      1           2   3                4
!  for the calculation sof the homogeneous propagators
!  the airy functions are defined in terms of their moduli and phases
!  for negative x these definitions are:
!      ai(-x) = m(x) cos[theta(x)]
!      bi(-x) = m(x) sin[theta(x)]
!      ai'(-x) = n(x) cos[phi(x)]
!      bi'(-x) = n(x) sin[phi(x)]
!  in other words
!          2              2        2
!      m(x)  = sqrt[ ai(x)  + bi(x)  ]
!          2               2         2
!      n(x)  = sqrt[ ai'(x)  + bi'(x)  ]
!      theta(x) = atan [ bi(x) / ai(x) ]
!      phi(x)   = atan [ bi'(x) / ai'(x) ]
!  for positive x the moduli and phases are defined by:
!      ai(x) = m(x) sinh[theta(x)]
!      bi(x) = m(x) cosh[theta(x)]
!      ai'(x) = n(x) sinh[phi(x)]
!      bi'(x) = n(x) cosh[phi(x)]
!  in other words
!          2              2        2
!      m(x)  = sqrt[ bi(x)  - ai(x)  ]
!          2               2         2
!      n(x)  = sqrt[ bi'(x)  - ai'(x)  ]
!      theta(x) = atanh [ ai(x) / bi(x) ]
!      phi(x)   = atanh [ ai'(x) / bi'(x) ]
!  here the the exponentially scaled airy functions
!  ai(x), ai'(x), bi(x), bi'(x) are:
!      ai(x)  = ai(x)  * exp[zeta]
!      ai'(x) = ai'(x) * exp[zeta]
!      bi(x)  = bi(x)  * exp[-zeta]
!      bi'(x) = bi'(x) * exp[-zeta]
!                          3/2
!      where zeta = (2/3) x
!  note that for positive x the phases are labeled chi and eta in
!  m. alexander and d. manolopoulos, "a stable linear reference potential
!  algorithm for solution ..."
!-----------------------------------------------------------------------------
!  for both x  and x  negative
!            1      2
!  (this corresponds to a channel which is classically open at both ends of th
!  interval)
!  we find:
!  y     = 1 / { m  m  sin[theta -theta ] }
!   2             1  2          2      1
!          n  sin[phi -theta ]
!           1        1      2
!  y    = ----------------------
!   1      m  sin[theta - theta ]
!           1          2       1
!          n  sin[phi -theta ]
!           2        2      1
!  y    = ----------------------
!   4      m  sin[theta - theta ]
!           2          2       1
!  here the subscripts 1 and 2 imply the moduli and phases evaluated at x = x
!                                                                            1
!  and x = x  , respectively
!           2
!-----------------------------------------------------------------------------
!  for both x  and x  positive
!            1      2
!  (this corresponds to a channel which is classically closed at both ends of
!  the interval)
!  we find:
!  1 / y  =  m  m  cosh[z -z ] { sinh[theta -theta ]
!       2     1  2       2  1              1      2
!                         + tanh[z -z ] sinh[theta +theta ] }
!                                 2  1            1      2
!                     3/2
!   where z  = (2/3) x     and similarly for z
!          1          1                       2
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           1              2    1          2  1            2    1
!  y    = --------------------------------------------------------
!   1      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           1              2      1          2  1            2      1
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           2              1    2          2  1            1    2
!  y     = --------------------------------------------------------
!   4      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           2              2      1          2  1            2      1
!-----------------------------------------------------------------------------
!  for x  positive and x  negative we find:
!       1               2
!  1 / y  = m  m  cosh[z ] cosh[theta ] { - cos[theta ] (1 + tanh[z ])
!       2    1  2       1            1               2             1
!                                   + tanh[theta ] sin[theta ] (1 - tanh[z ])
!                                               1           2             1
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       1           2            1             1           2             1
! y = ------------------------------------------------------------------------
!  1   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1           2            1               1           2             1
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       2         2            1               1         2             1
! y  = -----------------------------------------------------------------------
!  4   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2           2            1               1           2             1
!-----------------------------------------------------------------------------
!  for x  negative and x  positive we find:
!       1               2
! 1 / y  = m  m  cosh[z ] cosh[theta ] { cos[theta ] (1 + tanh[z ])
!      2    1  2       2            2             1             2
!                                  - tanh[theta ] sin[theta ] (1 - tanh[z ]) }
!                                              2           1             2
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       1         1            2               2         1             2
! y  = -----------------------------------------------------------------------
!  1   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1          1            2               2           1             2
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       2           1            2             2           1             2
! y  = -----------------------------------------------------------------------
!  4   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2          1            2               2           1             2
!-----------------------------------------------------------------------------
!  for the special case of a constant reference potential (hp=0)
!  then the propagators are:
!  for eignow .gt. 0 (the classically allowed region)
!    y1 = y4 = k cot (k width)
!    y2 = k / sin (k width)
!    where k = sqrt (eignow)
!  for eignow .lt. 0 (the classically forbidden region)
!    y1 = y4 = kap coth (kap width)
!    y2 = kap / sinh (kap width)
!
!    where kap = sqrt (-eignow)
!-----------------------------------------------------------------------------
!  this subroutine also calculates the diagona linhomogeneous log-derivative
!  propagators.  key equations are 9.13, 9.14, and 9.25 of the ph. d.
!  thesis of d. manolopoulos
!  the diagonal elements of the two inhomogeneous propagators are defined
!  in terms of the linearly intependent solutions psi+ and psi-, which
!  are
!  psi+ = [ - bi(x2)ai(x) + ai(x2)bi(x) ] y2 / w
!  and
!  psi- = [ - bi(x1)ai(x) + ai(x1)bi(x) ] y2 / w
!  where w is the wronskian (1/pi)
!  in the determination of these inhomogeneous propagators
!  the airy functions are defined as follows:
!  for negative x :
!      ai(-x) = ai(x) cos[th] + bi(x) sin[th]
!      bi(-x) = bi(x) cos[th] - ai(x) cos[th]
!                            3/2
!      where zeta = (2/3) |x|
!  for positive x :
!      ai(x) = exp(-zeta) ai(x)
!      bi(x) = exp(zeta) bi(x)

!  the integrals of the airy functions are defined as:

!  for a > 0
!      int[ai(x),{0,a}] = 1/3 - exp(-zeta)iai(a)
!      int[bi(x),{0,a}] = exp(zeta)ibi(a)
!  and for a < 0
!      int[ai(x),{a,0}] = int[ai(-x),{0,-a}]
!                       = 2/3 - q(a) cos(th) + p(a) sin(th)
!      int[bi(x),{a,0}] = int[bi(-x),{0,-b}] =
!                       = p(a) cos(th) + q(a) sin(th)
!  we further assume that the ground state wavefunction times the dipole
!  moment function can be expanded as phi(x) = phi0 + phi1 x
!  where phi1 = d[phi,rmid]/alpha and
!        phi0 = phi(rmid) - beta * d[phi,rmid]
!-----------------------------------------------------------------------------
!  for both x  and x  negative
!            1      2
!  (this corresponds to a channel which is classically open at both ends of th
!  interval)
!  we find:
!  y     = 1 / { m  m  sin[theta -theta ] }
!   2             1  2          2      1
!          n  sin[phi -theta ]
!           1        1      2
!  y    = ----------------------
!   1      m  sin[theta - theta ]
!           1          2       1
!          n  sin[phi -theta ]
!           2        2      1
!  y    = ----------------------
!   4      m  sin[theta - theta ]
!           2          2       1
!  here the subscripts 1 and 2 imply the moduli and phases evaluated at x = x
!                                                                            1
!  and x = x  , respectively
!           2
!-----------------------------------------------------------------------------
!  for both x  and x  positive
!            1      2
!  (this corresponds to a channel which is classically closed at both ends of
!  the interval)
!  we find:
!  1 / y  =  m  m  cosh[z -z ] { sinh[theta -theta ]
!       2     1  2       2  1              1      2
!                         + tanh[z -z ] sinh[theta +theta ] }
!                                 2  1            1      2
!                     3/2
!   where z  = (2/3) x     and similarly for z
!          1          1                       2
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           1              2    1          2  1            2    1
!  y    = --------------------------------------------------------
!   1      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           1              2      1          2  1            2      1
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           2              1    2          2  1            1    2
!  y     = --------------------------------------------------------
!   4      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           2              2      1          2  1            2      1
!-----------------------------------------------------------------------------
!  for x  positive and x  negative we find:
!       1               2
!  1 / y  = m  m  cosh[z ] cosh[theta ] { - cos[theta ] (1 + tanh[z ])
!       2    1  2       1            1               2             1
!                                   + tanh[theta ] sin[theta ] (1 - tanh[z ])
!                                               1           2             1
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       1           2            1             1           2             1
! y = ------------------------------------------------------------------------
!  1   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1           2            1               1           2             1
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       2         2            1               1         2             1
! y  = -----------------------------------------------------------------------
!  4   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2           2            1               1           2             1
!-----------------------------------------------------------------------------
!  for x  negative and x  positive we find:
!       1               2
! 1 / y  = m  m  cosh[z ] cosh[theta ] { cos[theta ] (1 + tanh[z ])
!      2    1  2       2            2             1             2
!                                  - tanh[theta ] sin[theta ] (1 - tanh[z ]) }
!                                              2           1             2
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       1         1            2               2         1             2
! y  = -----------------------------------------------------------------------
!  1   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1          1            2               2           1             2
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       2           1            2             2           1             2
! y  = -----------------------------------------------------------------------
!  4   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2          1            2               2           1             2
!-----------------------------------------------------------------------------
!  for the special case of a constant reference potential (hp=0)
!  then the propagators are:
!  for eignow .gt. 0 (the classically allowed region)
!    y1 = y4 = k cot (k width)
!    y2 = k / sin (k width)
!    where k = sqrt (eignow)
!  for eignow .lt. 0 (the classically forbidden region)
!    y1 = y4 = kap coth (kap width)
!    y2 = kap / sinh (kap width)
!
!    where kap = sqrt (-eignow)
!-----------------------------------------------------------------------------
!  written by:  millard alexander
!  current revision date (algorithm):  30-dec-1994
!-----------------------------------------------------------------------------
use mod_coqvec2, only: q => q2
use mod_phot, only: photof
use mod_hiutil, only: intairy
use mod_hivector, only: dset
use mod_hiamp, only: airymp
implicit double precision (a-h,o-z)
!      implicit none
double precision a, b, bfact, cs, cs1, cs2, csh, dalph2, dalpha, &
    darg, dbeta, dcay, delzet, denom, dhalf, dkap, dlzeta, &
    dmmod1, dmmod2, dnmod1, dnmod2, doneth, dphi1, dphi2, &
    dpi, droot, dslope, dthet1, dthet2, dtnhfm, &
    dtnhfp, dx1, dx2, dzeta1, dzeta2, emz1, emz2, &
    ez1, ez2, fact, oflow, one, rnow, scai1, scai2, scbi1, scbi2, &
    sn, sn1, sn2, snh, tnhfac, width, x1, x2, xairy1, xairy2, &
    xbiry1, xbiry2, zero
double precision eignow, gam1, gam2, hp, y1, y2, y4
integer i, nch
dimension eignow(1), hp(1), y1(1), y2(1), y4(1), gam1(1), gam2(1)
data     doneth,        dhalf &
  / 0.333333333333333d0, 0.5d0 /
data zero, one /0.d0, 1.d0/
data  dpi / 3.1415926535897932d0 /
!  the parameter oflow is the largest value of x for which exp(x)
!  does not cause a single precision overflow
!                                     n
!  a reasonable value is x = [ ln(2) 2 ] - 5, where n is the number of bits in
!  the characteristic of a floating point number
data oflow / 83.d0 /
if (.not. photof) then
 call dset(nch,zero,gam1,1)
 call dset(nch,zero,gam2,1)
endif
!     now determine propagators for all nch channels
do 10  i = 1, nch
  dslope = hp(i)
! activate next statement for constant reference potential
! force slope to equal zero, to force constant potential
!        dslope=0.d0
  darg = 1.e+10
  if (dslope .ne. 0.d0) &
    darg = log (abs(eignow(i))) - log (abs(dslope))
  if (darg .gt. 20.d0 .or. width .lt. 1.d-5) then
!  here if the relative slope in the wavevector matrix is less than 1.**(-20)
!  in magnitude, or sector width less than 1.e-5 bohr,
!  in which case the potential is assumed to be constant
    if (eignow(i) .gt. 0) then
!  here for classically allowed region (sines and cosines as reference
!  solutions)
      dcay = sqrt (eignow(i))
      darg = dcay * width
      sn=sin(darg)
      y1(i) = dcay / tan (darg)
      y4(i) = y1(i)
      y2(i) = dcay / sn
!  here for inhomogeneous propagators
      if (photof) then
        cs=cos(darg)
        b=rnow+width*dhalf
        a=rnow-width*dhalf
        denom=dcay*sn
        fact=(one-cs)*(q(i)-q(nch+i)*rnow)
        gam1(i)=(fact+(b-a*cs-sn/dcay)*q(nch+i))/denom
        gam2(i)=(fact+(a-b*cs+sn/dcay)*q(nch+i))/denom
      endif
    else
!  here for classically forbidden region (hyperbolic sines and cosines as
!  reference solutions)
      dkap = sqrt ( - eignow(i))
      darg = dkap * width
      snh=sinh(darg)
      y1(i) = dkap / tanh (darg)
      y4(i) = y1(i)
      y2(i) = dkap / snh
!  here for inhomogeneous propagators
      if (photof) then
        csh=cosh(darg)
        b=rnow+width*dhalf
        a=rnow-width*dhalf
        denom=dkap*snh
        fact=(-one+csh)*(q(i)-q(nch+i)*rnow)
        gam1(i)=(fact+(-b+a*csh+snh/dkap)*q(nch+i))/denom
        gam2(i)=(fact+(-a+b*csh-snh/dkap)*q(nch+i))/denom
      endif
    end if
  else
!  here if the relative slope in the wavevector matrix is greater than
!  1.**(-20) in magnitude, in which case a linear reference potential is used,
!  with airy functions as reference solutions
    droot = ( abs (dslope) ) ** doneth
    dalpha   = sign (droot, dslope)
    dbeta = - eignow(i) / dslope
    dx1 = dalpha * ( dbeta - width * dhalf)
    dx2 = dalpha * ( dbeta + width * dhalf)
    call airymp (dx1, dthet1, dphi1, dmmod1, dnmod1,scai1, scbi1, &
             dzeta1)
    call airymp (dx2, dthet2, dphi2, dmmod2, dnmod2,scai2, scbi2, &
             dzeta2)
    if (photof) then
! determine required airy integrals
      call intairy(dx1, xairy1, xbiry1)
      call intairy(dx2, xairy2, xbiry2)
! convert ground state wavefunction and its derivative from r as
! independent variable to x
      q(i)=q(i)-dbeta*q(nch+i)
      q(nch+i)=q(nch+i)
    endif

    x1 = dx1
    x2 = dx2
!-----------------------------------------------------------------------------
    if (x1 .gt. zero .and. x2 .gt. zero) then
!  here for both x  and x  positive
!                 1      2
      tnhfac = tanh(dzeta2 - dzeta1)
      bfact = sinh(dthet1 - dthet2) + &
              tnhfac * sinh(dthet1 + dthet2)
      dlzeta = dzeta2 - dzeta1
      y2(i) = zero
      if (abs(dlzeta) .le. oflow) then
        b = dmmod1 * dmmod2 * cosh(dzeta2 - dzeta1) * bfact
        y2(i) = 1. / b
      end if
      y1(i) = dnmod1 * (sinh(dthet2 - dphi1) &
            - tnhfac * sinh(dthet2 + dphi1) ) / (dmmod1 * bfact)
      y4(i) = dnmod2 * (sinh(dthet1 - dphi2) &
            + tnhfac * sinh(dthet1 + dphi2) ) / (dmmod2 * bfact)
      if (photof) then
        gam1(i)=-scbi2*xairy2-scai2*xbiry2 &
           +exp(dlzeta)*scbi2*xairy1+exp(-dlzeta)*scai2*xbiry1
        gam2(i)=-scbi1*xairy1-scai1*xbiry1 &
           +exp(dlzeta)*scai1*xbiry2+exp(-dlzeta)*scbi1*xairy2
      endif
!-----------------------------------------------------------------------------
    else if (x1 .le. zero .and. x2 .le. zero) then
!  here for both x  and x  negative
!                 1      2
      b =  dmmod1 * dmmod2 * sin(dthet2 - dthet1)
      y2(i) = 1. / b
      y1(i) = dnmod1 * sin(dphi1 - dthet2) &
            / (dmmod1 * sin(dthet2 - dthet1) )
      y4(i) = dnmod2 * sin(dphi2 - dthet1) &
            / (dmmod2 * sin(dthet2 - dthet1) )
      if (photof) then
        delzet=dzeta2-dzeta1
        cs=cos(delzet)
        sn=sin(delzet)
        gam1(i)=-scai2*xairy2+scbi2*xbiry2 &
               +cs*(scai2*xairy1-scbi2*xbiry1) &
               +sn*(scai2*xbiry1+scbi2*xairy1)
        gam2(i)=-scai1*xairy1+scbi1*xbiry1 &
               +cs*(scai1*xairy2-scbi1*xbiry2) &
               -sn*(scai1*xbiry2+scbi1*xairy2)
      endif
!-----------------------------------------------------------------------------
    else if (x1 .gt. zero .and. x2 .le. zero) then
!  here for x  positive and x  negative
!            1               2
      dtnhfp = 1 + tanh(dzeta1)
      dtnhfm = 1 - tanh(dzeta1)
      bfact = cosh(dthet1) * ( - cos(dthet2) * dtnhfp &
            + tanh(dthet1) * sin(dthet2) * dtnhfm)
      y2(i) = zero
      if (abs(dzeta1) .le. oflow) then
        y2(i) = cosh(dzeta1) * (dmmod1 * dmmod2 * bfact)
        y2(i) = one / y2(i)
      end if
      y1(i) = (dnmod1 * cosh(dphi1) * ( cos(dthet2) * dtnhfp &
            - tanh(dphi1) * sin(dthet2) * dtnhfm) ) &
            / (dmmod1 * bfact)
      y4(i) = (dnmod2 * cosh(dthet1) * ( - cos(dphi2) * dtnhfp &
            + tanh(dthet1) * sin(dphi2) * dtnhfm) ) &
            / (dmmod2 * bfact)
      if (photof) then
        cs2=cos(dzeta2)
        sn2=sin(dzeta2)
        ez1=exp(dzeta1)
        emz1=one/ez1
        gam1(i)=scbi2*(-cs2+xbiry2) &
               +scai2*(sn2-xairy2) &
                +emz1*xairy1*(scbi2*cs2-scai2*sn2) &
                +ez1*xbiry1*(scai2*cs2+scbi2*sn2)
        gam2(i)=-scai1*xbiry1-scbi1*xairy1 &
                +emz1*scai1*(xairy2*cs2-xbiry2*sn2) &
                +ez1*scbi1*(one-xbiry2*cs2-xairy2*sn2)
      endif
!-----------------------------------------------------------------------------
    else if (x2 .gt. zero .and. x1 .le. zero) then
!  here for x  positive and x  negative
!            2               1
      dtnhfp = 1 + tanh(dzeta2)
      dtnhfm = 1 - tanh(dzeta2)
      bfact = cosh(dthet2) * ( cos(dthet1) * dtnhfp &
            - tanh(dthet2) * sin(dthet1) * dtnhfm)
      y2(i) = zero
      if (abs(dzeta2) .le. oflow) then
        y2(i) =  cosh(dzeta2) * (dmmod1 * dmmod2 * bfact)
        y2(i) = one / y2(i)
      end if
      y4(i) = (dnmod2 * cosh(dphi2) * ( cos(dthet1) * dtnhfp &
            - tanh(dphi2) * sin(dthet1) * dtnhfm) ) &
            / (dmmod2 * bfact)
      y1(i) = (dnmod1 * cosh(dthet2) * ( - cos(dphi1) * dtnhfp &
            + tanh(dthet2) * sin(dphi1) * dtnhfm) ) &
            / (dmmod1 * bfact)
      if (photof) then
        ez2=exp(dzeta2)
        emz2=one/ez2
        cs1=cos(dzeta1)
        sn1=sin(dzeta1)
        gam1(i)=-scai2*xbiry2-scbi2*xairy2 &
                +emz2*scai2*(xairy1*cs1-xbiry1*sn1) &
                +ez2*scbi2*(one-xbiry1*cs1-xairy1*sn1)
! bug corrected here 4/14/94
        gam2(i)=scbi1*(-cs1+xbiry1) &
               +scai1*(sn1-xairy1) &
                +emz2*xairy2*(scbi1*cs1-scai1*sn1) &
                +ez2*xbiry2*(scai1*cs1+scbi1*sn1)
! here is the old, incorrect code`
!             gam1(i)=scbi1*(-cs1+xbiry1)
!    :               +scai1*(sn1-xairy1)
!    :                +emz2*xairy2*(scbi1*cs1-scai1*sn1)
!    :                +ez2*xbiry2*(scai1*cs1+scbi1*sn1)
      endif
    end if
!-----------------------------------------------------------------------------
    y1(i) = dalpha * y1(i)
    y4(i) = dalpha * y4(i)
    y2(i) = dalpha * y2(i) / dpi
    if (photof) then
        dalph2=dalpha*dalpha
        gam1(i)=q(i)*gam1(i)*y2(i)*dpi &
                +q(nch+i)*(y1(i)-y2(i))/dalpha
        gam2(i)=q(i)*gam2(i)*y2(i)*dpi &
                +q(nch+i)*(y4(i)-y2(i))/dalpha
        gam1(i)=gam1(i)/dalph2
        gam2(i)=gam2(i)/dalph2
    endif
!  at this point the y1, y2, and y4 propagators correspond identically to
!  eqs. (38)-(44) of m. alexander and d. manolopoulos, "a stable linear
!  reference potential algorithm for solution ..."
  end if
10 continue
return
end
! -----------------------------------------------------------------------
subroutine steppr (vecnow, vecnew, tmat, nmax, n)
!  determine matrix to transform log-deriv matrix into new interval
!  see eq. (22) of m.h. alexander, "hybrid quantum scattering algorithms ..."
! --------------------------------------------------------------------------
!  variables in call list:
!    vecnow:     on entry: matrix of eigenvectors of wavevector matrix in
!                current interval
!                on return: matrix of eigenvectors of wavevector matrix in
!                new interval - this is the matrix tn in eq. (22) of
!                m.h. alexander, "hybrid quantum scattering algorithms ..."
!    vecnew:     on entry:  contains matrix of eigenvectors in next interval
!    tmat:       on return: contains transformation matrix pn in eq. (22)
!    n:          number of channels
!    nmax:       maximum row dimension of matrices
!  subroutines called:
!     rgmmul:    generalized matrix multiply, called here to evaluate
!                a.b-transpose
! --------------------------------------------------------------------------
use mod_hivector, only: matmov
use mod_hiblas, only: dgemm
implicit double precision (a-h,o-z)
!      real vecnow, vecnew, tmat
integer n, nmax, isw
!  matrices of maximum row dimension nmax, stored in packed column form
dimension vecnow(1), vecnew(1), tmat(1)
data isw / 0/
#if defined(HIB_NONE)
call mxma (vecnew, 1, nmax, vecnow, nmax, 1, tmat, 1, nmax, &
            n, n, n)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
call dgemm('n','t',n,n,n,1.d0,vecnew,nmax,vecnow,nmax, &
           0d0,tmat,nmax)
#endif
!  restore eigenvectors
call matmov (vecnew, vecnow, n, n, nmax, nmax)
return
end
! -----------------------------------------------------------------------
function turn(e)
! current revision date: 23-sept-87
use constants, only: econv
use mod_hipot, only: pot
implicit none
real(8), intent(in) :: e
real(8) :: ee, r, dr, vv0
real(8) :: turn
ee = e/econv
r = 3.0d0
dr = 0.5d0
10 r = r+dr
call pot(vv0,r)
if(vv0-ee) 20,50,30
20 if(dr.lt.0) goto 10
goto 40
30 if(dr.gt.0) goto 10
40 dr = -dr*0.5d0
if(abs(dr).gt.0.01d0) goto 10
50 turn = r
return
end

! -----------------------------------------------------------------------
subroutine wavevc (w, eignow, rnow, nch, nmax, v2)
!  this subroutine first sets up the wavevector matrix at rnow
!  then diagonalizes this matrix
!  written by:  millard alexander
!  current revision date: 14-dec-2007
! ----------------------------------------------------------------
!  variables in call list:
!  w:           matrix of maximum row dimension nmax used to store
!               wavevector matrix
!  eignow:      on return:  array containing eigenvalues of wavevector matrix
!  scr1, scr2:  scratch vectors of dimension at least nch
!  rnow:        value of interparticle separation at which wavevector matrix
!               is to be evaluated
!  nch:         number of channels
!  nmax:        maximum number of channels
!  subroutines called:
!     potmat:         determines wavevector matrix
!     tred1,tqlrat:   eispack routines to obtain eigenvalues of real,
!                     matrix
!     dsyevr:         latest lapack eigenvalue routine
! ----------------------------------------------------------------
use mod_ancou, only: ancou_type
use mod_hiblas, only: dscal, dcopy, dsyevr_wrapper
implicit double precision (a-h,o-z)
real(8), intent(out) :: w(nmax*nmax)
real(8), intent(out) :: eignow(nch)
real(8), intent(in) :: rnow
integer, intent(in) :: nch
integer, intent(in) :: nmax
type(ancou_type), intent(in) :: v2
!      real rnow, xmin1
!      real eignow, w
integer, parameter :: ldz = 1
integer icol, ierr, ipt, nmaxm1, nmaxp1, nrow
real(8), dimension(ldz, nch):: vecnow_unused   ! this is the z array that dsyevr wants, even if it's not used when jobz = 'N'
!     external potmat, tred1, tqlrat
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
real(8) :: scr1(nch)
real(8) :: scr2(nch)
#endif
!  local arrays (for lapack dsyevr)
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
dimension isuppz(2*nch),iwork(10*nch),work(57*nch)
#endif

! ------------------------------------------------------------------
data xmin1 / -1.d0/
nmaxp1 = nmax + 1
nmaxm1 = nmax - 1
call potmat (w, rnow, nch, nmax, v2)
!  since potmat returns negative of lower triangle of w(r) matrix (eq.(3) of
!  m.h. alexander, "hybrid quantum scattering algorithms ..."),
!  next loop changes its sign
ipt = 1
do icol = 1, nch
!  nrow is the number of (diagonal plus subdiagonal) elements in column icol
!  ipt points to the diagonal element in column icol for a matrix stored in
!  packed column form
  nrow = nch - icol + 1
  call dscal (nrow, xmin1, w(ipt), 1)
  ipt = ipt + nmaxp1
end do
!  next loop fills in upper triangle of w
if (nch .gt. 1) then
  ipt = 2
  do icol = 1, nch -1
!  ipt points to the first subdiagonal element in column icol
!  nrow is the number of subdiagonal elements in column icol
    nrow = nch - icol
    call dcopy (nrow, w(ipt), 1, w(ipt + nmaxm1), nmax)
    ipt = ipt + nmaxp1
  end do
end if
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
lwork=57*nch
liwork=10*nch
abstol=1.e-16
lsup=2*nch
vl = 0.0
vu = 0.0
call dsyevr_wrapper('N','A','L',nch,w,nmax,vl,vu,il,iu,abstol,m, &
   eignow,vecnow_unused,ldz,isuppz,work,lwork,iwork,liwork,ierr)

if (ierr .ne. 0) then
  write (6, 115) ierr
  write (9, 115) ierr
115   format (' *** IERR =',i3,' IN WAVEVC/DSYEVR;  ABORT ***')
  write (9, 120) (eignow (i), i=1, nch)
120   format (' EIGENVALUES ARE:',/,8(1pe16.8) )
  call exit
end if
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
!  transform w to tridiagonal form
!  eignow, scr1 and scr2 are used as scratch vectors here
call tred1 (nmax, nch, w, eignow, scr1, scr2)
!  get eigenvalues of tridiagonal matrix
call tqlrat (nch, eignow, scr2, ierr)
if (ierr .ne. 0) then
  write (9, 130) ierr
  write (6, 130) ierr
130   format &
    (' *** TQLRAT IERR =', i3, ' .N.E. 0 IN WAVEVC; ABORT ***')
  call exit
end if
#endif
return
end

subroutine airprp (z, &
   xf, rend, drnow, en, &
   tolai, rincr, eshift, nch, nmax, itwo, iprint, twoen, noprin, v2)
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
use mod_ancou, only: ancou_type
use mod_hiba10_22p, only: energ22
use mod_par, only: par_iprint=>iprint
use mod_wave, only: irec, ifil, nchwfu, iendwv, get_wfu_airy_rec_length
use mod_selb, only: ibasty
use mod_phot, only: photof, wavefn, writs
use mod_himatrix, only: transp
#if (defined(HIB_UNIX) || defined(HIB_MAC)) && !defined(HIB_UNIX_IBM)
use mod_himatrix, only: mxma
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
use mod_himatrix, only: syminv
#endif
use mod_hivector, only: dset, vadd, vmul, matmov
use mod_hiblas, only: dscal, dcopy, daxpy_wrapper

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
type(ancou_type), intent(in) :: v2
logical :: airy_prop_completed
integer i, icol, ierr, ipt, izero, kstep, maxstp, &
        ncol, npt, nskip

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
real(8), allocatable :: w(:)
real(8), allocatable :: tmat(:)
real(8), allocatable :: vecnow(:)
real(8), allocatable :: vecnew(:)
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
integer(8) :: lrairy ! length of an airy record in bytes
! ----------------------------------------------------------------------------
! Save variables for subsequent energies
save spcmn, spcmx, rmin, maxstp

UNUSED_DUMMY(en)

allocate(w(nch*nmax))
allocate(tmat(nch*nmax))
allocate(vecnow(nch*nmax))
allocate(vecnew(nch*nmax))


if (.not.twoen) itwo = -1
if (itwo .le. 0) then
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
  call wavevc (w, eigold, rmin, nch, nmax, v2)
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
               rnow, drnow, xlarge, nch, nmax, v2)
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
    if (par_iprint .ge. 2) write (9, 55)
  55   format ('   ALSO ADIABATIC ENERGIES (HARTREE)')
  end if
end if
airy_prop_completed = .false.
if (itwo .ge. 0) then
  !  write or read relevant information
  call outmat (tmat, eigold, hp, eshift, drnow, rnow, &
               nch, nmax, itwo)
end if
!  start airy propagation
! ----------------------------------------------------------------------------
do kstep = 1, maxstp
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
    do i = 1, nphoto
      call vmul(y2,1,tmat(jnd:),1,cc,1,nch)
      call vadd(1,gam2(ind),1,cc,1,nch)
      call dcopy(nch,gam2(ind),1,q(ind),1)
      ind=ind+nch
      jnd=jnd+nmax
    end do
    ! q now contains gamma2(0,rnext) in local basis
  endif

  !                            -1
  !  evaluate  - y  ( y  + z )    y
  !               2    1    n      2
  !  in the next loop evaluate the full, rather than lower triangle
  !  changed for photodissociation calculation, old commands kept with
  !  c in first column
  npt = 1
  do i = 1, nch
    ncol = nch
    !        ncol = nch - i + 1
    fact = y2(i)
    call dscal (ncol, fact, z(npt), 1)
    npt = npt + nmax
    !        npt = npt + nskip
  end do
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
    !     The following information will not be written if writs set to F
    if (writs) then
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
      lrairy = get_wfu_airy_rec_length(nchwfu, 0)
    else
      lrairy = get_wfu_airy_rec_length(nchwfu, 1)
    end if
    !
    write (ifil, err=950) 'ENDWFUR', char(mod(irec, 256))
    iendwv = iendwv + lrairy

  end if
  !
  do i = 1, nch
    fact = - y2(i)
    call dscal (i, fact, z(i), nmax)
  end do
  !  add on  y
  !           4
  npt = 1
  call daxpy_wrapper (nch, one, y4, 1, z, nskip)
  !  now fill in upper half of z matrix
  ipt = 2
  if (nch .gt. 1) then
    do icol = 1, nch - 1
      !       ncol is number of subdiagonal elements in column icol
      !       ipt points to the first subdiagonal element in column icol
      !       (ipt + nmax - 1) points to the first superdiagonal element in row icol
      ncol = nch - icol
      call dcopy (ncol, z(ipt), 1, z(ipt + nmax -1), nmax)
      ipt = ipt + nskip
    end do
  end if
  if (itwo .le. 0) then
    !  obligatory write of step information if deviations from linear
    !  potential are unusually large
    !  this is only done if tolai .lt. 1, in which case the largest correction
    !  is used to estimate the next step
    if (tolai .lt. 1.d0) then
      cmax = max (cdiag, coff)
      if (cmax .gt. (5.d0* tolai)) then
        write (9,125)
        write (6,125)
125       format &
        (' ** ESTIMATED CORRECTIONS LARGER THAN 5*TOLAI IN AIRPRP')
        if (kstep .eq. 1) then
          write (9, 130)
          write (6, 130)
130         format ('    THE INITIAL VALUE OF DRNOW (SPAC*FSTFAC) IS', &
                  ' PROBABLY TOO LARGE')
        else
          write (9, 140)
          write (6, 140)
140         format &
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
150     format (i6, 4e10.3)
    end if
  end if
  !     get set for next step
  if (airy_prop_completed) then
    exit
  else
    if (itwo .le. 0) then
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
      else ! (tolai .lt. 1.)
        !  if tolai .ge. 1, then
        !  minimum step size is first interval width
        if (kstep .eq. 1) spcmn = drnow
        !  and next step size is tolai * present step size
        ! only if rnext > powr
        if (rnext .gt. rincr) drnow = tolai * drnow
      end if ! (tolai .lt. 1.)
      !  drnow is step size in next interval
      rlast = rnext
      rnext = rnext + drnow
      if (rnext .ge. rend) then
        airy_prop_completed = .true.
        rnext = rend
        drnow = rnext - rlast
      end if
      rnew = rlast + 0.5d0 * drnow
      if (kstep .gt. 1 .and. .not. airy_prop_completed) then
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
                   rnew, drnow, xlarge, nch, nmax, v2)
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
      if (itwo .lt. 0) then
        cycle ! next kstep
      end if
      if (airy_prop_completed) rnow = - rnow
    end if  ! (itwo .le. 0)
    ! 180
    !  write or read relevant information
    call outmat (tmat, eigold, hp, eshift, drnow, rnow, &
                 nch, nmax, itwo)
    if (itwo .ne. 0) then
      !  negative rnow is cue for last step in second energy calculation
      if (rnow .le. 0.d0) then
        rnow = - rnow
        airy_prop_completed = .true.
      end if
    end if
  end if  !  (airy_prop_completed)
  !     go back to start new step
end do
if (.not. airy_prop_completed) then
  !  the following statement is reached only if the integration has
  !  not reached the asymptotic region in maxstp steps
  write (9,210) maxstp, rnext
  210 format (' *** AIRY PROPAGATION NOT FINISHED IN', i4, &
          ' STEPS:  R-FIN SET TO', f8.4,' ***',/)
  xf = rnext
end if
if (itwo .ge. 0) then
  call outmat (vecnow, eigold, hp, eshift, drnow, xf, nch, &
               nmax, itwo)
end if
!  transform log-deriv matrix into free basis.  transformation matrix is
!  just vecnow-transpose; see eq.(24) of m.h. alexander, "hybrid quantum
!  scattering algorithms ..."
call transp (vecnow, nch, nmax)
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
  do i=1,nphoto
    call dcopy(nch,w(jnd),1,q(ind),1)
    ind=ind+nch
    jnd=jnd+nmax
  end do
endif
call dtrans (z, vecnow, w, hp, xlarge, nch, nmax, izero)

! if 2s-2p collisions, restore asymptotic case (e) energies
if (ibasty .eq. 10) call energ22

if (.not. noprin) then
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
end if

deallocate(w)
deallocate(tmat)
deallocate(vecnow)
deallocate(vecnew)

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
use mod_himatrix, only: mxma
use mod_hivector, only: matmov
use mod_hiblas, only: dscal, dcopy, daxpy_wrapper
use mod_hipot, only: ground
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
use mod_hiblas, only: dcopy
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
!  maxmgv:      find maximum (absolute value) element in a vector
! -----------------------------------------------------------------------
use mod_hivector, only: maxmgv
use mod_hiblas, only: dcopy, dgemm
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
use mod_cojhld, only: jout1 => jhold ! jout1(1)
use mod_coehld, only: jout2 => eholdint ! jout2(1)
use mod_coinhl, only: jlev => inhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_coz, only: sreal1 => z_as_vec ! sreal1(1)
use mod_cow, only: sreal2 => w_as_vec ! sreal2(1)
use mod_cozmat, only: simag1 => zmat_as_vec ! simag1(1)
use mod_hismat, only: sread, rdhead
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_hiutil, only: gennam
use mod_hitypes, only: bqs_type
use constants, only: zero
use mod_hiiolib1, only: openf
use mod_hiiolib1, only: closf
implicit none
character*(*), intent(in) :: fname1
character*(*), intent(in) :: fname2
integer, intent(in) :: ienerg
integer, intent(in) :: iprint
real(8), intent(out) :: acc
real(8), intent(out) :: accmx
real(8), intent(in) :: thrs
integer, intent(out) :: imx
integer, intent(out) :: jmx
integer, intent(out) :: ityp
type(bqs_type) :: row_bqs
type(bqs_type) :: pack1
type(bqs_type) :: pack2
real(8) :: erabs, ered1, ered2, ermabs, ermrel, errel
integer :: i, idif, ierr, ij, im
integer :: jfinal, jfirst, jlpar1, jlpar2, jm, jtot1, jtot2, jtotd
integer :: length, lenx1, lenx2
integer :: n, ncol, nlevel, nlevop, nnout1, nnout2, nopen1, nopen2, nu1, nu2, nud, numax, numin
real(8) :: rmu1, rmu2
character*20 cdate1,cdate2
character*40 xnam1,xnam2
character*48 potnam1,potnam2,label1,label2

logical existf,csflg1,csflg2,flghf1,flghf2,flgsu1,flgsu2
logical twoml1,twoml2,nucr1,nucr2

!
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
                  row_bqs, pack1, &
                  1, mmax, nopen1, ierr)
if(ierr.eq.-1) goto 200
if(ierr.lt.-1) then
  write(6,35) xnam1
35   format(' ERROR READING FILE ',(a))
  goto 200
end if
nopen2 = 0
call sread (0,sreal2, simag2, jtot2, jlpar2, nu2, &
                  row_bqs, pack2, &
                  2, mmax, nopen2, ierr)
if(ierr.eq.-1) goto 200
if(ierr.lt.-1) then
  write(6,35) xnam2
  goto 200
end if
if(nu1.ne.nu2) idif=idif+1
if(jtot1.ne.jtot2) idif=idif+1
if(jlpar1.ne.jlpar2) idif=idif+1
if(pack1%length /= pack2%length) idif=idif+1
do 60 i=1,pack1%length
if(pack1%jq(i) /= pack2%jq(i)) idif=idif+1
if(pack1%lq(i) /= pack2%lq(i)) idif=idif+1
60 if(pack1%inq(i) /= pack2%inq(i)) idif=idif+1
if(idif.ne.0) then
  write(6,70) jtot1,jtot2
70   format(/' PARAMETERS NOT EQUAL FOR JTOT1=',i3,'  JTOT2=',i3)
  goto 200
end if
ncol=pack1%length
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
  call compar(sreal1,sreal2,mmax,pack1%length,nopen1, &
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
  call compar(simag1,simag2,mmax,pack1%length,nopen1, &
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
  call compt2(sreal1,simag1,sreal2,simag2,mmax,pack1%length,nopen1, &
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
end module mod_hibrid1
