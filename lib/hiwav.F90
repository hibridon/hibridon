#include "assert.h"
module mod_wave  ! mod_wav is a modern replacement of common block cowave
  use mod_assert, only: fassert
  implicit none
  character, parameter :: wfu_format_version = char(3)  ! the version number of the wfu file format (graffy note on 13/09/2022 : unable to find documentation on the 6th byte of the wfu file, I interpreted it as a format version; its value was 2)
  integer(8), parameter :: ipos2_location = 9  ! location of the ipos2 field relative to the start of the wfu file

  type, public :: wfu_file_type
    integer :: irec    ! record number of last written g(a,b) matrix
    integer :: ifil    ! file unit used for wfu file (is it always FUNIT_WFU ?)
    integer :: nchwfu  ! number of channels in wfu file
    integer(8) :: ipos2
    integer(8) :: ipos3
    integer :: nrlogd  ! number of log derivatives
    integer(8) :: iendwv  ! number of bytes in wfu file (integer(8) to allow for wfu files bigger than 2gb)
    integer :: inflev  ! 0 or 1
  contains
   
  end type wfu_file_type

contains
function get_wfu_rec1_length(nchwfu)
  integer, intent(in) :: nchwfu  ! number of channels in wfu file
  integer(8) :: get_wfu_rec1_length
  !     The three variables below are used to determine the (machine
  !     dependent) size of the built-in types
  character :: char_t
  integer :: int_t
  real(8) :: dble_t

  !     4 chars '\128'      'W'     'F'     'U'
  ! +   4 chars '<inflev>' '<wfu_version>'    '\0'    '\0'
  ! +   8 bytes (ipos2)
  ! +   8 bytes (ipos3)
  ! +   4 bytes (nrlogd)
  ! +  20 bytes (cdate)
  ! +  48 bytes (label)
  ! +  48 bytes (potnam)
  ! -------------------
  ! = 144 bytes
  integer, parameter :: num_chars_in_r1 = 144 ! last element of header is potnam
  integer, parameter :: num_ints_in_r1 = 4 ! jtot, jlpar, nu, nch
  integer, parameter :: num_flags_in_r1 = 3 ! csflag, flaghf, photof
  integer, parameter :: num_quadchars_in_r1 = 3 ! 0000, 'ENDW', 'FUR\1'
  integer, parameter :: num_ints_like_in_r1 = num_ints_in_r1 + num_flags_in_r1 + num_quadchars_in_r1
  integer, parameter :: num_doubles_in_r1 = 4 ! ered, rmu, rstart, rendld
  integer, parameter :: num_ints_per_channel = 3 ! jq(channel_index), lq(channel_index), inq(channel_index)
  integer, parameter :: num_doubles_per_channel = 1 ! eint(channel_index)
  get_wfu_rec1_length = num_chars_in_r1 * sizeof(char_t) &
       + (num_ints_like_in_r1 + num_ints_per_channel * nchwfu) * sizeof(int_t) &
       + (num_doubles_in_r1 + num_doubles_per_channel * nchwfu) * sizeof(dble_t)
end function

!     Length (in bytes) for each record written by the LOGD propagator
function get_wfu_logd_rec_length(nchwfu, inflev)
  integer, intent(in) :: nchwfu  ! number of channels in wfu file
  integer, intent(in) :: inflev  ! 0 or 1
  integer(8) :: get_wfu_logd_rec_length
  !     The three variables below are used to determine the (machine
  !     dependent) size of the built-in types
  character :: char_t
  integer :: int_t
  real(8) :: dble_t
  integer, parameter :: char_size = int(sizeof(char_t), kind(int_t))
  integer, parameter :: dbl_size = int(sizeof(dble_t), kind(int_t))

  get_wfu_logd_rec_length = 0
    

  if (inflev .eq. 0) then
    get_wfu_logd_rec_length = get_wfu_logd_rec_length + dbl_size ! r-h
    get_wfu_logd_rec_length = get_wfu_logd_rec_length + dbl_size ! r
    get_wfu_logd_rec_length = get_wfu_logd_rec_length + nchwfu * dbl_size ! w (vector)
    get_wfu_logd_rec_length = get_wfu_logd_rec_length + (nchwfu ** 2) * dbl_size ! w (matrix)

    get_wfu_logd_rec_length = get_wfu_logd_rec_length + 8 * char_size  ! 'ENDWFUR' + <irec>
  else if (inflev .eq. 1) then
  else
    stop 'unexpected value for inflev'
  end if




end function

!     Length (in bytes) for each record written by the Airy propagator
function get_wfu_airy_rec_length(nchwfu, inflev)
  integer, intent(in) :: nchwfu  ! number of channels in wfu file
  integer, intent(in) :: inflev  ! 0 or 1
  integer(8) :: get_wfu_airy_rec_length
  !     The three variables below are used to determine the (machine
  !     dependent) size of the built-in types
  character :: char_t
  integer :: int_t
  real(8) :: dble_t
  integer, parameter :: char_size = int(sizeof(char_t), kind(int_t))
  integer, parameter :: dbl_size = int(sizeof(dble_t), kind(int_t))

  get_wfu_airy_rec_length = 0
  get_wfu_airy_rec_length = get_wfu_airy_rec_length + dbl_size  ! -rlast
  get_wfu_airy_rec_length = get_wfu_airy_rec_length + dbl_size  ! drnow
  get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! eigold

  if (inflev .eq. 0) then
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + (nchwfu ** 2) * dbl_size  ! z(nchwfu, nchwfu) (real part)
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + (nchwfu ** 2) * dbl_size  ! z(nchwfu, nchwfu) (imaginary part)
    ! propagators
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! y1
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! y2
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! y4
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! gam1
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! muab
  else if (inflev == 1) then
  else
    stop 'unexpected value for inflev'
  end if

  get_wfu_airy_rec_length = get_wfu_airy_rec_length + 8 * char_size  ! 'ENDWFUR' + <irec>

end function

!     ------------------------------------------------------------------
function iwavsk(wfu_file, irecr)
!     ------------------------------------------------------------------
!     Function to return offset of wfu file for recrod #irec (stream IO)
!
!     Written by: Qianli Ma
!     Latest revision: 20-apr-2012
!
!     This function needs nchwfu, ipos2 and ipos3 from the mod_wave
!     module.  These variables are set by waverd.
!
!     The stream IO counterpart for `dbrr(,,,irec)` is `read
!     (,pos=wavesk(irec))...
!     ------------------------------------------------------------------
implicit none
type(wfu_file_type), intent(in) :: wfu_file
integer, intent(in) :: irecr
integer(8) :: iwavsk  ! 1-based offset in the wfu file, in bytes

integer(8) :: lr1  ! length of record 1 in bytes
integer(8) :: lrlogd  ! length of a logd record
integer(8) :: lrairy  ! length of an airy record
!
if (irecr .le. 0) then
   iwavsk = -1
   goto 100
end if
if (irecr .eq. 1) then
   iwavsk = 1
   goto 100
end if
if (irecr .eq. 2) then
   iwavsk = wfu_file%ipos2
   goto 100
end if
if (irecr .eq. 3) then
   iwavsk = wfu_file%ipos3
   goto 100
end if
!     Length for record 1 (at the beginning of the file)
lr1 = get_wfu_rec1_length(wfu_file%nchwfu)
!     Length for each block written by the Airy and LOGDpropagator
lrairy = get_wfu_airy_rec_length(wfu_file%nchwfu, wfu_file%inflev)
lrlogd = get_wfu_logd_rec_length(wfu_file%nchwfu, wfu_file%inflev)
!
if ((irecr - 3) .le. wfu_file%nrlogd) then
!     within the logd range of the file
   iwavsk = lr1 + (irecr - 4) * lrlogd + 1
   goto 100
else
!     airy range of the file
   iwavsk = lr1 + wfu_file%nrlogd * lrlogd &
        + (irecr - 4 - wfu_file%nrlogd) * lrairy + 1
   goto 100
end if
!     should never reach here if function called properly
write (0, *) '*** OOPS! ERROR SEEKING WFU FILE. ABORT.'
call exit()
100 continue
ASSERT(iwavsk > 0)
end
!
! -------------------------------------------------------------------------
subroutine wavewr(jtot,jlpar,nu,nch,rstart,rendld, bqs, wfu_file)
! -------------------------------------------------------
!  subroutine to write initial header information on wavefunction file
!  (file jobname.WFU, logical unit 22), unit is opened in subroutine openfi
!     written by:  millard alexander
!     latest revision:  11-dec-2011
!     common blocks amat and bmat are used to store real and
!     imaginary parts of asymptotic wavefunction (only used in
!     read of wavefunction from saved file)
!
!     Major revision: 16-mar-2012 by Q. Ma
!     Use stream I/O for smaller file size and better compatibility
!
!     current revision: 20-apr-2012 by q. ma
!
! -------------------------------------------------------
#define AMAT_AS_VEC_METHOD_DISTINCT 1
#define AMAT_AS_VEC_METHOD_POINTER 2
#define AMAT_AS_VEC_METHOD_NOVEC 3
#define AMAT_AS_VEC_METHOD AMAT_AS_VEC_METHOD_DISTINCT
use mod_coeint, only: eint
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
use, intrinsic :: ISO_C_BINDING
use mod_coamat, only: amat ! amat(25)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
use mod_coamat, only: amat ! amat(25)
#endif
use mod_par, only: csflag, flaghf, wrsmat, photof
use funit
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_ered, only: ered, rmu
use mod_hiutil, only: dater
use mod_hitypes, only: bqs_type
implicit none
integer, intent(in) :: jtot
integer, intent(in) :: jlpar
integer, intent(in) :: nu
integer, intent(in) :: nch
real(8), intent(in) :: rstart
real(8), intent(in) :: rendld
type(bqs_type), intent(in) :: bqs
type(wfu_file_type), allocatable, intent(out) :: wfu_file
character*20 :: cdate

integer :: i
integer(8) :: end_of_rec1_pos
!
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
real, pointer :: amat_as_vec(:)
#endif
ASSERT( bqs%length == nch )
ASSERT(.not. allocated(wfu_file))
allocate(wfu_file)
wfu_file%ifil = FUNIT_WFU

wfu_file%nchwfu = nch
wfu_file%ipos2 = -1
wfu_file%ipos3 = -1
wfu_file%nrlogd = 0
!     Mark the position of the EOF of the WFU file in order to by pass
!     (likely) a bug in the intel compiler that INQUIRE does not return
!     the proper offset
wfu_file%iendwv = 1
!     Write magic number
write (wfu_file%ifil, err=950) char(128), 'WFU'

if (wrsmat) then
   write (wfu_file%ifil, err=950) char(0), wfu_format_version, char(0), char(0)
else
   write (wfu_file%ifil, err=950) char(1), wfu_format_version, char(0), char(0)
end if
!
write (wfu_file%ifil, err=950) wfu_file%ipos2, wfu_file%ipos3, wfu_file%nrlogd
call dater(cdate)
write (wfu_file%ifil, err=950) cdate, label, potnam
!     Four zero-bytes for alignment / C struct compatibility
write (wfu_file%ifil, err=950) char(0), char(0), char(0), char(0)
!
write (wfu_file%ifil, err=950) jtot, jlpar, nu, nch, csflag, flaghf, photof
write (wfu_file%ifil, err=950) ered, rmu, rstart, rendld
!
write (wfu_file%ifil, err=950) (bqs%jq(i), i=1, nch), (bqs%lq(i), i=1, nch), &
     (bqs%inq(i), i=1, nch)
write (wfu_file%ifil, err=950) (eint(i), i=1, nch)
!
write (wfu_file%ifil, err=950) 'ENDWFUR', char(1)
inquire(wfu_file%ifil, pos=end_of_rec1_pos)
ASSERT(end_of_rec1_pos == (get_wfu_rec1_length(wfu_file%nchwfu) + 1))
!
wfu_file%iendwv = wfu_file%iendwv + get_wfu_rec1_length(wfu_file%nchwfu)
wfu_file%irec=3  ! irec=2 and irec=3 are reserved and their position in the file are stored in ipos2 and ipos3
return

950 write (0, *) '*** ERROR WRITING WFU FILE. ABORT.'
call exit
return

end subroutine wavewr
!
!     ------------------------------------------------------------------
!     reads header file for wavefunction (wfu file)
subroutine waverd(jtot,jlpar,nu,nch,npts,nopen,nphoto,jflux, &
     rstart,rendld,rinf, rbesself, bqs, wfu_file)
use mod_coeint, only: eint
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_DISTINCT)
use mod_coamat, only: amat => psir ! amat(25) psir(nopen, nopen)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
use, intrinsic :: ISO_C_BINDING
use mod_coamat, only: amat ! amat(25)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
use mod_coamat, only: amat ! amat(25)
#endif
use mod_cobmat, only: bmat => psii ! bmat(25), here bmat is used as a vector 
use mod_cotq1, only: dpsir ! dpsir(25)
use mod_cotq2, only: dpsii ! dpsii(25)
use mod_coisc1, only: isc1 ! isc1(25)
use mod_cosc1, only: pk => sc1 ! sc1(10)  ! pk (asymptotic wavevectors)
use mod_cow, only: w => w_as_vec ! w(25)
use mod_cozmat, only: zmat => zmat_as_vec ! zmat(25)
use mod_par, only: csflag, flaghf, photof
use funit
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_ered, only: ered, rmu
use mod_hivector, only: dset
use mod_hitypes, only: rbesself_type, bqs_type
use constants, only: zero
implicit none
integer, intent(out) :: jtot
integer, intent(out) :: jlpar
integer, intent(out) :: nu
integer, intent(out) :: nch
integer, intent(out) :: npts
integer, intent(out) :: nopen
integer, intent(out) :: nphoto
integer, intent(in) :: jflux
real(8), intent(out) :: rstart
real(8), intent(out) :: rendld
real(8), intent(out) :: rinf
type(rbesself_type), intent(out) :: rbesself
type(bqs_type), intent(out) :: bqs
type(wfu_file_type), allocatable, intent(out) :: wfu_file

character*48 :: oldlab, oldpot
character*20 :: olddat

character :: csize8(8), csize4(4)
integer :: i
integer :: nopsq
integer :: nrecs
! integer, parameter :: izero=0
!
ASSERT(.not. allocated(wfu_file))
allocate(wfu_file)
wfu_file%ifil = FUNIT_WFU ! the wfu file is expected to be open using this unit
!     Read the magic number (from the start of the file)
read (wfu_file%ifil, pos=1, end=900, err=950) csize8
wfu_file%inflev = ichar(csize8(5))
if (csize8(6) /= wfu_format_version) then
  write (0,'(a,i3,a,i3,a)') '*** UNHANDLED VERSION OF WFU FORMAT : ', ichar(csize8(6)), ' (THIS VERSION OF HIBRIDON ONLY HANDLES WFU FORMAT VERSION ',  ichar(wfu_format_version),'). ABORT.'
  call exit()
end if
!
read (wfu_file%ifil, end=900, err=950) wfu_file%ipos2, wfu_file%ipos3, wfu_file%nrlogd
!
read (wfu_file%ifil, end=900, err=950) olddat, oldlab, oldpot
label = oldlab
potnam = oldpot
!     Read four zero bytes
read (wfu_file%ifil, end=900, err=950) csize4
!
read (wfu_file%ifil, end=900, err=950) jtot, jlpar, nu, nch, csflag, &
     flaghf, photof
!     nchwfu is used in locating the position for records
wfu_file%nchwfu = nch
read (wfu_file%ifil, end=900, err=950) ered, rmu, rstart, rendld
write (6, 245) olddat
if (jflux .ne. 0) write (FUNIT_FLX, 245) olddat
if (jflux .eq. 0) write (2, 245) olddat
245 format('    FROM CALCULATION ON: ',(a))
if (jflux .ne. 0) write (FUNIT_FLX, 250) oldlab
if (jflux .eq. 0) write (2, 250) oldlab
write (6, 250) oldlab
250 format('    INITIAL JOB LABEL: ', (a))
if (jflux .ne. 0) write (FUNIT_FLX, 251) oldpot
if (jflux .eq. 0) write (2, 251) oldpot
write (6, 251) oldpot
251 format('    INITIAL POT NAME: ', (a))
!
!     Read in channel labels
call bqs%init(nch)
read (wfu_file%ifil, end=900, err=950) (bqs%jq(i), i=1, nch), &
     (bqs%lq(i), i=1, nch), (bqs%inq(i), i=1, nch), &
     (eint(i), i=1, nch)
bqs%length = nch
!
! start reading in information from record 2 here
read (wfu_file%ifil, end=900, err=950, pos=iwavsk(wfu_file, 2)) nrecs, nopen, &
     nphoto, rinf
npts = nrecs - 3
! read in wavevectors, bessel functions j, j', n, n'
! first initialize to zero for all channels
call rbesself%init(nch)
call dset(nch,zero,pk,1)  ! pk (asymptotic wavevectors)
call dset(nch,zero,rbesself%fj,1)  ! fj
call dset(nch,zero,rbesself%fpj,1)  ! fpj
call dset(nch,zero,rbesself%fn,1)  ! fn
call dset(nch,zero,rbesself%fpn,1)  ! fpn
read (wfu_file%ifil, end=900, err=950) (pk(i), i=1, nopen), &
     (rbesself%fj(i), i=1, nopen), (rbesself%fpj(i), i=1, nopen), &
     (rbesself%fn(i), i=1, nopen), (rbesself%fpn(i), i=1, nopen)
rbesself%length = nopen
nopsq = nopen ** 2
! read in sreal and simag, store in w and zmat
read (wfu_file%ifil, end=900, err=950) (w(i), i=1, nopsq), &
     (zmat(i), i=1, nopsq)
if (photof) then
! read in number of initial photodissociation states
!        call dbri(mphoto,1,ifil,REC_LAST_USED)
!        nphoto=mphoto
! read in real part of photodissociation amplitude
! overlay sreal which is not needed for photodissociation problem
   read (wfu_file%ifil, end=900, err=950) (w(i), i=1, nphoto * nopen)
! read in imaginary part of photodissociation amplitude
! overlay simag which is not needed for photodissociation problem
   read (wfu_file%ifil, end=900, err=950) (zmat(i), i=1, nphoto * nopen)
endif
! read in channel packing list and real and imaginary parts
! of scattering wavefunction and derivative
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_DISTINCT)
read (wfu_file%ifil, end=900, err=950, pos=iwavsk(wfu_file, 3)) &
     (isc1(i), i=1, nopen), (amat(i), i=1, nopsq), &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
! amat_as_vec is a view of the matrix amat(nopen, nopen) as a vector(nopen*nopen)
call C_F_POINTER (C_LOC(amat), amat_as_vec, [nopsq])
read (wfu_file%ifil, end=900, err=950, pos=iwavsk(wfu_file, 3)) &
     (isc1(i), i=1, nopen), (amat_as_vec(i), i=1, nopsq), &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
read (wfu_file%ifil, end=900, err=950, pos=iwavsk(wfu_file, 3)) &
     (isc1(i), i=1, nopen), amat, &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
wfu_file%irec = 3
return
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE. ABORT.'
call exit
return
end

subroutine write_airy_record(wfu_file, rlast, drnow, nch, eigold, writs, nmax, z, vecnow, y1, y2, y4, gam1, lmuab)
implicit none
type(wfu_file_type), allocatable, intent(inout) :: wfu_file
real(8), intent(in) :: rlast
real(8), intent(in) :: drnow
integer, intent(in) :: nch
real(8), intent(in) :: eigold(nch)
logical, intent(in) :: writs
integer, intent(in) :: nmax
real(8), intent(in) :: z(nch*nmax)  ! this is G(n-1,n) in the local basis
real(8), intent(in) :: vecnow(nch*nmax)  ! vecnow is transformation from free basis into local basis in first interval
real(8), intent(in) :: y1(nch)
real(8), intent(in) :: y2(nch)
real(8), intent(in) :: y4(nch)
real(8), intent(in) :: gam1(nch)
real(8), intent(in) :: lmuab(nch)

integer(8) :: lrairy ! length of an airy record in bytes
integer :: i, ich, icol

  ASSERT(allocated(wfu_file))
  wfu_file%irec = wfu_file%irec + 1
  write (wfu_file%ifil, err=950) -rlast, drnow
  !     Adiabatic energies
  write (wfu_file%ifil, err=950) (eigold(i), i=1, nch)
  !     The following information will not be written if writs set to F
  if (writs) then
    icol = 1
    do ich = 1, nch
       write (wfu_file%ifil, err=950) (z(icol - 1 + i), i=1, nch)
       icol = icol + nmax
    end do
    icol = 1
    do ich = 1, nch
       write (wfu_file%ifil, err=950) (vecnow(icol - 1 + i), i=1, nch)
       icol = icol + nmax
    end do
    !
    write (wfu_file%ifil, err=950) (y1(i), i=1, nch), (y2(i), i=1, nch), &
         (y4(i), i=1, nch), (gam1(i), i=1, nch), &
         (lmuab(i), i=1, nch)
    lrairy = get_wfu_airy_rec_length(wfu_file%nchwfu, 0)
  else
    lrairy = get_wfu_airy_rec_length(wfu_file%nchwfu, 1)
  end if
  !
  write (wfu_file%ifil, err=950) 'ENDWFUR', char(mod(wfu_file%irec, 256))
  wfu_file%iendwv = wfu_file%iendwv + lrairy
  return

950 write (0, *) ' *** ERROR WRITING WFU FILE (AIRY). ABORT.'
  call exit()
end subroutine write_airy_record

subroutine write_logd_record(wfu_file, r, h, w, nch, nmax)
implicit none
type(wfu_file_type), allocatable, intent(inout) :: wfu_file
real(8), intent(in) :: r
real(8), intent(in) :: h
real(8), intent(in) :: w(nch*nmax)  ! contains the matrix g(a,m)g(m,b)=g(a,b)

integer, intent(in) :: nch
integer, intent(in) :: nmax

integer :: i, ich, icol

  wfu_file%irec = wfu_file%irec + 1
!     nrlogd is the number of LOGD records - used to seek the wfu file
  wfu_file%nrlogd = wfu_file%nrlogd + 1
  write (wfu_file%ifil, err=950) r - h, r, (w(i), i=1, nch)
  icol = 1
  do ich = 1, nch
    write (wfu_file%ifil, err=950) (w(icol - 1 + i), i=1, nch)
    icol = icol + nmax
  end do
  write (wfu_file%ifil, err=950) 'ENDWFUR', char(mod(wfu_file%irec, 256))
  wfu_file%iendwv = wfu_file%iendwv + get_wfu_logd_rec_length(wfu_file%nchwfu, 0)

  return
950 write (0, *) ' *** ERROR WRITING WFU FILE (write_logd_record). ABORT'
  call exit()

end subroutine write_logd_record


subroutine write_record2_header(wfu_file, nopen, nphoto, r, pk, fj, fpj, fn, fpn)
implicit none
type(wfu_file_type), allocatable, intent(inout) :: wfu_file
integer, intent(in) :: nopen
integer, intent(in) :: nphoto
real(8), intent(in) :: r
real(8), intent(in) :: pk(nopen)  ! wavevectors for open channels
real(8), intent(in) :: fj(nopen)  ! ricatti-bessel functions
real(8), intent(in) :: fpj(nopen)  ! ricatti-bessel functions
real(8), intent(in) :: fn(nopen)  ! ricatti-bessel functions
real(8), intent(in) :: fpn(nopen)  ! ricatti-bessel functions

integer int_t  ! used to determine the (machine dependent) size of built-in types
double precision dble_t    ! used to determine the (machine dependent) size of built-in types

integer :: i
  ! save in record 2 of direct access file:
  ! 1. number of open channels,
  ! 2. total number of records on wavefunction
  ! 3. asymptotic interparticle separation
  ! 4. wavevectors of open channels
  ! 5. bessel functions and derivatives for open channels
  !
  !     As of the 12.1.3 version of ifort, a big hole is likely to be
  !     created here in the wfu file.  Most likely to be a bug of ifort.
  !
  !     If the bug is fixed in the future, use the following line and
  !     remove ALL references to iendwv
  !$$$         inquire (ifil, pos=ipos2)
  !
  wfu_file%ipos2 = wfu_file%iendwv
  write (wfu_file%ifil, err=950, pos=ipos2_location) wfu_file%ipos2, wfu_file%ipos3, wfu_file%nrlogd
  write (wfu_file%ifil, err=950, pos=wfu_file%ipos2) wfu_file%irec, nopen, nphoto, &
       r, (pk(i), i=1, nopen), (fj(i), i=1, nopen), &
       (fpj(i), i=1, nopen), (fn(i), i=1, nopen), &
       (fpn(i), i=1, nopen)
  wfu_file%iendwv = wfu_file%iendwv + 3 * int(sizeof(int_t), kind(int_t)) + &
       (5 * nopen + 1) * int(sizeof(dble_t), kind(int_t))
  return

  950 write (0, *) '*** ERROR WRITING WFU FILE (write_record2_header). ABORT'
  call exit()
end subroutine write_record2_header

subroutine write_record2_scat_data(wfu_file, nopen, nmax, sr, si)
implicit none
type(wfu_file_type), allocatable, intent(inout) :: wfu_file
integer, intent(in) :: nopen
integer, intent(in) :: nmax
real(8), intent(in) :: sr(nmax, nmax)  ! s-matrix real part
real(8), intent(in) :: si(nmax, nmax)  ! s-matrix imaginary part

integer :: i, icol

integer char_t  ! used to determine the (machine dependent) size of built-in types
double precision dble_t    ! used to determine the (machine dependent) size of built-in types

  ! save real and imaginary part of s-matrix in record 2 of direct access file
  do icol=1, nopen
     write (wfu_file%ifil, err=950) (sr(i, icol), i=1, nopen)
  end do
  do icol=1, nopen
     write (wfu_file%ifil, err=950) (si(i, icol), i=1, nopen)
  end do
  write (wfu_file%ifil, err=950) 'ENDWFUR', char(2)
  wfu_file%iendwv = wfu_file%iendwv + 8 * sizeof(char_t) + (2 * nopen ** 2) * sizeof(dble_t)
  return

  950 write (0, *) '*** ERROR WRITING WFU FILE (write_record2_scat_data). ABORT'
  call exit()
end subroutine write_record2_scat_data


subroutine write_record2_photo_data(wfu_file, nphoto, nopen, nmax, sr, si)
implicit none
type(wfu_file_type), allocatable, intent(inout) :: wfu_file
integer, intent(in) :: nphoto
integer, intent(in) :: nopen
integer, intent(in) :: nmax
real(8), intent(in) :: sr(nmax, nmax)  ! s-matrix real part
real(8), intent(in) :: si(nmax, nmax)  ! s-matrix imaginary part

integer :: jrow, jcol

integer char_t  ! used to determine the (machine dependent) size of built-in types
double precision dble_t    ! used to determine the (machine dependent) size of built-in types

  ! call dbwi(nphoto,1,ifil,REC_LAST_USED)
  do jrow = 1, nphoto
    write (wfu_file%ifil, err=950) (sr(jrow, jcol), jcol=1, nopen)
  end do
  ! NB there were 2*nphoto!
  do jrow = 1, nphoto
    write (wfu_file%ifil, err=950) (si(jrow, jcol), jcol=1, nopen)
  end do

  write (wfu_file%ifil, err=950) 'ENDWFUR', char(2)
  wfu_file%iendwv = wfu_file%iendwv + 8 * sizeof(char_t) &
        + (2 * nopen ** 2) * sizeof(dble_t)
  return

  950 write (0, *) '*** ERROR WRITING WFU FILE (write_record2_photo_data). ABORT'
  call exit()
end subroutine write_record2_photo_data


end module mod_wave