#include "assert.h"
module mod_wave  ! mod_wav is a modern replacement of common block cowave
  implicit none
  integer :: irec    ! record number of last written g(a,b) matrix
  integer :: ifil    ! file unit used for wfu file
  integer :: nchwfu  ! number of channels in wfu file
  integer(8) :: ipos2
  integer(8) :: ipos3
  integer :: nrlogd  ! number of log derivatives
  integer(8) :: iendwv  ! number of bytes in wfu file (integer(8) to allow for wfu files bigger than 2gb)
  integer :: inflev  ! 0 or 1
  character, parameter :: wfu_format_version = char(3)  ! the version number of the wfu file format (graffy note on 13/09/2022 : unable to find documentation on the 6th byte of the wfu file, I interpreted it as a format version; its value was 2)
  integer(8), parameter :: ipos2_location = 9  ! locatioon of the ipos2 field relative to the start of the wfu file
  contains
function get_wfu_rec1_length(nchwfu)
  integer, intent(in) :: nchwfu  ! number of channels in wfu file
  integer :: get_wfu_rec1_length
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
  integer, parameter :: int_size = int(sizeof(int_t), kind(int_t))
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
  integer, parameter :: int_size = int(sizeof(int_t), kind(int_t))
  integer, parameter :: dbl_size = int(sizeof(dble_t), kind(int_t))

  get_wfu_airy_rec_length = 0
  get_wfu_airy_rec_length = get_wfu_airy_rec_length + dbl_size  ! -rlast
  get_wfu_airy_rec_length = get_wfu_airy_rec_length + dbl_size  ! drnow
  get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! eigold

  if (inflev .eq. 0) then
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + (nchwfu ** 2) * dbl_size  ! z(nchwfu, nchwfu)
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + (nchwfu ** 2) * dbl_size  ! vecnow(nchwfu, nchwfu)
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! y1
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! y2
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! y4
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! gam1
    get_wfu_airy_rec_length = get_wfu_airy_rec_length + nchwfu * dbl_size  ! sc10
  else if (inflev == 1) then
  else
    stop 'unexpected value for inflev'
  end if

  get_wfu_airy_rec_length = get_wfu_airy_rec_length + 8 * char_size  ! 'ENDWFUR' + <irec>

end function

end module mod_wave