!     This file contains various modules that replace common blocks in
!     hibridon.
!
!     Physical constants
module constants
implicit none
real(8), parameter :: pi=dacos(-1d0)
!     Below are constants used in Hibridon 4.4
real(8), parameter :: econv=219474.6d0, xmconv=5.485930d-4, &
     ang2c=0.280002846d0
!
!     All physical constants below are adapted from 2010 values from
!     http://physics.nist.gov/cuu/Constants/index.html
!
!     econv is the conversion from hartree to cm-1
!     The value of econv is the rydberg constant in m-1 devided by 50
!$$$      real(8), parameter :: econv=10973731.568539/50d0
!     xmconv converts mass in atom units (electron mass) to atom mass
!     unit (1/12 of the mass of 12C)
!$$$      real(8), parameter :: xmconv=5.4857990946d-4
!     ang2c converts bohr^2 to angstrom^2
!     The value of ang2 is (10^10 Bohr Radius)^2
!$$$      real(8), parameter :: ang2c=0.52917721092d0**2
end module constants

! file units used in hibridon
module funit
implicit none

   enum, bind( C )
   enumerator ::  &
      FUNIT_CS             =  1, &   ! cross sections
      FUNIT_SAV            =  3, &   ! <job-name>.sav
      FUNIT_STDOUT         =  6, &   ! standard output
      FUNIT_INP            =  8, &   ! input file (*.inp)
      FUNIT_OUT            =  9, &   ! output file (*.out)
      FUNIT_TRANS_MAT      = 10, &   ! temporary storage for transformation matrices (t<unit><pid>.tmp)
      FUNIT_QUAD_MAT       = 11, &   ! temporary storage for quadrature matrices (t<unit><pid>.tmp)
      FUNIT_CHANNEL_PARAMS = 12, &   ! temporary storage for channel parameters (t<unit><pid>.tmp)
      FUNIT_WFU            = 22, &   ! <job-name>.wfu
      FUNIT_PCS_START      = 25, &   ! <job-name><energy_index>.pcs files (partial cross sections; one for each energy)
      FUNIT_APCS_START     = 35, &   ! accumulation of cs partial cross sections at each projection index; one for each energy. (t<unit><pid>.tmp)
      FUNIT_SMT_START      = 45, &   ! <job-name><energy_index>.smt storage of selected s-matrix elements; one for each energy
      FUNIT_ICS_START      = 70      ! <job-name><energy_index>.ics files (integral cross sections; one for each energy)
   end enum


end module funit

! system independent parameters
module mod_si_params
  !  ncode is the number of bcod's
  !  bcod stores hibridon's commands
  !  fcod stores logical flags (length = lcode)
  integer, parameter :: ncode = 38
  integer, parameter :: lcode = 28
  !  iicode is the number of integer pcod's
  !  ircode is the number of real pcod's
  integer, parameter :: iicode = 10
  integer, parameter :: ircode = 9
  integer, parameter :: icode = iicode+ircode

contains
  subroutine set_param_names(boundc, param_names, param_names_size)
    !  subroutine to change param_names's for bound state or scattering
    use mod_hiparcst, only: IPAR_COUNT
    use rpar_enum
    implicit none
    logical, intent(in) :: boundc
    integer, intent(in) :: param_names_size  ! size of param_names array
    character*8, intent(out) :: param_names(param_names_size)  ! array containing the name of each parameter (old name: pcod)
    if (boundc) then
      param_names(IPAR_COUNT + RPAR_BOUND_R1)      = 'R1'
      param_names(IPAR_COUNT + RPAR_BOUND_R2)      = 'R2'
      param_names(IPAR_COUNT + RPAR_BOUND_C)       = 'C' 
      param_names(IPAR_COUNT + RPAR_BOUND_SPAC)    = 'SPAC' 
      param_names(IPAR_COUNT + RPAR_BOUND_DELR)    = 'DELR' 
      param_names(IPAR_COUNT + RPAR_BOUND_HSIMP)   = 'HSIMP' 
      param_names(IPAR_COUNT + RPAR_BOUND_EIGMIN)  = 'EIGMIN' 
      param_names(IPAR_COUNT + RPAR_BOUND_TOLAI)   = 'TOLAI' 
    else
      param_names(IPAR_COUNT + RPAR_SCAT_FSTFAC)  = 'FSTFAC'
      param_names(IPAR_COUNT + RPAR_SCAT_RINCR)   = 'RINCR'
      param_names(IPAR_COUNT + RPAR_SCAT_RCUT)    = 'RCUT'
      param_names(IPAR_COUNT + RPAR_SCAT_RENDAI)  = 'RENDAI'
      param_names(IPAR_COUNT + RPAR_SCAT_RENDLD)  = 'RENDLD'
      param_names(IPAR_COUNT + RPAR_SCAT_RSTART)  = 'RSTART'
      param_names(IPAR_COUNT + RPAR_SCAT_SPAC)    = 'SPAC'
      param_names(IPAR_COUNT + RPAR_SCAT_TOLAI)   = 'TOLAI'
    endif
    return
  end

end module mod_si_params

module mod_hinput_state
  use mod_hiparcst, only: LPAR_COUNT

  logical :: batch
  data batch /.false./

  !  lindx is pointer from fcod order to order in common block colpar
  ! graffy: colpar and fcod both contain the 28 logical parameters but in a different order, thus requiring a remapping through lindx. Why not simply having the same order, by making fcod match colpar?
  integer :: lindx(LPAR_COUNT)

end module mod_hinput_state
